import os
import re
import sys
import yaml
import copy
import argparse
import itertools
import time
import subprocess
import multiprocessing as mp


class YMLIncludeLoader(yaml.SafeLoader): 
    """Custom yaml loading to support including config files. Use `!include (file)` to insert content of `file` at that position."""
    
    def __init__(self, stream):
        super(YMLIncludeLoader, self).__init__(stream)
        self._root = os.path.split(stream.name)[0]

    def include(self, node):
        filename = os.path.join(self._root, self.construct_scalar(node))
        with open(filename, 'r') as f:
            return yaml.load(f, YMLIncludeLoader)

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def runBamboo(idx,cmd,queue):
    #cmd = 'python3 -u test_run.py'
    process = subprocess.Popen(cmd.split(),universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    # Poll process for new output until finished
    while True:
        nextline = process.stdout.readline()
        if nextline == '' and process.poll() is not None:
            break
        sys.stdout.write(nextline)
        sys.stdout.flush()
        queue.put_nowait([idx,nextline])
    process.communicate()
    exitCode = process.returncode
    queue.put_nowait([idx,'exitcode',exitCode])


class BambooLauncher:
    def __init__(self,mode,module,config,args,output,cores):
        self.mode   = mode
        self.module = module
        self.config = config
        self.args   = args
        self.output = output

        cmds = self.produceCommands()

        if self.mode == 'debug':
            print ('Printing commands [{}]'.format(len(cmds)))
            for cmd in cmds:
                print(cmd)
        elif self.mode == 'driver' or self.mode  == 'finalize':
            if which('bambooRun') is None:
                raise RuntimeError('bamboo venv not setup, will stop here')
            if cores == -1:
                cores = len(cmds)
            if cores > mp.cpu_count():
                print ("You asked for {} processes but the machine has only {} cpus, will reduce to latter".format(cores,mp.cpu_count()))
                cores = mp.cpu_count()

            self.runningLoop(cmds,cores)
        else:
            raise RuntimeError(f'Mode {self.mode} not understood')

    def runningLoop(self,cmds,cores):
        queue = mp.Queue()
        processes = [mp.Process(target=runBamboo, args=(idx,cmd,queue)) for idx,cmd in enumerate(cmds)]
        N = len(processes)
        launched  = [False] * N
        completed = [False] * N
        running   = [False] * N
        exitcodes = [None] * N
        slurm_ids = [None] * N
        sbatch_cmds = []
        while not all(completed):
            # Check queue #
            while not queue.empty():
                res = queue.get_nowait()
                if len(res) == 3:
                    idx,status,exitcode = res
                    assert status == 'exitcode'
                    exitcodes[idx] = exitcode
                elif len(res) == 2:
                    idx,line = res
                    if self.mode == 'driver':
                        # When a jobs has submitted, start new job
                        if 'job ID is' in line:
                            slurm_ids[idx] = int(re.findall(r'\d+',line)[0])
                            completed[idx] = True
                            cores += 1  
                    if self.mode == 'finalize':
                        # When a finalize spits out a sbatch command line, catch it
                        if 'sbatch' in line:
                            quoted = re.compile("'[^']*'")
                            sbatch_cmd = quoted.findall(line)[0].replace("'","")
                            sbatch_cmds.append(sbatch_cmd)
                            completed[idx] = True
                else:
                    print ('Something fishy ... '+res)

            # Check progress #
            for i in range(N):
                if not launched[i] and not running[i] and sum(running) < cores:
                    processes[i].start()
                    running[i] = True
                    launched[i] = True
                if running[i] and processes[i].exitcode is not None:
                    running[i] = False
                    completed[i] = True
            print ('Processes : {pending:2d} pending - {running:2d} running - {completed:2d} completed [{successes:2d} successes / {failures:2d} failures]'.format(
                    pending   = N-sum(launched),
                    running   = sum(running),
                    completed = sum(completed),
                    successes = sum([ec == 0 for ec in exitcodes if ec is not None]),
                    failures  = sum([ec != 0 for ec in exitcodes if ec is not None])))


            time.sleep(3)


        if self.mode == 'driver':
            if any([slurm_id is None for slurm_id in slurm_ids]):
                print ('Some processes did not return a jobid')
            if any([slurm_id is not None for slurm_id in slurm_ids]):
                print ('Submitted job ID : ',sorted(slurm_ids))
        if self.mode == 'finalize' and len(sbatch_cmds)>0:
            print ('Not all jobs have succeeded, see below for list of commands to resubmit')
            for sbatch_cmd in sbatch_cmds:
                print (sbatch_cmd.replace('sbatch','sbatch --licenses=cms_storage:3 '))
    
        for p in processes:
            p.join()

    def produceCommands(self):
        baseCmd = 'bambooRun {mode} -m {module} {config} {args} -o {output}'

        # Mode #
        if self.mode == 'driver':
            mode = '--distributed=driver'
        elif self.mode == 'finalize':
            mode = '--distributed=finalize'
        else:
            mode = ''

        # Yaml #
        if not isinstance(self.config,dict):
            configs = [{'key':None,'form':None,'arg':self.config}]
        else:
            if len(self.config.keys()) != 1:
                raise RuntimeError("Currently only one key can be used for the config entry")
            key = list(self.config.keys())[0]
            configs = []
            for entry in self.config[key]:
                if len(entry.keys()) != 1:
                    raise RuntimeError("Currently only one key can be used for each entry in the config")
                configs.append({'key':key,'form':list(entry.keys())[0],'arg':entry[list(entry.keys())[0]]})

        # Additionnal arguments #
        fixedArgs = []
        toVaryArgs = {}
        if not isinstance(self.args,list):
            raise RuntimeError('args items needs to be a list (currently {})'.format(self.args))
        for arg in self.args:
            if isinstance(arg,str):
                fixedArgs.append(arg)
            elif isinstance(arg,dict):
                if len(arg.keys()) != 1:
                    raise RuntimeError('Currently only one key is possible per argument')
                key = list(arg.keys())[0]
                toVaryArgs[key] = arg[key]
            else:
                raise RuntimeError("Could not understand arg {}".format(arg))

        toVaryIdx = {key:[i for i in range(len(val))] for key,val in toVaryArgs.items()}
        variedIdx = []
        for combIdx in list(itertools.product(*toVaryIdx.values())):
            variedIdx.append({key:val for key,val in zip(toVaryIdx.keys(),combIdx)})

        variedArgs = []
        for indices in variedIdx:
            entry = {}
            for key,idx in indices.items():
                arg = toVaryArgs[key][idx]
                if len(arg.keys()) != 1:
                    raise RuntimeError("Currently only one key can be used per entry in {}".format(key))
                form = list(arg.keys())[0]
                entry[key] = {'form':form,'arg':arg[form]}
            variedArgs.append(entry)

            
        # Full command #
        cmds = []
        for config in configs:
            # Yaml #
            for variedArg in variedArgs:
                #output = copy.copy(self.output)
                output_format = {}
                # config output format #
                if config['key'] is not None:
                    output_format[config['key']] = config['form']
                # Additional args #
                cmdArgs = " ".join(fixedArgs)+" "
                cmdArgs += " ".join([val['arg'] for val in variedArg.values()])+" "
                for key,val in variedArg.items():
                    output_format[key] = val['form']

                output_path = self.output.format(**output_format)
                if self.mode == 'driver' and os.path.exists(output_path):
                    print ('Driver mode and path {} already exists, will pass command'.format(output_path))
                if self.mode == 'finalize' and not os.path.exists(output_path):
                    print ('Finalize mode and path {} does not exist, will pass command'.format(output_path))

                # Line command #
                cmd = baseCmd.format(mode       = mode,
                                     module     = self.module,
                                     config     = config['arg'],
                                     args       = cmdArgs,
                                     output     = output_path)
                cmds.append(cmd)
            
        return cmds
            
        



YMLIncludeLoader.add_constructor('!include', YMLIncludeLoader.include)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Bamboo launcher')
    parser.add_argument('--yaml', action='store', required=True, type=str,
                        help='Yaml containing parameters')
    parser.add_argument('--mode', action='store', required=True, type=str,
                        help='Mode for launcher : driver | finalize | debug')
    parser.add_argument('-j', action='store', required=False, type=int, default=1,
                        help='Number of commands to run in parallel (default = 1), using -1 will spawn all the commands')
    args = parser.parse_args()

    if args.yaml is None:
        raise RuntimeError("Must provide the YAML file")
    if not os.path.isfile(args.yaml):
        raise RuntimeError("YAML file {} is not a valid file".format(args.yaml))
        
    with open(args.yaml,'r') as handle:
        f = yaml.load(handle,Loader=YMLIncludeLoader)
    instance = BambooLauncher(**f,mode=args.mode,cores=args.j)