#! /bin/env python

import copy
import os
import datetime
import sys
import glob
import logging

# Slurm configuration
from CP3SlurmUtils.Configuration import Configuration
from CP3SlurmUtils.SubmitWorker import SubmitWorker
from CP3SlurmUtils.Exceptions import CP3SlurmUtilsException

# Personal files #
import parameters

def submit_on_slurm(name,args,debug=False):
    # Check arguments #
    GPU = args.find("--GPU") != -1
    output = args.find("--output") != -1

    config = Configuration()
    config.sbatch_partition = parameters.partition
    config.sbatch_qos = parameters.QOS
    config.sbatch_chdir = parameters.main_path
    config.sbatch_time = parameters.time
    config.sbatch_memPerCPU= parameters.mem
    #config.sbatch_additionalOptions = ['-n '+str(parameters.tasks)]
    config.sbatch_additionalOptions = ['-c '+str(parameters.cpus)]#'--exclude=mb-sab[001-005,007-021,081-084,087-088,090,101-103],mb-opt[015-018,021,024-025,031,042,051-052,054,056-064,067-079,111,114-116],mb-ivy[201-208,211-212,214-217,219,220-222,224-227],mb-wes[001-002,003,005-019,021-051,053-055,057-074,076-086,251-252]']
    if parameters.partition == 'cp3-gpu':
        config.sbatch_additionalOptions += ['--gres gpu:1','--export=NONE']
    if parameters.partition == 'gpu':
        config.sbatch_additionalOptions += ['--gres=gpu:TeslaV100:'+str(parameters.gpus),'--export=NONE']
        
    config.inputSandboxContent = []
    config.useJobArray = True
    config.inputParamsNames = []
    config.inputParams = []
    if output:
        config.inputParamsNames += ["--verbose"]
        config.inputParams += [[""]]
    if not output:
        config.inputParamsNames += ['scan','task']
        if parameters.crossvalidation:
            config.inputParamsNames += ['modelId']

    config.payload = ""

    if parameters.partition == 'cp3-gpu':
        config.payload += "export PYTHONPATH=/python3/lib/python3.6/site-packages/:$PYTHONPATH\n" # GPU tf
        config.payload += "export PYTHONPATH=/root6/lib:$PYTHONPATH\n" # ROOT
        config.payload += "module load cp3\n" # needed on gpu to load slurm_utils
        config.payload += "module load python/python36_sl7_gcc73\n" 
        config.payload += "module load slurm/slurm_utils\n"
    if parameters.partition == 'gpu':
        config.payload += "module load root/6.12.04-sl7_gcc73 \n"
        config.payload += "module load root_numpy \n"
        config.payload += "module load TensorFlow \n"
        
    config.payload += "python3 {script} "
    if not output:
        config.payload += "--scan ${{scan}} --task ${{task}} "
    if parameters.crossvalidation:
        config.payload += "--modelId ${{modelId}}"
    config.payload += args

    timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    out_dir = parameters.main_path

    slurm_config = copy.deepcopy(config)
    slurm_working_dir = os.path.join(out_dir,'slurm',name+'_'+timestamp)

    slurm_config.batchScriptsDir = os.path.join(slurm_working_dir, 'scripts')
    slurm_config.inputSandboxDir = slurm_config.batchScriptsDir
    slurm_config.stageoutDir = os.path.join(slurm_working_dir, 'output')
    slurm_config.stageoutLogsDir = os.path.join(slurm_working_dir, 'logs')
    slurm_config.stageoutFiles = ["*.csv","*.zip","*.png"]

    slurm_config.payload = config.payload.format(script=os.path.join(out_dir,"HHMachineLearning.py"))

    if not output:
        for f in glob.glob(os.path.join(parameters.main_path,'split',name,'*.pkl')):
            task = os.path.basename(f)
            if parameters.crossvalidation:
                for N in range(parameters.N_models):
                    slurm_config.inputParams.append([name,task,N])
            else:
                slurm_config.inputParams.append([name,task])

    # Submit job!

    logging.info("Submitting job...")
    if not debug:
        submitWorker = SubmitWorker(slurm_config, submit=True, yes=True, debug=False, quiet=False)
        submitWorker()
        logging.info("Done")
    else:
        logging.info("Number of jobs : %d"%len(slurm_config.inputParams))
        logging.info(slurm_config.payload)
        logging.info(slurm_config.inputParamsNames)
        for inputParam in slurm_config.inputParams:
            logging.info(inputParam)
        logging.info('... don\'t worry, jobs not sent')

