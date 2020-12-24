import sys
import glob
import os
import math
import logging
import json
import pickle
import copy
import collections
import random
import yaml
import tensorflow as tf
#import enlighten

from prettytable import PrettyTable
import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from itertools import chain
from root_numpy import root2array, rec2array

import ROOT

import parameters
from import_tree import LoopOverTrees
from generate_mask import GenerateSampleMasks, GenerateSliceIndices


class DataGenerator(tf.keras.utils.Sequence):
    def __init__(self,path,inputs,outputs,inputsLBN=None,weight=None,cut='',batch_size=32,state_set='',model_idx=None):
        self.path       = path                          # Path to root file : can be single file, list or dir (in which case will take all files inside)
        self.inputs     = inputs                        # List of strings of the variables as inputs
        self.inputsLBN  = inputsLBN                     # List of inputs of the LBN
        self.outputs    = outputs                       # List of strings of the variables as outputs
        self.cut        = cut                           # Branch expression to use the cut
        self.weight     = weight                        # string for the branch containing the weight
        self.batch_size = batch_size                    # Batch size
        self.model_idx  = model_idx                     # model Idx for cross validation
        self.variables  = self.inputs + self.outputs    # List of all variables to be taken from root trees
        if self.inputsLBN is not None:
            self.variables += [inp for inp in self.inputsLBN if inp not in self.inputsLBN]

        if isinstance(self.path,str):
            if self.path.endswith('.yml'):
                with open (parameters.config,'r') as f:
                    sampleConfig = yaml.load(f)

                self.input_dir = sampleConfig["sampleDir"]
                self.sample_dict = sampleConfig["sampleDict"]

                strSelect = [f'{cat}_{channel}_{node}' for channel in parameters.channels for cat in parameters.categories for node in parameters.nodes]
                self.list_files = [os.path.join(self.input_dir,sample) for key in strSelect  for era in parameters.eras for sample in self.sample_dict[era][key]]

            elif os.path.isdir(path):
                self.list_files = glob.glob(path+'/*.root') 
            elif os.path.isfile(path):
                self.list_files = [path]                    
            else:
                raise RuntimeError("path '%s' is not a dir nor a file "%path)
        elif isinstance(self.path,list):
            self.list_files = path
        else:
            raise RuntimeError('path argument not understood') 
            
        self.state_set = state_set 
        if self.state_set not in ['training', 'validation', 'evaluation', 'output']:
            raise RuntimeError("Generator state %s not understood"%self.state_set)
        if parameters.crossvalidation:
            set_indices = GenerateSliceIndices(self.model_idx)
            if self.state_set == 'output': self.set_idx = set_indices[0]
            if self.state_set == 'evaluation': self.set_idx = set_indices[1]
            if self.state_set == 'training' or self.state_set == 'validation': self.set_idx = set_indices[2]
            self.cut += ' && (' + ' | '.join(['{}%{}=={}'.format(parameters.splitbranch,parameters.N_slices,idx) for idx in self.set_idx])+')'
        else:
            self.maskSamples = GenerateSampleMasks(self.list_files, parameters.suffix)

        self.xsec_dict = dict()
        self.event_weight_sum_dict = dict()
        for era in parameters.eras:
            with open(parameters.xsec_json.format(era=era),'r') as handle:
                self.xsec_dict[era] = json.load(handle)
            with open(parameters.event_weight_sum_json.format(era=era),'r') as handle:
                self.event_weight_sum_dict[era] = json.load(handle)

        if (len(self.list_files)>self.batch_size):
            logging.warning("More files than requested batch size, might be errors")
        logging.info("Starting importation for %s set"%self.state_set)

        self.get_indices()
        self.get_fractions()
        self.n          = 0
        self.max        = self.__len__() # Must be after get_indices because that's where self.n_batches is defined

    def get_indices(self):
        self.batch_sample = dict() 
        self.indices = dict()
        self.masks = dict()
        self.n_tot = 0
        #pbar = enlighten.Counter(total=len(self.list_files), desc='Indices', unit='File')
        logging.info('Starting indices importation')
        for i,f in enumerate(self.list_files):
            rootFile = ROOT.TFile(f)
            tree_exists = rootFile.GetListOfKeys().Contains(parameters.tree_name)
            if not tree_exists:
                continue 
            n = rootFile.Get(parameters.tree_name).GetEntries()
            rootFile.Close()
            if self.cut != '':
                indices_slice = np.zeros((n,1))
                incn = 100000
                p = 0
                while p < n:
                    indices_slice[p:p+incn] = rec2array(root2array(f,parameters.tree_name,branches=[self.cut],start=p,stop=p+incn))
                    p += incn
                self.indices[f],_ = np.where(indices_slice==1)
            else:
                self.indices[f] = np.arange(n)
            if not parameters.crossvalidation:
                mask = self.maskSamples[f]
                if self.indices[f].shape[0] != mask.shape[0]:
                    raise RuntimeError('Sample %s : size %d different from mask size %d'%(f,len(self.indices[f]),mask.shape[0]))
                if self.state_set == 'training' or self.state_set == 'validation':
                    self.masks[f] = mask == 0
                if self.state_set == 'evaluation':
                    self.masks[f] = mask == 1
                if self.state_set == 'output':
                    self.masks[f] = mask == 2
            else:
                self.masks[f] = np.full((self.indices[f].shape[0],), True, dtype=bool)
            if self.state_set == 'training':
                self.masks[f] = np.logical_and(self.masks[f],self.indices[f]%10!=0)
            if self.state_set == 'validation':
                self.masks[f] = np.logical_and(self.masks[f],self.indices[f]%10==0)
            self.indices[f] = list(self.indices[f])
            self.masks[f] = list(self.masks[f])
            assert len(self.masks[f]) == len(self.indices[f])
            self.n_tot += sum(self.masks[f])
            #pbar.update()
        logging.info("Number of events in %s set: %d"%(self.state_set,self.n_tot))

        if self.n_tot<self.batch_size:
            raise RuntimeError("Fewer events than required batch size for generator")

        self.n_batches = self.n_tot//self.batch_size
        logging.info("Will use %d batches of %d events (%d events will be lost for truncation purposes)"%(self.n_batches,self.batch_size,self.n_tot%self.batch_size))

        total_in_batch = 0
        for f in self.indices.keys():
            #size_in_batch = math.ceil((len(self.indices[f])/self.n_tot)*self.batch_size)
            if self.state_set == 'training' or self.state_set == 'validation': 
                size_in_batch = max(math.ceil((sum(self.masks[f])/self.n_tot)*self.batch_size),1)
            else:
                size_in_batch = self.batch_size
            self.batch_sample[f] = size_in_batch 
            total_in_batch += size_in_batch

    def get_fractions(self):
        self.indices_per_batch = []
        self.masks_per_batch = []
        pointers = {k:0 for k in self.indices.keys()}
        keys = list(self.indices.keys())
        tag_count = [{}]*self.n_batches
        era_count = [{}]*self.n_batches
        for i in range(self.n_batches):
            filled = False
            indices_sample = collections.defaultdict(list)
            masks_sample = collections.defaultdict(list)
            counter = 0
            while not filled:
                if self.state_set == 'training' or self.state_set == 'validation':
                    random.shuffle(keys) 
                for key in keys:
                    rem = self.batch_size-counter
                    inc = min(rem,self.batch_sample[key])
                    if pointers[key] < len(self.indices[key]):
                        indices_slice = self.indices[key][pointers[key]:min(pointers[key]+inc,len(self.indices[key]))]
                        masks_slice = self.masks[key][pointers[key]:min(pointers[key]+inc,len(self.masks[key]))]
                        indices_sample[key].extend(indices_slice)
                        masks_sample[key].extend(masks_slice)
                        pointers[key] += inc
                        counter += len([i for i,m in zip(indices_slice,masks_slice) if m])
                    if counter == self.batch_size:
                        filled = True
                        break
            tag_count[i] = {node:0 for node in parameters.nodes}
            era_count[i] = {era:0 for era in parameters.eras}
            for samplepath in indices_sample.keys():
                sample = samplepath.replace(self.input_dir+'/','') 
                keySelect = None
                eraSelect = None
                for era,keyDict in self.sample_dict.items():
                    for key,list_samples in keyDict.items():
                        if sample in list_samples:
                            keySelect = key
                            eraSelect = era
                            break
                    if keySelect is not None:
                        break
                if keySelect is None:
                    raise RuntimeError('Could not find sample %s in yaml file'%keySelect)
                tag = max([node for node in parameters.nodes if node in keySelect], key=len)
                cont = len([i for i,m in zip(indices_sample[samplepath],masks_sample[samplepath]) if m])
                tag_count[i][tag] += cont
                era_count[i][eraSelect] += cont

            for s, ind in indices_sample.items(): 
                indices_sample[s] = (min(ind),max(ind))

            self.indices_per_batch.append(dict(indices_sample))
            self.masks_per_batch.append(dict(masks_sample))
        tag_count_all = {k:sum([tag[k] for tag in tag_count]) for k in parameters.nodes}
        era_count_all = {k:sum([era[k] for era in era_count]) for k in parameters.eras}
    
        pt_tag = PrettyTable(["Batch"]+parameters.nodes) 
        for i in range(len(tag_count)):
            pt_tag.add_row(['Batch %d'%i]+["%d [%3.2f%%]"%(tag_count[i][node],tag_count[i][node]*100/self.batch_size) for node in parameters.nodes])  
        pt_tag.add_row(['Total']+["%d [%3.2f%%]"%(tag_count_all[node],tag_count_all[node]*100/self.n_tot) for node in parameters.nodes])  
        logging.info("Node content per batch")
        for line in pt_tag.get_string().split('\n'):
            logging.info(line)

        pt_era = PrettyTable(["Batch"]+parameters.eras) 
        for i in range(len(era_count)):
            pt_era.add_row(['Batch %d'%i]+["%d [%3.2f%%]"%(era_count[i][era],era_count[i][era]*100/self.batch_size) for era in parameters.eras])  
        pt_era.add_row(['Total']+["%d [%3.2f%%]"%(era_count_all[era],era_count_all[era]*100/self.n_tot)for era in parameters.eras])  
        logging.info("Era content per batch")
        for line in pt_era.get_string().split('\n'):
            logging.info(line)

    def __getitem__(self,index,additional_columns=False): # gets the batch for the supplied index
        # return a tuple (numpy array of image, numpy array of labels) or None at epoch end
        logging.debug("-"*80)
        logging.debug("New batch importation = index %d"%index)

        samples = []
        masks = []
        start = []
        stop = []
        eras = []
        tags = []

        for samplepath,ind in self.indices_per_batch[index].items():
            sample = samplepath.replace(self.input_dir,'') 
            keySelect = None
            eraSelect = None
            for era,keyDict in self.sample_dict.items():
                for key,list_samples in keyDict.items():
                    if sample in list_samples:
                        keySelect = key
                        eraSelect = era
                        break
                if keySelect is not None:
                    break
            if eraSelect is None:
                raise RuntimeError('Could not find the era in sampleDict for sample %s'%sample)
            if keySelect is None:
                raise RuntimeError('Could not find the key in sampleDict for sample %s'%sample)
            eras.append(eraSelect)
            tmp_nodes = [node for node in parameters.nodes if node in keySelect]
            tags.append(max([node for node in parameters.nodes if node in keySelect], key=len)) # Find longest node match (eg if one node name is included in another)
            samples.append(sample)
            masks.append(self.masks_per_batch[index][samplepath])
            start.append(ind[0])
            stop.append(ind[1]+1) # Python stop and start 

        # Import data #
        data = LoopOverTrees(input_dir             = self.input_dir,
                             variables             = self.variables,
                             weight                = self.weight,
                             list_sample           = samples,
                             start                 = start,
                             stop                  = stop,
                             xsec_dict             = self.xsec_dict,
                             event_weight_sum_dict = self.event_weight_sum_dict,
                             lumi_dict             = parameters.lumidict,
                             eras                  = eras,
                             tree_name             = parameters.tree_name,
                             additional_columns    = {'tag':tags,'era':eras},
                             cut                   = self.cut)
        mask = [m for m in chain.from_iterable(masks)]
        data = data[mask]
        data = data.sample(frac=1).reset_index(drop=True) # Randomize
        assert data.shape[0] == self.batch_size

        # weight Equalization #
        weight_per_tag = {tag:data[data['tag']==tag]['event_weight'].sum() for tag in pd.unique(data['tag'])}
        weight_scale = data['tag'].apply(lambda row: weight_per_tag[row])
        data['learning_weight'] = data['event_weight']*1e3/weight_scale

        inputs  = [var.replace('$','') for var in self.inputs]
        outputs = [var.replace('$','') for var in self.outputs]

        # Add target #
        label_encoder = LabelEncoder()
        onehot_encoder = OneHotEncoder(sparse=False)
        label_encoder.fit(outputs)
        # From strings to labels #
        integers = label_encoder.transform(data['tag']).reshape(-1, 1)
        # From labels to strings #
        onehotobj = onehot_encoder.fit(np.arange(len(outputs)).reshape(-1, 1))
        onehot = onehotobj.transform(integers)
        # From arrays to pd DF #
        cat = pd.DataFrame(onehot,columns=label_encoder.classes_,index=data.index)
        # Add to full #
        data = pd.concat([data,cat],axis=1)

        if not additional_columns:
            if self.inputsLBN is not None:
                data_input = (data[inputs].values,data[self.inputsLBN].values.reshape(-1,4,len(self.inputsLBN)//4))
            else:
                data_input = data[inputs].values
            data_output = data[outputs].values
            data_weight = data["learning_weight"].values
            return data_input,data_output,data_weight
        else:
            addcol = [col for col in data.columns if col not in inputs+outputs+["learning_weight"]]
            if self.inputsLBN is not None:
                return data[inputs+[inp for inp in self.inputsLBN if inp not in inputs]],data[outputs],data["learning_weight"],data[addcol]
            else:
                return data[inputs],data[outputs],data["learning_weight"],data[addcol]
             

    def __len__(self): # gets the number of batches
        # return the number of batches in this epoch (do not change in the middle of an epoch)
        return self.n_batches
    def on_epoch_end(self): # performs auto shuffle if enabled
        # Do what we need to do between epochs
        pass
        #self.get_fractions() # Change the batches after each epoch 
            # TODO : not sure we need it

    def __next__(self):
        if self.n >= self.max:
           self.n = 0
        result = self.__getitem__(self.n)
        self.n += 1
        return result
