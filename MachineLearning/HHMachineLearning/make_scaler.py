import os
import logging
import pickle
import glob
import enlighten
import numpy as np
import pandas as pd
from sklearn import preprocessing
from ROOT import TFile, TTree
from import_tree import Tree2Pandas
from root_numpy import tree2array, rec2array

import parameters

def MakeScaler(data=None,list_inputs=[],generator=False,batch=5000,list_samples=None,additional_columns={}):
    logging.info('Starting computation for the scaler')
    # Generate scaler #
    scaler = preprocessing.StandardScaler()
    if not os.path.exists(parameters.scaler_path):
        # Not generator #
        if data is not None:
            scaler.fit(data[list_inputs])
        # For generator #
        if generator:
            if list_samples is None:
                raise RuntimeError("Generator mask asked, you need to provide a sample list")
            logging.info("Computing mean")
            # Mean Loop #
            mean = np.zeros(len(list_inputs))
            Ntot = 0
            pbar = enlighten.Counter(total=len(list_samples), desc='Mean', unit='File')
            for f in list_samples:
                pbar.update()
                if not os.path.exists(f):
                    continue
                file_handle = TFile.Open(f)
                if not file_handle.GetListOfKeys().Contains(parameters.tree_name):
                    continue
                tree = file_handle.Get(parameters.tree_name)
                N = tree.GetEntries()
                Ntot += N
                file_handle.Close()
                logging.debug("Opening file %s (%d entries)"%(f,N))
                # Loop over batches #
                for i in range(0, N, batch):
                    array = Tree2Pandas(f,list_inputs,start=i,stop=i+batch,additional_columns=additional_columns,tree_name=parameters.tree_name)[[inp.replace('$','') for inp in list_inputs]].astype(np.float32).values
                    mean += np.sum(array,axis=0)
            mean /= Ntot
            
            # Var Loop #
            logging.info("Computing std")
            std = np.zeros(len(list_inputs))
            pbar = enlighten.Counter(total=len(list_samples), desc='Std', unit='File')
            for f in list_samples:
                pbar.update()
                if not os.path.exists(f):
                    continue
                file_handle = TFile.Open(f)
                if not file_handle.GetListOfKeys().Contains(parameters.tree_name):
                    continue
                tree = file_handle.Get(parameters.tree_name)
                N = tree.GetEntries()
                file_handle.Close()
                logging.debug("Opening file %s (%d entries)"%(f,N))
                # Loop over batches #
                for i in range(0, N, batch):
                    array = Tree2Pandas(f,list_inputs,start=i,stop=i+batch,additional_columns=additional_columns,tree_name=parameters.tree_name)[[inp.replace('$','') for inp in list_inputs]].astype(np.float32).values
                    std += np.sum(np.square(array-mean),axis=0)
            std = np.sqrt(std/Ntot)
            # Set manually #
            scaler.mean_ = mean
            scaler.scale_ = std

        # Disable preprocess on onehot variables #
        scaler.mean_[parameters.mask_op] = 0.
        scaler.scale_[parameters.mask_op] = 1.

        # Safe checks #
        scaler.mean_[np.isnan(scaler.mean_)] = 0.
        scaler.scale_[np.isnan(scaler.scale_)] = 1.
        scaler.scale_[scaler.scale_ == 0.] = 1.
        scaler.var_ =  scaler.scale_**2

        # Save #
        with open(parameters.scaler_path, 'wb') as handle:
            pickle.dump(scaler, handle)
        logging.info('Scaler %s has been created'%parameters.scaler_name)
    # If exists, will import it #
    else:
        with open(parameters.scaler_path, 'rb') as handle:
            scaler = pickle.load(handle)
        logging.info('Scaler %s has been imported'%parameters.scaler_name)
    # Test the scaler #
    if data is not None:
        try:
            y = scaler.transform(data[list_inputs])
            # Compute mean and var for inputs not in onehot encoding #
            mean_scale = np.mean(y[:,[not m for m in parameters.mask_op]])
            var_scale  = np.var(y[:,[not m for m in parameters.mask_op]])
        except ValueError:
            logging.warning("Problem with the scaler '%s' you imported, has the data changed since it was generated ?"%parameters.scaler_name)
        if abs(mean_scale)>0.01 or abs((var_scale-1)/var_scale)>0.1: # Check that scaling is correct to 1%
            logging.warning("Something is wrong with scaler '%s' (mean = %0.6f, var = %0.6f), maybe you loaded an incorrect scaler"%(parameters.scaler_name,mean_scale,var_scale))

