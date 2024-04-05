import os
import awkward as ak
import uproot
from collections import Mapping


class Data(Mapping):
    """ Cache for loading data from tree """
    def __init__(self,files,treename):
        self.trees = []
        self.files = files if isinstance(files,(list,tuple)) else [files]
        for f in self.files:
            if isinstance(f,uproot.reading.ReadOnlyDirectory):
                self.trees.append(f[treename])
            elif isinstance(f,str):
                assert os.path.isfile(f), f"File {f} does not exist"
                self.trees.append(uproot.open(f)[treename])
            else:
                raise ValueError

        self.data = {}

    def __getitem__(self,key):
        for f,t in zip(self.files,self.trees):
            if key not in t.keys():
                raise KeyError(f'Key {key} is not present in tree {t} of file {f}')
        if key not in self.data.keys():
            self.data[key] = ak.concatenate([tree[key].array() for tree in self.trees],axis=0)
        return self.data[key]

    def __len__(self):
        return len(self.data)

    def __iter__(self):
        return self.data.__iter__()

    @property
    def branches(self):
        return list(set([k for t in self.trees for k in t.keys()]))


    def keys(self):
        return self.data.keys()

