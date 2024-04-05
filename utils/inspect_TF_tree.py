import os
import sys
import math
import numpy as np
import uproot
import itertools
import awkward as ak
import ROOT

from inspect_utils import Data


#def get_data(f):
#     = uproot.open(f)
#    return {k.split(';')[0]:data[k] for k in data.keys() if data[k].classname=='TTree'}

f = sys.argv[-1]
if not os.path.isfile(f):
    raise RuntimeError(f'{f} is not a file')
F = uproot.open(f)

datas = {k.split(';')[0]:Data(F,k) for k in F.keys() if F[k].classname=='TTree'}
#trees = get_trees(f)
#arrays = {k:tree.arrays() for k,tree in trees.items()}


def print_gen_tree(flavour,row):
    data = datas[f'gen_{flavour}']
    print (f'Printing {data[f"n_{flavour}_gen"][row]} gen {flavour} : row = {row}, event = {data["event"][row]}')
    for i in range(data[f'n_{flavour}_gen'][row]):
        if f"{flavour}_{i}_idx" not in data.branches:
            print ('\tMissing following particles in tree!!!')
            break
        print (f'    idx = {data[f"{flavour}_{i}_idx"][row]:3d}, E = {data[f"{flavour}_{i}_E"][row]:8.2f}, pT = {data[f"{flavour}_{i}_pt"][row]:8.2f}, eta = {data[f"{flavour}_{i}_eta"][row]:+6.2f}, phi = {data[f"{flavour}_{i}_phi"][row]:+5.2f}, M = {data[f"{flavour}_{i}_M"][row]:6.2f}, pdg = {data[f"{flavour}_{i}_pdgId"][row]:3d}, clean = {data[f"{flavour}_{i}_clean"][row]:1d}',end='')
        if flavour in ['e','m','t']:
            print ()
        else:
            print (f', resolved = {data[f"{flavour}_{i}_resolved"][row]:1d}')

def print_match_tree(flavour,row):
    data = datas[f'match_{flavour}']
    print (f'Flavour {flavour} tree : row = {row}, event = {data["event"][row]}, content :')
    for br in [f'n_{flavour}_gen_match',f'n_{flavour}_reco_match',f'n_{flavour}_reco_final',f'n_{flavour}_reco_central',f'n_{flavour}_reco_forward']:
        if br in data.branches:
            print (f'  - {br:20s} : {data[br][row]}')

    for i in range(data[f'n_{flavour}_reco_match'][row]):
        if f"{flavour}_{i}_preFSR_idx" not in data.branches:
            print ('\tMissing following particles in tree!!!')
            break

        s = f'  Object {i} : '
        for t in ['clean','unique','unique_count','resolved','resolved_count','central','forward','DR','relE','check_FSR','sum_FSR']:
            br = f'{flavour}_{i}_match_{t}'
            if br in data.branches:
                s += f'{t} {data[br][row]:4.3f}  '
        print (s)
        print (f'    Pre FSR  : idx = {data[f"{flavour}_{i}_preFSR_idx"][row]:3d}, E = {data[f"{flavour}_{i}_preFSR_E"][row]:8.2f}, pT = {data[f"{flavour}_{i}_preFSR_pt"][row]:8.2f}, eta = {data[f"{flavour}_{i}_preFSR_eta"][row]:+6.2f}, phi = {data[f"{flavour}_{i}_preFSR_phi"][row]:+5.2f}, M = {data[f"{flavour}_{i}_preFSR_M"][row]:6.2f}, pdg = {data[f"{flavour}_{i}_preFSR_pdgId"][row]:3d}')
        print (f'    Post FSR : idx = {data[f"{flavour}_{i}_postFSR_idx"][row]:3d}, E = {data[f"{flavour}_{i}_postFSR_E"][row]:8.2f}, pT = {data[f"{flavour}_{i}_postFSR_pt"][row]:8.2f}, eta = {data[f"{flavour}_{i}_postFSR_eta"][row]:+6.2f}, phi = {data[f"{flavour}_{i}_postFSR_phi"][row]:+5.2f}, M = {data[f"{flavour}_{i}_postFSR_M"][row]:6.2f}, pdg = {data[f"{flavour}_{i}_postFSR_pdgId"][row]:3d}')

        if flavour in ['b','l','g']:
            print (f'    Gen jet  : idx = {data[f"{flavour}_{i}_genjet_idx"][row]:3d}, E = {data[f"{flavour}_{i}_genjet_E"][row]:8.2f}, pT = {data[f"{flavour}_{i}_genjet_pt"][row]:8.2f}, eta = {data[f"{flavour}_{i}_genjet_eta"][row]:+6.2f}, phi = {data[f"{flavour}_{i}_genjet_phi"][row]:+5.2f}, M = {data[f"{flavour}_{i}_genjet_M"][row]:6.2f}')
        print (f'    Reco     : idx = {data[f"{flavour}_{i}_reco_idx"][row]:3d}, E = {data[f"{flavour}_{i}_reco_E"][row]:8.2f}, pT = {data[f"{flavour}_{i}_reco_pt"][row]:8.2f}, eta = {data[f"{flavour}_{i}_reco_eta"][row]:+6.2f}, phi = {data[f"{flavour}_{i}_reco_phi"][row]:+5.2f}, M = {data[f"{flavour}_{i}_reco_M"][row]:6.2f}')




def print_fatjet_match_tree(row):
    data = datas[f'match_fatjet']
    print (f'Flavour fatjet tree : row = {row}, event = {data["event"][row]}, content :')
    for br in ['n_fatjets','n_fatjets_match','n_fatjets_match_final','n_fatjets_match_central']:
        if br in data.branches:
            print (f'  - {br:20s} : {data[br][row]}')

    for i in range(data[f'n_fatjets_match'][row]):
        if f"fat_{i}_clean" not in data.branches:
            print ('\tMissing following particles in tree!!!')
            break

        s = f'  Object {i} : '
        for t in ['clean','two_prongs','prong_count','check_FSR','sub1_sum_FSR','sub2_sum_FSR']:
            br = f'fat_{i}_{t}'
            if br in data.branches:
                s += f'{t} {data[br][row]:4.3f}  '
        print (s)
        print ('    Subjet 1')
        print (f'\tPre FSR  : idx = {data[f"fat_{i}_sub1_preFSR_idx"][row]:3d}, E = {data[f"fat_{i}_sub1_preFSR_E"][row]:8.2f}, pT = {data[f"fat_{i}_sub1_preFSR_pt"][row]:8.2f}, eta = {data[f"fat_{i}_sub1_preFSR_eta"][row]:+6.2f}, phi = {data[f"fat_{i}_sub1_preFSR_phi"][row]:+5.2f}, M = {data[f"fat_{i}_sub1_preFSR_M"][row]:6.2f}, pdg = {data[f"fat_{i}_sub1_preFSR_pdgId"][row]:3d}')
        print (f'\tPost FSR : idx = {data[f"fat_{i}_sub1_postFSR_idx"][row]:3d}, E = {data[f"fat_{i}_sub1_postFSR_E"][row]:8.2f}, pT = {data[f"fat_{i}_sub1_postFSR_pt"][row]:8.2f}, eta = {data[f"fat_{i}_sub1_postFSR_eta"][row]:+6.2f}, phi = {data[f"fat_{i}_sub1_postFSR_phi"][row]:+5.2f}, M = {data[f"fat_{i}_sub1_postFSR_M"][row]:6.2f}, pdg = {data[f"fat_{i}_sub1_postFSR_pdgId"][row]:3d}')
        print (f'\tSubjet   : idx = {data[f"fat_{i}_sub1_subjet_idx"][row]:3d}, E = {data[f"fat_{i}_sub1_subjet_E"][row]:8.2f}, pT = {data[f"fat_{i}_sub1_subjet_pt"][row]:8.2f}, eta = {data[f"fat_{i}_sub1_subjet_eta"][row]:+6.2f}, phi = {data[f"fat_{i}_sub1_subjet_phi"][row]:+5.2f}, M = {data[f"fat_{i}_sub1_subjet_M"][row]:6.2f}')

        print ('    Subjet 2')
        print (f'\tPre FSR  : idx = {data[f"fat_{i}_sub2_preFSR_idx"][row]:3d}, E = {data[f"fat_{i}_sub2_preFSR_E"][row]:8.2f}, pT = {data[f"fat_{i}_sub2_preFSR_pt"][row]:8.2f}, eta = {data[f"fat_{i}_sub2_preFSR_eta"][row]:+6.2f}, phi = {data[f"fat_{i}_sub2_preFSR_phi"][row]:+5.2f}, M = {data[f"fat_{i}_sub2_preFSR_M"][row]:6.2f}, pdg = {data[f"fat_{i}_sub2_preFSR_pdgId"][row]:3d}')
        print (f'\tPost FSR : idx = {data[f"fat_{i}_sub2_postFSR_idx"][row]:3d}, E = {data[f"fat_{i}_sub2_postFSR_E"][row]:8.2f}, pT = {data[f"fat_{i}_sub2_postFSR_pt"][row]:8.2f}, eta = {data[f"fat_{i}_sub2_postFSR_eta"][row]:+6.2f}, phi = {data[f"fat_{i}_sub2_postFSR_phi"][row]:+5.2f}, M = {data[f"fat_{i}_sub2_postFSR_M"][row]:6.2f}, pdg = {data[f"fat_{i}_sub2_postFSR_pdgId"][row]:3d}')
        print (f'\tSubjet   : idx = {data[f"fat_{i}_sub2_subjet_idx"][row]:3d}, E = {data[f"fat_{i}_sub2_subjet_E"][row]:8.2f}, pT = {data[f"fat_{i}_sub2_subjet_pt"][row]:8.2f}, eta = {data[f"fat_{i}_sub2_subjet_eta"][row]:+6.2f}, phi = {data[f"fat_{i}_sub2_subjet_phi"][row]:+5.2f}, M = {data[f"fat_{i}_sub2_subjet_M"][row]:6.2f}')

        print ('    Fatjet')
        print (f'\t         : idx = {data[f"fat_{i}_idx"][row]:3d}, E = {data[f"fat_{i}_E"][row]:8.2f}, pT = {data[f"fat_{i}_pt"][row]:8.2f}, eta = {data[f"fat_{i}_eta"][row]:+6.2f}, phi = {data[f"fat_{i}_phi"][row]:+5.2f}, M = {data[f"fat_{i}_M"][row]:6.2f}')


def print_content(row):
    print_gen_tree('e',row)
    print_gen_tree('m',row)
    print_gen_tree('t',row)
    print_gen_tree('d',row)
    print_gen_tree('u',row)
    print_gen_tree('s',row)
    print_gen_tree('c',row)
    print_gen_tree('b',row)
    print_gen_tree('g',row)

    print_match_tree('e',row)
    print_match_tree('m',row)
    print_match_tree('t',row)
    print_match_tree('b',row)
    print_match_tree('l',row)
    print_match_tree('g',row)

    print_fatjet_match_tree(row)


def make_LV(pt,eta,phi,m):
    lv = ROOT.Math.PtEtaPhiMVector()
    lv.SetPt(pt)
    lv.SetEta(eta)
    lv.SetPhi(phi)
    lv.SetM(m)
    return lv

