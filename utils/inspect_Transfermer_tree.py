import os
import sys
import math
import numpy as np
import uproot
import itertools
import awkward as ak
import ROOT

from inspect_utils import Data


f = sys.argv[-1]
if not os.path.isfile(f):
    raise RuntimeError(f'{f} is not a file')
F = uproot.open(f)

datas = {k.split(';')[0]:Data(F,k) for k in F.keys() if F[k].classname=='TTree'}


def print_reco_DL_tree(row):
    data = datas['reco_DL']
    print (f'Printing reco DL content : row = {row}, event = {data["event"][row]}')
    s = 'Flags :'
    for flag in ['flag_SR','flag_ElEl','flag_MuMu','flag_ElMu','flag_resolved','flag_boosted']:
        s += f' {flag.replace("flag_","")} = {data[flag][row]}'
    print (s)
    s = 'N     :'
    for n_obj in ['n_e','n_m','n_AK4','n_AK4B','n_AK8','n_AK8B']:
        s += f' {n_obj.replace("n_","")} = {data[n_obj][row]}'
    print (s)
    print ('Electrons')
    for i in range(1,data['n_e'][row]+1):
        print (f'    idx = {data[f"e{i}_idx"][row]:3.0f}, E = {data[f"e{i}_E"][row]:8.2f}, pT = {data[f"e{i}_pt"][row]:8.2f}, eta = {data[f"e{i}_eta"][row]:+6.2f}, phi = {data[f"e{i}_phi"][row]:+5.2f}, M = {data[f"e{i}_mass"][row]:6.2f}, pdg = {data[f"e{i}_pdgId"][row]:+3.0f}')
    print ('Muons')
    for i in range(1,data['n_m'][row]+1):
        print (f'    idx = {data[f"m{i}_idx"][row]:3.0f}, E = {data[f"m{i}_E"][row]:8.2f}, pT = {data[f"m{i}_pt"][row]:8.2f}, eta = {data[f"m{i}_eta"][row]:+6.2f}, phi = {data[f"m{i}_phi"][row]:+5.2f}, M = {data[f"m{i}_mass"][row]:6.2f}, pdg = {data[f"m{i}_pdgId"][row]:+3.0f}')
    print (f'AK4 jets : N = {data["n_AK4"][row]}')
    for i in range(1,data["n_AK4"][row]+1):
        if data[f"j{i}_idx"][row] >= 0:
            print (f'    idx = {data[f"j{i}_idx"][row]:3.0f}, E = {data[f"j{i}_E"][row]:8.2f}, pT = {data[f"j{i}_pt"][row]:8.2f}, eta = {data[f"j{i}_eta"][row]:+6.2f}, phi = {data[f"j{i}_phi"][row]:+5.2f}, M = {data[f"j{i}_mass"][row]:6.2f}, btag = {data[f"j{i}_btag"][row]:4.3f} (btagged = {data[f"j{i}_btagged"][row]:1.0f})')
    print (f'VBF jets : N = {data["n_VBF"][row]}')
    for i in range(1,data["n_VBF"][row]+1):
        if data[f"VBF{i}_idx"][row] >= 0:
            print (f'    idx = {data[f"VBF{i}_idx"][row]:3.0f}, E = {data[f"VBF{i}_E"][row]:8.2f}, pT = {data[f"VBF{i}_pt"][row]:8.2f}, eta = {data[f"VBF{i}_eta"][row]:+6.2f}, phi = {data[f"VBF{i}_phi"][row]:+5.2f}, M = {data[f"VBF{i}_mass"][row]:6.2f}, selected = {data[f"VBF{i}_sel"][row]:1.0f}')

    print (f'AK8B jets : N = {data["n_AK8B"][row]}')
    if data[f"fatj_idx"][row] >= 0:
        print (f'    idx = {data[f"fatj_idx"][row]:3.0f}, E = {data[f"fatj_E"][row]:8.2f}, pT = {data[f"fatj_pt"][row]:8.2f}, eta = {data[f"fatj_eta"][row]:+6.2f}, phi = {data[f"fatj_phi"][row]:+5.2f}, M = {data[f"fatj_mass"][row]:6.2f}')
        print ('AK8 sub jets')
        for i in [1,2]:
            print (f'    idx = {data[f"fatj_sub{i}_idx"][row]:3.0f}, E = {data[f"fatj_sub{i}_E"][row]:8.2f}, pT = {data[f"fatj_sub{i}_pt"][row]:8.2f}, eta = {data[f"fatj_sub{i}_eta"][row]:+6.2f}, phi = {data[f"fatj_sub{i}_phi"][row]:+5.2f}, M = {data[f"fatj_sub{i}_mass"][row]:6.2f}, btag = {data[f"fatj_sub{i}_btag"][row]:4.3f} (btagged = {data[f"fatj_sub{i}_btagged"][row]:1.0f})')

def print_gen_sample(tree,row):
    if tree == 'gen_TT':
        particles = [
            'top',
            'antitop',
            'bottom',
            'antibottom',
            'W_plus_from_top',
            'W_minus_from_antitop',
            'lep_plus',
            'lep_minus',
            'neutrino',
            'antineutrino',
            'quark_up',
            'quark_down',
            'antiquark_up',
            'antiquark_down',
        ]
    elif tree == 'gen_ST':
        particles = [
            'top',
            'antitop',
            'bottom',
            'antibottom',
            'W_plus_from_top',
            'W_minus_from_antitop',
            'lep_plus_from_top',
            'neutrino_from_top',
            'quark_up_from_top',
            'antiquark_down_from_top',
            'lep_minus_from_top',
            'antineutrino_from_top',
            'quark_down_from_top',
            'antiquark_up_from_top',
            'W_minus_prompt',
            'lep_minus_from_prompt_W',
            'antineutrino_from_prompt_W',
            'quark_down_from_prompt_W',
            'antiquark_up_from_prompt_W',
            'W_plus_prompt',
            'lep_plus_from_prompt_W',
            'neutrino_from_prompt_W',
            'quark_up_from_prompt_W',
            'antiquark_down_from_prompt_W',
        ]
    elif tree == 'gen_DY':
        particles = [
            'Z',
            'lep_from_Z',
            'antilep_from_Z',
            'quark_from_Z',
            'antiquark_from_Z',
            'lep_from_nonres',
            'antilep_from_nonres',
            'quark_from_nonres',
            'antiquark_from_nonres',
        ]
    elif tree == 'gen_ZZ':
        particles = [
            'Z1',
            'Z2',
            'lep_from_Z1',
            'antilep_from_Z1',
            'quark_from_Z1',
            'antiquark_from_Z1',
            'lep_from_Z2',
            'antilep_from_Z2',
            'quark_from_Z2',
            'antiquark_from_Z2',
            'lep_from_nonres',
            'antilep_from_nonres',
            'quark_from_nonres',
            'antiquark_from_nonres',
        ]
    elif tree == 'gen_ZH':
        particles = [
            'H',
            'Z',
            'bottom',
            'antibottom',
            'lep_plus',
            'lep_minus',
        ]
    elif tree == 'gen_HH':
        particles = [
            'radion',
            'graviton',
            'H1',
            'H2',
            'VBF_quark1',
            'VBF_quark2',
            'VBF_quark3',
            'VBF_quark4',
            'bottom',
            'antibottom',
            'W_plus',
            'W_minus',
            'Z1',
            'Z2',
            'lep_plus_from_W',
            'lep_minus_from_W',
            'neutrino_from_W',
            'antineutrino_from_W',
            'quark_up_from_W',
            'quark_down_from_W',
            'antiquark_up_from_W',
            'antiquark_down_from_W',
            'lep_plus_from_Z',
            'lep_minus_from_Z',
            'neutrino_from_Z',
            'antineutrino_from_Z',
            'quark_up_from_Z',
            'quark_down_from_Z',
            'antiquark_up_from_Z',
            'antiquark_down_from_Z',
        ]
    else:
        raise RuntimeError(f'Tree {tree} not found')

    data = datas[tree]
    print (f'Printing {tree} content : row = {row}, event = {data["event"][row]}, valid = {data["flag_valid"][row]}')
    max_len = max([len(particle) for particle in particles ])
    for particle in particles:
        if f"{particle}_idx" in data.branches and data[f"{particle}_idx"][row] >= 0:
            print (f'    {particle:{max_len+1}s}: idx = {data[f"{particle}_idx"][row]:3.0f}, E = {data[f"{particle}_E"][row]:8.2f}, pT = {data[f"{particle}_pt"][row]:8.2f}, eta = {data[f"{particle}_eta"][row]:+6.2f}, phi = {data[f"{particle}_phi"][row]:+5.2f}, M = {data[f"{particle}_mass"][row]:6.2f}, pdgId = {data[f"{particle}_pdgId"][row]:+3.0f}, sum_E = {data[f"{particle}_sum_E"][row]:8.2f}')
    print (f'ISR : N = {data["n_ISR"][row]}')
    for i in range(1,data['n_ISR'][row]+1):
        if f"ISR_{i}_idx" in data.branches:
            print (f'    idx = {data[f"ISR_{i}_idx"][row]:3.0f}, E = {data[f"ISR_{i}_E"][row]:8.2f}, pT = {data[f"ISR_{i}_pt"][row]:8.2f}, eta = {data[f"ISR_{i}_eta"][row]:+6.2f}, phi = {data[f"ISR_{i}_phi"][row]:+5.2f}, M = {data[f"ISR_{i}_mass"][row]:6.2f}, pdgId = {data[f"ISR_{i}_pdgId"][row]:+3.0f}, parent = {data[f"ISR_{i}_parent"][row]:3.0f}')


def print_gen_content(row):
    for key in datas.keys():
        if 'gen' in key:
            print_gen_sample(key,row)

def print_all_content(row):
    print_reco_DL_tree(row)
    print_gen_content(row)
