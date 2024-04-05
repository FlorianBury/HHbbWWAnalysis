import sys
import math
import numpy as np
import uproot
import itertools
import awkward as ak
import ROOT
import argparse

from IPython import embed
from inspect_utils import Data

parser = argparse.ArgumentParser(description='Inspect nanoaod content')
parser.add_argument('-f','--files',nargs='*',type=str,required=True)
parser.add_argument('-t','--tree',type=str,required=True)
args = parser.parse_args()
data = Data(args.files,args.tree)

def recurse(event,idx):
    while True:
        print (f'Current particle : pdg id = {data["GenPart_pdgId"][event][idx]}, statusflag  = {data["GenPart_statusFlags"][event][idx]:016b}')
        idx = data['GenPart_genPartIdxMother'][event][idx]
        if idx < 0:
            print ('idx = 0')
            break

def get_LV(pt,eta,phi,M):
    lv = ROOT.Math.PtEtaPhiMVector()
    lv.SetPt(pt)
    lv.SetEta(eta)
    lv.SetPhi(phi)
    lv.SetM(M)
    return lv

def convert_PtEtaPhiM(pt,eta,phi,M):
    lv = get_LV(pt,eta,phi,M)
    return (lv.Px(),lv.Py(),lv.Pz(),lv.E())

status_list = [
 "isPrompt",
 "isDecayedLeptonHadron",
 "isTauDecayProduct",
 "isPromptTauDecayProduct",
 "isDirectTauDecayProduct",
 "isDirectPromptTauDecayProduct",
 "isDirectHadronDecayProduct",
 "isHardProcess",
 "fromHardProcess",
 "isHardProcessTauDecayProduct",
 "isDirectHardProcessTauDecayProduct",
 "fromHardProcessBeforeFSR",
 "isFirstCopy",
 "isLastCopy",
 "isLastCopyBeforeFSR",
]

def format_status(status):
    s = []
    for i,bit in enumerate(reversed(f'{status:015b}')):
        if int(bit):
            s.append(status_list[i])
    return s

def event_to_row(event):
    if event is not None:
        row = ak.where(data["event"] == event)[0]
        if len(row) == 0:
            raise RuntimeError(f'Event number {event} not found in tree')
        assert len(row) == 1
        row = row[0]
    return row


def print_gen_content(event):
    pdgids = data["GenPart_pdgId"][event]
    mother_idx = data["GenPart_genPartIdxMother"][event]
    pt = data["GenPart_pt"][event]
    eta = data["GenPart_eta"][event]
    phi = data["GenPart_phi"][event]
    m = data["GenPart_mass"][event]
    E = [convert_PtEtaPhiM(pt[i],eta[i],phi[i],m[i])[3] for i in range(len(pt))]
    status = data["GenPart_statusFlags"][event]
    print ('event',event)
    for i in range(data['nGenPart'][event]):
        print (f'idx = {i:3d} : pdgid = {pdgids[i]:+5d}, mother idx = {mother_idx[i]:3d}, E = {E[i]:8.2f}, pT = {pt[i]:8.2f}, eta = {eta[i]:+6.2f}, phi = {phi[i]:+5.2f}, M = {m[i]:6.2f}, {format_status(status[i])}')

def print_ak4_jets(event):
    print (f'N ak4 jets = {data["nJet"][event]}')
    print (f'pt     {" ".join([f"{v:8.2f}" for v in data["Jet_pt"][event]])}')
    print (f'eta    {" ".join([f"{v:8.2f}" for v in data["Jet_eta"][event]])}')
    print (f'phi    {" ".join([f"{v:8.2f}" for v in data["Jet_phi"][event]])}')
    print (f'parton {" ".join([f"{v:8.2f}" for v in data["Jet_partonFlavour"][event]])}')
    print (f'hadron {" ".join([f"{v:8.2f}" for v in data["Jet_hadronFlavour"][event]])}')
    print (f'genjet {" ".join([f"{v:8.2f}" for v in data["Jet_genJetIdx"][event]])}')

def print_ak4_genjets(event):
    print (f'N ak4 gen jets = {data["nGenJet"][event]}')
    print (f'pt     {" ".join([f"{v:8.2f}" for v in data["GenJet_pt"][event]])}')
    print (f'eta    {" ".join([f"{v:8.2f}" for v in data["GenJet_eta"][event]])}')
    print (f'phi    {" ".join([f"{v:8.2f}" for v in data["GenJet_phi"][event]])}')
    print (f'parton {" ".join([f"{v:8.2f}" for v in data["GenJet_partonFlavour"][event]])}')
    print (f'hadron {" ".join([f"{v:8.2f}" for v in data["GenJet_hadronFlavour"][event]])}')

def print_ak8_jets(event):
    print (f'N ak8 jets = {data["nFatJet"][event]}')
    print (f'pt     {" ".join([f"{v:8.2f}" for v in data["FatJet_pt"][event]])}')
    print (f'eta    {" ".join([f"{v:8.2f}" for v in data["FatJet_eta"][event]])}')
    print (f'phi    {" ".join([f"{v:8.2f}" for v in data["FatJet_phi"][event]])}')
    print (f'hadron {" ".join([f"{v:8.2f}" for v in data["FatJet_hadronFlavour"][event]])}')
    print (f'nB     {" ".join([f"{v:8.2f}" for v in data["FatJet_nBHadrons"][event]])}')
    print (f'nC     {" ".join([f"{v:8.2f}" for v in data["FatJet_nCHadrons"][event]])}')
    print (f'genjet {" ".join([f"{v:8.2f}" for v in data["FatJet_genJetAK8Idx"][event]])}')
    print (f'sub1   {" ".join([f"{v:8.2f}" for v in data["FatJet_subJetIdx1"][event]])}')
    print (f'sub2   {" ".join([f"{v:8.2f}" for v in data["FatJet_subJetIdx2"][event]])}')

def print_ak8_genjets(event):
    print (f'N ak8 gen jets = {data["nGenJetAK8"][event]}')
    print (f'pt     {" ".join([f"{v:8.2f}" for v in data["GenJetAK8_pt"][event]])}')
    print (f'eta    {" ".join([f"{v:8.2f}" for v in data["GenJetAK8_eta"][event]])}')
    print (f'phi    {" ".join([f"{v:8.2f}" for v in data["GenJetAK8_phi"][event]])}')
    print (f'parton {" ".join([f"{v:8.2f}" for v in data["GenJetAK8_partonFlavour"][event]])}')
    print (f'hadron {" ".join([f"{v:8.2f}" for v in data["GenJetAK8_hadronFlavour"][event]])}')

def print_sub_jets(event):
    print (f'N sub jets = {data["nSubJet"][event]}')
    print (f'pt     {" ".join([f"{v:8.2f}" for v in data["SubJet_pt"][event]])}')
    print (f'eta    {" ".join([f"{v:8.2f}" for v in data["SubJet_eta"][event]])}')
    print (f'phi    {" ".join([f"{v:8.2f}" for v in data["SubJet_phi"][event]])}')
    print (f'nB     {" ".join([f"{v:8.2f}" for v in data["SubJet_nBHadrons"][event]])}')
    print (f'nC     {" ".join([f"{v:8.2f}" for v in data["SubJet_nCHadrons"][event]])}')


def print_sub_genjets(event):
    print (f'N ak8 subgenjets = {data["nSubGenJetAK8"][event]}')
    print (f'pt     {" ".join([f"{v:8.2f}" for v in data["SubGenJetAK8_pt"][event]])}')
    print (f'eta    {" ".join([f"{v:8.2f}" for v in data["SubGenJetAK8_eta"][event]])}')
    print (f'phi    {" ".join([f"{v:8.2f}" for v in data["SubGenJetAK8_phi"][event]])}')



def print_muons(event):
    print (f'N muons = {data["nMuon"][event]}')
    print (f'pt      {" ".join([f"{v:8.2f}" for v in data["Muon_pt"][event]])}')
    print (f'eta     {" ".join([f"{v:8.2f}" for v in data["Muon_eta"][event]])}')
    print (f'phi     {" ".join([f"{v:8.2f}" for v in data["Muon_phi"][event]])}')
    print (f'flav    {" ".join([f"{v:8.2f}" for v in data["Muon_genPartFlav"][event]])}')
    print (f'gen idx {" ".join([f"{v:8.2f}" for v in data["Muon_genPartIdx"][event]])}')

def print_electrons(event):
    print (f'N electrons = {data["nElectron"][event]}')
    print (f'pt      {" ".join([f"{v:8.2f}" for v in data["Electron_pt"][event]])}')
    print (f'eta     {" ".join([f"{v:8.2f}" for v in data["Electron_eta"][event]])}')
    print (f'phi     {" ".join([f"{v:8.2f}" for v in data["Electron_phi"][event]])}')
    print (f'flav    {" ".join([f"{v:8.2f}" for v in data["Electron_genPartFlav"][event]])}')
    print (f'jet idx {" ".join([f"{v:8.2f}" for v in data["Electron_jetIdx"][event]])}')
    print (f'gen idx {" ".join([f"{v:8.2f}" for v in data["Electron_genPartIdx"][event]])}')
    print (f'gam idx {" ".join([f"{v:8.2f}" for v in data["Electron_photonIdx"][event]])}')


def print_taus(event):
    print (f'N taus = {data["nTau"][event]}')
    print (f'pt      {" ".join([f"{v:8.2f}" for v in data["Tau_pt"][event]])}')
    print (f'eta     {" ".join([f"{v:8.2f}" for v in data["Tau_eta"][event]])}')
    print (f'phi     {" ".join([f"{v:8.2f}" for v in data["Tau_phi"][event]])}')
    print (f'flav    {" ".join([f"{v:8.2f}" for v in data["Tau_genPartFlav"][event]])}')
    print (f'jet idx {" ".join([f"{v:8.2f}" for v in data["Tau_jetIdx"][event]])}')
    print (f'gen idx {" ".join([f"{v:8.2f}" for v in data["Tau_genPartIdx"][event]])}')



def print_dressed_leptons(event):
    print (f'N leptons = {data["nGenDressedLepton"][event]}')
    print (f'pt      {" ".join([f"{v:8.2f}" for v in data["GenDressedLepton_pt"][event]])}')
    print (f'eta     {" ".join([f"{v:8.2f}" for v in data["GenDressedLepton_eta"][event]])}')
    print (f'phi     {" ".join([f"{v:8.2f}" for v in data["GenDressedLepton_phi"][event]])}')
    print (f'pdgId   {" ".join([f"{v:8.2f}" for v in data["GenDressedLepton_pdgId"][event]])}')
    print (f'hasTau  {" ".join([f"{v:8.2f}" for v in data["GenDressedLepton_hasTauAnc"][event]])}')

def print_reco_content(idx):
    print ('\nLeptons')
    print_electrons(idx)
    print_muons(idx)
    print_taus(idx)
    print_dressed_leptons(idx)
    print ('\nAK4 jets')
    print_ak4_jets(idx)
    print_ak4_genjets(idx)
    print ('\nAK8 jets')
    print_ak8_jets(idx)
    print_ak8_genjets(idx)
    print ('\nSubjets')
    print_sub_jets(idx)
    print_sub_genjets(idx)

def print_all_content(idx):
    print_gen_content(idx)
    print_reco_content(idx)




def find_b_data():
    pdg = data["GenPart_pdgId"]
    for i in range(len(pdg)):
        if ak.sum(abs(pdg[i])==5):
            print (i)

def find_tau_data():
    pdg = data["GenPart_pdgId"]
    for i in range(len(pdg)):
        if ak.sum(abs(pdg[i])==15):
            print (i)

def find_preFSR_gluon(event,idx):
    gluon_idx = idx
    mother_idx = idx
    while True:
        mother_idx = data["GenPart_genPartIdxMother"][mother_idx]
        if data["GenPart_genPartIdxMother"][mother_idx] == -1:
            break
        else:
            gluon_idx = mother_idx
    return gluon_idx

def gluon_is_ISR(event,idx):
    mother_idxs = data["GenPart_genPartIdxMother"][event]
    pdgIds = data["GenPart_pdgId"][event]
    while True:
        idx = mother_idxs[idx]
        if mother_idxs[idx] == -1:
            return True
        if pdgIds[idx] != 21:
            return False
    raise RuntimeError

def deltaR(e1,p1,e2,p2):
    deta = e1-e2
    dphi = p1-p2
    if dphi > math.pi:
        dphi -= 2*math.pi
    elif dphi <= -math.pi:
        dphi += 2*math.pi
    return math.sqrt(deta**2+dphi**2)


def print_ISR_gluons_DR(event):
    statuses = data["GenPart_statusFlags"][event]
    pdgIds = data["GenPart_pdgId"][event]
    pts  = data["GenPart_pt"][event]
    etas = data["GenPart_eta"][event]
    phis = data["GenPart_phi"][event]
    mother_idxs = data["GenPart_genPartIdxMother"][event]
    gluon_idxs = []
    for idx in range(data['nGenPart'][event]):
        if pdgIds[idx] == 21 and statuses[idx] & ( 0x1 << 13 ):
            if gluon_is_ISR(event,idx):
                if pts[idx] >= 10:
                    print (f'idx {idx} is last copy gluon from ISR (pt = {pts[idx]:6.2f}, eta = {etas[idx]:+5.2f}, phi = {phis[idx]:+5.2f}, mother idx = {mother_idxs[idx]:3d})')
                    gluon_idxs.append(idx)
    for i1,i2 in itertools.combinations(gluon_idxs,2):
        print (f'Gluon idx {i1} (pt = {pts[i1]:8.2f}, eta = {etas[i1]:+5.2f}, phi = {phis[i1]:+5.2f}), Gluon idx {i2} (pt = {pts[i2]:8.2f}, eta = {etas[i2]:+5.2f}, phi = {phis[i2]:+5.2f}) : DR = {deltaR(etas[i1],phis[i1],etas[i2],phis[i2]):3.2f}')

def get_DR_pair(event,cont1,idx1,cont2,idx2):
    pt1  = data[f"{cont1}_pt"][event][idx1]
    eta1 = data[f"{cont1}_eta"][event][idx1]
    phi1 = data[f"{cont1}_phi"][event][idx1]
    pt2  = data[f"{cont2}_pt"][event][idx2]
    eta2 = data[f"{cont2}_eta"][event][idx2]
    phi2 = data[f"{cont2}_phi"][event][idx2]
    print (f'Object 1 ({cont1+")":15s}: idx = {idx1:3d} (pt = {pt1:8.2f}, eta = {eta1:+5.2f}, phi = {phi1:+5.2f})')
    print (f'Object 2 ({cont2+")":15s}: idx = {idx2:3d} (pt = {pt2:8.2f}, eta = {eta2:+5.2f}, phi = {phi2:+5.2f})')
    print (f'DR = {deltaR(eta1,phi1,eta2,phi2):.3f}')

def scan_AK8(N=math.inf):
    DRs = []
    for i in range(len(data['event'])):
        if i >= N:
            break
        leptons_etaphis = [
            (data["GenDressedLepton_eta"][i][j],data["GenDressedLepton_phi"][i][j])
                for j in range(data["nGenDressedLepton"][i])
        ]
        jets_etaphis = [
            (data["Jet_eta"][i][j],data["Jet_phi"][i][j])
                for j in range(data["nJet"][i])
                if data["Jet_genJetIdx"][i][j] >= 0
        ]
        print (f'Event {i}')
        for j in range(data['nFatJet'][i]):
            print (f'\tFatJet {j}')
            if data['FatJet_subJetIdx1'][i][j] == -1:
                continue
            if data['FatJet_subJetIdx2'][i][j] == -1:
                continue
            idx1 = data['FatJet_subJetIdx1'][i][j]
            idx2 = data['FatJet_subJetIdx2'][i][j]
            eta1 = data['SubJet_eta'][i][idx1]
            eta2 = data['SubJet_eta'][i][idx2]
            phi1 = data['SubJet_phi'][i][idx1]
            phi2 = data['SubJet_phi'][i][idx2]
            DR = deltaR(eta1,phi1,eta2,phi2)
            DRs.append(DR)
            print (f'\t\tSubjet DR = {DR:.3f}')

            isolated = True
            for eta_lep,phi_lep in leptons_etaphis:
                if deltaR(eta_lep,phi_lep,eta1,phi1) < 0.4:
                    isolated = False
                if deltaR(eta_lep,phi_lep,eta2,phi2) < 0.4:
                    isolated = False

            print (f'\t\tIsolation = {isolated}')

            match_sub1 = False
            match_sub2 = False
            for eta_jet,phi_jet in jets_etaphis:
                if deltaR(eta_jet,phi_jet,eta1,phi1) < 0.4:
                    match_sub1 = True
                    print (f'\t\t\tMatch sub 1 DR = {deltaR(eta_jet,phi_jet,eta1,phi1)}')
                if deltaR(eta_jet,phi_jet,eta2,phi2) < 0.4:
                    match_sub2 = True
                    print (f'\t\t\tMatch sub 2 DR = {deltaR(eta_jet,phi_jet,eta2,phi2)}')


            print (f'\t\tMatch sub1 = {match_sub1}')
            print (f'\t\tMatch sub2 = {match_sub2}')
            if isolated and (not match_sub1 or not match_sub2):
                embed()
    return DRs

embed()
