import os.path
import sys
import bamboo
from copy import deepcopy
from bamboo.analysismodules import NanoAODModule, NanoAODHistoModule, DataDrivenBackgroundHistogramsModule
from bamboo import treefunctions as op
from bamboo import treeoperations as to
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection
from bamboo.plots import Plot, EquidistantBinning, SummedPlot, Selection, Skim

sys.path.append(os.path.dirname(os.path.abspath(__file__))) # Add scripts in this directory

from BaseHHtobbWW import BaseNanoHHtobbWW
from selectionDef import *

def addGenParticle(varDict,container,name,idx=None,map_sum_E=None):
    if idx is None:
        varDict[f'n_{name}']      = op.static_cast("UInt_t",op.rng_len(container))
        idx = 0
    varDict[f'{name}_idx']    = op.switch(op.rng_len(container)>idx, op.static_cast("Int_t",container[idx].idx), op.c_float(-9999.))
    varDict[f'{name}_pdgId']  = op.switch(op.rng_len(container)>idx, container[idx].pdgId, op.c_float(-9999.))
    varDict[f'{name}_pt']     = op.switch(op.rng_len(container)>idx, container[idx].pt, op.c_float(-9999.))
    varDict[f'{name}_eta']    = op.switch(op.rng_len(container)>idx, container[idx].eta, op.c_float(-9999.))
    varDict[f'{name}_phi']    = op.switch(op.rng_len(container)>idx, container[idx].phi, op.c_float(-9999.))
    varDict[f'{name}_mass']   = op.switch(op.rng_len(container)>idx, container[idx].mass, op.c_float(-9999.))
    varDict[f'{name}_Px']     = op.switch(op.rng_len(container)>idx, container[idx].p4.Px(), op.c_float(-9999.))
    varDict[f'{name}_Py']     = op.switch(op.rng_len(container)>idx, container[idx].p4.Py(), op.c_float(-9999.))
    varDict[f'{name}_Pz']     = op.switch(op.rng_len(container)>idx, container[idx].p4.Pz(), op.c_float(-9999.))
    varDict[f'{name}_E']      = op.switch(op.rng_len(container)>idx, container[idx].p4.E(), op.c_float(-9999.))
    if map_sum_E is not None:
        varDict[f'{name}_sum_E'] = op.switch(op.rng_len(container)>idx, map_sum_E[container[idx].idx], op.c_float(-9999.))

def addISRParticles(varDict,container,rng):
    for i in range(rng):
        varDict[f'n_ISR']          = op.static_cast("UInt_t",op.rng_len(container))
        varDict[f'ISR_{i+1}_idx']    = op.switch(op.rng_len(container)>i, op.static_cast("Int_t",container[i].idx), op.c_float(-9999.))
        varDict[f'ISR_{i+1}_parent'] = op.switch(op.rng_len(container)>i, op.static_cast("Int_t",container[i].parent.idx), op.c_float(-9999.))
        varDict[f'ISR_{i+1}_pdgId']  = op.switch(op.rng_len(container)>i, container[i].pdgId, op.c_float(-9999.))
        varDict[f'ISR_{i+1}_pt']     = op.switch(op.rng_len(container)>i, container[i].pt, op.c_float(-9999.))
        varDict[f'ISR_{i+1}_eta']    = op.switch(op.rng_len(container)>i, container[i].eta, op.c_float(-9999.))
        varDict[f'ISR_{i+1}_phi']    = op.switch(op.rng_len(container)>i, container[i].phi, op.c_float(-9999.))
        varDict[f'ISR_{i+1}_mass']   = op.switch(op.rng_len(container)>i, container[i].mass, op.c_float(-9999.))
        varDict[f'ISR_{i+1}_Px']     = op.switch(op.rng_len(container)>i, container[i].p4.Px(), op.c_float(-9999.))
        varDict[f'ISR_{i+1}_Py']     = op.switch(op.rng_len(container)>i, container[i].p4.Py(), op.c_float(-9999.))
        varDict[f'ISR_{i+1}_Pz']     = op.switch(op.rng_len(container)>i, container[i].p4.Pz(), op.c_float(-9999.))
        varDict[f'ISR_{i+1}_E']      = op.switch(op.rng_len(container)>i, container[i].p4.E(), op.c_float(-9999.))

def addRecoElectron(varDict,container,i):
    varDict[f'e{i+1}_idx']    = op.switch(op.rng_len(container)>i, op.static_cast("Int_t",container[i].idx), op.c_float(-9999.))
    varDict[f'e{i+1}_pdgId']  = op.switch(op.rng_len(container)>i, container[i].pdgId, op.c_float(-9999.))
    varDict[f'e{i+1}_charge'] = op.switch(op.rng_len(container)>i, container[i].charge, op.c_float(-9999.))
    varDict[f'e{i+1}_pt']     = op.switch(op.rng_len(container)>i, container[i].pt, op.c_float(-9999.))
    varDict[f'e{i+1}_eta']    = op.switch(op.rng_len(container)>i, container[i].eta, op.c_float(-9999.))
    varDict[f'e{i+1}_phi']    = op.switch(op.rng_len(container)>i, container[i].phi, op.c_float(-9999.))
    varDict[f'e{i+1}_mass']   = op.switch(op.rng_len(container)>i, container[i].mass, op.c_float(-9999.))
    varDict[f'e{i+1}_Px']     = op.switch(op.rng_len(container)>i, container[i].p4.Px(), op.c_float(-9999.))
    varDict[f'e{i+1}_Py']     = op.switch(op.rng_len(container)>i, container[i].p4.Py(), op.c_float(-9999.))
    varDict[f'e{i+1}_Pz']     = op.switch(op.rng_len(container)>i, container[i].p4.Pz(), op.c_float(-9999.))
    varDict[f'e{i+1}_E']      = op.switch(op.rng_len(container)>i, container[i].p4.E(), op.c_float(-9999.))

def addRecoMuon(varDict,container,i):
    varDict[f'm{i+1}_idx']    = op.switch(op.rng_len(container)>i, op.static_cast("Int_t",container[i].idx), op.c_float(-9999.))
    varDict[f'm{i+1}_pdgId']  = op.switch(op.rng_len(container)>i, container[i].pdgId, op.c_float(-9999.))
    varDict[f'm{i+1}_charge'] = op.switch(op.rng_len(container)>i, container[i].charge, op.c_float(-9999.))
    varDict[f'm{i+1}_pt']     = op.switch(op.rng_len(container)>i, container[i].pt, op.c_float(-9999.))
    varDict[f'm{i+1}_eta']    = op.switch(op.rng_len(container)>i, container[i].eta, op.c_float(-9999.))
    varDict[f'm{i+1}_phi']    = op.switch(op.rng_len(container)>i, container[i].phi, op.c_float(-9999.))
    varDict[f'm{i+1}_mass']   = op.switch(op.rng_len(container)>i, container[i].mass, op.c_float(-9999.))
    varDict[f'm{i+1}_Px']     = op.switch(op.rng_len(container)>i, container[i].p4.Px(), op.c_float(-9999.))
    varDict[f'm{i+1}_Py']     = op.switch(op.rng_len(container)>i, container[i].p4.Py(), op.c_float(-9999.))
    varDict[f'm{i+1}_Pz']     = op.switch(op.rng_len(container)>i, container[i].p4.Pz(), op.c_float(-9999.))
    varDict[f'm{i+1}_E']      = op.switch(op.rng_len(container)>i, container[i].p4.E(), op.c_float(-9999.))

def addRecoJet(varDict,container,i,lambda_btag):
    varDict[f'j{i+1}_idx']    = op.switch(op.rng_len(container)>i, op.static_cast("Int_t",container[i].idx), op.c_float(-9999.))
    varDict[f'j{i+1}_btag']   = op.switch(op.rng_len(container)>i, container[i].btagDeepFlavB, op.c_float(-9999.))
    varDict[f'j{i+1}_btagged']= op.switch(op.rng_len(container)>i, lambda_btag(container[i]), op.c_float(-9999.))
    varDict[f'j{i+1}_pt']     = op.switch(op.rng_len(container)>i, container[i].pt, op.c_float(-9999.))
    varDict[f'j{i+1}_eta']    = op.switch(op.rng_len(container)>i, container[i].eta, op.c_float(-9999.))
    varDict[f'j{i+1}_phi']    = op.switch(op.rng_len(container)>i, container[i].phi, op.c_float(-9999.))
    varDict[f'j{i+1}_mass']   = op.switch(op.rng_len(container)>i, container[i].mass, op.c_float(-9999.))
    varDict[f'j{i+1}_Px']     = op.switch(op.rng_len(container)>i, container[i].p4.Px(), op.c_float(-9999.))
    varDict[f'j{i+1}_Py']     = op.switch(op.rng_len(container)>i, container[i].p4.Py(), op.c_float(-9999.))
    varDict[f'j{i+1}_Pz']     = op.switch(op.rng_len(container)>i, container[i].p4.Pz(), op.c_float(-9999.))
    varDict[f'j{i+1}_E']      = op.switch(op.rng_len(container)>i, container[i].p4.E(), op.c_float(-9999.))

def addVBFJet(varDict,container,i,map_sel):
    varDict[f'VBF{i+1}_idx']    = op.switch(op.rng_len(container)>i, op.static_cast("Int_t",container[i].idx), op.c_float(-9999.))
    varDict[f'VBF{i+1}_sel']    = op.switch(op.rng_len(container)>i, map_sel[container[i].idx], op.c_float(-9999.))
    varDict[f'VBF{i+1}_pt']     = op.switch(op.rng_len(container)>i, container[i].pt, op.c_float(-9999.))
    varDict[f'VBF{i+1}_eta']    = op.switch(op.rng_len(container)>i, container[i].eta, op.c_float(-9999.))
    varDict[f'VBF{i+1}_phi']    = op.switch(op.rng_len(container)>i, container[i].phi, op.c_float(-9999.))
    varDict[f'VBF{i+1}_mass']   = op.switch(op.rng_len(container)>i, container[i].mass, op.c_float(-9999.))
    varDict[f'VBF{i+1}_Px']     = op.switch(op.rng_len(container)>i, container[i].p4.Px(), op.c_float(-9999.))
    varDict[f'VBF{i+1}_Py']     = op.switch(op.rng_len(container)>i, container[i].p4.Py(), op.c_float(-9999.))
    varDict[f'VBF{i+1}_Pz']     = op.switch(op.rng_len(container)>i, container[i].p4.Pz(), op.c_float(-9999.))
    varDict[f'VBF{i+1}_E']      = op.switch(op.rng_len(container)>i, container[i].p4.E(), op.c_float(-9999.))

def addMET(varDict,met):
    varDict[f'met_pt']     = met.p4.Pt()
    varDict[f'met_eta']    = met.p4.Eta()
    varDict[f'met_phi']    = met.p4.Phi()
    varDict[f'met_mass']   = met.p4.M()
    varDict[f'met_Px']     = met.p4.Px()
    varDict[f'met_Py']     = met.p4.Py()
    varDict[f'met_Pz']     = met.p4.Pz()
    varDict[f'met_E']      = met.p4.E()


def addRecoFatjet(varDict,container,lambda_btag):
    varDict[f'fatj_idx']    = op.switch(op.rng_len(container)>0, op.static_cast("Int_t",container[0].idx), op.c_float(-9999.))
    varDict[f'fatj_pt']     = op.switch(op.rng_len(container)>0, container[0].pt, op.c_float(-9999.))
    varDict[f'fatj_eta']    = op.switch(op.rng_len(container)>0, container[0].eta, op.c_float(-9999.))
    varDict[f'fatj_phi']    = op.switch(op.rng_len(container)>0, container[0].phi, op.c_float(-9999.))
    varDict[f'fatj_mass']   = op.switch(op.rng_len(container)>0, container[0].mass, op.c_float(-9999.))
    varDict[f'fatj_Px']     = op.switch(op.rng_len(container)>0, container[0].p4.Px(), op.c_float(-9999.))
    varDict[f'fatj_Py']     = op.switch(op.rng_len(container)>0, container[0].p4.Py(), op.c_float(-9999.))
    varDict[f'fatj_Pz']     = op.switch(op.rng_len(container)>0, container[0].p4.Pz(), op.c_float(-9999.))
    varDict[f'fatj_E']      = op.switch(op.rng_len(container)>0, container[0].p4.E(), op.c_float(-9999.))

    varDict[f'fatj_sub1_idx']    = op.switch(op.rng_len(container)>0, op.static_cast("Int_t",container[0].subJet1.idx), op.c_float(-9999.))
    varDict[f'fatj_sub1_btag']   = op.switch(op.rng_len(container)>0, container[0].subJet1.btagDeepB, op.c_float(-9999.))
    varDict[f'fatj_sub1_btagged']= op.switch(op.rng_len(container)>0, lambda_btag(container[0].subJet1), op.c_float(-9999.))
    varDict[f'fatj_sub1_pt']     = op.switch(op.rng_len(container)>0, container[0].subJet1.pt, op.c_float(-9999.))
    varDict[f'fatj_sub1_eta']    = op.switch(op.rng_len(container)>0, container[0].subJet1.eta, op.c_float(-9999.))
    varDict[f'fatj_sub1_phi']    = op.switch(op.rng_len(container)>0, container[0].subJet1.phi, op.c_float(-9999.))
    varDict[f'fatj_sub1_mass']   = op.switch(op.rng_len(container)>0, container[0].subJet1.mass, op.c_float(-9999.))
    varDict[f'fatj_sub1_Px']     = op.switch(op.rng_len(container)>0, container[0].subJet1.p4.Px(), op.c_float(-9999.))
    varDict[f'fatj_sub1_Py']     = op.switch(op.rng_len(container)>0, container[0].subJet1.p4.Py(), op.c_float(-9999.))
    varDict[f'fatj_sub1_Pz']     = op.switch(op.rng_len(container)>0, container[0].subJet1.p4.Pz(), op.c_float(-9999.))
    varDict[f'fatj_sub1_E']      = op.switch(op.rng_len(container)>0, container[0].subJet1.p4.E(), op.c_float(-9999.))

    varDict[f'fatj_sub2_idx']    = op.switch(op.rng_len(container)>0, op.static_cast("Int_t",container[0].subJet2.idx), op.c_float(-9999.))
    varDict[f'fatj_sub2_btag']   = op.switch(op.rng_len(container)>0, container[0].subJet2.btagDeepB, op.c_float(-9999.))
    varDict[f'fatj_sub2_btagged']= op.switch(op.rng_len(container)>0, lambda_btag(container[0].subJet2), op.c_float(-9999.))
    varDict[f'fatj_sub2_pt']     = op.switch(op.rng_len(container)>0, container[0].subJet2.pt, op.c_float(-9999.))
    varDict[f'fatj_sub2_eta']    = op.switch(op.rng_len(container)>0, container[0].subJet2.eta, op.c_float(-9999.))
    varDict[f'fatj_sub2_phi']    = op.switch(op.rng_len(container)>0, container[0].subJet2.phi, op.c_float(-9999.))
    varDict[f'fatj_sub2_mass']   = op.switch(op.rng_len(container)>0, container[0].subJet2.mass, op.c_float(-9999.))
    varDict[f'fatj_sub2_Px']     = op.switch(op.rng_len(container)>0, container[0].subJet2.p4.Px(), op.c_float(-9999.))
    varDict[f'fatj_sub2_Py']     = op.switch(op.rng_len(container)>0, container[0].subJet2.p4.Py(), op.c_float(-9999.))
    varDict[f'fatj_sub2_Pz']     = op.switch(op.rng_len(container)>0, container[0].subJet2.p4.Pz(), op.c_float(-9999.))
    varDict[f'fatj_sub2_E']      = op.switch(op.rng_len(container)>0, container[0].subJet2.p4.E(), op.c_float(-9999.))

def getGenParticle(tree,pdgIds=[],parent=None,parent_pdgIds=[],from_hardscattering=False,flags=[(0x1 << 0), (0x1 << 12)]):
    # Make pdgIds into a list #
    if not isinstance(pdgIds,(list,tuple)):
        pdgIds = [pdgIds]
    # statusflags #
    if len(flags) == 0:
        lambda_status = lambda p : op.c_bool(True)
    else:
        lambda_status = lambda p : op.AND(*[p.statusFlags & flag for flag in flags])
    # Lambda for parent in ancestors #
    if parent is None:
        lambda_check_parents = lambda a: op.c_bool(True)
    else:
        if isinstance(parent,bamboo.treeproxies.SelectionProxy):
            lambda_check_parents = lambda a: op.rng_any(
                parent,
                lambda par : a.idx == par.idx
            )
        else:
            lambda_check_parents = lambda a: a.idx == parent.idx
    # Lambda for parent pdgIds as direct parent #
    if isinstance(parent_pdgIds,int):
        lambda_check_parent_pdgIds = lambda a: p.parent.pdgId == parent_pdgIds
    elif isinstance(parent_pdgIds,(list,tuple)):
        if len(parent_pdgIds) == 0:
            lambda_check_parent_pdgIds = lambda p: op.c_bool(True)
        else:
            lambda_check_parent_pdgIds = lambda p: op.OR(*[p.parent.pdgId == pdgId for pdgId in parent_pdgIds])
    else:
        raise RuntimeError(f'Type of parent pdgId {type(parent_pdgIds)} not expected')
    # Hard scattering lambda #
    if from_hardscattering:
        lambda_scatter = lambda p: op.AND(
            op.abs(p.parent.eta) > 10000,
        )
    else:
        lambda_scatter = lambda p: op.abs(p.parent.eta) < 10000
    # return select #
    return op.select(
        tree.GenPart,
        lambda p : op.AND(
            op.OR(
                *[p.pdgId == pdgId for pdgId in pdgIds]
            ),
            lambda_status(p),
            lambda_scatter(p),
            op.rng_any(
                p.ancestors,
                lambda_check_parents,
            ),
            lambda_check_parent_pdgIds(p),
        )
    )

def getGenISR(tree,pdgIds=[]):
    return op.select(
        tree.GenPart,
        lambda p : op.AND(
            # PDG ID
            op.OR(
                *[p.pdgId == pdgId for pdgId in pdgIds]
            ),
            p.pt >= 10,
            # Origin
            op.OR(
                # Directly from the beam
                op.AND(
                    p.parent.parent.idx < 0,
                    op.abs(p.parent.eta) > 10000,
                ),
                # ISR
                op.AND(
                    p.parent.idx < 0,
                    op.abs(p.eta) < 6.,
                )
            )
        )
    )





class TransfermerNanoAOD(BaseNanoHHtobbWW,DataDrivenBackgroundHistogramsModule):
    def definePlots(self, tree, noSel, sample=None, sampleCfg=None):
        noSel = super().prepareObjects(tree, noSel, sample, sampleCfg, 'DL')
        noSel = self.beforeJetselection(noSel)
        era = sampleCfg['era']
        plots = []

        # protection against data #
        if not self.is_MC:
            return plots

        # Pre compute sum of descendance #
        map_sum_E = op.map(
            tree.GenPart,
            lambda p: op.rng_sum(
                tree.GenPart,
                lambda d: op.switch(
                    # if particle is direct descendant of pre-FSR, include energy, otherwise a 0
                    d.parent.idx == p.idx,
                    d.p4.E(),
                    op.c_float(0.)
                ),
            )
        )
        def get_lambda_valid(pObjs):
            return op.NOT(
                op.OR(
                    *[
                        op.rng_any(
                            pObj,
                            lambda p : p.p4.E() <= 0.99 * map_sum_E[p.idx]
                        )
                        for pObj in pObjs
                    ]
                )
            )

        ##########################################################
        #                       Reco tree                        #
        ##########################################################
        varDL = {
            'event': None,
        }
        # Make selections and record the flags #
        lepSelObjs = makeDoubleLeptonSelection(self,noSel,use_dd=False,fake_selection=False)
        #from IPython import embed; embed()
        channels = ['ElEl','MuMu','ElMu']
        for channel, selObj in zip(channels,lepSelObjs):
            varDL[f'flag_{channel}'] = op.switch(op.AND(*selObj.sel._cuts), op.c_int(1), op.c_int(0))

        resCut = op.AND(
            op.rng_len(self.ak4Jets)>=2,
            op.rng_len(self.ak4BJets)>=1,
            op.rng_len(self.ak8BJets)==0,
        )
        booCut = op.AND(
            op.rng_len(self.ak8Jets)>=1,
            op.rng_len(self.ak8BJets)>=1,
        )
        SRCut = op.AND(
            op.OR(
                op.AND(*lepSelObjs[0].sel._cuts),
                op.AND(*lepSelObjs[1].sel._cuts),
                op.AND(*lepSelObjs[2].sel._cuts),
            ),
            op.OR(
                resCut,
                booCut,
            )
        )

        varDL[f'flag_resolved'] = op.switch(resCut,op.c_int(1),op.c_int(0))
        varDL[f'flag_boosted'] = op.switch(booCut,op.c_int(1),op.c_int(0))
        varDL[f'flag_SR'] = op.switch(SRCut,op.c_int(1),op.c_int(0))

        # Counting #
        varDL["n_e"]    = op.static_cast("UInt_t",op.rng_len(self.electronsTightSel))
        varDL["n_m"]    = op.static_cast("UInt_t",op.rng_len(self.muonsTightSel))
        varDL["n_AK4"]  = op.static_cast("UInt_t",op.rng_len(self.ak4Jets))
        varDL["n_AK4B"] = op.static_cast("UInt_t",op.rng_len(self.ak4BJets))
        varDL["n_AK8"]  = op.static_cast("UInt_t",op.rng_len(self.ak8Jets))
        varDL["n_AK8B"] = op.static_cast("UInt_t",op.rng_len(self.ak8BJets))

        # Leptons #
        for i in range(4):
            addRecoElectron(varDL,self.electronsTightSel,i)
            addRecoMuon(varDL,self.muonsTightSel,i)

        # Jets #
        map_clean_jets = op.map(
            tree.Jet,
            lambda j: op.NOT(
                op.OR(
                    op.rng_any(
                        self.electronsTightSel,
                        lambda e : op.deltaR(e.p4,j.p4) < 0.4,
                    ),
                    op.rng_any(
                        self.muonsTightSel,
                        lambda m : op.deltaR(m.p4,j.p4) < 0.4,
                    ),
                    op.rng_any(
                        self.ak8BJets,
                        lambda f : op.deltaR(f.p4,j.p4) < 1.2,
                    ),
                )
            )
        )
        map_clean_fatjets = op.map(
            tree.Jet,
            lambda j: op.NOT(
                op.OR(
                    op.rng_any(
                        self.electronsTightSel,
                        lambda e : op.deltaR(e.p4,j.p4) < 0.4,
                    ),
                    op.rng_any(
                        self.muonsTightSel,
                        lambda m : op.deltaR(m.p4,j.p4) < 0.4,
                    ),
                )
            )
        )

        jets = op.select(self.ak4JetsByBtagScore,lambda j : map_clean_jets[j.idx])
        fatjets = op.select(self.ak8BJets,lambda j : map_clean_fatjets[j.idx])
        for i in range(15):
            addRecoJet(varDL,jets,i,self.lambda_ak4Btag)

        # Fatjets #
        addRecoFatjet(varDL,self.ak8BJets,self.lambda_subjetBtag)

        # Select VBF jets (remove AK4 double counting) #
        VBF = op.select(
            self.VBFJets,
            lambda j: op.AND(
                map_clean_jets[j.idx],
                op.NOT(
                    op.rng_any(
                        jets,
                        lambda ak4: j.idx == ak4.idx,
                    )
                )
            )
        )
        varDL["n_VBF"] = op.static_cast("UInt_t",op.rng_len(VBF))
        map_sel_VBF = op.map(
            tree.Jet,
            lambda j : op.multiSwitch(
                # Boosted
                (
                    op.AND(
                        booCut,
                        op.rng_len(self.VBFJetPairsBoosted)>0,
                    ),
                    op.OR(
                        j.idx == self.VBFJetPairsBoosted[0][0].idx,
                        j.idx == self.VBFJetPairsBoosted[0][1].idx,
                    ),
                ),
                # Resolved
                (
                    op.AND(
                        resCut,
                        op.rng_len(self.VBFJetPairsResolved)>0,
                    ),
                    op.OR(
                        j.idx == self.VBFJetPairsResolved[0][0].idx,
                        j.idx == self.VBFJetPairsResolved[0][1].idx,
                    ),
                ),
                # Default
                op.c_bool(False),
            )
        )
        for i in range(8):
            addVBFJet(varDL,VBF,i,map_sel_VBF)

        # MET #
        addMET(varDL,self.corrMET)


        # Save into plots #
        plots.append(Skim(f'reco_DL',varDL,noSel))


        ##########################################################
        #                      ttbar sample                      #
        ##########################################################
        if sample in ['TTTo2L2Nu','TTToSemiLeptonic','TTToHadronic']:
            top = getGenParticle(tree,6,from_hardscattering=True)
            antitop = getGenParticle(tree,-6,from_hardscattering=True)
            bottom = getGenParticle(tree,5,top)
            antibottom = getGenParticle(tree,-5,antitop)
            W_plus_from_top = getGenParticle(tree,+24,top)
            W_minus_from_antitop = getGenParticle(tree,-24,antitop)
            lep_plus = getGenParticle(tree,[-11,-13,-15],W_plus_from_top)
            lep_minus = getGenParticle(tree,[+11,+13,+15],W_minus_from_antitop)
            neutrino = getGenParticle(tree,[+12,+14,+16],W_plus_from_top)
            antineutrino = getGenParticle(tree,[-12,-14,-16],W_minus_from_antitop)
            quark_up = getGenParticle(tree,[+2,+4],W_plus_from_top)
            quark_down = getGenParticle(tree,[+1,+3],W_minus_from_antitop)
            antiquark_up = getGenParticle(tree,[-2,-4],W_minus_from_antitop)
            antiquark_down = getGenParticle(tree,[-1,-3],W_plus_from_top)

            particles = {
                'top': top,
                'antitop': antitop,
                'bottom': bottom,
                'antibottom': antibottom,
                'W_plus_from_top': W_plus_from_top,
                'W_minus_from_antitop': W_minus_from_antitop,
                'lep_plus': lep_plus,
                'lep_minus': lep_minus,
                'neutrino': neutrino,
                'antineutrino': antineutrino,
                'quark_up': quark_up,
                'quark_down': quark_down,
                'antiquark_up': antiquark_up,
                'antiquark_down': antiquark_down,
            }

            genTTVar = {'event': None}
            for pName,pObj in particles.items():
                addGenParticle(genTTVar,pObj,pName,map_sum_E=map_sum_E)
            ISR = getGenISR(tree,[+1,+2,+3,+4,+5,-1,-2,-3,-4,-5,21])
            addISRParticles(genTTVar,ISR,15)
            genTTVar['flag_valid'] = op.static_cast("Int_t",get_lambda_valid(particles.values()))
            plots.append(Skim(f'gen_TT',genTTVar,noSel))

        ##########################################################
        #                  single top sample                     #
        ##########################################################
        if sample == 'ST_tW_top_5f' or sample == 'ST_tW_antitop_5f':
            genSTVar = {'event': None}

            # ST_tW_top_5f #
            top = getGenParticle(tree,6,from_hardscattering=True)
            bottom = getGenParticle(tree,+5,top)
            W_plus_from_top = getGenParticle(tree,+24,top)
            lep_plus_from_top = getGenParticle(tree,[-11,-13,-15],W_plus_from_top)
            neutrino_from_top = getGenParticle(tree,[+12,+14,+16],W_plus_from_top)
            quark_up_from_top = getGenParticle(tree,[+2,+4],W_plus_from_top)
            antiquark_down_from_top = getGenParticle(tree,[-1,-3],W_plus_from_top)

            W_minus_prompt = getGenParticle(tree,-24,from_hardscattering=True)
            lep_minus_from_prompt_W = getGenParticle(tree,[+11,+13,+15],W_minus_prompt)
            antineutrino_from_prompt_W = getGenParticle(tree,[-12,-14,-16],W_minus_prompt)
            quark_down_from_prompt_W = getGenParticle(tree,[+1,+3],W_minus_prompt)
            antiquark_up_from_prompt_W = getGenParticle(tree,[-2,-4],W_minus_prompt)

            # ST_tW_antitop_5f #
            antitop = getGenParticle(tree,-6,from_hardscattering=True)
            antibottom = getGenParticle(tree,-5,antitop)
            W_minus_from_antitop = getGenParticle(tree,-24,antitop)
            lep_minus_from_top = getGenParticle(tree,[+11,+13,+15],W_minus_from_antitop)
            antineutrino_from_top = getGenParticle(tree,[-12,-14,-16],W_minus_from_antitop)
            quark_down_from_top = getGenParticle(tree,[+1,+3],W_minus_from_antitop)
            antiquark_up_from_top = getGenParticle(tree,[-2,-4],W_minus_from_antitop)

            W_plus_prompt = getGenParticle(tree,+24,from_hardscattering=True)
            lep_plus_from_prompt_W = getGenParticle(tree,[-11,-13,-15],W_plus_prompt)
            neutrino_from_prompt_W = getGenParticle(tree,[+12,+14,+16],W_plus_prompt)
            quark_up_from_prompt_W = getGenParticle(tree,[+2,+4],W_plus_prompt)
            antiquark_down_from_prompt_W = getGenParticle(tree,[-1,-3],W_plus_prompt)

            # Fill tree #
            particles = {
                # ST_tW_top_5f
                'top': top,
                'bottom': bottom,
                'W_plus_from_top': W_plus_from_top,
                'lep_plus_from_top': lep_plus_from_top,
                'neutrino_from_top': neutrino_from_top,
                'quark_up_from_top': quark_up_from_top,
                'antiquark_down_from_top': antiquark_down_from_top,
                'W_minus_prompt': W_minus_prompt,
                'lep_minus_from_prompt_W': lep_minus_from_prompt_W,
                'antineutrino_from_prompt_W': antineutrino_from_prompt_W,
                'quark_down_from_prompt_W': quark_down_from_prompt_W,
                'antiquark_up_from_prompt_W': antiquark_up_from_prompt_W,
                # ST_tW_antitop_5f
                'antitop': antitop,
                'antibottom': antibottom,
                'W_minus_from_antitop': W_minus_from_antitop,
                'lep_minus_from_top': lep_minus_from_top,
                'antineutrino_from_top': antineutrino_from_top,
                'quark_down_from_top': quark_down_from_top,
                'antiquark_up_from_top': antiquark_up_from_top,
                'W_plus_prompt': W_plus_prompt,
                'lep_plus_from_prompt_W': lep_plus_from_prompt_W,
                'neutrino_from_prompt_W': neutrino_from_prompt_W,
                'quark_up_from_prompt_W': quark_up_from_prompt_W,
                'antiquark_down_from_prompt_W': antiquark_down_from_prompt_W,
            }
            genSTVar = {'event': None}
            for pName,pObj in particles.items():
                addGenParticle(genSTVar,pObj,pName,map_sum_E=map_sum_E)
            ISR = getGenISR(tree,[+1,+2,+3,+4,+5,-1,-2,-3,-4,-5,21])
            addISRParticles(genSTVar,ISR,15)
            genSTVar['flag_valid'] = op.static_cast("Int_t",get_lambda_valid(particles.values()))
            plots.append(Skim(f'gen_ST',genSTVar,noSel))
            #if sample in ['ST_tW_top_5f','ST_tW_antitop_5f','ST_tchannel_antitop_4f','ST_tchannel_top_4f','ST_schannel_4f','ST_tWll_5f']:

        ##########################################################
        #                        Drell-Yan                       #
        ##########################################################
        if sample.startswith('DY'):
            Z = getGenParticle(tree,+23,from_hardscattering=True)
            lep_from_Z = getGenParticle(tree,[11,12,13,14,15,16],Z)
            antilep_from_Z = getGenParticle(tree,[-11,-12,-13,-14,-15,-16],Z)
            quark_from_Z = getGenParticle(tree,[1,2,3,4,5],Z)
            antiquark_from_Z = getGenParticle(tree,[-1,-2,-3,-4,-5],Z)

            lep_from_nonres = getGenParticle(tree,[11,12,13,14,15,16],from_hardscattering=True)
            antilep_from_nonres = getGenParticle(tree,[-11,-12,-13,-14,-15,-16],from_hardscattering=True)
            quark_from_nonres = getGenParticle(tree,[1,2,3,4,5],from_hardscattering=True)
            antiquark_from_nonres = getGenParticle(tree,[-1,-2,-3,-4,-5],from_hardscattering=True)

            genDYVar = {'event': None}
            # Fill tree #
            particles = {
                'Z': Z,
                'lep_from_Z': lep_from_Z,
                'antilep_from_Z': antilep_from_Z,
                'quark_from_Z': quark_from_Z,
                'antiquark_from_Z': antiquark_from_Z,
                'lep_from_nonres': lep_from_nonres,
                'antilep_from_nonres': antilep_from_nonres,
                'quark_from_nonres': quark_from_nonres,
                'antiquark_from_nonres': antiquark_from_nonres,
            }
            genDYVar = {'event': None}
            for pName,pObj in particles.items():
                addGenParticle(genDYVar,pObj,pName,map_sum_E=map_sum_E)
            ISR = getGenISR(tree,[+1,+2,+3,+4,+5,-1,-2,-3,-4,-5,21])
            addISRParticles(genDYVar,ISR,15)
            genDYVar['flag_valid'] = op.static_cast("Int_t",get_lambda_valid(particles.values()))
            plots.append(Skim(f'gen_DY',genDYVar,noSel))

        ##########################################################
        #                    diboson sample                      #
        ##########################################################
        if sample.startswith("ZZTo"):
            Zs = getGenParticle(tree,+23,from_hardscattering=True)

            lep_from_Z         = getGenParticle(tree,[11,12,13,14,15,16],Zs)
            antilep_from_Z     = getGenParticle(tree,[-11,-12,-13,-14,-15,-16],Zs)
            quark_from_Z       = getGenParticle(tree,[1,2,3,4,5],Zs)
            antiquark_from_Z   = getGenParticle(tree,[-1,-2,-3,-4,-5],Zs)

            # Want to split the p=leptons/quarks into Z1 and Z2
            # If there is a p from Z -> at least one Z
            # So for the p from Z1, just need to check the first one
            # For p from Z2, need to make sure there are two before checking it
            # if N(Z)<2, then no way it can come from Z2 so just return false
            def split_from_Z(cont_from_Z):
                cont_from_Z1 = op.select(
                    cont_from_Z,
                    lambda p : op.rng_any(
                        p.ancestors,
                        lambda a: a.idx == Zs[0].idx
                    )
                )
                cont_from_Z2 = op.select(
                    cont_from_Z,
                    lambda p : op.switch(
                        op.rng_len(Zs) == 2,
                        op.rng_any(
                            p.ancestors,
                            lambda a: a.idx == Zs[1].idx
                        ),
                        op.c_bool(False),
                    )
                )
                return cont_from_Z1,cont_from_Z2

            lep_from_Z1,lep_from_Z2             = split_from_Z(lep_from_Z)
            antilep_from_Z1,antilep_from_Z2     = split_from_Z(antilep_from_Z)
            quark_from_Z1,quark_from_Z2         = split_from_Z(quark_from_Z)
            antiquark_from_Z1,antiquark_from_Z2 = split_from_Z(antiquark_from_Z)

            lep_from_nonres         = getGenParticle(tree,[11,12,13,14,15,16],from_hardscattering=True)
            antilep_from_nonres     = getGenParticle(tree,[-11,-12,-13,-14,-15,-16],from_hardscattering=True)
            quark_from_nonres       = getGenParticle(tree,[1,2,3,4,5],from_hardscattering=True)
            antiquark_from_nonres   = getGenParticle(tree,[-1,-2,-3,-4,-5],from_hardscattering=True)

            # Fill tree #
            particles = {
                'Zs' : Zs,
                'lep_from_Z1': lep_from_Z1,
                'antilep_from_Z1': antilep_from_Z1,
                'quark_from_Z1': quark_from_Z1,
                'antiquark_from_Z1': antiquark_from_Z1,
                'lep_from_Z2': lep_from_Z2,
                'antilep_from_Z2': antilep_from_Z2,
                'quark_from_Z2': quark_from_Z2,
                'antiquark_from_Z2': antiquark_from_Z2,
                'lep_from_nonres': lep_from_nonres,
                'antilep_from_nonres': antilep_from_nonres,
                'quark_from_nonres': quark_from_nonres,
                'antiquark_from_nonres': antiquark_from_nonres,
            }
            genZZVar = {'event': None}
            for pName,pObj in particles.items():
                if pName == 'Zs':
                    addGenParticle(genZZVar,pObj,'Z1',idx=0,map_sum_E=map_sum_E)
                    addGenParticle(genZZVar,pObj,'Z2',idx=1,map_sum_E=map_sum_E)
                else:
                    addGenParticle(genZZVar,pObj,pName,map_sum_E=map_sum_E)
            ISR = getGenISR(tree,[+1,+2,+3,+4,+5,-1,-2,-3,-4,-5,21])
            addISRParticles(genZZVar,ISR,15)
            genZZVar['flag_valid'] = op.static_cast("Int_t",get_lambda_valid(particles.values()))
            plots.append(Skim(f'gen_ZZ',genZZVar,noSel))


        ##########################################################
        #                      single Higgs                      #
        ##########################################################
        if sample.startswith("ZH_HToBB_ZToLL"):
            H = getGenParticle(tree,+25,from_hardscattering=True)
            Z = getGenParticle(tree,+23,from_hardscattering=True)

            bottom = getGenParticle(tree,+5,H)
            antibottom = getGenParticle(tree,-5,H)

            lep_plus = getGenParticle(tree,[-11,-13,-15],Z)
            lep_minus = getGenParticle(tree,[+11,+13,+15],Z)

            # Fill tree #
            particles = {
                'H': H,
                'Z': Z,
                'bottom': bottom,
                'antibottom': antibottom,
                'lep_plus': lep_plus,
                'lep_minus': lep_minus,
            }
            genZHVar = {'event': None}
            for pName,pObj in particles.items():
                addGenParticle(genZHVar,pObj,pName,map_sum_E=map_sum_E)
            ISR = getGenISR(tree,[+1,+2,+3,+4,+5,-1,-2,-3,-4,-5,21])
            addISRParticles(genZHVar,ISR,15)
            genZHVar['flag_valid'] = op.static_cast("Int_t",get_lambda_valid(particles.values()))
            plots.append(Skim(f'gen_ZH',genZHVar,noSel))



        ##########################################################
        #                           HH                           #
        ##########################################################
        if 'HH' in sample:
            genHHVar = {'event': None}
            particles = {}
            if sample.startswith("GluGluToHH") or sample.startswith('VBFHH'):
                higgs = getGenParticle(tree,+25,from_hardscattering=True)
            elif sample.startswith('GluGluToBulkGravitonToHH'):
                graviton = getGenParticle(tree,+39,from_hardscattering=True)
                higgs = getGenParticle(tree,+25,parent=graviton)
                particles['graviton'] = graviton
            elif sample.startswith('GluGluToRadionToHH'):
                radion = getGenParticle(tree,+35,from_hardscattering=True)
                higgs = getGenParticle(tree,+25,radion)
                particles['radion'] = radion
            else:
                raise RuntimeError(f'HH sample name {sample} not implemented')
            bottom = getGenParticle(tree,+5,higgs)
            antibottom = getGenParticle(tree,-5,higgs)

            particles = {
                **particles,
                'higgs': higgs,
                'bottom': bottom,
                'antibottom': antibottom,
            }

            if 'HHTo2B2Tau' in sample:
                tau = getGenParticle(tree,[+15],higgs)
                antitau = getGenParticle(tree,[-15],higgs)
                particles = {
                    **particles,
                    'tau': tau,
                    'antitau': antitau,
                }
            elif 'HHTo2B2VTo2L2Nu' in sample or 'HHTo2B2WToLNu2J' in sample or 'HHTo2B2VLNu2J' in sample:
                W_plus = getGenParticle(tree,+24,parent=higgs)
                W_minus = getGenParticle(tree,-24,parent=higgs)
                lep_plus_from_W = getGenParticle(tree,[-11,-13,-15],parent=W_plus)
                lep_minus_from_W = getGenParticle(tree,[+11,+13,+15],parent=W_minus)
                neutrino_from_W = getGenParticle(tree,[+12,+14,+16],parent=W_plus)
                antineutrino_from_W = getGenParticle(tree,[-12,-14,-16],parent=W_minus)
                quark_up_from_W = getGenParticle(tree,[+2,+4],parent=W_plus)
                quark_down_from_W = getGenParticle(tree,[+1,+3],parent=W_minus)
                antiquark_up_from_W = getGenParticle(tree,[-2,-4],parent=W_minus)
                antiquark_down_from_W = getGenParticle(tree,[-1,-3],parent=W_plus)

                Zs = getGenParticle(tree,+23,parent=higgs)
                lep_plus_from_Z = getGenParticle(tree,[-11,-13,-15],parent=Zs)
                lep_minus_from_Z = getGenParticle(tree,[+11,+13,+15],parent=Zs)
                neutrino_from_Z = getGenParticle(tree,[+12,+14,+16],parent=Zs)
                antineutrino_from_Z = getGenParticle(tree,[-12,-14,-16],parent=Zs)
                quark_up_from_Z = getGenParticle(tree,[+2,+4],parent=Zs)
                quark_down_from_Z = getGenParticle(tree,[+1,+3],parent=Zs)
                antiquark_up_from_Z = getGenParticle(tree,[-2,-4],parent=Zs)
                antiquark_down_from_Z = getGenParticle(tree,[-1,-3],parent=Zs)

                particles = {
                    **particles,
                    'W_plus': W_plus,
                    'W_minus': W_minus,
                    'Zs': Zs,
                    'lep_plus_from_W': lep_plus_from_W,
                    'lep_minus_from_W': lep_minus_from_W,
                    'neutrino_from_W': neutrino_from_W,
                    'antineutrino_from_W': antineutrino_from_W,
                    'quark_up_from_W': quark_up_from_W,
                    'quark_down_from_W': quark_down_from_W,
                    'antiquark_up_from_W': antiquark_up_from_W,
                    'antiquark_down_from_W': antiquark_down_from_W,
                    'lep_plus_from_Z': lep_plus_from_Z,
                    'lep_minus_from_Z': lep_minus_from_Z,
                    'neutrino_from_Z': neutrino_from_Z,
                    'antineutrino_from_Z': antineutrino_from_Z,
                    'quark_up_from_Z': quark_up_from_Z,
                    'quark_down_from_Z': quark_down_from_Z,
                    'antiquark_up_from_Z': antiquark_up_from_Z,
                    'antiquark_down_from_Z': antiquark_down_from_Z,
                }
            else:
                raise RuntimeError(f'HH sample name {sample} not implemented')

            ISR = getGenISR(tree,[+1,+2,+3,+4,+5,-1,-2,-3,-4,-5,21])
            if sample.startswith('VBFHH'):
                VBF_quark = getGenParticle(
                    tree,
                    pdgIds = [+1,+2,+3,+4,+5,-1,-2,-3,-4,-5],
                    flags = [(0x1 << 0),(0x1 << 7)], # isPrompt and isHardProcess
                    parent_pdgIds = [+1,+2,+3,+4,+5,-1,-2,-3,-4,-5],
                    from_hardscattering=True,
                )
                particles = {
                    **particles,
                    'VBF_quark': VBF_quark,
                }
                ISR = op.select(
                    ISR,
                    lambda qi: op.NOT(
                        op.rng_any(
                            VBF_quark,
                            lambda qv: qi.idx == qv.idx
                        )
                    )
                )
            addISRParticles(genHHVar,ISR,15)
            for pName,pObj in particles.items():
                if pName == 'higgs':
                    addGenParticle(genHHVar,pObj,'H1',idx=0,map_sum_E=map_sum_E)
                    addGenParticle(genHHVar,pObj,'H2',idx=1,map_sum_E=map_sum_E)
                elif pName == 'Zs':
                    addGenParticle(genHHVar,pObj,'Z1',idx=0,map_sum_E=map_sum_E)
                    addGenParticle(genHHVar,pObj,'Z2',idx=1,map_sum_E=map_sum_E)
                elif pName == 'VBF_quark':
                    addGenParticle(genHHVar,pObj,'VBF_quark1',idx=0,map_sum_E=map_sum_E)
                    addGenParticle(genHHVar,pObj,'VBF_quark2',idx=1,map_sum_E=map_sum_E)
                    addGenParticle(genHHVar,pObj,'VBF_quark3',idx=2,map_sum_E=map_sum_E)
                    addGenParticle(genHHVar,pObj,'VBF_quark4',idx=3,map_sum_E=map_sum_E)
                else:
                    addGenParticle(genHHVar,pObj,pName,map_sum_E=map_sum_E)
            genHHVar['flag_valid'] = op.static_cast("Int_t",get_lambda_valid(particles.values()))
            plots.append(Skim(f'gen_HH',genHHVar,noSel))



        return plots


    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        pass # No need to do plotting
