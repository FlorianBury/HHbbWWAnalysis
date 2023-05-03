import os
import sys
from operator import mul
from functools import reduce

from bamboo.analysismodules import SkimmerModule
from bamboo import treefunctions as op
from bamboo.analysisutils import makePileupWeight

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)))) # Add scripts in this directory
from BaseHHtobbWW import BaseNanoHHtobbWW
from selectionDef import *

#===============================================================================================#
#                                 SkimmerHHtobbWW                                               #
#===============================================================================================#
class SkimmerNanoHHtobbWWDL(BaseNanoHHtobbWW,SkimmerModule):
    """ Plotter module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
        super(SkimmerNanoHHtobbWWDL, self).__init__(args)

    def initialize(self):
        super(SkimmerNanoHHtobbWWDL, self).initialize(True) # avoids doing the pseudo-data for skimmer


    def defineSkimSelection(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(SkimmerNanoHHtobbWWDL,self).prepareObjects(t, noSel, sample, sampleCfg, "DL", forSkimmer=True)
            # For the Skimmer, SF must not use defineOnFirstUse -> segmentation fault

        era = sampleCfg['era'] 

        # Initialize varsToKeep dict #
        varsToKeep = dict()  

        #---------------------------------------------------------------------------------------# 
        #                                     Selections                                        #
        #---------------------------------------------------------------------------------------#
        if not self.inclusive_sel:
            #----- Check arguments -----#
            jet_level = ["Ak4","Ak8","Resolved0Btag","Resolved1Btag","Resolved2Btag","Boosted0Btag","Boosted1Btag"] # Only one must be in args

            if [boolean for (level,boolean) in self.args.__dict__.items() if level in jet_level].count(True) != 1:
                raise RuntimeError("Only one of the jet arguments must be used, check --help")
            if self.args.Channel not in ["ElEl","MuMu","ElMu"]:
                raise RuntimeError("Channel must be either 'ElEl', 'MuMu' or 'ElMu'")

            #----- Lepton selection -----#
            # Args are passed within the self #
            ElElSelObj,MuMuSelObj,ElMuSelObj = makeDoubleLeptonSelection(self,noSel,use_dd=False,fake_selection=self.args.FakeCR)
            if self.args.Channel == "ElEl":
                selObj = ElElSelObj
                dilepton = self.ElElFakeSel[0]
            if self.args.Channel == "MuMu":
                selObj = MuMuSelObj
                dilepton = self.MuMuFakeSel[0]
            if self.args.Channel == "ElMu":
                selObj = ElMuSelObj
                dilepton = self.ElMuFakeSel[0]

            #----- HME -----#
            if self.args.analysis == 'res':
                if self.args.Resolved1Btag or self.args.Resolved2Btag:
                    HME,HME_eff = self.computeResolvedHMEAfterLeptonSelections(
                                        sel   = selObj.sel,
                                        l1    = dilepton[0],
                                        l2    = dilepton[1],
                                        bjets = self.ak4JetsByBtagScore,
                                        met   = self.corrMET)
                elif self.args.Boosted1Btag:
                    HME,HME_eff = self.computeBoostedHMEAfterLeptonSelections(
                                        sel     = selObj.sel,
                                        l1      = dilepton[0],
                                        l2      = dilepton[1],
                                        fatjets = self.ak8Jets,
                                        met     = self.corrMET)
                else:
                    raise RuntimeError("Wrong category for resonant HME computations")

            #----- Apply jet corrections -----#
            ElElSelObj.sel = self.beforeJetselection(ElElSelObj.sel,'ElEl')
            MuMuSelObj.sel = self.beforeJetselection(MuMuSelObj.sel,'MuMu')
            ElMuSelObj.sel = self.beforeJetselection(ElMuSelObj.sel,'ElMu')

            #----- Jet selection -----#
            # Since the selections in one line, we can use the non copy option of the selection to modify the selection object internally
            if any([self.args.__dict__[item] for item in ["Ak4","Resolved0Btag","Resolved1Btag","Resolved2Btag"]]):
                makeAtLeastTwoAk4JetSelection(self,selObj,use_dd=False) 
            if any([self.args.__dict__[item] for item in ["Ak8","Boosted0Btag","Boosted1Btag"]]):
               makeAtLeastOneAk8JetSelection(self,selObj,use_dd=False) 
            if self.args.Resolved0Btag:
                makeExclusiveResolvedNoBtagSelection(self,selObj,use_dd=False)
            if self.args.Resolved1Btag:
                makeExclusiveResolvedOneBtagSelection(self,selObj,use_dd=False,dy_selection=self.args.DYCR)
            if self.args.Resolved2Btag:
                makeExclusiveResolvedTwoBtagsSelection(self,selObj,use_dd=False,dy_selection=self.args.DYCR)
            if self.args.Boosted0Btag:
                makeInclusiveBoostedNoBtagSelection(self,selObj,use_dd=False)
            if self.args.Boosted1Btag:
                makeInclusiveBoostedOneBtagSelection(self,selObj,use_dd=False,dy_selection=self.args.DYCR)
        else:
            noSel = self.beforeJetselection(noSel)

        def getVariation(sf, variation):
            from bamboo.treeoperations import adaptArg
            sf = adaptArg(sf)
            toChange = [sf]
            clNds = []
            sf_v = sf.clone(select=toChange.__contains__,selClones=clNds)
            assert clNds
            for nd in clNds:
                nd.changeVariation(variation)
            return sf_v.result 
                

        #---------------------------------------------------------------------------------------# 
        #                                 Synchronization tree                                  #
        #---------------------------------------------------------------------------------------#
        if self.args.Synchronization:
            if self.args.analysis == 'res':
                raise RuntimeError("This part of the Skimmer is not planned for resonant")
            # Event variables #
            varsToKeep["event"]             = None # Already in tree
            varsToKeep["run"]               = None # Already in tree 
            varsToKeep["ls"]                = t.luminosityBlock
            varsToKeep["n_presel_mu"]       = op.static_cast("UInt_t",op.rng_len(self.muonsPreSel))
            varsToKeep["n_fakeablesel_mu"]  = op.static_cast("UInt_t",op.rng_len(self.muonsFakeSel))
            varsToKeep["n_mvasel_mu"]       = op.static_cast("UInt_t",op.rng_len(self.muonsTightSel))
            varsToKeep["n_presel_ele"]      = op.static_cast("UInt_t",op.rng_len(self.electronsPreSel))
            varsToKeep["n_fakeablesel_ele"] = op.static_cast("UInt_t",op.rng_len(self.electronsFakeSel))
            varsToKeep["n_mvasel_ele"]      = op.static_cast("UInt_t",op.rng_len(self.electronsTightSel))
            varsToKeep["n_presel_ak4Jet"]   = op.static_cast("UInt_t",op.rng_len(self.ak4Jets))    
            varsToKeep["n_presel_ak8Jet"]   = op.static_cast("UInt_t",op.rng_len(self.ak8Jets))    
            varsToKeep["n_presel_ak4JetVBF"]= op.static_cast("UInt_t",op.rng_len(self.VBFJetsPreSel))
            varsToKeep["n_presel_ak4JetVBF_postLepClean"] = op.static_cast("UInt_t",op.rng_len(self.VBFJets))
            if self.args.Resolved0Btag or self.args.Resolved1Btag or self.args.Resolved2Btag:
                varsToKeep["n_presel_ak4JetVBF_postJetClean"] = op.static_cast("UInt_t",op.rng_len(self.VBFJetsResolved))
                varsToKeep["n_presel_ak4JetVBFpairs"] = op.static_cast("UInt_t",op.rng_len(self.VBFJetPairsResolved)>0)
            if self.args.Boosted0Btag or self.args.Boosted1Btag:
                varsToKeep["n_presel_ak4JetVBF_postJetClean"] = op.static_cast("UInt_t",op.rng_len(self.VBFJetsBoosted))
                varsToKeep["n_presel_ak4JetVBFpairs"] = op.static_cast("UInt_t",op.rng_len(self.VBFJetPairsBoosted)>0)
            varsToKeep["n_medium_ak4BJet"]  = op.static_cast("UInt_t",op.rng_len(self.ak4BJets))    
            varsToKeep["n_medium_ak8BJet"]  = op.static_cast("UInt_t",op.rng_len(self.ak8BJets))    
            varsToKeep["is_SR"]             = op.static_cast("UInt_t",op.OR(op.rng_len(self.ElElTightSel)>=1,
                                                                            op.rng_len(self.MuMuTightSel)>=1,
                                                                            op.rng_len(self.ElMuTightSel)>=1))
            varsToKeep["is_FR"]             = op.c_bool(self.args.FakeCR)

            if self.args.Channel == 'ElEl':
                varsToKeep["is_ee"] = op.c_bool(True)
                varsToKeep["is_mm"] = op.c_bool(False)
                varsToKeep["is_em"] = op.c_bool(False)
            if self.args.Channel == 'MuMu':
                varsToKeep["is_ee"] = op.c_bool(False)
                varsToKeep["is_mm"] = op.c_bool(True)
                varsToKeep["is_em"] = op.c_bool(False)
            if self.args.Channel == 'ElMu':
                varsToKeep["is_ee"] = op.c_bool(False)
                varsToKeep["is_mm"] = op.c_bool(False)
                varsToKeep["is_em"] = op.c_bool(True)

            varsToKeep["is_resolved"]       = op.static_cast("UInt_t",op.AND(op.rng_len(self.ak4BJets)>=1,op.rng_len(self.ak8BJets)==0))
            varsToKeep["is_boosted"]        = op.static_cast("UInt_t",op.AND(op.rng_len(self.ak8BJets)>0))

            varsToKeep['resolved_tag']      = op.static_cast("UInt_t",op.AND(op.rng_len(self.ak4BJets)>=1,op.rng_len(self.ak8BJets)==0))
            varsToKeep['boosted_tag']       = op.static_cast("UInt_t",op.AND(op.rng_len(self.ak8BJets)>0))

            # Triggers #
            varsToKeep["triggers_SingleElectron"]   = op.OR(*self.triggersPerPrimaryDataset['SingleElectron'])
            varsToKeep["triggers_SingleMuon"]       = op.OR(*self.triggersPerPrimaryDataset['SingleMuon'])
            varsToKeep["triggers_DoubleElectron"]   = op.OR(*self.triggersPerPrimaryDataset['DoubleEGamma'])
            varsToKeep["triggers_DoubleMuon"]       = op.OR(*self.triggersPerPrimaryDataset['DoubleMuon'])
            varsToKeep["triggers_MuonElectron"]     = op.OR(*self.triggersPerPrimaryDataset['MuonEG'])

            # Muons #
            for i in range(1,3): # 2 leading muons
                varsToKeep["mu{}_pt".format(i)]                    = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].pt, op.c_float(-9999., "float"))
                varsToKeep["mu{}_eta".format(i)]                   = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].eta, op.c_float(-9999.))
                varsToKeep["mu{}_phi".format(i)]                   = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].phi, op.c_float(-9999.))
                varsToKeep["mu{}_E".format(i)]                     = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].p4.E(), op.c_float(-9999., "float"))
                varsToKeep["mu{}_charge".format(i)]                = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].charge, op.c_int(-9999.))
                varsToKeep["mu{}_conept".format(i)]                = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muon_conept[self.muonsPreSel[i-1].idx], op.c_float(-9999.))
                varsToKeep["mu{}_miniRelIso".format(i)]            = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].miniPFRelIso_all, op.c_float(-9999.))
                varsToKeep["mu{}_PFRelIso04".format(i)]            = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].pfRelIso04_all, op.c_float(-9999.))
                varsToKeep["mu{}_jetNDauChargedMVASel".format(i)]  = op.c_float(-9999.)
                varsToKeep["mu{}_jetPtRel".format(i)]              = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].jetPtRelv2, op.c_float(-9999.))
                varsToKeep["mu{}_jetRelIso".format(i)]             = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].jetRelIso, op.c_float(-9999.))
                if self.inclusive_sel:
                    varsToKeep["mu{}_jetDeepJet".format(i)]        = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].jet.btagDeepFlavB, op.c_float(-9999.))
                else:
                    varsToKeep["mu{}_jetDeepJet".format(i)]        = op.c_float(-9999.)
                varsToKeep["mu{}_sip3D".format(i)]                 = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].sip3d, op.c_float(-9999.))
                varsToKeep["mu{}_dxy".format(i)]                   = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].dxy, op.c_float(-9999.))
                varsToKeep["mu{}_dxyAbs".format(i)]                = op.switch(op.rng_len(self.muonsPreSel) >= i, op.abs(self.muonsPreSel[i-1].dxy), op.c_float(-9999.))
                varsToKeep["mu{}_dz".format(i)]                    = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].dz, op.c_float(-9999.))
                varsToKeep["mu{}_segmentCompatibility".format(i)]  = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].segmentComp, op.c_float(-9999.))
                varsToKeep["mu{}_leptonMVA".format(i)]             = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].mvaTTH, op.c_float(-9999.))
                varsToKeep["mu{}_mediumID".format(i)]              = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].mediumId, op.c_float(-9999.,"Bool_t"))
                varsToKeep["mu{}_dpt_div_pt".format(i)]            = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].tunepRelPt, op.c_float(-9999.))  # Not sure
                varsToKeep["mu{}_isfakeablesel".format(i)]         = op.switch(op.rng_len(self.muonsPreSel) >= i, op.switch(self.lambda_muonFakeSel(self.muonsPreSel[i-1]), op.c_int(1), op.c_int(0)), op.c_int(-9999))
                varsToKeep["mu{}_ismvasel".format(i)]              = op.switch(op.rng_len(self.muonsPreSel) >= i, op.switch(op.AND(self.lambda_muonTightSel(self.muonsPreSel[i-1]), self.lambda_muonFakeSel(self.muonsPreSel[i-1])), op.c_int(1), op.c_int(0)), op.c_int(-9999)) # mvasel encompasses fakeablesel
                if self.is_MC:
                    varsToKeep["mu{}_isGenMatched".format(i)]          = op.switch(op.rng_len(self.muonsPreSel) >= i, op.switch(self.lambda_is_matched(self.muonsPreSel[i-1]), op.c_int(1), op.c_int(0)), op.c_int(-9999))
                    varsToKeep["mu{}_genPartFlav".format(i)]           = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].genPartFlav, op.c_int(-9999))
                    varsToKeep["mu{}_looseSF".format(i)]               = op.switch(op.rng_len(self.muonsPreSel) >= i, reduce(mul,self.lambda_MuonLooseSF(self.muonsPreSel[i-1])), op.c_int(-9999))
                    varsToKeep["mu{}_tightSF".format(i)]               = op.switch(op.rng_len(self.muonsPreSel) >= i, reduce(mul,self.lambda_MuonTightSF(self.muonsPreSel[i-1])), op.c_int(-9999))
                varsToKeep["mu{}_FR".format(i)]                    = op.switch(op.rng_len(self.muonsPreSel) >= i, self.lambda_FR_mu(self.muonsPreSel[i-1]), op.c_int(-9999))
                varsToKeep["mu{}_FRcorr".format(i)]                = op.switch(op.rng_len(self.muonsPreSel) >= i, self.lambda_FRcorr_mu(self.muonsPreSel[i-1]), op.c_int(-9999))
                varsToKeep["mu{}_FF".format(i)]                    = op.switch(op.rng_len(self.muonsPreSel) >= i, self.lambda_FF_mu(self.muonsPreSel[i-1]), op.c_int(-9999))
                #for syst in self.lambda_FR_mu(self.muonsPreSel[i-1]).op.varMap.keys():
                #    varsToKeep["mu{}_FR_{}".format(i,syst)] = op.switch(op.rng_len(self.muonsPreSel) >= i, self.lambda_FR_mu(self.muonsPreSel[i-1]).op.varMap[syst].result,op.c_int(-9999))
                #for idx,muonFR in enumerate(self.muonFRList):
                #    varsToKeep["mu{}_FR{}".format(i,idx)] = op.switch(op.rng_len(self.muonsPreSel) >= i, muonFR(self.muonsPreSel[i-1]),op.c_int(-9999))
                #    varsToKeep["mu{}_FR{}up".format(i,idx)] = op.switch(op.rng_len(self.muonsPreSel) >= i, getVariation(muonFR(self.muonsPreSel[i-1]),muonFR._systName+'up'), op.c_int(-9999))
                #    varsToKeep["mu{}_FR{}down".format(i,idx)] = op.switch(op.rng_len(self.muonsPreSel) >= i, getVariation(muonFR(self.muonsPreSel[i-1]),muonFR._systName+'down'), op.c_int(-9999))

            
            # Electrons #
            for i in range(1,3): # 2 leading electrons 
                varsToKeep["ele{}_pt".format(i)]                    = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].pt, op.c_float(-9999.))
                varsToKeep["ele{}_eta".format(i)]                   = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].eta, op.c_float(-9999.))
                varsToKeep["ele{}_phi".format(i)]                   = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].phi, op.c_float(-9999.))
                varsToKeep["ele{}_E".format(i)]                     = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].p4.E(), op.c_float(-9999.,))
                varsToKeep["ele{}_charge".format(i)]                = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].charge, op.c_int(-9999.))
                varsToKeep["ele{}_conept".format(i)]                = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electron_conept[self.electronsPreSel[i-1].idx], op.c_float(-9999.))
                varsToKeep["ele{}_miniRelIso".format(i)]            = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].miniPFRelIso_all, op.c_float(-9999.))
                varsToKeep["ele{}_PFRelIso03".format(i)]            = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].pfRelIso03_all, op.c_float(-9999.)) # Iso03, Iso04 not in NanoAOD
                varsToKeep["ele{}_jetNDauChargedMVASel".format(i)]  = op.c_float(-9999.)
                varsToKeep["ele{}_jetPtRel".format(i)]              = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].jetPtRelv2, op.c_float(-9999.))
                varsToKeep["ele{}_jetRelIso".format(i)]             = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].jetRelIso, op.c_float(-9999.))
                if self.inclusive_sel:
                    varsToKeep["mu{}_jetDeepJet".format(i)]        = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].jet.btagDeepFlavB, op.c_float(-9999.))
                else:
                    varsToKeep["mu{}_jetDeepJet".format(i)]        = op.c_float(-9999.)
                varsToKeep["ele{}_sip3D".format(i)]                 = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].sip3d, op.c_float(-9999.))
                varsToKeep["ele{}_dxy".format(i)]                   = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].dxy, op.c_float(-9999.))
                varsToKeep["ele{}_dxyAbs".format(i)]                = op.switch(op.rng_len(self.electronsPreSel) >= i, op.abs(self.electronsPreSel[i-1].dxy), op.c_float(-9999.))
                varsToKeep["ele{}_dz".format(i)]                    = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].dz, op.c_float(-9999.))
                varsToKeep["ele{}_ntMVAeleID".format(i)]            = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].mvaFall17V2noIso, op.c_float(-9999.))
                varsToKeep["ele{}_leptonMVA".format(i)]             = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].mvaTTH, op.c_float(-9999.))
                varsToKeep["ele{}_passesConversionVeto".format(i)]  = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].convVeto, op.c_float(-9999.,"Bool_t"))
                varsToKeep["ele{}_nMissingHits".format(i)]          = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].lostHits, op.c_float(-9999.,"UChar_t"))
                varsToKeep["ele{}_sigmaEtaEta".format(i)]           = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].sieie, op.c_float(-9999.))
                varsToKeep["ele{}_HoE".format(i)]                   = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].hoe, op.c_float(-9999.))
                varsToKeep["ele{}_OoEminusOoP".format(i)]           = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].eInvMinusPInv, op.c_float(-9999.))
                varsToKeep["ele{}_isfakeablesel".format(i)]         = op.switch(op.rng_len(self.electronsPreSel) >= i, op.switch(self.lambda_electronFakeSel(self.electronsPreSel[i-1]), op.c_int(1), op.c_int(0)), op.c_int(-9999))
                varsToKeep["ele{}_ismvasel".format(i)]              = op.switch(op.rng_len(self.electronsPreSel) >= i, op.switch(op.AND(self.lambda_electronTightSel(self.electronsPreSel[i-1]), self.lambda_electronFakeSel(self.electronsPreSel[i-1])), op.c_int(1), op.c_int(0)), op.c_int(-9999)) # mvasel encompasses fakeablesel
                if self.is_MC:
                    varsToKeep["ele{}_isGenMatched".format(i)]          = op.switch(op.rng_len(self.electronsPreSel) >= i, op.switch(self.lambda_is_matched(self.electronsPreSel[i-1]), op.c_int(1), op.c_int(0)), op.c_int(-9999))
                    varsToKeep["ele{}_genPartFlav".format(i)]           = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].genPartFlav, op.c_int(-9999))
                    varsToKeep["ele{}_looseSF".format(i)]               = op.switch(op.rng_len(self.electronsPreSel) >= i, reduce(mul,self.lambda_ElectronLooseSF(self.electronsPreSel[i-1])), op.c_int(-9999))
                    varsToKeep["ele{}_tightSF".format(i)]               = op.switch(op.rng_len(self.electronsPreSel) >= i, reduce(mul,self.lambda_ElectronTightSF(self.electronsPreSel[i-1])), op.c_int(-9999))
                varsToKeep["ele{}_deltaEtaSC".format(i)]            = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].deltaEtaSC, op.c_int(-9999))
                varsToKeep["ele{}_FR".format(i)]                    = op.switch(op.rng_len(self.electronsPreSel) >= i, self.lambda_FR_el(self.electronsPreSel[i-1]), op.c_int(-9999))
                varsToKeep["ele{}_FRcorr".format(i)]                = op.switch(op.rng_len(self.electronsPreSel) >= i, self.lambda_FRcorr_el(self.electronsPreSel[i-1]), op.c_int(-9999))
                varsToKeep["ele{}_FF".format(i)]                    = op.switch(op.rng_len(self.electronsPreSel) >= i, self.lambda_FF_el(self.electronsPreSel[i-1]), op.c_int(-9999))
               # for syst in self.lambda_FR_el(self.electronsPreSel[i-1]).op.varMap.keys():
               #     varsToKeep["ele{}_FR_{}".format(i,syst)] = op.switch(op.rng_len(self.electronsPreSel) >= i, self.lambda_FR_el(self.electronsPreSel[i-1]).op.varMap[syst].result,op.c_int(-9999))

            # AK4 Jets #
            for i in range(1,5): # 4 leading jets 
                varsToKeep["ak4Jet{}_pt".format(i)]                 = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].pt, op.c_float(-9999.))
                varsToKeep["ak4Jet{}_eta".format(i)]                = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].eta, op.c_float(-9999.))
                varsToKeep["ak4Jet{}_phi".format(i)]                = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].phi, op.c_float(-9999.))
                varsToKeep["ak4Jet{}_E".format(i)]                  = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].p4.E(), op.c_float(-9999.))
                varsToKeep["ak4Jet{}_CSV".format(i)]                = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].btagDeepFlavB, op.c_float(-9999.))
                if self.is_MC:
                    varsToKeep["ak4Jet{}_hadronFlavour".format(i)]      = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].hadronFlavour, op.c_float(-9999.))
                    varsToKeep["ak4Jet{}_btagSF".format(i)]             = op.switch(op.rng_len(self.ak4Jets) >= i, self.DeepJetDiscReshapingSF(self.ak4Jets[i-1]), op.c_float(-9999.))
                    varsToKeep["ak4Jet{}_puid_eff".format(i)]           = op.switch(op.rng_len(self.ak4Jets) >= i, self.jetpuid_mc_eff(self.ak4Jets[i-1]), op.c_float(-9999.))
                    varsToKeep["ak4Jet{}_puid_sfeff".format(i)]         = op.switch(op.rng_len(self.ak4Jets) >= i, self.jetpuid_sf_eff(self.ak4Jets[i-1]), op.c_float(-9999.))
                    varsToKeep["ak4Jet{}_puid_mis".format(i)]           = op.switch(op.rng_len(self.ak4Jets) >= i, self.jetpuid_mc_mis(self.ak4Jets[i-1]), op.c_float(-9999.))
                    varsToKeep["ak4Jet{}_puid_sfmis".format(i)]         = op.switch(op.rng_len(self.ak4Jets) >= i, self.jetpuid_sf_mis(self.ak4Jets[i-1]), op.c_float(-9999.))

            # VBF Jets #
            for i in range(1,6): # 5 leading jets
                varsToKeep["ak4JetVBF{}_pt".format(i)]              = op.switch(op.rng_len(self.VBFJetsPreSel) >= i, self.VBFJetsPreSel[i-1].pt, op.c_float(-9999.))
                varsToKeep["ak4JetVBF{}_eta".format(i)]             = op.switch(op.rng_len(self.VBFJetsPreSel) >= i, self.VBFJetsPreSel[i-1].eta, op.c_float(-9999.))
                varsToKeep["ak4JetVBF{}_phi".format(i)]             = op.switch(op.rng_len(self.VBFJetsPreSel) >= i, self.VBFJetsPreSel[i-1].phi, op.c_float(-9999.))
                varsToKeep["ak4JetVBF{}_E".format(i)]               = op.switch(op.rng_len(self.VBFJetsPreSel) >= i, self.VBFJetsPreSel[i-1].p4.E(), op.c_float(-9999.))
                varsToKeep["ak4JetVBF{}_CSV".format(i)]             = op.switch(op.rng_len(self.VBFJetsPreSel) >= i, self.VBFJetsPreSel[i-1].btagDeepFlavB, op.c_float(-9999.))
                if self.is_MC:
                    varsToKeep["ak4JetVBF{}_btagSF".format(i)]          = op.switch(op.rng_len(self.VBFJetsPreSel) >= i, self.DeepJetDiscReshapingSF(self.VBFJetsPreSel[i-1]), op.c_float(-9999.))

            if self.args.Resolved0Btag or self.args.Resolved1Btag or self.args.Resolved2Btag:
                VBFJetPair = self.VBFJetPairsResolved
            if self.args.Boosted0Btag or self.args.Boosted1Btag:
                VBFJetPair = self.VBFJetPairsBoosted
            
            if self.args.Resolved0Btag or self.args.Resolved1Btag or self.args.Resolved2Btag or self.args.Boosted0Btag or self.args.Boosted1Btag:
                varsToKeep["ak4JetVBFPair1_pt"]              = op.switch(op.rng_len(VBFJetPair) >= 1, VBFJetPair[0][0].pt, op.c_float(-9999.))
                varsToKeep["ak4JetVBFPair1_eta"]             = op.switch(op.rng_len(VBFJetPair) >= 1, VBFJetPair[0][0].eta, op.c_float(-9999.))
                varsToKeep["ak4JetVBFPair1_phi"]             = op.switch(op.rng_len(VBFJetPair) >= 1, VBFJetPair[0][0].phi, op.c_float(-9999.))
                varsToKeep["ak4JetVBFPair1_E"]               = op.switch(op.rng_len(VBFJetPair) >= 1, VBFJetPair[0][0].p4.E(), op.c_float(-9999.))
                varsToKeep["ak4JetVBFPair1_CSV"]             = op.switch(op.rng_len(VBFJetPair) >= 1, VBFJetPair[0][0].btagDeepFlavB, op.c_float(-9999.))
                if self.is_MC:
                    varsToKeep["ak4JetVBFPair1_btagSF"]          = op.switch(op.rng_len(VBFJetPair) >= 1, self.DeepJetDiscReshapingSF(VBFJetPair[0][0]), op.c_float(-9999.))
                varsToKeep["ak4JetVBFPair2_pt"]              = op.switch(op.rng_len(VBFJetPair) >= 1, VBFJetPair[0][1].pt, op.c_float(-9999.))
                varsToKeep["ak4JetVBFPair2_eta"]             = op.switch(op.rng_len(VBFJetPair) >= 1, VBFJetPair[0][1].eta, op.c_float(-9999.))
                varsToKeep["ak4JetVBFPair2_phi"]             = op.switch(op.rng_len(VBFJetPair) >= 1, VBFJetPair[0][1].phi, op.c_float(-9999.))
                varsToKeep["ak4JetVBFPair2_E"]               = op.switch(op.rng_len(VBFJetPair) >= 1, VBFJetPair[0][1].p4.E(), op.c_float(-9999.))
                varsToKeep["ak4JetVBFPair2_CSV"]             = op.switch(op.rng_len(VBFJetPair) >= 1, VBFJetPair[0][1].btagDeepFlavB, op.c_float(-9999.))
                if self.is_MC:
                    varsToKeep["ak4JetVBFPair2_btagSF"]          = op.switch(op.rng_len(VBFJetPair) >= 1, self.DeepJetDiscReshapingSF(VBFJetPair[0][1]), op.c_float(-9999.))


            # AK8 Jets #
            for i in range(1,3): # 2 leading fatjets 
                varsToKeep["ak8Jet{}_pt".format(i)]                 = op.switch(op.rng_len(self.ak8Jets) >= i, self.ak8Jets[i-1].pt, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_eta".format(i)]                = op.switch(op.rng_len(self.ak8Jets) >= i, self.ak8Jets[i-1].eta, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_phi".format(i)]                = op.switch(op.rng_len(self.ak8Jets) >= i, self.ak8Jets[i-1].phi, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_E".format(i)]                  = op.switch(op.rng_len(self.ak8Jets) >= i, self.ak8Jets[i-1].p4.E(), op.c_float(-9999.))
                varsToKeep["ak8Jet{}_msoftdrop".format(i)]          = op.switch(op.rng_len(self.ak8Jets) >= i, self.ak8Jets[i-1].msoftdrop, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_tau1".format(i)]               = op.switch(op.rng_len(self.ak8Jets) >= i, self.ak8Jets[i-1].tau1, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_tau2".format(i)]               = op.switch(op.rng_len(self.ak8Jets) >= i, self.ak8Jets[i-1].tau2, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_subjet0_pt".format(i)]         = op.switch(op.rng_len(self.ak8Jets) >= i, self.ak8Jets[i-1].subJet1.pt, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_subjet0_eta".format(i)]        = op.switch(op.rng_len(self.ak8Jets) >= i, self.ak8Jets[i-1].subJet1.eta, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_subjet0_phi".format(i)]        = op.switch(op.rng_len(self.ak8Jets) >= i, self.ak8Jets[i-1].subJet1.phi, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_subjet0_CSV".format(i)]        = op.switch(op.rng_len(self.ak8Jets) >= i, self.ak8Jets[i-1].subJet1.btagDeepB, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_subjet1_pt".format(i)]         = op.switch(op.rng_len(self.ak8Jets) >= i, self.ak8Jets[i-1].subJet2.pt, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_subjet1_eta".format(i)]        = op.switch(op.rng_len(self.ak8Jets) >= i, self.ak8Jets[i-1].subJet2.eta, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_subjet1_phi".format(i)]        = op.switch(op.rng_len(self.ak8Jets) >= i, self.ak8Jets[i-1].subJet2.phi, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_subjet1_CSV".format(i)]        = op.switch(op.rng_len(self.ak8Jets) >= i, self.ak8Jets[i-1].subJet2.btagDeepB, op.c_float(-9999.))

            # MET #
             
            varsToKeep["PFMET"]    = self.corrMET.pt
            varsToKeep["PFMETphi"] = self.corrMET.phi
            varsToKeep["met1_E"]   = self.corrMET.p4.E()
            varsToKeep["met1_pt"]  = self.corrMET.pt
            varsToKeep["met1_eta"] = self.corrMET.eta
            varsToKeep["met1_phi"] = self.corrMET.phi

            # VBF pair #
            if self.inclusive_sel:
                varsToKeep["vbf_m_jj"]    = op.c_float(-9999.)
                varsToKeep["vbf_dEta_jj"] = op.c_float(-9999.)
            else:
                if self.args.Resolved0Btag or self.args.Resolved1Btag or self.args.Resolved2Btag:
                    varsToKeep["vbf_m_jj"]    = op.switch(op.rng_len(self.VBFJetPairsResolved) >= 1, op.invariant_mass(self.VBFJetPairsResolved[0][0].p4,self.VBFJetPairsResolved[0][1].p4) , op.c_float(-9999.))
                    varsToKeep["vbf_pair_mass"] = op.switch(op.rng_len(self.VBFJetPairsResolved) >= 1, op.invariant_mass(self.VBFJetPairsResolved[0][0].p4,self.VBFJetPairsResolved[0][1].p4) , op.c_float(-9999.))
                    varsToKeep["vbf_dEta_jj"] = op.switch(op.rng_len(self.VBFJetPairsResolved) >= 1, op.abs(self.VBFJetPairsResolved[0][0].eta-self.VBFJetPairsResolved[0][1].eta), op.c_float(-9999.))
                    varsToKeep["vbf_pairs_absdeltaeta"] = op.switch(op.rng_len(self.VBFJetPairsResolved) >= 1, op.abs(self.VBFJetPairsResolved[0][0].eta-self.VBFJetPairsResolved[0][1].eta), op.c_float(-9999.))
                if self.args.Boosted0Btag or self.args.Boosted1Btag:
                    varsToKeep["vbf_m_jj"]    = op.switch(op.rng_len(self.VBFJetPairsBoosted) >= 1, op.invariant_mass(self.VBFJetPairsBoosted[0][0].p4,self.VBFJetPairsBoosted[0][1].p4) , op.c_float(-9999.))
                    varsToKeep["vbf_pair_mass"] = op.switch(op.rng_len(self.VBFJetPairsBoosted) >= 1, op.invariant_mass(self.VBFJetPairsBoosted[0][0].p4,self.VBFJetPairsBoosted[0][1].p4) , op.c_float(-9999.))
                    varsToKeep["vbf_dEta_jj"] = op.switch(op.rng_len(self.VBFJetPairsBoosted) >= 1, op.abs(self.VBFJetPairsBoosted[0][0].eta-self.VBFJetPairsBoosted[0][1].eta), op.c_float(-9999.))
                    varsToKeep["vbf_pairs_absdeltaeta"] = op.switch(op.rng_len(self.VBFJetPairsBoosted) >= 1, op.abs(self.VBFJetPairsBoosted[0][0].eta-self.VBFJetPairsBoosted[0][1].eta), op.c_float(-9999.))

            # SF #
            if self.is_MC:
                electronMuon_cont = op.combine((self.electronsFakeSel, self.muonsFakeSel))
                #varsToKeep["trigger_SF"] = op.multiSwitch(
                #        (op.AND(op.rng_len(self.electronsTightSel)==1,op.rng_len(self.muonsTightSel)==0) , self.ttH_singleElectron_trigSF(self.electronsTightSel[0])),
                #        (op.AND(op.rng_len(self.electronsTightSel)==0,op.rng_len(self.muonsTightSel)==1) , self.ttH_singleMuon_trigSF(self.muonsTightSel[0])),
                #        (op.AND(op.rng_len(self.electronsTightSel)>=2,op.rng_len(self.muonsTightSel)==0) , self.lambda_ttH_doubleElectron_trigSF(self.electronsTightSel)),
                #        (op.AND(op.rng_len(self.electronsTightSel)==0,op.rng_len(self.muonsTightSel)>=2) , self.lambda_ttH_doubleMuon_trigSF(self.muonsTightSel)),
                #        (op.AND(op.rng_len(self.electronsTightSel)>=1,op.rng_len(self.muonsTightSel)>=1) , self.lambda_ttH_electronMuon_trigSF(electronMuon_cont[0])),
                #         op.c_float(1.))
                if self.args.Channel == 'ElEl':
                    varsToKeep["weight_trigger"] = op.switch(op.rng_len(self.ElElFakeSel) > 0, self.lambda_ttH_doubleElectron_trigSF(self.ElElFakeSel[0]), op.c_float(0.))
                if self.args.Channel == 'MuMu':
                    varsToKeep["weight_trigger"] = op.switch(op.rng_len(self.MuMuFakeSel) > 0, self.lambda_ttH_doubleMuon_trigSF(self.MuMuFakeSel[0]), op.c_float(0.))
                if self.args.Channel == 'ElMu':
                    varsToKeep["weight_trigger"] = op.switch(op.rng_len(self.ElMuFakeSel) > 0, self.lambda_ttH_electronMuon_trigSF(self.ElMuFakeSel[0]), op.c_float(0.))

                if self.args.Channel == "ElEl":
                    varsToKeep["weight_electron_id_loose"] = op.switch(op.rng_len(self.ElElFakeSel) > 0,
                                                                       reduce(mul,self.lambda_ElectronLooseSF(self.ElElFakeSel[0][0]))* \
                                                                       reduce(mul,self.lambda_ElectronLooseSF(self.ElElFakeSel[0][1])), 
                                                                       op.c_float(1.))
                    varsToKeep["weight_muon_idiso_loose"] = op.c_float(1.)
                    varsToKeep["weight_electron_tth_tight"] = op.switch(op.rng_len(self.ElElFakeSel) > 0 ,
                                                                       op.switch(op.AND(self.lambda_electronTightSel(self.ElElFakeSel[0][0]), self.lambda_is_matched(self.ElElFakeSel[0][0])),
                                                                                 self.elTightMVA(self.ElElFakeSel[0][0]),
                                                                                 op.c_float(1.)) * \
                                                                       op.switch(op.AND(self.lambda_electronTightSel(self.ElElFakeSel[0][1]), self.lambda_is_matched(self.ElElFakeSel[0][1])),
                                                                                 self.elTightMVA(self.ElElFakeSel[0][1]),
                                                                                 op.c_float(1.)),
                                                                       op.c_float(1.))
                    varsToKeep["weight_electron_tth_relaxed"] = op.switch(op.rng_len(self.ElElFakeSel) > 0 ,
                                                                       op.switch(op.AND(self.lambda_electronTightSel(self.ElElFakeSel[0][0]), self.lambda_is_matched(self.ElElFakeSel[0][0])),
                                                                                 self.elRelaxedTightMVA(self.ElElFakeSel[0][0]),
                                                                                 op.c_float(1.)) * \
                                                                       op.switch(op.AND(self.lambda_electronTightSel(self.ElElFakeSel[0][1]), self.lambda_is_matched(self.ElElFakeSel[0][1])),
                                                                                 self.elRelaxedTightMVA(self.ElElFakeSel[0][1]),
                                                                                 op.c_float(1.)),
                                                                       op.c_float(1.))
                    varsToKeep["weight_electron_total_tight"] = op.switch(op.rng_len(self.ElElFakeSel) > 0 ,
                                                                       reduce(mul,self.lambda_ElectronTightSF(self.ElElFakeSel[0][0]))*\
                                                                       reduce(mul,self.lambda_ElectronTightSF(self.ElElFakeSel[0][1])), 
                                                                       op.c_float(1.))
                    varsToKeep["weight_muon_tth_tight"] = op.c_float(1.)
                    varsToKeep["weight_muon_tth_relaxed"] = op.c_float(1.)
                    varsToKeep["weight_muon_total_tight"] = op.c_float(1.)
                    varsToKeep["weight_lepton_total_SF"] = op.switch(op.rng_len(self.ElElFakeSel) > 0 ,
                                                               reduce(mul,self.lambda_ElectronLooseSF(self.ElElFakeSel[0][0]))* \
                                                               reduce(mul,self.lambda_ElectronLooseSF(self.ElElFakeSel[0][1]))* \
                                                               reduce(mul,self.lambda_ElectronTightSF(self.ElElFakeSel[0][0]))* \
                                                               reduce(mul,self.lambda_ElectronTightSF(self.ElElFakeSel[0][1])),
                                                               op.c_float(1.))
                if self.args.Channel == "MuMu":
                    varsToKeep["weight_electron_id_loose"] = op.c_float(1.)
                    varsToKeep["weight_muon_idiso_loose"] = op.switch(op.rng_len(self.MuMuFakeSel) > 0, 
                                                                       reduce(mul,self.lambda_MuonLooseSF(self.MuMuFakeSel[0][0]))* \
                                                                       reduce(mul,self.lambda_MuonLooseSF(self.MuMuFakeSel[0][1])), 
                                                                       op.c_float(1.))
                    varsToKeep["weight_electron_tth_tight"] = op.c_float(1.)
                    varsToKeep["weight_electron_tth_relaxed"] = op.c_float(1.)

                    varsToKeep["weight_electron_total_tight"] = op.c_float(1.)
                    varsToKeep["weight_muon_tth_tight"] = op.switch(op.rng_len(self.MuMuFakeSel) > 0 ,
                                                                       op.switch(op.AND(self.lambda_muonTightSel(self.MuMuFakeSel[0][0]), self.lambda_is_matched(self.MuMuFakeSel[0][0])),
                                                                                 self.muTightMVA(self.MuMuFakeSel[0][0]),
                                                                                 op.c_float(1.)) * \
                                                                       op.switch(op.AND(self.lambda_muonTightSel(self.MuMuFakeSel[0][1]), self.lambda_is_matched(self.MuMuFakeSel[0][1])),
                                                                                 self.muTightMVA(self.MuMuFakeSel[0][1]),
                                                                                 op.c_float(1.)),
                                                                       op.c_float(1.))
                    varsToKeep["weight_muon_tth_relaxed"] = op.switch(op.rng_len(self.MuMuFakeSel) > 0, 
                                                                       op.switch(op.AND(self.lambda_muonTightSel(self.MuMuFakeSel[0][0]), self.lambda_is_matched(self.MuMuFakeSel[0][0])),
                                                                                 self.muRelaxedTightMVA(self.MuMuFakeSel[0][0]),
                                                                                 op.c_float(1.)) * \
                                                                       op.switch(op.AND(self.lambda_muonTightSel(self.MuMuFakeSel[0][1]), self.lambda_is_matched(self.MuMuFakeSel[0][1])),
                                                                                 self.muRelaxedTightMVA(self.MuMuFakeSel[0][1]),
                                                                                 op.c_float(1.)),
                                                                       op.c_float(1.))
                    varsToKeep["weight_muon_total_tight"] = op.switch(op.rng_len(self.MuMuFakeSel) > 0 ,
                                                                       reduce(mul,self.lambda_MuonTightSF(self.MuMuFakeSel[0][0]))*\
                                                                       reduce(mul,self.lambda_MuonTightSF(self.MuMuFakeSel[0][1])), 
                                                                       op.c_float(1.))
                    varsToKeep["weight_lepton_total_SF"] = op.switch(op.rng_len(self.MuMuFakeSel) > 0,
                                                               reduce(mul,self.lambda_MuonLooseSF(self.MuMuFakeSel[0][0]))* \
                                                               reduce(mul,self.lambda_MuonLooseSF(self.MuMuFakeSel[0][1]))* \
                                                               reduce(mul,self.lambda_MuonTightSF(self.MuMuFakeSel[0][0]))* \
                                                               reduce(mul,self.lambda_MuonTightSF(self.MuMuFakeSel[0][1])),
                                                               op.c_float(1.))
                if self.args.Channel == "ElMu":
                    varsToKeep["weight_electron_id_loose"] = op.switch(op.rng_len(self.ElMuFakeSel) > 0, 
                                                                       reduce(mul,self.lambda_ElectronLooseSF(self.ElMuFakeSel[0][0])), 
                                                                       op.c_float(1.))
                    varsToKeep["weight_muon_idiso_loose"] = op.switch(op.rng_len(self.ElMuFakeSel) > 0,
                                                                       reduce(mul,self.lambda_MuonLooseSF(self.ElMuFakeSel[0][1])), 
                                                                       op.c_float(1.))
                    varsToKeep["weight_electron_tth_tight"] = op.switch(op.rng_len(self.ElMuFakeSel) > 0,
                                                                       op.switch(op.AND(self.lambda_electronTightSel(self.ElMuFakeSel[0][0]), self.lambda_is_matched(self.ElMuFakeSel[0][0])),
                                                                                 self.elTightMVA(self.ElMuFakeSel[0][0]),
                                                                                 op.c_float(1.)),
                                                                       op.c_float(1.))
                    varsToKeep["weight_electron_tth_relaxed"] = op.switch(op.rng_len(self.ElMuFakeSel) > 0, 
                                                                       op.switch(op.AND(self.lambda_electronTightSel(self.ElMuFakeSel[0][0]), self.lambda_is_matched(self.ElMuFakeSel[0][0])),
                                                                                 self.elRelaxedTightMVA(self.ElMuFakeSel[0][0]),
                                                                                 op.c_float(1.)),
                                                                       op.c_float(1.))
                    varsToKeep["weight_electron_total_tight"] = op.switch(op.rng_len(self.ElMuFakeSel) > 0, 
                                                                       reduce(mul,self.lambda_ElectronTightSF(self.ElMuFakeSel[0][0])), 
                                                                       op.c_float(1.))
                    varsToKeep["weight_muon_tth_tight"] = op.switch(op.rng_len(self.ElMuFakeSel) > 0, 
                                                                       op.switch(op.AND(self.lambda_muonTightSel(self.ElMuFakeSel[0][1]), self.lambda_is_matched(self.ElMuFakeSel[0][1])),
                                                                                 self.muTightMVA(self.ElMuFakeSel[0][1]),
                                                                                 op.c_float(1.)),
                                                                       op.c_float(1.))
                    varsToKeep["weight_muon_tth_relaxed"] = op.switch(op.rng_len(self.ElMuFakeSel) > 0, 
                                                                       op.switch(op.AND(self.lambda_muonTightSel(self.ElMuFakeSel[0][1]), self.lambda_is_matched(self.ElMuFakeSel[0][1])),
                                                                                 self.muRelaxedTightMVA(self.ElMuFakeSel[0][1]),
                                                                                 op.c_float(1.)),
                                                                       self.muRelaxedTightMVA(self.ElMuFakeSel[0][1]), 
                                                                       op.c_float(1.))
                    varsToKeep["weight_muon_total_tight"] = op.switch(op.rng_len(self.ElMuFakeSel) > 0, 
                                                                       reduce(mul,self.lambda_MuonTightSF(self.ElMuFakeSel[0][1])), 
                                                                       op.c_float(1.))
                    varsToKeep["weight_lepton_total_SF"] = op.switch(op.rng_len(self.ElMuFakeSel) > 0, 
                                                               reduce(mul,self.lambda_ElectronLooseSF(self.ElMuFakeSel[0][0]))* \
                                                               reduce(mul,self.lambda_MuonLooseSF(self.ElMuFakeSel[0][1]))* \
                                                               reduce(mul,self.lambda_ElectronTightSF(self.ElMuFakeSel[0][0]))* \
                                                               reduce(mul,self.lambda_MuonTightSF(self.ElMuFakeSel[0][1])),
                                                               op.c_float(1.))


                # L1 Prefire #
                if era in ["2016","2017"]:
                    #varsToKeep["L1prefire"] = self.L1Prefiring
                    varsToKeep["weight_l1_ecal_prefiring"] = self.L1Prefiring
                else:
                    #varsToKeep["L1prefire"] = op.c_float(-9999.)
                    varsToKeep["weight_l1_ecal_prefiring"] = op.c_float(-9999.)

            # Fake rate #
            if self.args.Channel == "ElEl":
                varsToKeep["weight_fakeRate"] = self.ElElFakeFactor(self.ElElFakeSel[0]) if self.args.FakeCR else op.c_float(1.)
                varsToKeep["weight_fake_electrons"] = op.abs(self.ElElFakeFactor(self.ElElFakeSel[0]))
                varsToKeep["weight_fake_muons"]     = op.c_float(1.)
                varsToKeep["weight_fake_two_non_tight"] = op.static_cast("Float_t",op.sign(self.ElElFakeFactor(self.ElElFakeSel[0])))
            if self.args.Channel == "MuMu":
                varsToKeep["weight_fakeRate"] = self.MuMuFakeFactor(self.MuMuFakeSel[0]) if self.args.FakeCR else op.c_float(1.)
                varsToKeep["weight_fake_electrons"] = op.c_float(1.)
                varsToKeep["weight_fake_muons"]     = op.abs(self.MuMuFakeFactor(self.MuMuFakeSel[0]))
                varsToKeep["weight_fake_two_non_tight"] = op.static_cast("Float_t",op.sign(self.MuMuFakeFactor(self.MuMuFakeSel[0])))
            if self.args.Channel == "ElMu":
                varsToKeep["weight_fakeRate"] = self.ElMuFakeFactor(self.ElMuFakeSel[0]) if self.args.FakeCR else op.c_float(1.)
                varsToKeep["weight_fake_electrons"] = op.switch(self.lambda_electronTightSel(self.ElMuFakeSel[0][0]),
                                                                op.c_float(1.),
                                                                self.lambda_FF_el(self.ElMuFakeSel[0][0]))
                varsToKeep["weight_fake_muons"]     = op.switch(self.lambda_muonTightSel(self.ElMuFakeSel[0][1]),
                                                                op.c_float(1.),
                                                                self.lambda_FF_mu(self.ElMuFakeSel[0][1]))
                varsToKeep["weight_fake_two_non_tight"] = op.static_cast("Float_t",op.sign(self.ElMuFakeFactor(self.ElMuFakeSel[0])))
            if self.is_MC:
                varsToKeep["weight_fake_is_mc"] = op.c_float(-1.)
            else:
                varsToKeep["weight_fake_is_mc"] = op.c_float(1.)

            # PU ID SF #
            if self.is_MC:
                varsToKeep["weight_jet_PUid"] = self.puid_reweighting
                varsToKeep["weight_jet_PUid_efficiency"] = self.puid_reweighting_efficiency
                varsToKeep["weight_jet_PUid_mistag"] = self.puid_reweighting_mistag

            # Btagging SF #
            if self.is_MC:
                #varsToKeep["btag_SF"] = self.btagAk4SF
                varsToKeep["weight_ak4BtagWeight"] = self.btagAk4SF
                if "BtagRatioWeight" in self.__dict__.keys():
                    #varsToKeep["btag_ratio_SF"] = self.BtagRatioWeight
                    varsToKeep["weight_ak4BtagNorm"] = self.BtagRatioWeight
                varsToKeep["weight_ak8BtagWeight"] = self.ak8BtagReweighting

            # PS weights #
            if self.is_MC:
                varsToKeep["weight_PSWeight_ISR"] = self.psISRSyst
                varsToKeep["weight_PSWeight_FSR"] = self.psFSRSyst

            # PDF weights #
            if self.is_MC:
                varsToKeep["weight_scaleWeight"] = self.scaleWeight
                #varsToKeep["weight_LHEScaleWeight_len"] = op.static_cast("UInt_t",op.rng_len(t.LHEScaleWeight)) if hasattr(t,'LHEScaleWeight') else op.c_float(-9999)
                #for i in range(0,10):
                #    varsToKeep["weight_LHEScaleWeight_{}".format(i)] = op.switch(op.rng_len(t.LHEScaleWeight)>i, t.LHEScaleWeight[i], op.c_float(-9999)) if hasattr(t,'LHEScaleWeight') else op.c_float(-9999)


            # ttbar PT reweighting #
            if self.is_MC:
                if "group" in sampleCfg and sampleCfg["group"] == 'ttbar':
                    varsToKeep["topPt_wgt"] = self.ttbar_weight(self.genTop[0],self.genAntitop[0])
            # GenVar #
            #if 'type' in sampleCfg.keys() and sampleCfg["type"] == "signal":
            if 'genh' in self.__dict__.keys():
                varsToKeep["nh"] = op.static_cast("UInt_t",op.rng_len(self.genh))
                varsToKeep["mHH"] = op.switch(op.rng_len(self.genh)==2, self.mHH, op.c_float(-9999))
                varsToKeep["cosHH"] = op.switch(op.rng_len(self.genh)==2, self.cosHH, op.c_float(-9999))
                varsToKeep["reweightLO"] = op.switch(op.rng_len(self.genh)==2, self.reweightLO(op.c_float(1.)), op.c_float(-9999))

           # Event Weight #
            if self.is_MC:
                #varsToKeep["MC_weight"] = t.genWeight
                varsToKeep["weight_gen_weight"] = t.genWeight
                #varsToKeep["PU_weight"] = self.PUWeight
                varsToKeep["weight_pileup"] = self.PUWeight
                varsToKeep["total_weight"] = noSel.weight if self.inclusive_sel else selObj.sel.weight

            # Selection weights #
            if self.is_MC and not self.inclusive_sel:
                currentSel = selObj.sel
                converter = {f'none'                                                                                                                 : 'initial',
                             f'HHMCWeight'                                                                                                           : None,
                             f'passMETFlags'                                                                                                         : 'MET_flags',
                             f'genWeight'                                                                                                            : 'gen_weight',
                             f'stitching'                                                                                                            : 'stitching',
                             f'SystOff'                                                                                                              : None,
                             f'PDFScaleWeights'                                                                                                      : 'scale_weights',
                             f'PSweights'                                                                                                            : 'PSWeight',
                             f'L1PreFiringRate'                                                                                                      : 'l1_ecal_prefiring',
                             f'puWeight'                                                                                                             : 'pileup',
                             f'Has2Fakeable{self.args.Channel}'                                                                                      : 'two_leptons',
                             f'Has2Fakeable{self.args.Channel}OS'                                                                                    : 'opposite_sign',
                             f'Has2Fakeable{self.args.Channel}OSWithTriggers'                                                                        : 'trigger_loose_lepton',
                             f'Has2Fakeable{self.args.Channel}OSWithTriggersPtCuts'                                                                  : 'pt_cuts',
                             f'Has2Fakeable{self.args.Channel}OSWithTriggersPtCutsPreMllCut'                                                         : 'mll_cut',
                             f'Has2Fakeable{self.args.Channel}OSWithTriggersPtCutsPreMllCutOutZ'                                                     : 'Z_cut',
                             f'Has2Fakeable{self.args.Channel}OSWithTriggersPtCutsPreMllCutOutZTightSelected'                                        : 'tight_lepton',
                             f'jetPUIDReweighting{self.args.Channel}'                                                                                : 'jet_PUid',
                             f'BtagAk4SF{self.args.Channel}'                                                                                         : 'btag_ak4',
                             f'BtagAk8SF{self.args.Channel}'                                                                                         : 'btag_ak8',
                             f'Has2Fakeable{self.args.Channel}OSWithTriggersPtCutsPreMllCutOutZTightSelectedTwoAk4Jets'                              : 'two_ak4',               
                             f'Has2Fakeable{self.args.Channel}OSWithTriggersPtCutsPreMllCutOutZTightSelectedTwoAk4JetsExclusiveResolvedOneBtag'      : 'resolved_1b',               
                             f'Has2Fakeable{self.args.Channel}OSWithTriggersPtCutsPreMllCutOutZTightSelectedTwoAk4JetsExclusiveResolvedTwoBtags'     : 'resolved_2b',               
                             f'Has2Fakeable{self.args.Channel}OSWithTriggersPtCutsPreMllCutOutZTightSelectedOneAk8Jet'                               : 'one_ak8',               
                             f'Has2Fakeable{self.args.Channel}OSWithTriggersPtCutsPreMllCutOutZTightSelectedOneAk8JetInclusiveBoostedOneBtag'        : 'boosted'}

                while hasattr(currentSel,'parent'):
                    parentSel = currentSel.parent
                    if not currentSel.name in converter.keys():
                        raise RuntimeError(f'Unknown {currentSel.name}, maybe add it to the list ?')
                    if converter[currentSel.name] is not None:
                        varsToKeep[f"weight_sel_acc_{converter[currentSel.name]}"] = op.static_cast('float',currentSel.weight)
                        if parentSel is None:
                            varsToKeep[f"weight_sel_{converter[currentSel.name]}"] = op.c_float(1.)
                        else:
                            varsToKeep[f"weight_sel_{converter[currentSel.name]}"] = op.static_cast('float',op.static_cast('float',currentSel.weight)/op.static_cast('float',parentSel.weight))
                    currentSel = parentSel

                for brName in converter.values():
                    if brName is None:
                        continue
                    if f"weight_sel_{brName}" not in varsToKeep.keys():
                        varsToKeep[f"weight_sel_{brName}"] = op.static_cast('float',op.c_float(0.))
                        varsToKeep[f"weight_sel_acc_{brName}"] = op.static_cast('float',op.c_float(0.))

            if not self.inclusive_sel:
                import mvaEvaluatorDL_nonres
                inputsLeps = mvaEvaluatorDL_nonres.returnLeptonsMVAInputs(
                                                self      = self,
                                                l1        = dilepton[0],
                                                l2        = dilepton[1],
                                                channel   = self.args.Channel)
                inputsJets = mvaEvaluatorDL_nonres.returnJetsMVAInputs(
                                                self      = self,
                                                jets      = self.ak4Jets)
                inputsMET = mvaEvaluatorDL_nonres.returnMETMVAInputs(
                                                self      = self,
                                                met       = self.corrMET)     
                inputsFatjet = mvaEvaluatorDL_nonres.returnFatjetMVAInputs(
                                                self      = self)
                inputsHL = mvaEvaluatorDL_nonres.returnHighLevelMVAInputs(
                                                self      = self,
                                                l1        = dilepton[0],
                                                l2        = dilepton[1],
                                                met       = self.corrMET,
                                                jets      = self.ak4Jets,
                                                bjets     = self.ak4JetsByBtagScore[:op.min(op.rng_len(self.ak4JetsByBtagScore),op.static_cast("std::size_t",op.c_int(2)))],
                                                electrons = self.electronsFakeSel,
                                                muons     = self.muonsFakeSel,
                                                channel   = self.args.Channel)

                inputsParam = mvaEvaluatorDL_nonres.returnParamMVAInputs(self)
                inputsEventNr = mvaEvaluatorDL_nonres.returnEventNrMVAInputs(self,t)

                inputsDict = {**inputsLeps,**inputsJets,**inputsMET,**inputsFatjet,**inputsHL,**inputsParam,**inputsEventNr}

                for (varname,_,_),var in inputsDict.items():
                    varsToKeep[varname] = var
           
                path_model = os.path.join(os.path.abspath(os.path.dirname(__file__)),'MachineLearning','ml-models','models','multi-classification','dnn','DL','12','model','model.pb')
                nodes = ['GGF','VBF','H', 'DY', 'ST', 'TT', 'TTVX', 'VVV', 'Rare']
                input_names = ["lep","jet","fat","met","hl","param","eventnr"]
                output_name = "Identity"
                if not os.path.exists(path_model):
                    raise RuntimeError('Could not find model file %s'%path_model)
                try:
                    DNN = op.mvaEvaluator(path_model,mvaType='Tensorflow',otherArgs=(input_names, output_name))
                except:
                    raise RuntimeError('Could not load model %s'%path_model)

                inputs = [op.array("double",*mvaEvaluatorDL_nonres.inputStaticCast(inputsLeps,"float")),
                          op.array("double",*mvaEvaluatorDL_nonres.inputStaticCast(inputsJets,"float")),
                          op.array("double",*mvaEvaluatorDL_nonres.inputStaticCast(inputsFatjet,"float")),
                          op.array("double",*mvaEvaluatorDL_nonres.inputStaticCast(inputsMET,"float")),
                          op.array("double",*mvaEvaluatorDL_nonres.inputStaticCast(inputsHL,"float")),
                          op.array("double",*mvaEvaluatorDL_nonres.inputStaticCast(inputsParam,"float")),
                          op.array("long",*mvaEvaluatorDL_nonres.inputStaticCast(inputsEventNr,"long"))]

                outputs = DNN(*inputs)
                for node, output in zip(nodes,outputs): 
                   varsToKeep[node] = output

            # Return #
           
            if self.inclusive_sel:
                return noSel, varsToKeep
            else:
                return selObj.sel, varsToKeep
                

        #---------------------------------------------------------------------------------------# 
        #                                    Selection tree                                     #
        #---------------------------------------------------------------------------------------#
        if not self.args.analysis == 'res':
            raise RuntimeError("This part of the Skimmer is only for resonant")
        import mvaEvaluatorDL_res

        if self.args.Channel is None:
            raise RuntimeError("You need to specify --Channel")
        if self.args.Channel == "ElEl": dilepton = self.ElElTightSel[0]
        if self.args.Channel == "MuMu": dilepton = self.MuMuTightSel[0]
        if self.args.Channel == "ElMu": dilepton = self.ElMuTightSel[0]

        varsToKeep["is_SR"] = op.static_cast("UInt_t",op.OR(op.rng_len(self.ElElTightSel)>0,
                                                            op.rng_len(self.MuMuTightSel)>0,
                                                            op.rng_len(self.ElMuTightSel)>0))
        varsToKeep['is_ee'] = op.static_cast("UInt_t",op.rng_len(self.ElElTightSel)>0)
        varsToKeep['is_mm'] = op.static_cast("UInt_t",op.rng_len(self.MuMuTightSel)>0)
        varsToKeep['is_em'] = op.static_cast("UInt_t",op.rng_len(self.ElMuTightSel)>0)

        l1 = dilepton[0]
        l2 = dilepton[1]

        # DNN inputs #
        inputsAll = mvaEvaluatorDL_res.returnResonantMVAInputs(
                                self      = self,
                                l1        = dilepton[0],
                                l2        = dilepton[1],
                                channel   = self.args.Channel,
                                jets      = self.ak4Jets,
                                fatjets   = self.ak8Jets,                                                                                                                                   
                                bjets     = self.ak4JetsByBtagScore[:op.min(op.rng_len(self.ak4JetsByBtagScore),op.static_cast("std::size_t",op.c_int(2)))],
                                met       = self.corrMET,
                                electrons = self.electronsTightSel,
                                muons     = self.muonsTightSel)

        for (varName,_,_), var in inputsAll.items():
            if not isinstance(var,list):
                varsToKeep[varName] = var

        # DNN outputs #
        dnn_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)),'MachineLearning','ResonantModels')
        path_model_HighMass = os.path.join(dnn_dir,'Resonant_HighMass_Final_512x4_w0p1.pb')
        path_model_LowMass = os.path.join(dnn_dir,'Resonant_LowMass_Final_512x4_w1.pb')
        input_names_HighMass = []
        input_names_LowMass = []
        with open(os.path.join(dnn_dir,'Resonant_HighMass_Final_512x4_w0p1_inputs.txt'),'r') as handle:
            for line in handle:
                input_names_HighMass.append(line.split()[0])
        with open(os.path.join(dnn_dir,'Resonant_LowMass_Final_512x4_w1_inputs.txt'),'r') as handle:
            for line in handle:
                input_names_LowMass.append(line.split()[0])
        output_name = "Identity"

        print ("DNN model : %s"%path_model_HighMass)
        print ("DNN model : %s"%path_model_LowMass)
        if not os.path.exists(path_model_HighMass):
            raise RuntimeError('Could not find model file %s'%path_model_HighMass)
        if not os.path.exists(path_model_LowMass):
            raise RuntimeError('Could not find model file %s'%path_model_LowMass)
        try:
            DNN_HighMass = op.mvaEvaluator(path_model_HighMass,mvaType='Tensorflow',otherArgs=(input_names_HighMass, output_name))
        except:
            raise RuntimeError('Could not load model %s'%path_model_HighMass)
        try:
            DNN_LowMass = op.mvaEvaluator(path_model_LowMass,mvaType='Tensorflow',otherArgs=(input_names_LowMass, output_name))
        except:
            raise RuntimeError('Could not load model %s'%path_model_LowMass)



        if self.args.mass is None:
            raise RuntimeError('--mass needs to be used')
        self.args.mass = self.args.mass[0]
        if self.args.mass <= 500:
            DNN = DNN_LowMass
            input_names = input_names_LowMass
        else:
            DNN = DNN_HighMass
            input_names = input_names_LowMass

        inputs = {inpName:val  for (inpName,_,_),val in inputsAll.items()}
        inputs['param'] = op.c_float(float(self.args.mass)) 
        inputsArr = []

        for inpName in input_names:
            if inpName not in inputs.keys():
                for key in inputs.keys():
                    if inpName == key.replace('$','').replace(' ','').replace('_',''):
                        inpName = key
            
            if inpName not in inputs.keys():
                raise RuntimeError(f"Input node {inpName} not found in the inputs in bamboo")
            inpVal = inputs[inpName]

            if inpName == "eventnr":
                inpType = "long"
            else:
                inpType = "float"
            if isinstance(inpVal,list):
                inpVal = [op.static_cast(inpType,inp) for inp in inpVal] 
                inputsArr.append(op.array(inpType,*inpVal))
            else:
                inpVal = op.static_cast(inpType,inpVal)
                inputsArr.append(op.array(inpType,inpVal))

        nodes = ['DY','GGF','H','Rare','ST','TT','TTVX','VVV']
        outputs = DNN(*inputsArr)

        for node, output in zip(nodes,outputs): 
           varsToKeep[node] = output
        #----- Additional -----#
        #if self.args.Resolved1Btag or self.args.Resolved2Btag:
        #    varsToKeep["b1_E"]  = self.ak4JetsByBtagScore[0].p4.E()
        #    varsToKeep["b1_Px"] = self.ak4JetsByBtagScore[0].p4.Px()
        #    varsToKeep["b1_Py"] = self.ak4JetsByBtagScore[0].p4.Py()
        #    varsToKeep["b1_Pz"] = self.ak4JetsByBtagScore[0].p4.Pz()
        #    varsToKeep["b2_E"]  = self.ak4JetsByBtagScore[1].p4.E()
        #    varsToKeep["b2_Px"] = self.ak4JetsByBtagScore[1].p4.Px()
        #    varsToKeep["b2_Py"] = self.ak4JetsByBtagScore[1].p4.Py()
        #    varsToKeep["b2_Pz"] = self.ak4JetsByBtagScore[1].p4.Pz()
        #if self.args.Boosted1Btag:
        #    varsToKeep["fatbjet_E"]  = self.ak8BJets[0].p4.E()
        #    varsToKeep["fatbjet_Px"] = self.ak8BJets[0].p4.Px()
        #    varsToKeep["fatbjet_Py"] = self.ak8BJets[0].p4.Py()
        #    varsToKeep["fatbjet_Pz"] = self.ak8BJets[0].p4.Pz()
        #    varsToKeep["fatbjet_subjet1_E"]  = self.ak8BJets[0].subJet1.p4.E()
        #    varsToKeep["fatbjet_subjet1_Px"] = self.ak8BJets[0].subJet1.p4.Px()
        #    varsToKeep["fatbjet_subjet1_Py"] = self.ak8BJets[0].subJet1.p4.Py()
        #    varsToKeep["fatbjet_subjet1_Pz"] = self.ak8BJets[0].subJet1.p4.Pz()
        #    varsToKeep["fatbjet_subjet2_E"]  = self.ak8BJets[0].subJet2.p4.E()
        #    varsToKeep["fatbjet_subjet2_Px"] = self.ak8BJets[0].subJet2.p4.Px()
        #    varsToKeep["fatbjet_subjet2_Py"] = self.ak8BJets[0].subJet2.p4.Py()
        #    varsToKeep["fatbjet_subjet2_Pz"] = self.ak8BJets[0].subJet2.p4.Pz()

        #    gen = op.select(t.GenPart,lambda g : op.AND(op.OR(op.abs(g.pdgId) == 1,
        #                                                      op.abs(g.pdgId) == 2,
        #                                                      op.abs(g.pdgId) == 3,
        #                                                      op.abs(g.pdgId) == 4,
        #                                                      op.abs(g.pdgId) == 5),
        #                                                g.statusFlags & ( 0x1 << 13),
        #                                                g.pt>=20,
        #                                                op.abs(g.eta)<2.4))
        #    
        #    ak8_sub1 = self.ak8BJets[0].subJet1
        #    ak8_sub2 = self.ak8BJets[0].subJet2
        #    gen_sub1 = op.sort(gen, lambda g : -op.deltaR(g.p4,ak8_sub1.p4))[0]
        #    gen_sub2 = op.sort(gen, lambda g : -op.deltaR(g.p4,ak8_sub2.p4))[0]
        #        
        #    varsToKeep["fatbjet_subjet1_Pt"]  = self.ak8BJets[0].subJet1.pt
        #    varsToKeep["fatbjet_subjet2_Pt"]  = self.ak8BJets[0].subJet2.pt
        #    varsToKeep["gen_subjet1_Pt"]  = gen_sub1.pt
        #    varsToKeep["gen_subjet2_Pt"]  = gen_sub2.pt
        #    varsToKeep["gen_subjet1_pdgId"]  = gen_sub1.pdgId
        #    varsToKeep["gen_subjet2_pdgId"]  = gen_sub2.pdgId



        #----- HME ----#
        varsToKeep["HME"] = HME
        varsToKeep["HME_eff"] = HME_eff

        #----- Additional variables -----#
        if self.is_MC:
            varsToKeep["MC_weight"] = t.genWeight
        varsToKeep['total_weight']  = selObj.sel.weight
        varsToKeep['era']           = op.c_int(int(self.era))
        varsToKeep["event"]         = None # Already in tree
        varsToKeep["run"]           = None # Already in tree 
        varsToKeep["ls"]            = t.luminosityBlock

        return selObj.sel, varsToKeep

    ### PostProcess ###
    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        super(SkimmerNanoHHtobbWWDL, self).postProcess(taskList, config, workdir, resultsdir, forSkimmer=True)

