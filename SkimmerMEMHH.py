import os
import sys
import copy
from operator import mul
from functools import reduce

from bamboo.analysismodules import HistogramsModule, DataDrivenBackgroundHistogramsModule
from bamboo import treefunctions as op
from bamboo.analysisutils import makePileupWeight
from bamboo.plots import Plot, Skim, EquidistantBinning

from IPython import embed

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)))) # Add scripts in this directory
from BaseHHtobbWW import BaseNanoHHtobbWW
from selectionDef import *
from highlevelLambdas import *

class SkimmerMEMHHDL(BaseNanoHHtobbWW,HistogramsModule):
    """ Plotter module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
        super(SkimmerMEMHHDL, self).__init__(args)

    def initialize(self):
        super(SkimmerMEMHHDL, self).initialize(True) # avoids doing the pseudo-data for skimmer


    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        noSel = super(SkimmerMEMHHDL,self).prepareObjects(t, noSel, sample, sampleCfg, "DL", forSkimmer=True)
            # For the Skimmer, SF must not use defineOnFirstUse -> segmentation fault

        # Initialize #
        era = sampleCfg['era']
        plots = []

        if self.inclusive_sel:
            raise RuntimeError("Inclusive analysis not possible")

        #---------------------------------------------------------------------------------------#
        #                                     Selections                                        #
        #---------------------------------------------------------------------------------------#
        #----- Check arguments -----#
        jet_levels = ["Resolved0Btag","Resolved1Btag","Resolved2Btag","Boosted0Btag","Boosted1Btag"]
        jet_selections = []
        for jet_level in jet_levels:
            assert hasattr(self.args,jet_level)
            if getattr(self.args,jet_level):
                jet_selections.append(jet_level)
        # If none specified, run on all #
        if len(jet_selections) == 0:
            jet_selections = jet_levels



        #----- Lepton selection -----#
        # Args are passed within the self #
        channels = ['ElEl','MuMu','ElMu']
        if self.args.Channel is not None:
            assert self.args.Channel in channels
            channels = [self.args.Channel]
        dileptons = [self.ElElFakeSel[0],self.MuMuFakeSel[0],self.ElMuFakeSel[0]]
        lepSelObjs = makeDoubleLeptonSelection(self,noSel,use_dd=False,fake_selection=self.args.FakeCR)

        # Lepton loop #
        for channel,lepSelObj,dilepton in zip(channels,lepSelObjs,dileptons):
            #----- Apply jet corrections -----#
            lepSelObj.sel = self.beforeJetselection(lepSelObj.sel,channel)

            #----- Apply jet basic selections -----#
            resSelObj = makeAtLeastTwoAk4JetSelection(self,lepSelObj,copy_sel=True)
            booSelObj = makeAtLeastOneAk8JetSelection(self,lepSelObj,copy_sel=True)

            # Jet loop #
            for jet_selection in jet_selections:

                #----- Make jet selection ----#
                if "Resolved" in jet_selection:
                    if jet_selection == 'Resolved0Btag':
                        selObj = makeExclusiveResolvedNoBtagSelection(self,resSelObj,use_dd=False,copy_sel=True)
                    if jet_selection == 'Resolved1Btag':
                        selObj = makeExclusiveResolvedOneBtagSelection(self,resSelObj,use_dd=False,copy_sel=True)
                    if jet_selection == 'Resolved2Btag':
                        selObj = makeExclusiveResolvedTwoBtagsSelection(self,resSelObj,use_dd=False,copy_sel=True)
                elif "Boosted" in jet_selection:
                    if jet_selection == "Boosted0Btag":
                        selObj = makeInclusiveBoostedNoBtagSelection(self,booSelObj,use_dd=False,copy_sel=True)
                    if jet_selection == "Boosted1Btag":
                        selObj = makeInclusiveBoostedOneBtagSelection(self,booSelObj,use_dd=False,copy_sel=True)
                else:
                    raise RuntimeError



                #----- Initialise variable dict -----#
                varsToKeep = dict()

                #---------------------------------------------------------------------------------------#
                #                                    Selection tree                                     #
                #---------------------------------------------------------------------------------------#

                #----- MET variables -----#
                MET = self.corrMET

                varsToKeep['met_pt']  = MET.pt
                varsToKeep['met_phi'] = MET.phi
                varsToKeep['met_E']   = MET.p4.E()
                varsToKeep['met_Px']  = MET.p4.Px()
                varsToKeep['met_Py']  = MET.p4.Py()
                varsToKeep['met_Pz']  = MET.p4.Pz()

                #----- Lepton variables -----#
                varsToKeep["is_SR"] = op.static_cast("UInt_t",op.OR(op.rng_len(self.ElElTightSel)>0,
                                                                    op.rng_len(self.MuMuTightSel)>0,
                                                                    op.rng_len(self.ElMuTightSel)>0))
                varsToKeep['is_ee'] = op.static_cast("UInt_t",op.rng_len(self.ElElTightSel)>0)
                varsToKeep['is_mm'] = op.static_cast("UInt_t",op.rng_len(self.MuMuTightSel)>0)
                varsToKeep['is_em'] = op.static_cast("UInt_t",op.rng_len(self.ElMuTightSel)>0)
                varsToKeep['resolved1b_tag'] = op.static_cast("UInt_t",op.AND(op.rng_len(self.ak4BJets)==1,op.rng_len(self.ak8BJets)==0))
                varsToKeep['resolved2b_tag'] = op.static_cast("UInt_t",op.AND(op.rng_len(self.ak4BJets)>=2,op.rng_len(self.ak8BJets)==0))
                varsToKeep['boosted1b_tag'] = op.static_cast("UInt_t",op.AND(op.rng_len(self.ak8BJets)>0))

                l1 = dilepton[0]
                l2 = dilepton[1]

                varsToKeep['l1_Px']     = l1.p4.Px()
                varsToKeep['l1_Py']     = l1.p4.Py()
                varsToKeep['l1_Pz']     = l1.p4.Pz()
                varsToKeep['l1_E']      = l1.p4.E()
                varsToKeep['l1_pt']     = l1.pt
                varsToKeep['l1_eta']    = l1.eta
                varsToKeep['l1_phi']    = l1.phi
                varsToKeep['l1_pdgId']  = l1.pdgId
                varsToKeep['l1_charge'] = l1.charge

                varsToKeep['l2_Px']     = l2.p4.Px()
                varsToKeep['l2_Py']     = l2.p4.Py()
                varsToKeep['l2_Pz']     = l2.p4.Pz()
                varsToKeep['l2_E']      = l2.p4.E()
                varsToKeep['l2_pt']     = l2.pt
                varsToKeep['l2_eta']    = l2.eta
                varsToKeep['l2_phi']    = l2.phi
                varsToKeep['l2_pdgId']  = l2.pdgId
                varsToKeep['l2_charge'] = l2.charge

                varsToKeep['m_ll'] = op.invariant_mass(l1.p4,l2.p4)

                ##----- Jet variables -----#
                jets = self.ak4JetsByBtagScore
                for idx in range(1,5):
                    varsToKeep[f'j{idx}_Px']  = op.switch(op.rng_len(jets)>=idx, jets[idx-1].p4.Px(), op.c_float(-9999))
                    varsToKeep[f'j{idx}_Py']  = op.switch(op.rng_len(jets)>=idx, jets[idx-1].p4.Py(), op.c_float(-9999))
                    varsToKeep[f'j{idx}_Pz']  = op.switch(op.rng_len(jets)>=idx, jets[idx-1].p4.Pz(), op.c_float(-9999))
                    varsToKeep[f'j{idx}_E']   = op.switch(op.rng_len(jets)>=idx, jets[idx-1].p4.E(), op.c_float(-9999))
                    varsToKeep[f'j{idx}_pt']  = op.switch(op.rng_len(jets)>=idx, jets[idx-1].pt, op.c_float(-9999))
                    varsToKeep[f'j{idx}_eta'] = op.switch(op.rng_len(jets)>=idx, jets[idx-1].eta, op.c_float(-9999))
                    varsToKeep[f'j{idx}_phi'] = op.switch(op.rng_len(jets)>=idx, jets[idx-1].phi, op.c_float(-9999))
                    varsToKeep[f'j{idx}_btag']= op.switch(op.rng_len(jets)>=idx, jets[idx-1].btagDeepFlavB, op.c_float(-9999))
                    varsToKeep[f'j{idx}_btagged']= op.switch(op.rng_len(jets)>=idx, op.switch(self.lambda_ak4Btag(jets[idx-1]),1,0), op.c_float(-9999))

                varsToKeep['n_ak4'] = op.static_cast("UInt_t",op.rng_len(self.ak4Jets))
                varsToKeep['n_ak4_btag'] = op.static_cast("UInt_t",op.rng_len(self.ak4BJets))
                if era == "2016":
                    varsToKeep['btag_threshold'] = op.c_float(0.3093)
                if era == "2017":
                    varsToKeep['btag_threshold'] = op.c_float(0.3033)
                if era == "2018":
                    varsToKeep['btag_threshold'] = op.c_float(0.2770)

                #----- Fatjet variables -----#
                fatjets = self.ak8BJets
                subJet1 = fatjets[0].subJet1
                subJet2 = fatjets[0].subJet2

                varsToKeep['fatj_sub1_Px']   = op.switch(op.rng_len(fatjets)>0, subJet1.p4.Px(), op.c_float(-9999))
                varsToKeep['fatj_sub1_Py']   = op.switch(op.rng_len(fatjets)>0, subJet1.p4.Py(), op.c_float(-9999))
                varsToKeep['fatj_sub1_Pz']   = op.switch(op.rng_len(fatjets)>0, subJet1.p4.Pz(), op.c_float(-9999))
                varsToKeep['fatj_sub1_E']    = op.switch(op.rng_len(fatjets)>0, subJet1.p4.E(), op.c_float(-9999))
                varsToKeep['fatj_sub1_pt']   = op.switch(op.rng_len(fatjets)>0, subJet1.pt, op.c_float(-9999))
                varsToKeep['fatj_sub1_eta']  = op.switch(op.rng_len(fatjets)>0, subJet1.eta, op.c_float(-9999))
                varsToKeep['fatj_sub1_phi']  = op.switch(op.rng_len(fatjets)>0, subJet1.phi, op.c_float(-9999))
                varsToKeep['fatj_sub1_btag'] = op.switch(op.rng_len(fatjets)>0, subJet1.btagDeepB, op.c_float(-9999))

                varsToKeep['fatj_sub2_Px']   = op.switch(op.rng_len(fatjets)>0, subJet2.p4.Px(), op.c_float(-9999))
                varsToKeep['fatj_sub2_Py']   = op.switch(op.rng_len(fatjets)>0, subJet2.p4.Py(), op.c_float(-9999))
                varsToKeep['fatj_sub2_Pz']   = op.switch(op.rng_len(fatjets)>0, subJet2.p4.Pz(), op.c_float(-9999))
                varsToKeep['fatj_sub2_E']    = op.switch(op.rng_len(fatjets)>0, subJet2.p4.E(), op.c_float(-9999))
                varsToKeep['fatj_sub2_pt']   = op.switch(op.rng_len(fatjets)>0, subJet2.pt, op.c_float(-9999))
                varsToKeep['fatj_sub2_eta']  = op.switch(op.rng_len(fatjets)>0, subJet2.eta, op.c_float(-9999))
                varsToKeep['fatj_sub2_phi']  = op.switch(op.rng_len(fatjets)>0, subJet2.phi, op.c_float(-9999))
                varsToKeep['fatj_sub2_btag'] = op.switch(op.rng_len(fatjets)>0, subJet2.btagDeepB, op.c_float(-9999))

                varsToKeep['fatj_Px']             = op.switch(op.rng_len(fatjets)>0, fatjets[0].p4.Px(), op.c_float(-9999))
                varsToKeep['fatj_Py']             = op.switch(op.rng_len(fatjets)>0, fatjets[0].p4.Py(), op.c_float(-9999))
                varsToKeep['fatj_Pz']             = op.switch(op.rng_len(fatjets)>0, fatjets[0].p4.Pz(), op.c_float(-9999))
                varsToKeep['fatj_E']              = op.switch(op.rng_len(fatjets)>0, fatjets[0].p4.E(), op.c_float(-9999))
                varsToKeep['fatj_pt']             = op.switch(op.rng_len(fatjets)>0, fatjets[0].pt, op.c_float(-9999))
                varsToKeep['fatj_eta']            = op.switch(op.rng_len(fatjets)>0, fatjets[0].eta, op.c_float(-9999))
                varsToKeep['fatj_phi']            = op.switch(op.rng_len(fatjets)>0, fatjets[0].phi, op.c_float(-9999))
                varsToKeep['fatj_softdropMass']   = op.switch(op.rng_len(fatjets)>0, fatjets[0].msoftdrop, op.c_float(-9999))
                varsToKeep['fatj_btagDeepB']      = op.switch(op.rng_len(fatjets)>0, fatjets[0].btagDeepB, op.c_float(-9999))
                varsToKeep['fatj_btagHbb']        = op.switch(op.rng_len(fatjets)>0, fatjets[0].btagHbb, op.c_float(-9999))

                varsToKeep['n_ak8'] = op.static_cast("UInt_t",op.rng_len(self.ak8Jets))
                varsToKeep['n_ak8_btag'] = op.static_cast("UInt_t",op.rng_len(self.ak8BJets))

                #----- VBF variables -----#
                if "Resolved" in jet_selection:
                    VBFpairs = self.VBFJetPairsResolved
                if "Boosted" in jet_selection:
                    VBFpairs = self.VBFJetPairsBoosted

                varsToKeep['has_VBF']     = op.static_cast("UInt_t",op.rng_len(VBFpairs)>0)

                varsToKeep['VBF_j1_Px']   = op.switch(op.rng_len(VBFpairs)>0, VBFpairs[0][0].p4.Px(), op.c_float(-9999))
                varsToKeep['VBF_j1_Py']   = op.switch(op.rng_len(VBFpairs)>0, VBFpairs[0][0].p4.Py(), op.c_float(-9999))
                varsToKeep['VBF_j1_Pz']   = op.switch(op.rng_len(VBFpairs)>0, VBFpairs[0][0].p4.Pz(), op.c_float(-9999))
                varsToKeep['VBF_j1_E']    = op.switch(op.rng_len(VBFpairs)>0, VBFpairs[0][0].p4.E(), op.c_float(-9999))
                varsToKeep['VBF_j1_pt']   = op.switch(op.rng_len(VBFpairs)>0, VBFpairs[0][0].pt, op.c_float(-9999))
                varsToKeep['VBF_j1_eta']  = op.switch(op.rng_len(VBFpairs)>0, VBFpairs[0][0].eta, op.c_float(-9999))
                varsToKeep['VBF_j1_phi']  = op.switch(op.rng_len(VBFpairs)>0, VBFpairs[0][0].phi, op.c_float(-9999))

                varsToKeep['VBF_j2_Px']   = op.switch(op.rng_len(VBFpairs)>0, VBFpairs[0][1].p4.Px(), op.c_float(-9999))
                varsToKeep['VBF_j2_Py']   = op.switch(op.rng_len(VBFpairs)>0, VBFpairs[0][1].p4.Py(), op.c_float(-9999))
                varsToKeep['VBF_j2_Pz']   = op.switch(op.rng_len(VBFpairs)>0, VBFpairs[0][1].p4.Pz(), op.c_float(-9999))
                varsToKeep['VBF_j2_E']    = op.switch(op.rng_len(VBFpairs)>0, VBFpairs[0][1].p4.E(), op.c_float(-9999))
                varsToKeep['VBF_j2_pt']   = op.switch(op.rng_len(VBFpairs)>0, VBFpairs[0][1].pt, op.c_float(-9999))
                varsToKeep['VBF_j2_eta']  = op.switch(op.rng_len(VBFpairs)>0, VBFpairs[0][1].eta, op.c_float(-9999))
                varsToKeep['VBF_j2_phi']  = op.switch(op.rng_len(VBFpairs)>0, VBFpairs[0][1].phi, op.c_float(-9999))


                #----- Additional variables -----#
                if self.is_MC:
                    varsToKeep["MC_weight"]     = t.genWeight
                varsToKeep['total_weight']      = selObj.sel.weight
                varsToKeep["event"]             = None # Already in tree
                varsToKeep["run"]               = None # Already in tree
                varsToKeep["ls"]                = t.luminosityBlock


                #----- Add to plots -----#
                if self.args.NoSystematics or self.args.Systematic == 'nominal':
                    plots.append(Skim(f"{channel}_{jet_selection}_nominal", varsToKeep, selObj.sel))
                else:
                    # Get all possible variations #
                    systematics = []
                    for name,var in varsToKeep.items():
                        if var is None:
                            continue
                        if name == 'total_weight':
                            continue # otherwise we catch also pure-weight systematics
                        for syst in op.getSystematicVariations(var):
                            if syst not in systematics:
                                systematics.append(syst)
                    # Check if valid systematic #
                    if self.args.Systematic not in systematics:
                        raise RuntimeError(f'Systematic {self.args.Systematic} not in list {systematics}')
                    # Produce up and down variations #
                    for name, var in varsToKeep.items():
                        if var is not None:
                            varsToKeep[name] = op.forSystematicVariation(var,self.args.Systematic)
                    plots.append(Skim(f"{channel}_{jet_selection}_{self.args.Systematic}", varsToKeep, selObj.sel))

                #----- Systematic variatons -----#
                #if not self.args.NoSystematics:
                #    #systDicts = []
                #    # Produce new skim for each systematic #
                #    for systematic in systematics:
                #        varsToKeepSyst = {}
                #        for name, var in varsToKeep.items():
                #            if var is None:
                #                varsToKeepSyst[name] = var
                #            else:
                #                name += '_'+systematic
                #                varsToKeepSyst[name] = op.forSystematicVariation(var,systematic)
                #        plots.append(Skim(f"{channel}_{jet_selection}_{systematic}", varsToKeepSyst, selObj.sel))
                #        break
        return plots

    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        super(SkimmerMEMHHDL, self).postProcess(taskList, config, workdir, resultsdir, forSkimmer=True)
        import pandas as pd
        from bamboo.root import gbl
        from bamboo.plots import Skim
        from bamboo.analysisutils import loadPlotIt
        plotItConfig, samples, plots, systematics, legend = loadPlotIt(config, [], eras=None, workdir=workdir, resultsdir=resultsdir, readCounters=self.readCounters)

        skims = [ap for ap in self.plotList if isinstance(ap, Skim)]
        for skim in skims:
            print (f'Starting parquet file production on skim {skim.name}')
            for smp in samples:
                print (f'\tLooking at sample {smp.name}')
                frames = []
                for cb in (smp.files if hasattr(smp, "files") else [smp]):
                    tree = cb.tFile.Get(skim.treeName)
                    print (f'\t\tLooking at file {cb.name}')
                    if not tree:
                        print( f"\t\t... KEY TTree {skim.treeName} does not exist, skipping this files")
                    else:
                        cols = gbl.ROOT.RDataFrame(cb.tFile.Get(skim.treeName)).AsNumpy()
                        cols["scale"] = cb.scale
                        cols["process"] = smp.name
                        cols["file"] = cb.name
                        frames.append(pd.DataFrame(cols))
                if len(frames) > 0:
                    df = pd.concat(frames)
                    pqoutname = os.path.join(resultsdir, f"{smp.name}.parquet")
                    df.to_parquet(pqoutname)
                    print(f"\t-> Dataframe for sample {smp.name} saved to {pqoutname}")

