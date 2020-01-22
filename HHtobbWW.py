import sys
import os
from itertools import chain

from bamboo.analysismodules import NanoAODHistoModule
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection
from bamboo.scalefactors import binningVariables_nano

from bamboo import treefunctions as op
from bamboo.plots import Plot, EquidistantBinning, SummedPlot
from bamboo import treefunctions as op

sys.path.append('/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/') # Add scripts in this directory -- TODO : make cleaner
from plotDef import makeDileptonPlots, makeJetsPlots, makeFatJetPlots
from scalefactorsbbWW import ScaleFactorsbbWW

class NanoHHTobbWW(NanoAODHistoModule):
    """ Example module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
            super(NanoHHTobbWW, self).__init__(args)
            self.calcToAdd += ["nJet", "nMuon"]
            # Set ratio plots #
            self.plotDefaults = {"show-ratio": True,
                                 "y-axis": "Events",
                                 #"log-y"  : "both",
                                 "ratio-y-axis-range" : [0.5,2],
                                 "ratio-y-axis" : 'Ratio Data/MC'}


    def prepareTree(self, tree, sample=None, sampleCfg=None, enableSystematics=None):
        # JEC's Recommendation for Full RunII: https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
        # JER : -----------------------------: https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
        # Get base aguments #
        tree,noSel,be,lumiArgs = super(NanoHHTobbWW,self).prepareTree(tree, sample=sample, sampleCfg=sampleCfg)
        era = sampleCfg['era']
        triggersPerPrimaryDataset = {}
        from bamboo.analysisutils import configureJets ,configureRochesterCorrection

        ## Check distributed option #
        isNotWorker = (self.args.distributed != "worker") 

        # Rochester and JEC corrections (depends on era) #     

        ############################################################################################
        # ERA 2016 #
        ############################################################################################
        if era == "2016":
            # Rochester corrections #
            configureRochesterCorrection(variProxy  = tree._Muon,
                                         paramsFile = os.path.join(os.path.dirname(__file__), "data", "RoccoR2016.txt"),
                                         isMC       = self.isMC(sample),
                                         backend    = be, 
                                         uName      = sample)

            # Trigger efficiencies #
            # tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_(DZ) 
            # tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_(DZ) 
            # -> can be present with, without and without both the DZ
            if self.isMC(sample) or "2016F" in sample or "2016G" in sample:# or "2016H" in sample:
                # Found in 2016F : both
                # Found in 2016G : both
                triggersPerPrimaryDataset = {
                    "SingleMuon" :  [ tree.HLT.IsoMu24],
                    "SingleEG"   :  [ tree.HLT.Ele27_WPTight_Gsf],
                    "DoubleMuon" :  [ tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL,
                                      tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ],
                    "DoubleEGamma": [ tree.HLT.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ],
                    "MuonEG":       [ tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,
                                      tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,
                                      tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,
                                      tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL]}
            elif "2016H" in sample:
                # Found in 2016H : has DZ but not without
                triggersPerPrimaryDataset = {
                    "SingleMuon" :  [ tree.HLT.IsoMu24],
                    "SingleEG"   :  [ tree.HLT.Ele27_WPTight_Gsf],
                    "DoubleMuon" :  [ tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL,
                                      tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ],
                    "DoubleEGamma": [ tree.HLT.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ],
                    "MuonEG":       [ tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,
                                      tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ]}

            elif "2016B" in sample or "2016C" in sample or "2016D" in sample or "2016E" in sample : 
                # Found in 2016B : only without DZ
                # Found in 2016C : only without DZ
                # Found in 2016D : only without DZ
                # Found in 2016E : only without DZ

                triggersPerPrimaryDataset = {
                    "SingleMuon" :  [ tree.HLT.IsoMu24],
                    "SingleEG"   :  [ tree.HLT.Ele27_WPTight_Gsf],
                    "DoubleMuon" :  [ tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL,
                                      tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ],
                    "DoubleEGamma": [ tree.HLT.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ],
                    "MuonEG":       [ tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, 
                                      tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL]} 

            # Jet treatment #
            #cachJEC_dir = '/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/cacheJEC'
            #if self.isMC(sample):   # if MC -> needs smearing
            #    configureJets(variProxy             = tree._Jet, 
            #                  jetType               = "AK4PFchs",
            #                  jec                   = "Summer16_07Aug2017_V20_MC",
            #                  smear                 = "Summer16_25nsV1_MC",
            #                  jesUncertaintySources = ["Total"],
            #                  mayWriteCache         = isNotWorker,
            #                  isMC                  = self.isMC(sample),
            #                  backend               = be, 
            #                  uName                 = sample,
            #                  cachedir              = cachJEC_dir)
            #else:                   # If data -> extract info from config 
            #    jecTag = None
            #    if "2016B" in sample or "2016C" in sample or "2016D" in sample:
            #        jecTag = "Summer16_07Aug2017BCD_V11_DATA"
            #    elif "2016E" in sample or "2016F" in sample:
            #        jecTag = "Summer16_07Aug2017EF_V11_DATA"
            #    elif "2016G" in sample or "2016H" in sample:
            #        jecTag = "Summer16_07Aug2017GH_V11_DATA"
            #        configureJets(tree,"Jet","AK4PFchs",
            #            jec="Summer16_07Aug2017GH_V11_DATA", mayWriteCache=isNotWorker,cachedir=cachJEC_dir)
            #    configureJets(variProxy             = tree._Jet, 
            #                  jetType               = "AK4PFchs",
            #                  jec                   = jecTag,
            #                  mayWriteCache         = isNotWorker,
            #                  cachedir              = cachJEC_dir)

        ############################################################################################
        # ERA 2017 #
        ############################################################################################
        #elif era == "2017":
        #    configureRochesterCorrection(tree._Muon.calc,os.path.join(os.path.dirname(__file__), "data", "RoccoR2017.txt"))
        #    triggersPerPrimaryDataset = {
        #        "DoubleMuon" : [ tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL,
        #                         tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ,
        #                         tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8,
        #                         tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 ],
        #        "DoubleEG"   : [ tree.HLT.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL ],
        #        # it's recommended to not use the DZ  version  for 2017 and 2018, it would be a needless efficiency loss
        #        #---> https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary
        #        "MuonEG"     : [ tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,
        #                         tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,
        #                         tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,
        #                         tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ ]
        #        }
        #    if self.isMC(sample):
        #        configureJets(tree, "Jet", "AK4PFchs",
        #            jec="Fall17_17Nov2017_V32_MC",
        #            smear="Fall17_V3_MC",
        #            jesUncertaintySources=["Total"], mayWriteCache=isNotWorker)
        #    else:
        #        if "2017B" in sample:
        #            configureJets(tree, "Jet", "AK4PFchs", jec="Fall17_17Nov2017B_V32_DATA", mayWriteCache=isNotWorker)
        #        elif "2017C" in sample:
        #            configureJets(tree, "Jet", "AK4PFchs", jec="Fall17_17Nov2017C_V32_DATA", mayWriteCache=isNotWorker)
        #        elif "2017D" in sample or "2017E" in sample:
        #            configureJets(tree, "Jet", "AK4PFchs", jec="Fall17_17Nov2017DE_V32_DATA", mayWriteCache=isNotWorker)
        #        elif "2017F" in sample:
        #            configureJets(tree, "Jet", "AK4PFchs", jec="Fall17_17Nov2017F_V32_DATA", mayWriteCache=isNotWorker)

        ############################################################################################
        # ERA 2018 #
        ############################################################################################
        #elif era == "2018":
        #    configureRochesterCorrection(tree._Muon.calc,os.path.join(os.path.dirname(__file__), "data", "RoccoR2018.txt"))
        #    triggersPerPrimaryDataset = {
        #        "DoubleMuon" : [ tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL,
        #                         tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ,
        #                         tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8,
        #                         tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 ],
        #        "EGamma"     : [ tree.HLT.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL ], 
        #        "MuonEG"     : [ tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,
        #                         tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,
        #                         tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,
        #                         tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ ]
        #        }
        #    if self.isMC(sample):
        #        configureJets(tree, "Jet", "AK4PFchs",
        #            jec="Autumn18_V8_MC",
        #            smear="Autumn18_V1_MC",
        #            jesUncertaintySources=["Total"], mayWriteCache=isNotWorker)
        #    else:
        #        if "2018A" in sample:
        #            configureJets(tree, "Jet", "AK4PFchs",
        #                jec="Autumn18_RunA_V8_DATA", mayWriteCache=isNotWorker)
        #        elif "2018B" in sample:
        #            configureJets(tree, "Jet", "AK4PFchs",
        #                jec="Autumn18_RunB_V8_DATA", mayWriteCache=isNotWorker)
        #        elif "2018C" in sample:
        #            configureJets(tree, "Jet", "AK4PFchs",
        #                jec="Autumn18_RunC_V8_DATA", mayWriteCache=isNotWorker)
        #        elif "2018D" in sample:
        #            configureJets(tree, "Jet", "AK4PFchs",
        #                jec="Autumn18_RunD_V8_DATA", mayWriteCache=isNotWorker)
        else:
            raise RuntimeError("Unknown era {0}".format(era))

        # Get weights #
        if self.isMC(sample):
            noSel = noSel.refine("genWeight", weight=tree.genWeight, cut=op.OR(*chain.from_iterable(triggersPerPrimaryDataset.values())))
        else:
            noSel = noSel.refine("withTrig", cut=makeMultiPrimaryDatasetTriggerSelection(sample, triggersPerPrimaryDataset))

        return tree,noSel,be,lumiArgs
        
    def definePlots(self, t, noSel, sample=None, sampleCfg=None): 
        # Some imports #
        from bamboo.analysisutils import forceDefine

        era = sampleCfg['era']

        # Initialize scalefactors class #
        SF = ScaleFactorsbbWW()

        #############################################################################
        ##########################    Pile-up    ####################################
        #############################################################################
        puWeightsFile = None
        if era == "2016":
            sfTag="94X"
            puWeightsFile = os.path.join(os.path.dirname(__file__), "data", "puweights2016.json")
        elif era == "2017":
            sfTag="94X"     
            puWeightsFile = os.path.join(os.path.dirname(__file__), "data", "puweights2017.json")
        elif era == "2018":
            sfTag="102X"
            puWeightsFile = os.path.join(os.path.dirname(__file__), "data", "puweights2018.json")
        if self.isMC(sample) and puWeightsFile is not None:
            from bamboo.analysisutils import makePileupWeight
            noSel = noSel.refine("puWeight", weight=makePileupWeight(puWeightsFile, t.Pileup_nTrueInt, systName="pileup"))
        isMC = self.isMC(sample)
        plots = []

        forceDefine(t._Muon.calcProd, noSel)
        #forceDefine(t._Jet.calcProd, noSel) # calculate once per event (for every event)

        #############################################################################
        ################################  Muons #####################################
        #############################################################################
        # Wp // 2016- 2017 -2018   // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Muon_Isolation
        muonsByPt = op.sort(t.Muon, lambda mu : -mu.p4.Pt())
        muons = op.select(muonsByPt, lambda mu : op.AND(mu.p4.Pt() > 15., op.abs(mu.p4.Eta()) < 2.4, mu.tightId, mu.pfRelIso04_all<0.15))
            # Subleading lepton pt cut is at 15 GeV so better start at that point
            # isolation : tight (Muon::PFIsoTight) cut value = 0.15 (ε~0.95)
      
        # Scalefactors #
        if self.isMC(sample):
            if era=="2016":
                muMediumIDSF = SF.get_scalefactor("lepton", ("muon_{0}_{1}".format(era, sfTag), "id_tight"), combine="weight", systName="muid")
                muMediumISOSF = SF.get_scalefactor("lepton", ("muon_{0}_{1}".format(era, sfTag), "iso_tight_id_medium"), combine="weight", systName="muiso")
                TrkIDSF = SF.get_scalefactor("lepton", ("muon_{0}_{1}".format(era, sfTag), "idtrk_highpt"), combine="weight")         # Need to ask what it is
                TrkISOSF = SF.get_scalefactor("lepton", ("muon_{0}_{1}".format(era, sfTag), "isotrk_loose_idtrk_highptidandipcut"), combine="weight") # Need to ask what it is
            else:
                muMediumIDSF = SF.get_scalefactor("lepton", ("muon_{0}_{1}".format(era, sfTag), "id_medium"), systName="muid")
                muMediumISOSF = SF.get_scalefactor("lepton", ("muon_{0}_{1}".format(era, sfTag), "iso_tight_id_medium"), systName="muiso") 

        #############################################################################
        #############################  Electrons  ###################################
        #############################################################################
        #Wp  // 2016: Electron_cutBased_Sum16==3  -> medium     // 2017 -2018  : Electron_cutBased ==3   --> medium ( Fall17_V2)
        # asking for electrons to be in the Barrel region with dz<1mm & dxy< 0.5mm   //   Endcap region dz<2mm & dxy< 0.5mm 
        electronsByPt = op.sort(t.Electron, lambda ele : -ele.p4.Pt())
        electrons = op.select(electronsByPt, lambda ele : op.AND(ele.p4.Pt() > 15., op.abs(ele.p4.Eta()) < 2.5 , ele.cutBased>=4 )) # //cut-based ID Fall17 V2 the recommended one from POG for the FullRunII
            # electron cut based ID = 0:fail, 1: veto, 2:loose, 3:medium, 4:tight (From root file -> Events tree)
            # Subleading lepton pt cut is at 15 GeV so better start at that point

        # Scalefactors #
        if self.isMC(sample):
            elMediumIDSF = SF.get_scalefactor("lepton", ("electron_{0}_{1}".format(era,sfTag), "id_tight"), systName="elid")

        #############################################################################
        #############################   Dilepton  ###################################
        #############################################################################
        # Combine dileptons #
        OsElEl = op.combine(electrons, N=2, pred=lambda el1,el2 : op.AND(el1.charge != el2.charge , el1.p4.Pt() > 25, el2.p4.Pt() > 15 ))
        OsMuMu = op.combine(muons, N=2, pred=lambda mu1,mu2 : op.AND(mu1.charge != mu2.charge , mu1.p4.Pt() > 25, mu2.p4.Pt() > 15 ))
        OsElMu = op.combine((electrons, muons), pred=lambda el,mu : op.AND(el.charge != mu.charge ,el.p4.Pt() > 25, mu.p4.Pt() > 15 ))
        OsMuEl = op.combine((muons, electrons), pred=lambda mu,el : op.AND(el.charge != mu.charge ,el.p4.Pt() > 15, mu.p4.Pt() > 25 ))

        # Plots Numbers of dilepton in each channel #
        plots.append(Plot.make1D("ElEl_channel",op.rng_len(OsElEl),noSel,EquidistantBinning(5,0,5.),title='Number of dilepton events in ElEl channel',xTitle='N_{dilepton} (ElEl channel)'))
        plots.append(Plot.make1D("MuMu_channel",op.rng_len(OsMuMu),noSel,EquidistantBinning(5,0,5.),title='Number of dilepton events in MuMu channel',xTitle='N_{dilepton} (MuMu channel)'))
        plots.append(Plot.make1D("ElMu_channel",op.rng_len(OsElMu),noSel,EquidistantBinning(5,0,5.),title='Number of dilepton events in ElMu channel',xTitle='N_{dilepton} (ElMu channel)'))
        plots.append(Plot.make1D("MuEl_channel",op.rng_len(OsMuEl),noSel,EquidistantBinning(5,0,5.),title='Number of dilepton events in MuEl channel',xTitle='N_{dilepton} (MuEl channel)'))

        # Scalefactors #
       # if self.isMC(sample):
       #     doubleEleTrigSF = SF.get_scalefactor("dilepton", ("doubleEleLeg_HHMoriond17_2016"), systName="eleltrig")     
       #     doubleMuTrigSF  = SF.get_scalefactor("dilepton", ("doubleMuLeg_HHMoriond17_2016"), systName="mumutrig")    
       #     elemuTrigSF     = SF.get_scalefactor("dilepton", ("elemuLeg_HHMoriond17_2016"), systName="elmutrig")
       #     mueleTrigSF     = SF.get_scalefactor("dilepton", ("mueleLeg_HHMoriond17_2016"), systName="mueltrig")

       # llSF =  {
       #     "ElEl" : (lambda ll : [ elMediumIDSF(ll[0][0]),                                                                     # First lepton SF
       #                             elMediumIDSF(ll[0][1]),                                                                     # Second lepton SF
       #                             doubleEleTrigSF(ll[0])]),                                                                   # Dilepton SF
       #     "MuMu" : (lambda ll : [ muMediumIDSF(ll[0][0]), muMediumISOSF(ll[0][0]), TrkIDSF(ll[0][0]), TrkISOSF(ll[0][0]),     # First lepton SF
       #                             muMediumIDSF(ll[0][1]), muMediumISOSF(ll[0][1]), TrkIDSF(ll[0][1]), TrkISOSF(ll[0][1]),     # Second lepton SF
       #                             doubleMuTrigSF(ll[0])]),                                                                    # Dilepton SF
       #     "ElMu" : (lambda ll : [ elMediumIDSF(ll[0][0]),                                                                     # First lepton SF
       #                             muMediumIDSF(ll[0][1]), muMediumISOSF(ll[0][1]), TrkIDSF(ll[0][1]), TrkISOSF(ll[0][1]),     # Second lepton SF
       #                             elemuTrigSF(ll[0])]),                                                                       # Dilepton SF
       #     "MuEl" : (lambda ll : [ muMediumIDSF(ll[0][0]), muMediumISOSF(ll[0][0]), TrkIDSF(ll[0][0]), TrkISOSF(ll[0][0]),     # First lepton SF
       #                             elMediumIDSF(ll[0][1]),                                                                     # Second lepton SF
       #                             mueleTrigSF(ll[0])]),                                                                       # Dilepton SF
       #         # ll is a proxy list of dileptons 
       #         # ll[0] is the first dilepton 
       #         # ll[0][0] is the first lepton and ll[0][1] the second in the dilepton
       #         }
        llSFApplied = {
            "ElEl": None, #llSF["ElEl"](OsElEl) if isMC else None,
            "MuMu": None, #llSF["MuMu"](OsMuMu) if isMC else None,
            "ElMu": None, #llSF["ElMu"](OsElMu) if isMC else None,
            "MuEl": None, #llSF["MuEl"](OsMuEl) if isMC else None,
                      }

        # Selection #
        hasOsElEl = noSel.refine("hasOsElEl",
                                 cut = [op.rng_len(OsElEl) >= 1,                        # Require at least one dilepton ElEl
                                        (op.rng_len(electrons)+op.rng_len(muons))<=2],  # Not more than two tight leptons
                                 weight = llSFApplied["ElEl"])
        hasOsMuMu = noSel.refine("hasOsMuMu",
                                 cut = [op.rng_len(OsMuMu) >= 1,                        # Require at least one dilepton MuMu
                                        (op.rng_len(electrons)+op.rng_len(muons))<=2],  # Not more than two tight leptons
                                 weight = llSFApplied["MuMu"])
        hasOsElMu = noSel.refine("hasOsElMu",
                                 cut = [op.rng_len(OsElMu) >= 1,                        # Require at least one dilepton ElMu
                                        (op.rng_len(electrons)+op.rng_len(muons))<=2],  # Not more than two tight leptons
                                 weight = llSFApplied["ElMu"])
        hasOsMuEl = noSel.refine("hasOsMuEl",
                                 cut = [op.rng_len(OsMuEl) >= 1,                        # Require at least one dilepton MuEl
                                        (op.rng_len(electrons)+op.rng_len(muons))<=2],  # Not more than two tight leptons
                                 weight = llSFApplied["MuEl"])

        hasOsCutChannelList = [
                             {'channel':'ElEl','sel':hasOsElEl,'dilepton':OsElEl[0],'suffix':'hasOsElEl'},
                             {'channel':'MuMu','sel':hasOsMuMu,'dilepton':OsMuMu[0],'suffix':'hasOsMuMu'},
                             {'channel':'ElMu','sel':hasOsElMu,'dilepton':OsElMu[0],'suffix':'hasOsElMu'},
                             {'channel':'MuEl','sel':hasOsMuEl,'dilepton':OsMuEl[0],'suffix':'hasOsMuEl'},
                           ]
        for channelDict in hasOsCutChannelList:
            plots.extend(makeDileptonPlots(self, **channelDict))

        # Dilepton Z peak exclusion (charge already done in previous selection) # 
        lambda_lowMllCut    = lambda dilep: op.invariant_mass(dilep[0].p4, dilep[1].p4)>12.
        lambda_outZ         = lambda dilep: op.NOT(op.in_range(80.,op.invariant_mass(dilep[0].p4, dilep[1].p4),100.))
        #lambda_mllLowerband = lambda dilep: op.in_range(12.,op.invariant_mass(dilep[0].p4, dilep[1].p4),80.)
        #lambda_mllUpperband = lambda dilep: op.invariant_mass(dilep[0].p4, dilep[1].p4)>100.
        #lambda_mllCut       = lambda dilep: op.OR(lambda_mllLowerband(dilep),lambda_mllUpperband(dilep))
        #lambda_OSdilep      = lambda dilep: dilep[0].charge*dilep[1].charge == -1 # Must have been done before

        hasOsElElLowMllCutOutZ = hasOsElEl.refine("hasOsElElLowMllCutOutZ",cut=[lambda_lowMllCut(OsElEl[0]),lambda_outZ(OsElEl[0])])
        hasOsMuMuLowMllCutOutZ = hasOsMuMu.refine("hasOsMuMuLowMllCutOutZ",cut=[lambda_lowMllCut(OsMuMu[0]),lambda_outZ(OsMuMu[0])])
        hasOsElMuLowMllCut     = hasOsElMu.refine("hasOsElMuLowMllCut",cut=[lambda_lowMllCut(OsElMu[0])]) # Z peak cut not needed because Opposite Flavour
        hasOsMuElLowMllCut     = hasOsMuEl.refine("hasOsMuElLowMllCut",cut=[lambda_lowMllCut(OsMuEl[0])]) # Z peak cut not needed because Opposite Flavour

        hasOsMllCutChannelList = [
                             {'channel':'ElEl','sel':hasOsElElLowMllCutOutZ,'dilepton':OsElEl[0],'suffix':'hasOsElElLowMllCutOutZ'},
                             {'channel':'MuMu','sel':hasOsMuMuLowMllCutOutZ,'dilepton':OsMuMu[0],'suffix':'hasOsMuMuLowMllCutOutZ'},
                             {'channel':'ElMu','sel':hasOsElMuLowMllCut,    'dilepton':OsElMu[0],'suffix':'hasOsElMuLowMllCut'},
                             {'channel':'MuEl','sel':hasOsMuElLowMllCut,    'dilepton':OsMuEl[0],'suffix':'hasOsMuElLowMllCut'},
                           ]

        for channelDict in hasOsMllCutChannelList:
            plots.extend(makeDileptonPlots(self, **channelDict))

        #############################################################################
        ################################  Jets  #####################################
        #############################################################################
        # select jets   // 2016 - 2017 - 2018   ( j.jetId &2) ->      tight jet ID
        jetsByPt = op.sort(t.Jet, lambda jet : -jet.p4.Pt())
        jetsSel = op.select(jetsByPt, lambda j : op.AND(j.p4.Pt() > 20., op.abs(j.p4.Eta())< 2.4, (j.jetId &2)))        # Jets = AK4 jets
        fatjetsByPt = op.sort(t.FatJet, lambda fatjet : -fatjet.p4.Pt())
        fatjetsSel = op.select(fatjetsByPt, lambda j : op.AND(j.p4.Pt() > 20., op.abs(j.p4.Eta())< 2.4, (j.jetId &2)))        # FatJets = AK8 jets

        # exclude from the jetsSel any jet that happens to include within its reconstruction cone a muon or an electron.
        jets = op.select(jetsSel, lambda j : op.AND(op.NOT(op.rng_any(electrons, lambda ele : op.deltaR(j.p4, ele.p4) < 0.3 )), op.NOT(op.rng_any(muons, lambda mu : op.deltaR(j.p4, mu.p4) < 0.3 ))))
        fatjets = op.select(fatjetsSel, lambda j : op.AND(op.NOT(op.rng_any(electrons, lambda ele : op.deltaR(j.p4, ele.p4) < 0.3 )), op.NOT(op.rng_any(muons, lambda mu : op.deltaR(j.p4, mu.p4) < 0.3 ))))
        # Boosted and resolved jets categories #
        # Boosted -> at least one AK8 jet (fatjet) with at least one subjet passing medium working point of DeepCSV (btagDeepB branch)
        # Resolved -> at least two Ak4 jets (jet) with at least one passing the medium working point of DeepJet (btagDeepFlavB branch)
        if era == "2016": # Must check that subJet exists before looking at the btag
            lambda_boosted  = lambda fatjet : op.OR(op.AND(fatjet.subJet1._idx.result != -1,fatjet.subJet1.btagDeepB > 0.6321), op.AND(fatjet.subJet2._idx.result != -1,fatjet.subJet2.btagDeepB > 0.6321))
            lambda_resolved = lambda jet    : jet.btagDeepFlavB > 0.3093
            lambda_notResolved = lambda jet : jet.btagDeepFlavB <= 0.3093

        elif era =="2017":
            lambda_boosted  = lambda fatjet : op.OR(op.AND(fatjet.subJet1._idx.result != -1,fatjet.subJet1.btagDeepB > 0.4941), op.AND(fatjet.subJet2._idx.result != -1,fatjet.subJet2.btagDeepB > 0.4941))
            lambda_resolved = lambda jet    : jet.btagDeepFlavB > 0.3033
            lambda_notResolved = lambda jet : jet.btagDeepFlavB <= 0.3033
        elif era == "2018":
            lambda_boosted  = lambda fatjet : op.OR(op.AND(fatjet.subJet1._idx.result != -1,fatjet.subJet1.btagDeepB > 0.4184), op.AND(fatjet.subJet2._idx.result != -1,fatjet.subJet2.btagDeepB > 0.4184))
            lambda_resolved = lambda jet    : jet.btagDeepFlavB > 0.2770
            lambda_notResolved = lambda jet : jet.btagDeepFlavB <= 0.2770

        # Select the bjets we want #
        bjetsBoosted        = op.select(fatjets, lambda_boosted)
        bjetsResolved       = op.select(jets, lambda_resolved)
        lightjetsResolved   = op.select(jets, lambda_notResolved) # To plot the case when 1 bjet + 1 lightjet

        # Scalefactors : 2017 and 2018 not yet present in the dict #
        DeepCSVMediumSFApplied = None
        DeepJetMediumSFApplied = None
        if self.isMC(sample):
            pass
            #DeepCSVTag_discriVar = {"BTagDiscri": lambda j : j.btagDeepB}
            #DeepJetTag_discriVar = {"BTagDiscri": lambda j : j.btagDeepFlavB}
            #DeepCSVMediumSF = SF.get_scalefactor("jet", ("subjet_btag_"+era+"_"+sfTag, "DeepCSV_medium"), additionalVariables=DeepCSVTag_discriVar, systName="deepcsv") # For BOOSTED (btag on subjet)
            #DeepJetMediumSF = SF.get_scalefactor("jet", ("btag_"+era+"_"+sfTag, "DeepJet_medium"), additionalVariables=DeepJetTag_discriVar, systName="deepjet") # For RESOLVED
            #DeepCSVMediumSFApplied = [DeepCSVMediumSF(bjetsBoosted[0].subJet1)] # Must be applied on subjets : need to check each time which one has been btagged
            #DeepJetMediumSFApplied = [DeepJetMediumSF(bjetsResolved[0])]

    
        # Define the boosted and Resolved (+exclusive) selections #
        hasBoostedJets = noSel.refine("hasBoostedJets",
                               cut=[op.rng_len(bjetsBoosted)>=1],
                               weight = DeepCSVMediumSFApplied)
        hasResolvedJets = noSel.refine("hasResolvedJets",
                               cut=[op.rng_len(jets)>=2,op.rng_len(bjetsResolved)>=1],
                               weight = DeepJetMediumSFApplied)
        hasExclusiveResolvedJets = noSel.refine("hasExclusiveResolved",
                               cut=[op.rng_len(jets)>=2, op.rng_len(bjetsResolved)>=1,op.rng_len(bjetsBoosted)==0],
                               weight = DeepJetMediumSFApplied)
        hasExclusiveBoostedJets = noSel.refine("hasExclusiveBoostedJets",
                               cut=[op.rng_len(bjetsBoosted)>=1,op.OR(op.rng_len(jets)<=1,op.rng_len(bjetsResolved)==0)],
                               weight = DeepCSVMediumSFApplied)

        hasNotBoostedJets = noSel.refine("hasNotBoostedJets",
                               cut=[op.rng_len(bjetsBoosted)==0])
        hasNotResolvedJets = noSel.refine("hasNotResolvedJets",
                               cut=[op.OR(op.rng_len(jets)<=1,op.rng_len(bjetsResolved)==0)])
        hasBoostedAndResolvedJets = noSel.refine("hasBoostedAndResolvedJets", 
                               cut=[op.rng_len(bjetsBoosted)>=1,op.rng_len(jets)>=2,op.rng_len(bjetsResolved)>=1])
        hasNotBoostedAndResolvedJets = noSel.refine("hasNotBoostedAndResolvedJets", 
                               cut=[op.OR(op.rng_len(bjetsBoosted)==0,op.rng_len(jets)<=1,op.rng_len(bjetsResolved)==0)])
        hasNotExclusiveResolvedJets = noSel.refine("hasNotExclusiveResolved",
                               cut=[op.OR(op.OR(op.rng_len(jets)<=1,op.rng_len(bjetsResolved)==0),op.AND(op.rng_len(bjetsBoosted)>=1,op.rng_len(jets)>=2,op.rng_len(bjetsResolved)>=1))])
        hasNotExclusiveBoostedJets = noSel.refine("hasNotExclusiveBoostedJets",
                               cut=[op.OR(op.rng_len(bjetsBoosted)==0,op.AND(op.rng_len(jets)>=2,op.rng_len(bjetsResolved)>=1))])
        # Note that these selection should be done the same (and SF) when combining the OS lepton selection #
        # Because so far there is no way to "concatenate" two refine's #

        # Counting events from different selections for debugging #
        # Passing Boosted selection #
        PassedBoosted = Plot.make1D("PassedBoosted",
                                    op.c_int(1),
                                    hasBoostedJets,
                                    EquidistantBinning(2,0.,2.),
                                    title='Passed Boosted',
                                    xTitle='Passed Boosted')
        FailedBoosted = Plot.make1D("FailedBoosted",
                                    op.c_int(0),
                                    hasNotBoostedJets,
                                    EquidistantBinning(2,0.,2.),
                                    title='Failed Boosted',
                                    xTitle='Failed Boosted')
        plots.append(SummedPlot("BoostedCase",
                                [FailedBoosted,PassedBoosted],
                                xTitle="Boosted selection"))

        # Passing Resolved selection #
        PassedResolved = Plot.make1D("PassedResolved",
                                     op.c_int(1),
                                     hasResolvedJets,
                                     EquidistantBinning(2,0.,2.),
                                     title='Passed Resolved',
                                     xTitle='Passed Resolved')
        FailedResolved = Plot.make1D("FailedResolved",
                                     op.c_int(0),
                                     hasNotResolvedJets,
                                     EquidistantBinning(2,0.,2.),
                                     title='Failed Resolved',
                                     xTitle='Failed Resolved')
        plots.append(SummedPlot("ResolvedCase",
                                [FailedResolved,PassedResolved],
                                xTitle="Resolved selection"))

        # Passing Exclusive Resolved (Resolved AND NOT Boosted) #
        PassedExclusiveResolved = Plot.make1D("PassedExclusiveResolved",
                                              op.c_int(1),
                                              hasExclusiveResolvedJets,
                                              EquidistantBinning(2,0.,2.),
                                              title='Passed Exclusive Resolved',
                                              xTitle='Passed Exclusive Resolved')
        FailedExclusiveResolved = Plot.make1D("FailedExclusiveResolved",
                                              op.c_int(0),
                                              hasNotExclusiveResolvedJets,
                                              EquidistantBinning(2,0.,2.),
                                              title='Failed Exclusive Resolved',
                                              xTitle='Failed Exclusive Resolved')
        plots.append(SummedPlot("ExclusiveResolvedCase",
                                [FailedExclusiveResolved,PassedExclusiveResolved],
                                xTitle="Exclusive Resolved selection"))

        # Passing Exclusive Boosted (Boosted AND NOT Resolved) #
        PassedExclusiveBoosted = Plot.make1D("PassedExclusiveBoosted",
                                             op.c_int(1),
                                             hasExclusiveBoostedJets,
                                             EquidistantBinning(2,0.,2.),
                                             title='Passed Exclusive Boosted',
                                             xTitle='Passed Exclusive Boosted')
        FailedExclusiveBoosted = Plot.make1D("FailedExclusiveBoosted",
                                             op.c_int(0),
                                             hasNotExclusiveBoostedJets,
                                             EquidistantBinning(2,0.,2.),
                                             title='Failed Exclusive Boosted',
                                             xTitle='Failed Exclusive Boosted')
        plots.append(SummedPlot("ExclusiveBoostedCase",
                                [FailedExclusiveBoosted,PassedExclusiveBoosted],
                                xTitle="Exclusive Boosted selection"))

        # Passing Boosted AND Resolved #
        PassedBoth = Plot.make1D("PassedBoth",
                                 op.c_int(1),
                                 hasBoostedAndResolvedJets,
                                 EquidistantBinning(2,0.,2.),
                                 title='Passed Both Boosted and Resolved',
                                 xTitle='Passed Boosted and Resolved')
        FailedBoth = Plot.make1D("FailedBoth", # Means failed the (Boosted AND Resolved) = either one or the other 
                                 op.c_int(0),
                                 hasNotBoostedAndResolvedJets,
                                 EquidistantBinning(2,0.,2.),
                                 title='Failed combination Boosted and Resolved',
                                 xTitle='Failed combination')
        plots.append(SummedPlot("BoostedAndResolvedCase",[FailedBoth,PassedBoth],xTitle="Boosted and Resolved selection"))

        # Count number of boosted and resolved jets #
        plots.append(Plot.make1D("NBoostedJets",
                     op.rng_len(bjetsBoosted),
                     hasBoostedJets,
                     EquidistantBinning(5,0.,5.),
                     title='Number of boosted jets in boosted case',
                     xTitle='N boosted bjets'))
        plots.append(Plot.make1D("NResolvedJets",
                     op.rng_len(bjetsResolved),
                     hasResolvedJets,
                     EquidistantBinning(5,0.,5.),
                     title='Number of resolved jets in resolved case',
                     xTitle='N resolved bjets'))

        # Plot number of subjets in the boosted fatjets #
        lambda_noSubjet  = lambda fatjet : op.AND(fatjet.subJet1._idx.result == -1, op.AND(fatjet.subJet2._idx.result == -1 )) 
        lambda_oneSubjet = lambda fatjet : op.AND(fatjet.subJet1._idx.result != -1, op.AND(fatjet.subJet2._idx.result == -1 ))
        lambda_twoSubjet = lambda fatjet : op.AND(fatjet.subJet1._idx.result != -1, op.AND(fatjet.subJet2._idx.result != -1 ))

        hasNoSubjet   = hasBoostedJets.refine("hasNoSubjet",
                       cut=[lambda_noSubjet(bjetsBoosted[0])])
        hasOneSubjet  = hasBoostedJets.refine("hasOneSubjet",
                       cut=[lambda_oneSubjet(bjetsBoosted[0])])
        hasTwoSubjet  = hasBoostedJets.refine("hasTwoSubjet",
                       cut=[lambda_twoSubjet(bjetsBoosted[0])])

        plot_hasNoSubjet = Plot.make1D("plot_hasNoSubjet", # Fill bin 0
                                       op.c_int(0),
                                       hasNoSubjet,
                                       EquidistantBinning(3,0.,3.),
                                       title='Boosted jet without subjet')
        plot_hasOneSubjet = Plot.make1D("plot_hasOneSubjet", # Fill bin 1
                                        op.c_int(1),
                                        hasOneSubjet,
                                        EquidistantBinning(3,0.,3.),
                                        title='Boosted jet with one subjet')
        plot_hasTwoSubjet = Plot.make1D("plot_hasTwoSubjet", # Fill bin 2
                                        op.c_int(2),
                                        hasTwoSubjet,
                                        EquidistantBinning(3,0.,3.),title='Boosted jet with two subjets')
        plots.append(SummedPlot("NumberOfSubjets",
                                [plot_hasNoSubjet,plot_hasOneSubjet,plot_hasTwoSubjet],
                                xTitle="Number of subjets in boosted jet"))

        # Plot jets quantities without the dilepton selections #
        
        JetNoChannelBoostedList = [
                             {'channel':'NoChannel','sel':hasBoostedJets,'fatjets':bjetsBoosted,'suffix':"hasBoostedJets"},
                             {'channel':'NoChannel','sel':hasExclusiveBoostedJets,'fatjets':bjetsBoosted,'suffix':"hasExclusiveBoostedJets"}
                           ]
        JetNoChannelResolvedList = [
                             {'channel':'NoChannel','sel':hasResolvedJets,'bjets':bjetsResolved,'lightjets':lightjetsResolved,'alljets':jets,'suffix':"hasResolvedJets"},
                             {'channel':'NoChannel','sel':hasExclusiveResolvedJets,'bjets':bjetsResolved,'lightjets':lightjetsResolved,'alljets':jets,'suffix':"hasExclusiveResolvedJets"}
                           ]
        for channelDict in JetNoChannelBoostedList:
            plots.extend(makeFatJetPlots(self, **channelDict))

        for channelDict in JetNoChannelResolvedList:
            plots.extend(makeJetsPlots(self, **channelDict))

        #############################################################################
        ##################### Jets + Dilepton combination ###########################
        #############################################################################

        ##### BOOSTED #####
        # Combine dilepton and Boosted selections #
        hasOsElElLowMllCutOutZBoostedJets = hasOsElElLowMllCutOutZ.refine("hasOsElElLowMllCutOutZBoostedJets", 
                                      cut=[op.rng_len(bjetsBoosted)>=1],
                                      weight = DeepCSVMediumSFApplied)
        hasOsMuMuLowMllCutOutZBoostedJets = hasOsMuMuLowMllCutOutZ.refine("hasOsMuMuLowMllCutOutZBoostedJets", 
                                      cut=[op.rng_len(bjetsBoosted)>=1],
                                      weight = DeepCSVMediumSFApplied)
        hasOsElMuLowMllCutBoostedJets = hasOsElMuLowMllCut.refine("hasOsElMuLowMllCutBoostedJets", 
                                      cut=[op.rng_len(bjetsBoosted)>=1],
                                      weight = DeepCSVMediumSFApplied)
        hasOsMuElLowMllCutBoostedJets = hasOsMuElLowMllCut.refine("hasOsMuElLowMllCutBoostedJets", 
                                      cut=[op.rng_len(bjetsBoosted)>=1],
                                      weight = DeepCSVMediumSFApplied)

        hasOsElElLowMllCutOutZExclusiveBoostedJets = hasOsElElLowMllCutOutZ.refine("hasOsElElLowMllCutOutZExclusiveBoostedJets", 
                                      cut=[op.rng_len(bjetsBoosted)>=1,op.OR(op.rng_len(jets)<=1,op.rng_len(bjetsResolved)==0)],
                                      weight = DeepCSVMediumSFApplied)
        hasOsMuMuLowMllCutOutZExclusiveBoostedJets = hasOsMuMuLowMllCutOutZ.refine("hasOsMuMuLowMllCutOutZExclusiveBoostedJets", 
                                      cut=[op.rng_len(bjetsBoosted)>=1,op.OR(op.rng_len(jets)<=1,op.rng_len(bjetsResolved)==0)],
                                      weight = DeepCSVMediumSFApplied)
        hasOsElMuLowMllCutExclusiveBoostedJets = hasOsElMuLowMllCut.refine("hasOsElMuLowMllCutExclusiveBoostedJets", 
                                      cut=[op.rng_len(bjetsBoosted)>=1,op.OR(op.rng_len(jets)<=1,op.rng_len(bjetsResolved)==0)],
                                      weight = DeepCSVMediumSFApplied)
        hasOsMuElLowMllCutExclusiveBoostedJets = hasOsMuElLowMllCut.refine("hasOsMuElLowMllCutExclusiveBoostedJets", 
                                      cut=[op.rng_len(bjetsBoosted)>=1,op.OR(op.rng_len(jets)<=1,op.rng_len(bjetsResolved)==0)],
                                      weight = DeepCSVMediumSFApplied)


        # Boosted + OS dilepton plots #
        hasOsMllCutBoostedChannelList = [
                             {'channel':'ElEl','sel':hasOsElElLowMllCutOutZBoostedJets,'dilepton':OsElEl[0],'suffix':'hasOsElElLowMllCutOutZBoostedJets'},
                             {'channel':'MuMu','sel':hasOsMuMuLowMllCutOutZBoostedJets,'dilepton':OsMuMu[0],'suffix':'hasOsMuMuLowMllCutOutZBoostedJets'},
                             {'channel':'ElMu','sel':hasOsElMuLowMllCutBoostedJets,    'dilepton':OsElMu[0],'suffix':'hasOsElMuLowMllCutBoostedJets'},
                             {'channel':'MuEl','sel':hasOsMuElLowMllCutBoostedJets,    'dilepton':OsMuEl[0],'suffix':'hasOsMuElLowMllCutBoostedJets'},

                             {'channel':'ElEl','sel':hasOsElElLowMllCutOutZExclusiveBoostedJets,'dilepton':OsElEl[0],'suffix':'hasOsElElLowMllCutOutZExclusiveBoostedJets'},
                             {'channel':'MuMu','sel':hasOsMuMuLowMllCutOutZExclusiveBoostedJets,'dilepton':OsMuMu[0],'suffix':'hasOsMuMuLowMllCutOutZExclusiveBoostedJets'},
                             {'channel':'ElMu','sel':hasOsElMuLowMllCutExclusiveBoostedJets,    'dilepton':OsElMu[0],'suffix':'hasOsElMuLowMllCutExclusiveBoostedJets'},
                             {'channel':'MuEl','sel':hasOsMuElLowMllCutExclusiveBoostedJets,    'dilepton':OsMuEl[0],'suffix':'hasOsMuElLowMllCutExclusiveBoostedJets'},
                           ]

        for channelDict in hasOsMllCutBoostedChannelList:
            # Dilepton plots #
            plots.extend(makeDileptonPlots(self, **channelDict))
            # Fatjet plots #
            plots.extend(makeFatJetPlots(self,
                                         sel       = channelDict['sel'],
                                         fatjets   = bjetsBoosted,
                                         suffix    = channelDict['suffix'],
                                         channel   = channelDict['channel']))
  
  
        ##### RESOLVED #####
        # Combine dilepton and Exclusive Resolved (Exclusive = NOT Boosted) selections #
        hasOsElElLowMllCutOutZResolvedJets = hasOsElElLowMllCutOutZ.refine("hasOsElElLowMllCutOutZResolvedJets", 
                                      cut=[op.rng_len(jets)>=2, op.rng_len(bjetsResolved)>=1],
                                      weight = DeepJetMediumSFApplied)
        hasOsMuMuLowMllCutOutZResolvedJets = hasOsMuMuLowMllCutOutZ.refine("hasOsMuMuLowMllCutOutZResolvedJets", 
                                      cut=[op.rng_len(jets)>=2, op.rng_len(bjetsResolved)>=1], 
                                      weight = DeepJetMediumSFApplied)
        hasOsElMuLowMllCutResolvedJets = hasOsElMuLowMllCut.refine("hasOsElMuLowMllCutResolvedJets", 
                                      cut=[op.rng_len(jets)>=2, op.rng_len(bjetsResolved)>=1],
                                      weight = DeepJetMediumSFApplied)
        hasOsMuElLowMllCutResolvedJets = hasOsMuElLowMllCut.refine("hasOsMuElLowMllCutResolvedJets", 
                                      cut=[op.rng_len(jets)>=2, op.rng_len(bjetsResolved)>=1],
                                      weight = DeepJetMediumSFApplied)

        hasOsElElLowMllCutOutZExclusiveResolvedJets = hasOsElElLowMllCutOutZ.refine("hasOsElElLowMllCutOutZExclusiveResolvedJets", 
                                      cut=[op.rng_len(jets)>=2, op.rng_len(bjetsResolved)>=1,op.rng_len(bjetsBoosted)==0],
                                      weight = DeepJetMediumSFApplied)
        hasOsMuMuLowMllCutOutZExclusiveResolvedJets = hasOsMuMuLowMllCutOutZ.refine("hasOsMuMuLowMllCutOutZExclusiveResolvedJets", 
                                      cut=[op.rng_len(jets)>=2, op.rng_len(bjetsResolved)>=1,op.rng_len(bjetsBoosted)==0], 
                                      weight = DeepJetMediumSFApplied)
        hasOsElMuLowMllCutExclusiveResolvedJets = hasOsElMuLowMllCut.refine("hasOsElMuLowMllCutExclusiveResolvedJets", 
                                      cut=[op.rng_len(jets)>=2, op.rng_len(bjetsResolved)>=1,op.rng_len(bjetsBoosted)==0],
                                      weight = DeepJetMediumSFApplied)
        hasOsMuElLowMllCutExclusiveResolvedJets = hasOsMuElLowMllCut.refine("hasOsMuElLowMllCutExclusiveResolvedJets", 
                                      cut=[op.rng_len(jets)>=2, op.rng_len(bjetsResolved)>=1,op.rng_len(bjetsBoosted)==0],
                                      weight = DeepJetMediumSFApplied)

        # ExclusiveResolved + OS dilepton plots #
        hasOsMllCutExclusiveResolvedChannelList = [
                             {'channel':'ElEl','sel':hasOsElElLowMllCutOutZResolvedJets,'dilepton':OsElEl[0],'suffix':'hasOsElElLowMllCutOutZResolvedJets'},
                             {'channel':'MuMu','sel':hasOsMuMuLowMllCutOutZResolvedJets,'dilepton':OsMuMu[0],'suffix':'hasOsMuMuLowMllCutOutZResolvedJets'},
                             {'channel':'ElMu','sel':hasOsElMuLowMllCutResolvedJets,    'dilepton':OsElMu[0],'suffix':'hasOsElMuLowMllCutResolvedJets'},
                             {'channel':'MuEl','sel':hasOsMuElLowMllCutResolvedJets,    'dilepton':OsMuEl[0],'suffix':'hasOsMuElLowMllCutResolvedJets'},

                             {'channel':'ElEl','sel':hasOsElElLowMllCutOutZExclusiveResolvedJets,'dilepton':OsElEl[0],'suffix':'hasOsElElLowMllCutOutZExclusiveResolvedJets'},
                             {'channel':'MuMu','sel':hasOsMuMuLowMllCutOutZExclusiveResolvedJets,'dilepton':OsMuMu[0],'suffix':'hasOsMuMuLowMllCutOutZExclusiveResolvedJets'},
                             {'channel':'ElMu','sel':hasOsElMuLowMllCutExclusiveResolvedJets,    'dilepton':OsElMu[0],'suffix':'hasOsElMuLowMllCutExclusiveResolvedJets'},
                             {'channel':'MuEl','sel':hasOsMuElLowMllCutExclusiveResolvedJets,    'dilepton':OsMuEl[0],'suffix':'hasOsMuElLowMllCutExclusiveResolvedJets'},
                           ]

        for channelDict in hasOsMllCutExclusiveResolvedChannelList:
            # Dilepton plots #
            plots.extend(makeDileptonPlots(self, **channelDict))
            # Jets plots #
            plots.extend(makeJetsPlots(self,
                                       sel         = channelDict['sel'],
                                       bjets       = bjetsResolved,
                                       lightjets   = lightjetsResolved,
                                       alljets     = jets,
                                       suffix      = channelDict['suffix'],
                                       channel     = channelDict['channel']))

       ## helper selection (OR) to make sure jet calculations are only done once
       #hasOSLL = noSel.refine("hasOSLL", cut=op.OR(*( hasOSLL_cmbRng(rng) for rng in osLLRng.values())))
       #forceDefine(t._Jet.calcProd, hasOSLL)
       #for varNm in t._Jet.available:
       #    forceDefine(t._Jet[varNm], hasOSLL)
        return plots

