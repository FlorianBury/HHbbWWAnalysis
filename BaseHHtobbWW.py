import os
import re
import sys
import json
import yaml
import copy
from itertools import chain
from functools import partial

import bamboo
from bamboo import treefunctions as op
from bamboo.analysismodules import NanoAODModule, NanoAODHistoModule, NanoAODSkimmerModule
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection
from bamboo.scalefactors import binningVariables_nano, BtagSF, get_scalefactor
from bamboo.plots import Plot, EquidistantBinning, Selection, SelectionWithDataDriven, CutFlowReport
from bamboo.analysisutils import forceDefine
from bamboo.root import loadLibrary, loadHeader

from METScripts import METFilter, METcorrection
from scalefactorsbbWW import ScaleFactorsbbWW
from btagHelper import makeBtagRatioReweighting
from triggers import returnTriggerRanges
from highlevelLambdas import highlevelLambdas
from DDHelper import DataDrivenPseudoData, DataDrivenLOReweighting
from selectionDef import SelectionObject

import logging
logger = logging.getLogger(__name__)

#===============================================================================================#
#                                  BaseHHtobbWW                                                 #
#===============================================================================================#
class BaseNanoHHtobbWW(NanoAODModule):
    """ Base module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
        super(BaseNanoHHtobbWW, self).__init__(args)
        # Set plots options #
        self.plotDefaults = {"show-ratio": True,
                             "y-axis-show-zero" : True,
                             #"normalized": True,
                             "y-axis": "Events",
                             "log-y"  : "both",
                             "ratio-y-axis-range" : [0.8,1.2],
                             "ratio-y-axis" : '#frac{Data}{MC}',
                             "sort-by-yields" : True}


    #-------------------------------------------------------------------------------------------#
    #                                       addArgs                                             #
    #-------------------------------------------------------------------------------------------#
    def addArgs(self,parser):
        super(BaseNanoHHtobbWW, self).addArgs(parser)

        parser.title = """

Arguments for the HH->bbWW analysis on bamboo framework

----- Argument groups -----
    * Lepton arguments *
        --NoZVeto
        --ZPeak
        --NoTauVeto

    * Jet arguments *
        --Ak4
        --Ak8
        --Resolved0Btag
        --Resolved1Btag
        --Resolved2Btag
        --TightResolved0b4j
        --TightResolved1b3j
        --TightResolved2b2j
        --SemiBoostedHbb
        --SemiBoostedWjj
        --TTBarCR
        --BtagReweightingOn
        --BtagReweightingOff

    * Skimmer arguments *
        --Synchronization
        --Channel

    * Plotter arguments *
        --OnlyYield

    * Technical *
        --backend
        --NoSystematics

----- Plotter mode -----
Every combination of lepton and jet arguments can be used, if none specified they will all be plotted

----- Skimmer mode -----
One lepton and and one jet argument must be specified in addition to the required channel
(this is because only one selection at a time can produce a ntuple)

----- Detailed explanation -----

                    """
        #----- Technical arguments -----#
        parser.add_argument("--backend", 
                            type=str, 
                            default="dataframe", 
                            help="Backend to use, 'dataframe' (default) or 'lazy'")
        parser.add_argument("-s",
                            "--subset", 
                            type        = str,
                            required    = False,
                            help="Subset of samples to be run over, keys defined in the 'samples' entry of the yaml analysis config")
        parser.add_argument("--NoSystematics", 
                            action      = "store_true",
                            default     = False,
                            help="Disable all systematic variations (default=False)")
        parser.add_argument("--analysis", 
                            type        = str,
                            default     = 'res',
                            help="Analysis type : res (default) | nonres ")
        parser.add_argument("--Events", 
                            nargs       = '+',
                            type        = int,
                            help="Cut on events (as list)")
        parser.add_argument("--HHReweighting", 
                            action      = 'append',
                            type        = str,
                            default     = None,
                            help="GGF LO->NLO reweighting benchmarks to use (can be several)")
        parser.add_argument("--mass", 
                            nargs       = '+',
                            action      = 'extend',
                            type        = float,
                            default     = None,
                            help="Mass to use for the parametric DNN (can be several)")
        parser.add_argument("--era", 
                            action      = 'store',
                            type        = int,
                            default     = None,
                            help="Era to be fed to the parametric DNN")
        parser.add_argument("--PrintYield", 
                            action      = "store_true",
                            default     = False,
                            help="Print yield to screen (for debugging)")



        #----- Lepton selection arguments -----#
        parser.add_argument("--NoZVeto", 
                            action      = "store_true",
                            default     = False,
                            help        = "Remove the cut of preselected leptons |M_ll-M_Z|>10 GeV")
        parser.add_argument("--ZPeak", 
                            action      = "store_true",
                            default     = False,
                            help        = "Select the Z peak at tight level |M_ll-M_Z|<10 GeV (must be used with --NoZVeto, only effective with --Tight)")
        parser.add_argument("--NoTauVeto", 
                            action      = "store_true",
                            default     = False,
                            help        = "Select the events do not have any tau overlapped with fakeable leptons")
        parser.add_argument("--FakeRateNonClosureMCFakes", 
                            action      = "store_true",
                            default     = False,
                            help        = "Use tight non-prompt lepton to estimate MC fakes [Only on TTHardronic events]")
        parser.add_argument("--FakeRateNonClosureMCClosure", 
                            action      = "store_true",
                            default     = False,
                            help        = "Use tight non-prompt lepton to estimate MC closure [Only on TTHardronic events]")



        #----- Jet selection arguments -----#
        parser.add_argument("--Ak4", 
                                action      = "store_true",
                                default     = False,
                                help        = "Produce the plots/skim for two Ak4 jets passing the selection criteria")
        parser.add_argument("--Ak8", 
                                action      = "store_true",
                                default     = False,
                                help        = "Produce the plots/skim for one Ak8 jet passing the selection criteria")
        parser.add_argument("--Resolved0Btag", 
                                action      = "store_true",
                                default     = False,
                                help        = "Produce the plots/skim for the exclusive resolved category with no btagged jet")
        parser.add_argument("--Resolved1Btag", 
                                action      = "store_true",
                                default     = False,
                                help        = "Produce the plots/skim for the exclusive resolved category with only one btagged jet")
        parser.add_argument("--Resolved2Btag", 
                                action      = "store_true",
                                default     = False,
                                help        = "Produce the plots/skim for the exclusive resolved category with two btagged jets")

        # SL Categories
        # Resolved
        parser.add_argument("--Res2b2Wj", 
                                action      = "store_true",
                                default     = False,
                                help        = "Produce the plots/skim for the exclusive 2b2Wj JPA category")
        parser.add_argument("--Res1b3Wj", 
                                action      = "store_true",
                                default     = False,
                                help        = "Produce the plots/skim for the exclusive 1b3Wj JPA category") # for basic reco only
        parser.add_argument("--Res2b1Wj", 
                                action      = "store_true",
                                default     = False,
                                help        = "Produce the plots/skim for the exclusive 2b1Wj JPA category")
        parser.add_argument("--Res2b0Wj", 
                                action      = "store_true",
                                default     = False,
                                help        = "Produce the plots/skim for the exclusive 2b0Wj JPA category")
        parser.add_argument("--Res1b2Wj", 
                                action      = "store_true",
                                default     = False,
                                help        = "Produce the plots/skim for the exclusive 1b2Wj JPA category")
        parser.add_argument("--Res1b1Wj", 
                                action      = "store_true",
                                default     = False,
                                help        = "Produce the plots/skim for the exclusive 1b1Wj JPA category")
        parser.add_argument("--Res1b0Wj", 
                                action      = "store_true",
                                default     = False,
                                help        = "Produce the plots/skim for the exclusive 1b0Wj JPA category")
        parser.add_argument("--Res0b", 
                                action      = "store_true",
                                default     = False,
                                help        = "Produce the plots/skim for the exclusive 0b JPA category")

        # Semi-Boosted
        parser.add_argument("--Hbb2Wj", 
                                action      = "store_true",
                                default     = False,
                                help        = "Produce the plots/skim for the exclusive Hbb2Wj JPA category")
        parser.add_argument("--Hbb1Wj", 
                                action      = "store_true",
                                default     = False,
                                help        = "Produce the plots/skim for the exclusive Hbb1Wj JPA category")
        parser.add_argument("--Hbb0Wj", 
                                action      = "store_true",
                                default     = False,
                                help        = "Produce the plots/skim for the exclusive Hbb0Wj JPA category")

        parser.add_argument("--Resolved", 
                                action      = "store_true",
                                default     = False,
                                help        = "Produce the plots/skim for all JPA resolved categories")
        parser.add_argument("--Boosted", 
                                action      = "store_true",
                                default     = False,
                                help        = "Produce the plots/skim for all JPA boosted categories")
        ################

        parser.add_argument("--Boosted0Btag", 
                                action      = "store_true",
                                default     = False,
                                help        = "Produce the plots/skim for the inclusive boosted category")
        parser.add_argument("--Boosted1Btag", 
                                action      = "store_true",
                                default     = False,
                                help        = "Produce the plots/skim for the inclusive boosted category")
        parser.add_argument("--TTBarCR", 
                                action      = "store_true",
                                default     = False,
                                help        = "Apply cut on Mbb for ttbar CR (only effective with --Resolved2Btag)")
        parser.add_argument("--BtagReweightingOn", 
                                action      = "store_true",
                                default     = False,
                                help        = "Btag ratio study : Btag SF applied (without the ratio), will only do the plots for reweighting (jets and leptons args are ignored)")
        parser.add_argument("--BtagReweightingOff", 
                                action      = "store_true",
                                default     = False,
                                help        = "Btag ratio study : Btag Sf not applied (without the ratio), will only do the plots for reweighting (jets and leptons args are ignored)")
        parser.add_argument("--DYStitchingPlots", 
                                action      = "store_true",
                                default     = False,
                                help        = "DY stitching studies : only produce LHE jet multiplicities (inclusive analysis, only on DY events, rest of plots ignored)")
        parser.add_argument("--WJetsStitchingPlots", 
                                action      = "store_true",
                                default     = False,
                                help        = "W+jets stitching studies : only produce LHE jet multiplicities (inclusive analysis, only on W+jets events, rest of plots ignored)")
        parser.add_argument("--NoStitching", 
                                action      = "store_true",
                                default     = False,
                                help        = "To not apply the stitching weights to DY and WJets samples")

        #----- Skimmer arguments -----#
        parser.add_argument("--Synchronization", 
                            action      = "store_true",
                            default     = False,
                            help        = "Produce the skims for the synchronization (without triggers, corrections of flags) if alone. If sync for specific selection, lepton, jet and channel arguments need to be used")
        parser.add_argument("--Channel", 
                            action      = "store",
                            type        = str,
                            help        = "Specify the channel for the tuple : ElEl, MuMu, ElMu")
        parser.add_argument("--FakeCR", 
                            action      = "store_true",
                            default     = False,
                            help        = "Use the Fake CR instead of the SR")
        parser.add_argument("--DYCR", 
                            action      = "store_true",
                            default     = False,
                            help        = "Use the DY CR instead of the SR")


        #----- Plotter arguments -----#
        parser.add_argument("--OnlyYield", 
                            action      = "store_true",
                            default     = False,
                            help        = "Only produce the yield plots")
        parser.add_argument("--Classifier", 
                            action      = "store",
                            type        = str,
                            help        = "BDT-SM | BDT-Rad900 | DNN | LBN")
        parser.add_argument("--WhadTagger", 
                            action      = "store",
                            type        = str,
                            help        = "BDT | simple")


    #-------------------------------------------------------------------------------------------#
    #                                      initialize                                           #
    #-------------------------------------------------------------------------------------------#
    def initialize(self,forSkimmer=False):
        # Include all the contributions from the subsets in the yaml #
        self._customizeAnalysisCfg(self.analysisConfig)

        # Add the LO reweighting to the datadriven parts #
        if self.args.HHReweighting is not None and 'all' in self.args.HHReweighting:
            self.args.HHReweighting = self.analysisConfig['benchmarks']['targets']
        if not forSkimmer and self.args.HHReweighting is not None:
            self.analysisConfig['datadriven'].update({benchmark:{'replaces': 'all', 'uses': 'all'} 
                                                        for benchmark in self.analysisConfig['benchmarks']['targets'] 
                                                            if benchmark in self.args.HHReweighting})
            newContrib = [benchmark for benchmark in self.analysisConfig['benchmarks']['targets'] if benchmark in self.args.HHReweighting]
            if self.args.datadriven is None:
                self.args.datadriven = newContrib
            else:
                self.args.datadriven += newContrib

        super(BaseNanoHHtobbWW, self).initialize()

        # PseudoData #
        if not forSkimmer:
            if "PseudoData" in self.datadrivenContributions:
                contrib = self.datadrivenContributions["PseudoData"]
                self.datadrivenContributions["PseudoData"] = DataDrivenPseudoData(contrib.name, contrib.config)

        # Include the datadriven reweighting #
        if self.args.HHReweighting is not None:
            self.datadrivenContributions.update({benchmark:DataDrivenLOReweighting(benchmark,self.analysisConfig['datadriven'][benchmark],substrs=self.analysisConfig['benchmarks']['uses']) 
                                                        for benchmark in self.analysisConfig['benchmarks']['targets'] 
                                                            if benchmark in self.args.HHReweighting}) 

        # Check analysis type #
        if self.args.analysis not in ['res','nonres']:
            raise RuntimeError("Analysis type '{}' not understood".format(self.args.analysis))


    #-------------------------------------------------------------------------------------------#
    #                                   customizeAnalysisCfg                                    #
    #-------------------------------------------------------------------------------------------#
    def _customizeAnalysisCfg(self,analysisCfg):
        samples = {} 
        if self.args.subset is None:
            return
        reqArgs = self.args.subset.split(',')
        foundArgs = set()
        subsets = []
        for item in analysisCfg['samples']:
            if not 'keys' in item.keys() or not 'config' in item.keys():
                continue
            keys = item['keys']
            if not isinstance(keys,list):
                keys = [keys]
            if all([key in reqArgs for key in keys]) or 'all' in reqArgs:
                foundArgs.update(keys)
                subsets.append(item['config'])
                with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'Yaml',item['config'])) as handle:
                    subsetDict = yaml.load(handle,yaml.SafeLoader)
                for sampleName,sampleCfg in subsetDict.items():
                    if self.args.analysis == 'res' and 'type' in sampleCfg.keys() and sampleCfg['type'] == 'signal':
                        if 'mass' not in sampleCfg.keys():
                            raise RuntimeError(f'Yaml config {item["config"]} sample {sampleName} does not have `mass` entry')
                        mass = float(re.findall('M-\d+',sampleName)[0].replace('M-',''))
                        if float(mass) != float(sampleCfg['mass']):
                            raise RuntimeError(f'Yaml config {item["config"]} sample {sampleName} have `mass` entry {sampleCfg["mass"]} but name tells {mass}')
                        if self.args.mass is not None and mass not in self.args.mass:
                            continue
                    samples[sampleName] = sampleCfg
        self.analysisConfig['samples'] = samples
        notFoundArgs = [arg for arg in reqArgs if arg not in foundArgs]
        if len(notFoundArgs)>0:
            raise RuntimeError('The following subsets have not been found in the keys of the analysis yaml file : '+', '.join(notFoundArgs))
        if len(subsets) > 0:
            print ("Imported following yaml subsets :")
            for subset in subsets:
                print ('... {}'.format(subset))

    #-------------------------------------------------------------------------------------------#
    #                                        counters                                           #
    #-------------------------------------------------------------------------------------------#
    def readCounters(self, resultsFile):
        counters = super(BaseNanoHHtobbWW, self).readCounters(resultsFile)
        # Corrections to the generated sum "
        if resultsFile.GetListOfKeys().FindObject('generated_sum_corrected'):
            sample = os.path.basename(resultsFile.GetName())
            print (f'Sample {sample} : genEventSumw correction from {counters["genEventSumw"]:.3f} to {resultsFile.Get("generated_sum_corrected").GetBinContent(1):.3f}')
            counters["genEventSumw"] = resultsFile.Get('generated_sum_corrected').GetBinContent(1)
        return counters

    def mergeCounters(self, outF, infileNames, sample=None):
        super(BaseNanoHHtobbWW, self).mergeCounters(outF, infileNames, sample)
        if outF.GetListOfKeys().FindObject('generated_sum_corrected'): # Main file 
            self.generated_sum_corrected = copy.deepcopy(outF.Get('generated_sum_corrected'))
        else: # All the additional files ("datadriven")
            if hasattr(self,'generated_sum_corrected'): 
                self.generated_sum_corrected.Write()

    #-------------------------------------------------------------------------------------------#
    #                                       prepareTree                                         #
    #-------------------------------------------------------------------------------------------#
    def prepareTree(self, tree, sample=None, sampleCfg=None):
        from bamboo.treedecorators import NanoAODDescription, nanoRochesterCalc, nanoJetMETCalc, nanoJetMETCalc_METFixEE2017, nanoFatJetCalc
        # JEC's Recommendation for Full RunII: https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
        # JER : -----------------------------: https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution

        # Get base aguments #
        era = sampleCfg['era']
        self.is_MC = self.isMC(sample) # no confusion between boolean is_MC and method isMC()
        metName = "METFixEE2017" if era == "2017" else "MET"
        tree,noSel,backend,lumiArgs = super(BaseNanoHHtobbWW,self).prepareTree(tree          = tree, 
                                                                               sample        = sample, 
                                                                               sampleCfg     = sampleCfg, 
                                                                               description   = NanoAODDescription.get(
                                                                                                 tag             = "v7", 
                                                                                                 year            = (era if era else "2016"),
                                                                                                 isMC            = self.is_MC,
                                                                                                 systVariations  = [ (nanoJetMETCalc_METFixEE2017 if era == "2017" else nanoJetMETCalc), nanoFatJetCalc]),
                                                                                                 # will do Jet and MET variations, and not the Rochester correction
                                                                               backend       = ("lazy" if self.args.onlypost else self.args.backend))
    
        # Plots in base that need to be propagated to the Plotters #
        self.base_plots = []


        #----- Helper classes declaration ----#
        if self.args.analysis == "res":
            loadLibrary(os.path.join(os.path.dirname(os.path.abspath(__file__)),'HMEStudy','build','libBambooHMEEvaluator'))
            loadHeader(os.path.join(os.path.dirname(os.path.abspath(__file__)),'HMEStudy','include','HME.h'))
            loadHeader(os.path.join(os.path.dirname(os.path.abspath(__file__)),'HMEStudy','include','Reader.h'))
            self.hmeEval = op.define("hme::HMEEvaluator", "hme::HMEEvaluator <<name>>{{}}; // for {sample}".format(sample=sample.replace('-','')), nameHint="bamboo_hmeEval{sample}".format(sample=sample.replace('-','')))

        #----- CutFlow report -----#
        self.yields = CutFlowReport("yields",printInLog=self.args.PrintYield,recursive=self.args.PrintYield)


        #----- Safeguards for signals -----#
        if "HH" in sample:
            noSel = noSel.refine('HHMCWeight',cut=[op.abs(tree.genWeight)<100])
            if self.args.PrintYield:
                self.yields.add(noSel)
            # Correct the gen event weight sum #
            self.base_plots.append(Plot.make1D("generated_sum_corrected",
                                               op.c_float(0.5),
                                               noSel,
                                               EquidistantBinning(1,0.,1.),
                                               weight=tree.genWeight, 
                                               autoSyst=False))

        #----- Genweight -----#
        if self.is_MC:
            noSel = noSel.refine("genWeight", weight=tree.genWeight)
            if self.args.PrintYield:
                self.yields.add(noSel)

        # Event cut #
        if self.args.Events:
            print ("Events to use only :")
            for e in self.args.Events:
                print ('... %d'%e)
            noSel = noSel.refine('eventcut',cut = [op.OR(*[tree.event == e for e in self.args.Events])])
            if self.args.PrintYield:
                self.yields.add(noSel)

        # Save some useful stuff in self #
        self.sample = sample
        self.sampleCfg = sampleCfg
        self.era = era

        # Check distributed option #
        isNotWorker = (self.args.distributed != "worker") 

        # Turn off systs #
        if self.args.NoSystematics:
            noSel = noSel.refine('SystOff',autoSyst=False)
            if self.args.PrintYield:
                self.yields.add(noSel)

        # Check era #
        if era != "2016" and era != "2017" and era != "2018":
            raise RuntimeError("Unknown era {0}".format(era))

        # Rochester and JEC corrections (depends on era) #     

        # Check if basic synchronization is required (no corrections and triggers) #
        self.inclusive_sel = ((self.args.Synchronization 
                               and not any([self.args.__dict__[key] for key in['Ak4', 'Ak8', 'Resolved0Btag', 'Resolved1Btag', 'Resolved2Btag', 'Boosted0Btag','Boosted1Btag',
                                                                               'Resolved','Boosted','Res2b2Wj','Res2b1Wj','Res2b0Wj','Res1b2Wj','Res1b1Wj','Res1b1Wj','Res0b',
                                                                               'Hbb2Wj','Hbb1Wj','Hbb0Wj']]) \
                               and (self.args.Channel is None or self.args.Channel=='None')) \
                              # No channel selection
                              # None is local mode, "None" in distributed mode
                              or self.args.BtagReweightingOn 
                              or self.args.BtagReweightingOff
                              or self.args.DYStitchingPlots
                              or self.args.WJetsStitchingPlots)
        # Inclusive plots
        # If no lepton, jet and channel selection : basic object selection (no trigger nor corrections)
        if self.inclusive_sel:
            print ("Inclusive analysis, no selections applied")

        #----- Theory uncertainties -----#
        if self.is_MC:
            # PDF scale weights #
            # Twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopSystematics#Factorization_and_renormalizatio 
            # len(LHEScaleWeight) == 9
            #   -> nominal = LHEScaleWeight[4]
            #   -> Fact     : up = LHEScaleWeight[5] and down = LHEScaleWeight[3]
            #   -> Renorm   : up = LHEScaleWeight[7] and down = LHEScaleWeight[1]
            #   -> Mixed    : up = LHEScaleWeight[8] and down = LHEScaleWeight[0]
            # len(LHEScaleWeight) == 8
            #   -> nominal = 1.
            #   -> Fact     : up = LHEScaleWeight[4] and down = LHEScaleWeight[3]
            #   -> Renorm   : up = LHEScaleWeight[6] and down = LHEScaleWeight[1]
            #   -> Mixed    : up = LHEScaleWeight[7] and down = LHEScaleWeight[0]
            # len(LHEScaleWeight) different (usually 44 ?!)
            #   -> nominal = up = down = 1.
            # Meaning :
            #        [' LHE scale variation weights (w_var / w_nominal)',         
            #         ' [0] is renscfact=0.5d0 facscfact=0.5d0 ',         
            #         ' [1] is renscfact=0.5d0 facscfact=1d0 ',         
            #         ' [2] is renscfact=0.5d0 facscfact=2d0 ',         
            #         ' [3] is renscfact=1d0 facscfact=0.5d0 ',         
            #         ' [4] is renscfact=1d0 facscfact=1d0 ',         
            #         ' [5] is renscfact=1d0 facscfact=2d0 ',         
            #         ' [6] is renscfact=2d0 facscfact=0.5d0 ',         
            #         ' [7] is renscfact=2d0 facscfact=1d0 ',         
            #         ' [8] is renscfact=2d0 facscfact=2d0 ']    
            # Clipping is done to avoid malicious files in ST samples
            basicScaleWeight = op.systematic(op.c_float(1.),
                                             name         = "ScaleWeight",
                                             #Factup       = op.c_float(1.),
                                             #Factdown     = op.c_float(1.),
                                             #Renormup     = op.c_float(1.),
                                             #Renormdown   = op.c_float(1.),
                                             #Mixedup      = op.c_float(1.),
                                             #Mixeddown    = op.c_float(1.),
                                             Envelopeup   = op.c_float(1.),
                                             Envelopedown = op.c_float(1.))
            if sample.startswith('ST'): # Dropped because bugs in the LHE scale weights
                self.scaleWeight = basicScaleWeight
            elif hasattr(tree,"LHEScaleWeight"): # If has tree -> find the values
                factor = 1.
                if sample.startswith('DYToLL_0J') or sample.startswith('DYToLL_1J'):
                    # Bug of factor 0.5, see https://hypernews.cern.ch/HyperNews/CMS/get/generators/4383.html?inline=-1 (only in the "8" weights case)
                    factor = 2.
                self.scaleWeight = op.multiSwitch((op.AND(op.rng_len(tree.LHEScaleWeight) == 9, tree.LHEScaleWeight[4] != 0.),
                                                                        op.systematic(op.c_float(1.),    #tree.LHEScaleWeight[4],
                                                                                      name          = "ScaleWeight",
                                                                                      #Factup        = op.min(op.c_float(10.),op.max(op.c_float(0.),tree.LHEScaleWeight[5]/tree.LHEScaleWeight[4])),
                                                                                      #Factdown      = op.min(op.c_float(10.),op.max(op.c_float(0.),tree.LHEScaleWeight[3]/tree.LHEScaleWeight[4])),
                                                                                      #Renormup      = op.min(op.c_float(10.),op.max(op.c_float(0.),tree.LHEScaleWeight[7]/tree.LHEScaleWeight[4])),
                                                                                      #Renormdown    = op.min(op.c_float(10.),op.max(op.c_float(0.),tree.LHEScaleWeight[1]/tree.LHEScaleWeight[4])),
                                                                                      #Mixedup       = op.min(op.c_float(10.),op.max(op.c_float(0.),tree.LHEScaleWeight[8]/tree.LHEScaleWeight[4])),
                                                                                      #Mixeddown     = op.min(op.c_float(10.),op.max(op.c_float(0.),tree.LHEScaleWeight[0]/tree.LHEScaleWeight[4])),
                                                                                      Envelopeup    = op.min(op.c_float(10.),
                                                                                                             op.max(op.max(op.max(tree.LHEScaleWeight[5]/tree.LHEScaleWeight[4], # Fact up 
                                                                                                                                  tree.LHEScaleWeight[7]/tree.LHEScaleWeight[4]), # Renorm up
                                                                                                                           tree.LHEScaleWeight[8]/tree.LHEScaleWeight[4]), # Mixed up
                                                                                                                    op.c_float(0.))),
                                                                                      Envelopedown  = op.max(op.c_float(0.),
                                                                                                             op.min(op.min(op.min(tree.LHEScaleWeight[3]/tree.LHEScaleWeight[4], # Fact down
                                                                                                                                  tree.LHEScaleWeight[1]/tree.LHEScaleWeight[4]), # Renorm down
                                                                                                                           tree.LHEScaleWeight[0]/tree.LHEScaleWeight[4]), # Mixed own 
                                                                                                                    op.c_float(10.))))),
                                                  (op.rng_len(tree.LHEScaleWeight) == 8, 
                                                              op.systematic(op.c_float(1.),
                                                                            name         = "ScaleWeight",
                                                                            #Factup       = factor * tree.LHEScaleWeight[4],
                                                                            #Factdown     = factor * tree.LHEScaleWeight[3],
                                                                            #Renormup     = factor * tree.LHEScaleWeight[6],
                                                                            #Renormdown   = factor * tree.LHEScaleWeight[1],
                                                                            #Mixedup      = factor * tree.LHEScaleWeight[7],
                                                                            #Mixeddown    = factor * tree.LHEScaleWeight[0],
                                                                            Envelopeup   = op.max(op.max(factor * tree.LHEScaleWeight[4],    # Fact up
                                                                                                         factor * tree.LHEScaleWeight[6]),    # Renorm up
                                                                                                  factor * tree.LHEScaleWeight[7]),   # Mixed up
                                                                            Envelopedown = op.min(op.min(factor * tree.LHEScaleWeight[3],    # Fact down
                                                                                                         factor * tree.LHEScaleWeight[1]),    # Renorm down
                                                                                                  factor * tree.LHEScaleWeight[0]))), # Mixed down
                                                  basicScaleWeight)

                                                                                                 
            else: # No tree -> use 1
                self.scaleWeight = basicScaleWeight
                                          
            noSel = noSel.refine("PDFScaleWeights", weight = [self.scaleWeight])
            if self.args.PrintYield:
                self.yields.add(noSel)

            # PDF #
            #noSel = noSel.refine("PDFWeights", weight = [tree.LHEPdfWeight])
            # TODO : add that next time 
            # https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/mc102X_doca.html#LHEPdfWeight
            # https://lhapdf.hepforge.org/pdfsets.html
            # -> range of ID
            # Branch LHEPdfWeight -> docstring -> get the id -> check the list above
            # we want "NNPDF31_nnlo_hessian_pdfas" (ID : 306000) -> if not use 1.
            

            # PS weights #
            self.psISRSyst = op.switch(op.rng_len(tree.PSWeight) == 4,
                                       op.systematic(op.c_float(1.), name="psISR", up=tree.PSWeight[2], down=tree.PSWeight[0]),
                                       op.systematic(op.c_float(1.), name="psISR", up=op.c_float(1.), down=op.c_float(1.)))
            self.psFSRSyst = op.switch(op.rng_len(tree.PSWeight) == 4,
                                       op.systematic(op.c_float(1.), name="psFSR", up=tree.PSWeight[3], down=tree.PSWeight[1]),
                                       op.systematic(op.c_float(1.), name="psFSR", up=op.c_float(1.), down=op.c_float(1.)))
            noSel = noSel.refine("PSweights", weight = [self.psISRSyst, self.psFSRSyst])
            if self.args.PrintYield:
                self.yields.add(noSel)

        #----- Triggers and Corrections -----#
        self.triggersPerPrimaryDataset = {}
        def addHLTPath(key,HLT):
            if key not in self.triggersPerPrimaryDataset.keys():
                self.triggersPerPrimaryDataset[key] = []
            try:
                self.triggersPerPrimaryDataset[key].append(getattr(tree.HLT,HLT))
            except AttributeError:
                print ("Could not find branch tree.HLT.%s, will omit it"%HLT)

        from bamboo.analysisutils import configureJets ,configureRochesterCorrection, configureType1MET 
        ############################################################################################
        # ERA 2016 #
        ############################################################################################
        if era == "2016":
#            # Rochester corrections #
#            configureRochesterCorrection(variProxy  = tree._Muon,
#                                         paramsFile = os.path.join(os.path.dirname(__file__), "data", "RoccoR2016.txt"),
#                                         isMC       = self.is_MC,
#                                         backend    = backend, 
#                                         uName      = sample)
 
            # SingleMuon #
            addHLTPath("SingleMuon","IsoMu22")
            addHLTPath("SingleMuon","IsoTkMu22")
            addHLTPath("SingleMuon","IsoMu22_eta2p1")
            addHLTPath("SingleMuon","IsoTkMu22_eta2p1")
            addHLTPath("SingleMuon","IsoMu24")
            addHLTPath("SingleMuon","IsoTkMu24")
            # SingleElectron #
            addHLTPath("SingleElectron","Ele27_WPTight_Gsf")
            addHLTPath("SingleElectron","Ele25_eta2p1_WPTight_Gsf")
            addHLTPath("SingleElectron","Ele27_eta2p1_WPLoose_Gsf")
            # DoubleMuon #
            addHLTPath("DoubleMuon","Mu17_TrkIsoVVL_Mu8_TrkIsoVVL")
            addHLTPath("DoubleMuon","Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ")
            addHLTPath("DoubleMuon","Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL")
            addHLTPath("DoubleMuon","Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ")
            # DoubleEGamma #
            addHLTPath("DoubleEGamma","Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ")
            # MuonEG #
            addHLTPath("MuonEG","Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL")
            addHLTPath("MuonEG","Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")
            addHLTPath("MuonEG","Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL")
            addHLTPath("MuonEG","Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ")

            # Links : 
            # JEC : https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC
            # JER (smear) : https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
            # JetMET treatment #
            if self.is_MC:   # if MC -> needs smearing
                configureJets(variProxy             = tree._Jet, 
                              jetType               = "AK4PFchs",
                              jec                   = "Summer16_07Aug2017_V11_MC",
                              smear                 = "Summer16_25nsV1_MC",
                              jesUncertaintySources = "Merged",
                              regroupTag            = "V2",
                              enableSystematics     = lambda v : not "jesTotal" in v,
                              mayWriteCache         = isNotWorker,
                              isMC                  = self.is_MC,
                              backend               = backend, 
                              uName                 = sample)
                configureJets(variProxy             = tree._FatJet, 
                              jetType               = "AK8PFPuppi", 
                              jec                   = "Summer16_07Aug2017_V11_MC", 
                              smear                 = "Summer16_25nsV1_MC", 
                              jesUncertaintySources = "Merged",
                              uncertaintiesFallbackJetType = "AK4PFchs",
                              mcYearForFatJets      = era if self.args.analysis == 'res' else None, 
                              regroupTag            = "V2",
                              enableSystematics     = lambda v : not "jesTotal" in v,
                              mayWriteCache         = isNotWorker, 
                              isMC                  = self.is_MC,
                              backend               = backend, 
                              uName                 = sample)
                configureType1MET(variProxy             = getattr(tree, f"_{metName}"),
                                  jec                   = "Summer16_07Aug2017_V11_MC",
                                  smear                 = "Summer16_25nsV1_MC",
                                  isT1Smear             = True,
                                  jesUncertaintySources = "Merged",
                                  regroupTag            = "V2",
                                  enableSystematics     = lambda v : not "jesTotal" in v,
                                  mayWriteCache         = isNotWorker,
                                  isMC                  = self.is_MC,
                                  backend               = backend,
                                  uName                 = sample)

            else:                   # If data -> extract info from config 
                jecTag = None
                if "2016B" in sample or "2016C" in sample or "2016D" in sample:
                    jecTag = "Summer16_07Aug2017BCD_V11_DATA"
                elif "2016E" in sample or "2016F" in sample:
                    jecTag = "Summer16_07Aug2017EF_V11_DATA"
                elif "2016G" in sample or "2016H" in sample:
                    jecTag = "Summer16_07Aug2017GH_V11_DATA"
                else:
                    raise RuntimeError("Could not find appropriate JEC tag for data")
                configureJets(variProxy             = tree._Jet, 
                              jetType               = "AK4PFchs",
                              jec                   = jecTag,
                              mayWriteCache         = isNotWorker,
                              isMC                  = self.is_MC,
                              backend               = backend, 
                              uName                 = sample)
                configureJets(variProxy             = tree._FatJet, 
                              jetType               = "AK8PFPuppi", 
                              jec                   = jecTag,
                              regroupTag            = "V2",
                              mayWriteCache         = isNotWorker, 
                              isMC                  = self.is_MC,
                              backend               = backend, 
                              uName                 = sample)
                configureType1MET(variProxy         = getattr(tree, f"_{metName}"),
                                  jec               = jecTag,
                                  mayWriteCache     = isNotWorker,
                                  isMC              = self.is_MC,
                                  backend           = backend, 
                                  uName             = sample)

        ############################################################################################
        # ERA 2017 #
        ############################################################################################
        elif era == "2017":
            #configureRochesterCorrection(variProxy  = tree._Muon,
            #                             paramsFile = os.path.join(os.path.dirname(__file__), "data", "RoccoR2017.txt"),
            #                             isMC       = self.is_MC,
            #                             backend    = backend, 
            #                             uName      = sample)

            # SingleMuon #
            addHLTPath("SingleMuon","IsoMu24")
            addHLTPath("SingleMuon","IsoMu27")
            # SingleElectron #
            addHLTPath("SingleElectron","Ele35_WPTight_Gsf")
            addHLTPath("SingleElectron","Ele32_WPTight_Gsf")
            # DoubleMuon #
            addHLTPath("DoubleMuon","Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8")
            addHLTPath("DoubleMuon","Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8")
            # DoubleEGamma #
            addHLTPath("DoubleEGamma","Ele23_Ele12_CaloIdL_TrackIdL_IsoVL")
            # MuonEG #
            addHLTPath("MuonEG","Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")
            addHLTPath("MuonEG","Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")
            addHLTPath("MuonEG","Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL")
            addHLTPath("MuonEG","Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ")
            
            # JetMET treatment #
            if self.is_MC:   # if MC -> needs smearing
                configureJets(variProxy             = tree._Jet, 
                              jetType               = "AK4PFchs",
                              jec                   = "Fall17_17Nov2017_V32_MC",
                              smear                 = "Fall17_V3b_MC",
                              jesUncertaintySources = "Merged",
                              regroupTag            = "V2",
                              enableSystematics     = lambda v : not "jesTotal" in v,
                              mayWriteCache         = isNotWorker,
                              isMC                  = self.is_MC,
                              backend               = backend, 
                              uName                 = sample)
                configureJets(variProxy             = tree._FatJet, 
                              jetType               = "AK8PFPuppi", 
                              jec                   = "Fall17_17Nov2017_V32_MC",
                              smear                 = "Fall17_V3b_MC",
                              jesUncertaintySources = "Merged",
                              regroupTag            = "V2",
                              enableSystematics     = lambda v : not "jesTotal" in v,
                              uncertaintiesFallbackJetType = "AK4PFchs",
                              mcYearForFatJets      = era if self.args.analysis == 'res' else None, 
                              mayWriteCache         = isNotWorker, 
                              isMC                  = self.is_MC,
                              backend               = backend, 
                              uName                 = sample)
                configureType1MET(variProxy             = getattr(tree, f"_{metName}"),
                                  jec                   = "Fall17_17Nov2017_V32_MC",
                                  smear                 = "Fall17_V3b_MC",
                                  isT1Smear             = True,
                                  jesUncertaintySources = "Merged",
                                  regroupTag            = "V2",
                                  enableSystematics     = lambda v : not "jesTotal" in v,
                                  mayWriteCache         = isNotWorker,
                                  isMC                  = self.is_MC,
                                  backend               = backend,
                                  uName                 = sample)

            else:                   # If data -> extract info from config 
                jecTag = None
                if "2017B" in sample:
                    jecTag = "Fall17_17Nov2017B_V32_DATA"
                elif "2017C" in sample:
                    jecTag = "Fall17_17Nov2017C_V32_DATA"
                elif "2017D" in sample or "2017E" in sample:
                    jecTag = "Fall17_17Nov2017DE_V32_DATA"
                elif "2017F" in sample:
                    jecTag = "Fall17_17Nov2017F_V32_DATA"
                else:
                    raise RuntimeError("Could not find appropriate JEC tag for data")
                configureJets(variProxy             = tree._Jet, 
                              jetType               = "AK4PFchs",
                              jec                   = jecTag,
                              mayWriteCache         = isNotWorker,
                              isMC                  = self.is_MC,
                              backend               = backend, 
                              uName                 = sample)
                configureJets(variProxy             = tree._FatJet, 
                              jetType               = "AK8PFPuppi", 
                              jec                   = jecTag,
                              mayWriteCache         = isNotWorker, 
                              isMC                  = self.is_MC,
                              backend               = backend, 
                              uName                 = sample)
                configureType1MET(variProxy         = getattr(tree, f"_{metName}"),
                                  jec               = jecTag,
                                  mayWriteCache     = isNotWorker,
                                  isMC              = self.is_MC,
                                  backend           = backend, 
                                  uName             = sample)

        ############################################################################################
        # ERA 2018 #
        ############################################################################################
        elif era == "2018":
            #configureRochesterCorrection(variProxy  = tree._Muon,
            #                             paramsFile = os.path.join(os.path.dirname(__file__), "data", "RoccoR2018.txt"),
            #                             isMC       = self.is_MC,
            #                             backend    = backend, 
            #                             uName      = sample)

            # SingleMuon #
            addHLTPath("SingleMuon","IsoMu24")
            addHLTPath("SingleMuon","IsoMu27")
            # SingleElectron #
            addHLTPath("SingleElectron","Ele32_WPTight_Gsf")
            # DoubleMuon #
            addHLTPath("DoubleMuon","Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8")
            # DoubleEGamma #
            addHLTPath("DoubleEGamma","Ele23_Ele12_CaloIdL_TrackIdL_IsoVL")
            # MuonEG #
            addHLTPath("MuonEG","Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")
            addHLTPath("MuonEG","Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")
            addHLTPath("MuonEG","Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL")

            # JetMET treatment #
            if self.is_MC:   # if MC -> needs smearing
                configureJets(variProxy             = tree._Jet, 
                              jetType               = "AK4PFchs",
                              jec                   = "Autumn18_V19_MC",
                              smear                 = "Autumn18_V7b_MC",
                              jesUncertaintySources = "Merged",
                              regroupTag            = "V2",
                              enableSystematics     = lambda v : not "jesTotal" in v,
                              addHEM2018Issue       = True,
                              mayWriteCache         = isNotWorker,
                              isMC                  = self.is_MC,
                              backend               = backend, 
                              uName                 = sample)
                configureJets(variProxy             = tree._FatJet, 
                              jetType               = "AK8PFPuppi", 
                              jec                   = "Autumn18_V19_MC",
                              smear                 = "Autumn18_V7b_MC",
                              jesUncertaintySources = "Merged",
                              regroupTag            = "V2",
                              enableSystematics     = lambda v : not "jesTotal" in v,
                              addHEM2018Issue       = True,
                              uncertaintiesFallbackJetType = "AK4PFchs",
                              mayWriteCache         = isNotWorker, 
                              isMC                  = self.is_MC,
                              backend               = backend, 
                              uName                 = sample)
                configureType1MET(variProxy             = getattr(tree, f"_{metName}"),
                                  jec                   = "Autumn18_V19_MC",
                                  smear                 = "Autumn18_V7b_MC",
                                  isT1Smear             = True,
                                  jesUncertaintySources = "Merged",
                                  regroupTag            = "V2",
                                  enableSystematics     = lambda v : not "jesTotal" in v,
                                  addHEM2018Issue       = True,
                                  mayWriteCache         = isNotWorker,
                                  isMC                  = self.is_MC,
                                  backend               = backend,
                                  uName                 = sample)

            else:                   # If data -> extract info from config 
                jecTag = None
                if "2018A" in sample:
                    jecTag = "Autumn18_RunA_V19_DATA"
                elif "2018B" in sample:
                    jecTag = "Autumn18_RunB_V19_DATA"
                elif "2018C" in sample:
                    jecTag = "Autumn18_RunC_V19_DATA"
                elif "2018D" in sample:
                    jecTag = "Autumn18_RunD_V19_DATA"
                else:
                    raise RuntimeError("Could not find appropriate JEC tag for data")
                configureJets(variProxy             = tree._Jet, 
                              jetType               = "AK4PFchs",
                              jec                   = jecTag,
                              mayWriteCache         = isNotWorker,
                              isMC                  = self.is_MC,
                              backend               = backend, 
                              uName                 = sample)
                configureJets(variProxy             = tree._FatJet, 
                              jetType               = "AK8PFPuppi", 
                              jec                   = jecTag,
                              mayWriteCache         = isNotWorker, 
                              isMC                  = self.is_MC,
                              backend               = backend, 
                              uName                 = sample)
                configureType1MET(variProxy         = getattr(tree, f"_{metName}"),
                                  jec               = jecTag,
                                  mayWriteCache     = isNotWorker,
                                  isMC              = self.is_MC,
                                  backend           = backend, 
                                  uName             = sample)

        return tree,noSel,backend,lumiArgs

    #-------------------------------------------------------------------------------------------#
    #                                     prepareObjects                                        #
    #-------------------------------------------------------------------------------------------#
    def prepareObjects(self, t, noSel, sample, sampleCfg, channel, forSkimmer=False):
        # Some imports #

        if channel not in ["DL","SL"]:
            raise RuntimeError('Channel %s not understood'%channel)

        era = sampleCfg['era']
        self.era = era
        self.tree = t


        ###########################################################################
        #                              Pseudo-data                                #
        ###########################################################################
        if not forSkimmer and "PseudoData" in self.datadrivenContributions and self.is_MC: 
            # Skimmer does not know about self.datadrivenContributions
            noSel = SelectionWithDataDriven.create(parent   = noSel,
                                                   name     = 'pseudodata',
                                                   ddSuffix = 'Pseudodata',
                                                   enable   = "PseudoData" in self.datadrivenContributions 
                                    and self.datadrivenContributions["PseudoData"].usesSample(self.sample, self.sampleCfg))

        ###########################################################################
        #                          Signal Reweighting                             #
        ###########################################################################
        self.signalReweightBenchmarks = {}
        if 'type' in sampleCfg.keys() and sampleCfg["type"] == "signal" \
                                      and not forSkimmer \
                                      and 'benchmarks' in self.analysisConfig.keys() \
                                      and any(useFile in sample for useFile in self.analysisConfig['benchmarks']['uses']) \
                                      and self.args.HHReweighting is not None:
            # Get gen level Higgs #
            self.genh = op.select(t.GenPart,lambda g : op.AND(g.pdgId==25, g.statusFlags & ( 0x1 << 13)))
            # Get gen level variables mHH and cos(theta*) #
            HH_p4 = self.genh[0].p4 + self.genh[1].p4 
            cm = HH_p4.BoostToCM() 
            boosted_h = op.extMethod("ROOT::Math::VectorUtil::boost", returnType=self.genh[0].p4._typeName)(self.genh[0].p4,cm)
            self.mHH = op.invariant_mass(self.genh[0].p4,self.genh[1].p4) 
            self.cosHH = op.abs(boosted_h.Pz()/boosted_h.P())

            for v in (self.mHH, self.cosHH):
                forceDefine(v, noSel)
            def funConst(x, v=None):
                return v
            signalReweightParams = {
                    'Eta': partial(funConst, v=self.mHH),
                    'Pt' : partial(funConst, v=self.cosHH)
            }

            if forSkimmer:
                # In the case of the skimmer, we cannot use the create from datadriven.
                # At the same time we do not want the many-to-many procedure for the DNN training
                # because it would mean repeating events with different weights.
                # -> We use the one-to-one method 
                json_file = os.path.join(os.path.abspath(os.path.dirname(__file__)),'data','ScaleFactors_GGF_LO','{}_to_BenchmarkSM_{}.json'.format(sample,era))
                if os.path.exists(json_file):
                    print("Found file {} -> Will apply for LO->NLO reweighting".format(json_file))
                    reweightLO = get_scalefactor("lepton", json_file, paramDefs=signalReweightParams)(None)
                    noSel = noSel.refine("LoToNLO",weight = self.reweightLO)
                    if self.args.PrintYield:
                        self.yields.add(noSel)

            else:
                # For the Plotter we want the many-to-many worfklow.
                # This means create as many files (with the create) 
                # as there are NLO weight files (+ original LO file)
                for benchmark in self.analysisConfig['benchmarks']['targets']:
                    if benchmark in self.args.HHReweighting:
                        json_file = os.path.join(os.path.abspath(os.path.dirname(__file__)),'data','ScaleFactors_GGF_LO','{}_to_{}_{}.json'.format(sample,benchmark,era))
                        if os.path.exists(json_file):
                            print("Found file {} -> Will apply for LO->NLO reweighting".format(json_file))
                            # Apply reweighting to create #
                            reweightLO = get_scalefactor("lepton", json_file, paramDefs=signalReweightParams)(None)
                            self.signalReweightBenchmarks[benchmark] = reweightLO
                            forceDefine(reweightLO, noSel)

        ###########################################################################
        #                           TTbar reweighting                             #
        ###########################################################################
        if "group" in sampleCfg and sampleCfg["group"] == 'ttbar': 
            print ('Applied TT top reweighting')
            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#Use_case_3_ttbar_MC_is_used_to_m -> do not use when ttbar is background
            # Correct : https://indico.cern.ch/event/904971/contributions/3857701/attachments/2036949/3410728/TopPt_20.05.12.pdf
            #       -> Weight formula : slide 2
            #       -> top/antitop SF : slide 12 bottom left 
            # Get tops #
            self.genTop = op.select(t.GenPart,lambda g : op.AND(g.pdgId==6, g.statusFlags & ( 0x1 << 13)))
            self.genAntitop = op.select(t.GenPart,lambda g : op.AND(g.pdgId==-6, g.statusFlags & ( 0x1 << 13)))
                # statusFlags==13 : isLastCopy
                # Pdgid == 6 : top
            # Lambda to compute weight if there is a ttbar #
            self.t_SF = lambda t : op.exp(-2.02274e-01 + 1.09734e-04*t.pt + -1.30088e-07*t.pt**2 + (5.83494e+01/(t.pt+1.96252e+02)))
            self.ttbar_weight = lambda t,tbar : op.sqrt(self.t_SF(t)*self.t_SF(tbar))
            self.ttbar_sys = op.systematic(self.ttbar_weight(self.genTop[0],self.genAntitop[0]),
                                           name = "ttbarweightsyst",
                                           up   = op.pow(self.ttbar_weight(self.genTop[0],self.genAntitop[0]),2), # Up = apply twice
                                           down = op.c_float(1))                                                  # Down = not apply
            # Apply correction to TT #
            noSel = noSel.refine("ttbarWeight",weight=self.ttbar_sys)
            if self.args.PrintYield:
                self.yields.add(noSel)

        ###########################################################################
        #                               Stitching                                 #
        ###########################################################################
        if "group" in sampleCfg and (sampleCfg["group"] == 'DY' or sampleCfg["group"] == 'Wjets') and not self.args.NoStitching:
            stitch_file = os.path.abspath(os.path.join(os.path.dirname(__file__),'data','Stitching','stitching_weights_{}_{}.json'.format(sampleCfg["group"],era)))
            if not os.path.exists(stitch_file):
                raise RuntimeError("Could not find stitching file %s"%stitch_file)
            with open(stitch_file) as handle:
                dict_weights = json.load(handle)
            if sample in dict_weights.keys():
                dict_weights = dict_weights[sample]
                if isinstance(dict_weights[list(dict_weights.keys())[0]],float): # Only binned in jet multiplicity
                    stitch_op = op.multiSwitch(*[(t.LHE.Njets==int(njets),op.c_float(weight)) for njets,weight in dict_weights.items()], op.c_float(1.))
                elif isinstance(dict_weights[list(dict_weights.keys())[0]],list): # Binned in jet mult + HT bins
                    stitch_op = op.multiSwitch(*[(op.AND(t.LHE.Njets==int(njets),op.in_range(weights['low'],t.LHE.HT,weights['up'])),op.c_float(weights['value'])) 
                                                  for njets,listBin in dict_weights.items() for weights in listBin], op.c_float(1.))
                else:
                    raise RuntimeError("Stitching weight format not understood")
                noSel = noSel.refine("DYStitching",weight = stitch_op)
                if self.args.PrintYield:
                    self.yields.add(noSel)
                print ('Applied stitching')
                
        ###########################################################################
        #                               tH samples                                #
        ###########################################################################
        if sample.startswith('TH'):
            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SingleTopHiggsGeneration13TeV
            # https://twiki.cern.ch/twiki/pub/CMS/SingleTopHiggsGeneration13TeV/reweight_encondig.txt
            # 
            noSel = noSel.refine("tHWeight",weight=t.LHEReweightingWeight[11])
            print ('Applied tH LHE weights')
            if self.args.PrintYield:
                self.yields.add(noSel)
            # Record value for generated-sum
            self.base_plots.append(Plot.make1D("generated_sum_corrected",
                                               op.c_float(0.5),
                                               noSel,
                                               EquidistantBinning(1,0.,1.),
                                               weight=t.LHEReweightingWeight[11]*t.genWeight, 
                                               autoSyst=False))

        #############################################################################
        #                            Pre-firing rates                               #
        #############################################################################
        if era in ["2016","2017"] and self.is_MC and hasattr(t,'L1PreFiringWeight_Nom'):
            self.L1Prefiring = op.systematic(t.L1PreFiringWeight_Nom,
                                             name = "L1PreFiring",
                                             up   = t.L1PreFiringWeight_Up,
                                             down = t.L1PreFiringWeight_Dn)
            
            noSel = noSel.refine("L1PreFiringRate", weight = self.L1Prefiring)
            if self.args.PrintYield:
                self.yields.add(noSel)

        #############################################################################
        #                             Pile-up                                       #
        #############################################################################
        # Get MC PU weight file #
        puWeightsFile = None
        if self.is_MC:
            if 'related-sample' in sampleCfg.keys():
                puWeightsFile = os.path.join(os.path.dirname(__file__), "data", "pileup",f'{sampleCfg["related-sample"]}_{era}.json')
                nameHint = f'puweightFromFile{sampleCfg["related-sample"]}'.replace('-','_')
            else:
                puWeightsFile = os.path.join(os.path.dirname(__file__), "data", "pileup",f'{sample}_{era}.json')
                nameHint = f'puweightFromFile{sample}'.replace('-','_')
            if not os.path.exists(puWeightsFile):
                raise RuntimeError("Could not find pileup file %s"%puWeightsFile)
            from bamboo.analysisutils import makePileupWeight
            self.PUWeight = makePileupWeight(puWeightsFile, t.Pileup_nTrueInt, systName="pileup",nameHint=nameHint)
            noSel = noSel.refine("puWeight", weight = self.PUWeight)
            if self.args.PrintYield:
                self.yields.add(noSel)

        #############################################################################
        #                                 MET                                       #
        #############################################################################
        # MET filter #
        if not self.inclusive_sel:
            noSel = noSel.refine("passMETFlags", cut=METFilter(t.Flag, era, self.is_MC) )
            if self.args.PrintYield:
                self.yields.add(noSel)

        # MET corrections #
        self.rawMET = t.MET if era != "2017" else t.METFixEE2017
        self.corrMET = METcorrection(self.rawMET,t.PV,sample,era,self.is_MC) 


        #############################################################################
        #                      Lepton Lambdas Variables                             #
        #############################################################################
        # Associated jet Btagging #
        self.lambda_hasAssociatedJet = lambda lep : lep.jet.idx != -1
        if era == "2016": 
            self.lambda_lepton_associatedJetLessThanMediumBtag = lambda lep : op.OR(op.NOT(self.lambda_hasAssociatedJet(lep)),
                                                                        lep.jet.btagDeepFlavB <= 0.3093)
            self.lambda_lepton_associatedJetLessThanTightBtag  = lambda lep : op.OR(op.NOT(self.lambda_hasAssociatedJet(lep)),
                                                                        lep.jet.btagDeepFlavB <= 0.7221)
        elif era =="2017":
            self.lambda_lepton_associatedJetLessThanMediumBtag = lambda lep : op.OR(op.NOT(self.lambda_hasAssociatedJet(lep)),
                                                                        lep.jet.btagDeepFlavB <= 0.3033)
            self.lambda_lepton_associatedJetLessThanTightBtag  = lambda lep : op.OR(op.NOT(self.lambda_hasAssociatedJet(lep)),
                                                                        lep.jet.btagDeepFlavB <= 0.7489)
        elif era == "2018":
            self.lambda_lepton_associatedJetLessThanMediumBtag = lambda lep : op.OR(op.NOT(self.lambda_hasAssociatedJet(lep)),
                                                                        lep.jet.btagDeepFlavB <= 0.2770)
            self.lambda_lepton_associatedJetLessThanTightBtag  = lambda lep : op.OR(op.NOT(self.lambda_hasAssociatedJet(lep)),
                                                                        lep.jet.btagDeepFlavB <= 0.7264)
                     # If no associated jet, isolated lepton : cool !  

        # Cone pt #
        self.lambda_conept_electron = lambda lep : op.multiSwitch((op.AND(op.abs(lep.pdgId)!=11 , op.abs(lep.pdgId)!=13) , op.static_cast("Float_t",lep.pt)),
                                                                  # if (abs(lep.pdgId)!=11 and abs(lep.pdgId)!=13): return lep.pt : anything that is not muon or electron
                                                                  (op.AND(op.abs(lep.pdgId)==11 , lep.mvaTTH > 0.30) , op.static_cast("Float_t",lep.pt)),
                                                                  # if electron, check above MVA 
                                                                   op.static_cast("Float_t",0.9*lep.pt*(1.+lep.jetRelIso)))
                                                                  # else: return 0.90 * lep.pt / jetPtRatio where jetPtRatio = 1./(Electron_jetRelIso + 1.)
        self.lambda_conept_muon = lambda lep : op.multiSwitch((op.AND(op.abs(lep.pdgId)!=11 , op.abs(lep.pdgId)!=13) , op.static_cast("Float_t",lep.pt)),
                                                               # if (abs(lep.pdgId)!=11 and abs(lep.pdgId)!=13): return lep.pt : anything that is not muon or electron
                                                              (op.AND(op.abs(lep.pdgId)==13 , lep.mediumId ,lep.mvaTTH > 0.50) , op.static_cast("Float_t",lep.pt)),
                                                               # if muon, check that passes medium and above MVA
                                                               op.static_cast("Float_t",0.9*lep.pt*(1.+lep.jetRelIso)))
                                                               # else: return 0.90 * lep.pt / lep.jetPtRatiov2

        self.electron_conept = op.map(t.Electron, self.lambda_conept_electron)
        self.muon_conept = op.map(t.Muon, self.lambda_conept_muon)

        # Btag interpolation #
                    # https://indico.cern.ch/event/812025/contributions/3475878/attachments/1867083/3070589/gp-fr-run2b.pdf (slide 7)
        self.lambda_muon_x = lambda mu : op.min(op.max(0.,(0.9*mu.pt*(1+mu.jetRelIso))-20.)/(45.-20.), 1.)
                    # x = min(max(0, jet_pt-PT_min)/(PT_max-PT_min), 1) where jet_pt = 0.9*PT_muon*(1+MuonJetRelIso), PT_min=25, PT_max=40
        if era == "2016": 
            self.lambda_muon_btagInterpolation = lambda mu : self.lambda_muon_x(mu)*0.0614 + (1-self.lambda_muon_x(mu))*0.3093
        elif era =="2017":
            self.lambda_muon_btagInterpolation = lambda mu : self.lambda_muon_x(mu)*0.0521 + (1-self.lambda_muon_x(mu))*0.3033
        elif era == "2018":
            self.lambda_muon_btagInterpolation = lambda mu : self.lambda_muon_x(mu)*0.0494 + (1-self.lambda_muon_x(mu))*0.2770
            # return x*WP_loose+(1-x)*WP_medium
        #self.lambda_muon_deepJetInterpIfMvaFailed = lambda mu : mu.jet.btagDeepFlavB < self.lambda_muon_btagInterpolation(mu)
        self.lambda_muon_deepJetInterpIfMvaFailed = lambda mu : op.OR(op.NOT(self.lambda_hasAssociatedJet(mu)),     # If no associated jet, isolated lepton : cool !
                                                                      mu.jet.btagDeepFlavB < self.lambda_muon_btagInterpolation(mu))

        # Dilepton lambdas #
        self.lambda_leptonOS  = lambda l1,l2 : l1.charge != l2.charge
        if self.is_MC:
            self.lambda_is_matched = lambda lep : op.OR(lep.genPartFlav==1,  # Prompt muon or electron
                                                        lep.genPartFlav==15) # From tau decay
                                                        #lep.genPartFlav==22) # From photon conversion (only available for electrons)
        else:
            self.lambda_is_matched = lambda lep : op.c_bool(True)
        
        # Tight and fake selections 
        # Tight Dilepton : must also be Gen matched if MC #
        self.lambda_dilepton_matched = lambda dilep : op.AND(self.lambda_is_matched(dilep[0]),self.lambda_is_matched(dilep[1]))
        self.lambda_tightpair_ElEl = lambda dilep : op.AND(self.lambda_dilepton_matched(dilep),
                                                           self.lambda_electronTightSel(dilep[0]),
                                                           self.lambda_electronTightSel(dilep[1]))
        self.lambda_tightpair_MuMu = lambda dilep : op.AND(self.lambda_dilepton_matched(dilep),
                                                           self.lambda_muonTightSel(dilep[0]),
                                                           self.lambda_muonTightSel(dilep[1]))
        self.lambda_tightpair_ElMu = lambda dilep : op.AND(self.lambda_dilepton_matched(dilep),
                                                           self.lambda_electronTightSel(dilep[0]),
                                                           self.lambda_muonTightSel(dilep[1]))
             
        # Fake Extrapolation dilepton #               
        self.lambda_fakepair_ElEl = lambda dilep : op.AND(self.lambda_dilepton_matched(dilep),
                                                          op.NOT(op.AND(self.lambda_electronTightSel(dilep[0]),
                                                                        self.lambda_electronTightSel(dilep[1]))))
        self.lambda_fakepair_MuMu = lambda dilep : op.AND(self.lambda_dilepton_matched(dilep),
                                                          op.NOT(op.AND(self.lambda_muonTightSel(dilep[0]),
                                                                        self.lambda_muonTightSel(dilep[1]))))
        self.lambda_fakepair_ElMu = lambda dilep : op.AND(self.lambda_dilepton_matched(dilep),
                                                          op.NOT(op.AND(self.lambda_electronTightSel(dilep[0]),
                                                                        self.lambda_muonTightSel(dilep[1]))))
                

        #############################################################################
        #                                 Muons                                     #
        #############################################################################
        muonsByPt = op.sort(t.Muon, lambda mu : -self.muon_conept[mu.idx]) # Ordering done by conept
        # Preselection #
        self.lambda_muonPreSel = lambda mu : op.AND(mu.pt >= 5.,
                                                    op.abs(mu.eta) <= 2.4,
                                                    op.abs(mu.dxy) <= 0.05,
                                                    op.abs(mu.dz) <= 0.1,
                                                    mu.miniPFRelIso_all <= 0.4, # mini PF relative isolation, total (with scaled rho*EA PU corrections)
                                                    mu.sip3d <= 8,
                                                    mu.looseId)

        self.muonsPreSel = op.select(muonsByPt, self.lambda_muonPreSel)
        # Fakeable selection #
        self.lambda_muonFakeSel = lambda mu : op.AND(self.muon_conept[mu.idx] >= op.c_float(10.),
                                                     self.lambda_lepton_associatedJetLessThanMediumBtag(mu),
                                                     op.OR(mu.mvaTTH >= 0.50, op.AND(mu.jetRelIso<0.8 , self.lambda_muon_deepJetInterpIfMvaFailed(mu))))
                                                    # If mvaTTH < 0.50 : jetRelIso <0.8 and < deepJet medium with interpolation
        self.muonsFakeSel = op.select(self.muonsPreSel, self.lambda_muonFakeSel)
        # Tight selection #
        self.lambda_muonTightSel = lambda mu : op.AND(mu.mvaTTH >= 0.50, # Lepton MVA id from ttH
                                                      mu.mediumId)

        self.muonsTightSel = op.select(self.muonsFakeSel, self.lambda_muonTightSel)

        # Cone-pt 4-vector lambda #
        self.getMuonConeP4 = lambda lep : op.construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >", (self.muon_conept[lep.idx], 
                                                                                                                      lep.eta, lep.phi, lep.mass))

        #############################################################################
        #                              Electrons                                    #
        #############################################################################
        electronsByPt = op.sort(t.Electron, lambda ele : -self.electron_conept[ele.idx]) # Ordering done by conept
        self.lambda_electronPreSel = lambda ele : op.AND(ele.pt >= 7.,
                                                         op.abs(ele.eta) <= 2.5,
                                                         op.abs(ele.dxy) <= 0.05,
                                                         op.abs(ele.dz) <= 0.1,
                                                         ele.miniPFRelIso_all <= 0.4, # mini PF relative isolation, total (with scaled rho*EA PU corrections)

                                                         ele.sip3d <= 8,
                                                         ele.mvaFall17V2noIso_WPL, 
                                                         ele.lostHits <=1)    # number of missing inner hits

        self.electronsPreSelInclu = op.select(electronsByPt, self.lambda_electronPreSel) # can include a muon in cone
        self.lambda_cleanElectron = lambda ele : op.NOT(op.rng_any(self.muonsPreSel, lambda mu : op.deltaR(mu.p4, ele.p4) <= 0.3 ))
            # No overlap between electron and muon in cone of DR<=0.3
        self.electronsPreSel = op.select(self.electronsPreSelInclu, self.lambda_cleanElectron)
        # Fakeable selection #
        self.lambda_electronFakeSel = lambda ele : op.AND(self.electron_conept[ele.idx] >= 10,
                                                          op.OR(
                                                                  op.AND(op.abs(ele.eta+ele.deltaEtaSC)<=1.479, ele.sieie<=0.011), 
                                                                  op.AND(op.abs(ele.eta+ele.deltaEtaSC)>1.479, ele.sieie<=0.030)),
                                                          ele.hoe <= 0.10,
                                                          ele.eInvMinusPInv >= -0.04,
                                                          op.OR(ele.mvaTTH >= 0.30, op.AND(ele.jetRelIso<0.7, ele.mvaFall17V2noIso_WP90)), # Lepton MVA id from ttH
                                                              # If mvaTTH < 0.30 : jetRelIso <0.7 and Fall17V2noIso_WP90
                                                          op.switch(ele.mvaTTH < 0.30, 
                                                                    self.lambda_lepton_associatedJetLessThanTightBtag(ele),
                                                                    self.lambda_lepton_associatedJetLessThanMediumBtag(ele)),
                                                          ele.lostHits == 0,    # number of missing inner hits
                                                          ele.convVeto)        # Passes conversion veto

        self.electronsFakeSel = op.select(self.electronsPreSel, self.lambda_electronFakeSel)
        # Tight selection #
        self.lambda_electronTightSel = lambda ele : ele.mvaTTH >= 0.30
        self.electronsTightSel = op.select(self.electronsFakeSel, self.lambda_electronTightSel)

        # Cone-pt 4-vector lambda #
        self.getElectronConeP4 = lambda lep : op.construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >", (self.electron_conept[lep.idx], 
                                                                                                                          lep.eta, lep.phi, lep.mass))


        #############################################################################
        #                               Dileptons                                   #
        #############################################################################

        #----- Needed for both SL and DL -----#
        # Preselected dilepton -> Mll>12 cut #
        self.ElElDileptonPreSel = op.combine(self.electronsPreSelInclu, N=2)
        self.MuMuDileptonPreSel = op.combine(self.muonsPreSel, N=2)
        self.ElMuDileptonPreSel = op.combine((self.electronsPreSelInclu, self.muonsPreSel))
            # Need pairs without sign for one of the selection
            # We use electrons before muon cleaning : electron-muon pair that has low mass are typically close to each other
        # OS Preselected dilepton -> Z Veto #
        self.OSElElDileptonPreSel = op.combine(self.electronsPreSel, N=2, pred=self.lambda_leptonOS)
        self.OSMuMuDileptonPreSel = op.combine(self.muonsPreSel, N=2, pred=self.lambda_leptonOS)
        self.OSElMuDileptonPreSel = op.combine((self.electronsPreSel, self.muonsPreSel), pred=self.lambda_leptonOS) 

        # Dilepton for selection #
        if channel == "DL":
            self.ElElFakeSel = op.combine(self.electronsFakeSel, N=2)
            self.MuMuFakeSel = op.combine(self.muonsFakeSel, N=2)
            self.ElMuFakeSel = op.combine((self.electronsFakeSel,self.muonsFakeSel))

            self.ElElTightSel = op.combine(self.electronsTightSel, N=2)
            self.MuMuTightSel = op.combine(self.muonsTightSel, N=2)
            self.ElMuTightSel = op.combine((self.electronsTightSel,self.muonsTightSel))

        #############################################################################
        #                               Triggers                                    #
        #############################################################################
          
        #----- Trigger mask for data -----#
        if not self.is_MC:
            if era == "2018": # For 2018 the electron samples are merged, need to build according dict for primary dataset
                trigger_per_pd = { pd : trig for pd, trig in self.triggersPerPrimaryDataset.items() if pd not in ("SingleElectron", "DoubleEGamma") } 
                trigger_per_pd["EGamma"] = self.triggersPerPrimaryDataset["SingleElectron"] + self.triggersPerPrimaryDataset["DoubleEGamma"]
                noSel = noSel.refine("DataMaskTriggers", cut=[makeMultiPrimaryDatasetTriggerSelection(sample, trigger_per_pd)])
            else:
                noSel = noSel.refine("DataMaskTriggers", cut=[makeMultiPrimaryDatasetTriggerSelection(sample, self.triggersPerPrimaryDataset)])
            if self.args.PrintYield:
                self.yields.add(noSel)
            # makeMultiPrimaryDatasetTriggerSelection must be done first, to make sure an event is not passed several times from several datasets
            # Then the trigger selection is based on fakeable leptons ->  done at lepton selection (before jet forceDefine)

        ##############################################################################
        #                                  Tau                                       #
        ##############################################################################
        self.lambda_tauSel = lambda ta : op.AND(ta.pt > 20.,
                                                op.abs(ta.p4.Eta()) < 2.3,
                                                op.abs(ta.dxy) <= 1000.0,
                                                op.abs(ta.dz) <= 0.2,
                                                ta.idDecayModeNewDMs,
                                                op.OR(ta.decayMode == 0,
                                                      ta.decayMode == 1,
                                                      ta.decayMode == 2,
                                                      ta.decayMode == 10,
                                                      ta.decayMode == 11),
                                                (ta.idDeepTau2017v2p1VSjet >> 4 & 0x1) == 1,
                                                (ta.idDeepTau2017v2p1VSe >> 0 & 0x1)   == 1,
                                                (ta.idDeepTau2017v2p1VSmu >> 0 & 0x1)  == 1
                                            )
        self.tauSel = op.select (t.Tau, self.lambda_tauSel)
        # Cleaning #
        self.lambda_tauClean = lambda ta : op.AND(op.NOT(op.rng_any(self.electronsFakeSel, lambda el : op.deltaR(ta.p4, el.p4) <= 0.3)), 
                                                  op.NOT(op.rng_any(self.muonsFakeSel, lambda mu : op.deltaR(ta.p4, mu.p4) <= 0.3)))
        self.tauCleanSel = op.select(self.tauSel, self.lambda_tauClean)

        #############################################################################
        #                                AK4 Jets                                   #
        #############################################################################
        self.ak4JetsByPt = op.sort(t.Jet, lambda jet : -jet.pt)
        # Preselection #
        self.lambda_ak4JetsPreSel = lambda j : op.AND(j.jetId & 1 if era == "2016" else j.jetId & 2, # Jet ID flags bit1 is loose, bit2 is tight, bit3 is tightLepVeto
                                                      j.pt >= 25.,
                                                      op.abs(j.eta) <= 2.4)
        self.lambda_jetPUID = lambda j : op.OR(((j.puId >> 2) & 1) ,j.pt > 50.) # Jet PU ID bit1 is loose (only to be applied to jets with pt<50)

        self.ak4JetsPreSel        = op.select(self.ak4JetsByPt, lambda j : op.AND(self.lambda_ak4JetsPreSel(j),self.lambda_jetPUID(j)))

        # Cleaning #
        if channel == 'SL':
            def returnLambdaCleaningWithRespectToLeadingLepton(DR):
               return lambda j : op.multiSwitch(
                        (op.AND(op.rng_len(self.electronsFakeSel) >= 1, op.rng_len(self.muonsFakeSel) == 0), op.deltaR(j.p4, self.electronsFakeSel[0].p4)>=DR),
                        (op.AND(op.rng_len(self.electronsFakeSel) == 0, op.rng_len(self.muonsFakeSel) >= 1), op.deltaR(j.p4, self.muonsFakeSel[0].p4)>=DR),
                        (op.AND(op.rng_len(self.muonsFakeSel) >= 1, op.rng_len(self.electronsFakeSel) >= 1),
                               op.switch(self.electron_conept[self.electronsFakeSel[0].idx] >= self.muon_conept[self.muonsFakeSel[0].idx],
                                         op.deltaR(j.p4, self.electronsFakeSel[0].p4)>=DR,
                                         op.deltaR(j.p4, self.muonsFakeSel[0].p4)>=DR)),
                        op.c_bool(True))
            self.lambda_cleanAk4Jets = returnLambdaCleaningWithRespectToLeadingLepton(0.4)

        # remove jets within cone of DR<0.4 of leading lept
        if channel == 'DL':
            def returnLambdaCleaningWithRespectToLeadingLeptons(DR):
                return lambda j : op.multiSwitch(
                      (op.AND(op.rng_len(self.electronsFakeSel) >= 2,op.rng_len(self.muonsFakeSel) == 0), 
                          # Only electrons 
                          op.AND(op.deltaR(j.p4, self.electronsFakeSel[0].p4)>=DR, op.deltaR(j.p4, self.electronsFakeSel[1].p4)>=DR)),
                      (op.AND(op.rng_len(self.electronsFakeSel) == 0,op.rng_len(self.muonsFakeSel) >= 2), 
                          # Only muons  
                          op.AND(op.deltaR(j.p4, self.muonsFakeSel[0].p4)>=DR, op.deltaR(j.p4, self.muonsFakeSel[1].p4)>=DR)),
                      (op.AND(op.rng_len(self.electronsFakeSel) == 1,op.rng_len(self.muonsFakeSel) == 1),
                          # One electron + one muon
                          op.AND(op.deltaR(j.p4, self.electronsFakeSel[0].p4)>=DR, op.deltaR(j.p4, self.muonsFakeSel[0].p4)>=DR)),
                      (op.AND(op.rng_len(self.electronsFakeSel) >= 1,op.rng_len(self.muonsFakeSel) >= 1),
                          # At least one electron + at least one muon
                       op.switch(self.electron_conept[self.electronsFakeSel[0].idx] > self.muon_conept[self.muonsFakeSel[0].idx],
                                 # Electron is leading #
                                 op.switch(op.rng_len(self.electronsFakeSel) == 1,
                                           op.AND(op.deltaR(j.p4, self.electronsFakeSel[0].p4)>=DR, op.deltaR(j.p4, self.muonsFakeSel[0].p4)>=DR),
                                           op.switch(self.electron_conept[self.electronsFakeSel[1].idx] > self.muon_conept[self.muonsFakeSel[0].idx],
                                                     op.AND(op.deltaR(j.p4, self.electronsFakeSel[0].p4)>=DR, op.deltaR(j.p4, self.electronsFakeSel[1].p4)>=DR),
                                                     op.AND(op.deltaR(j.p4, self.electronsFakeSel[0].p4)>=DR, op.deltaR(j.p4, self.muonsFakeSel[0].p4)>=DR))),
                                 # Muon is leading #
                                 op.switch(op.rng_len(self.muonsFakeSel) == 1,
                                           op.AND(op.deltaR(j.p4, self.muonsFakeSel[0].p4)>=DR, op.deltaR(j.p4, self.electronsFakeSel[0].p4)>=DR),
                                           op.switch(self.muon_conept[self.muonsFakeSel[1].idx] > self.electron_conept[self.electronsFakeSel[0].idx],
                                                     op.AND(op.deltaR(j.p4, self.muonsFakeSel[0].p4)>=DR, op.deltaR(j.p4, self.muonsFakeSel[1].p4)>=DR),
                                                     op.AND(op.deltaR(j.p4, self.muonsFakeSel[0].p4)>=DR, op.deltaR(j.p4, self.electronsFakeSel[0].p4)>=DR))))),
                       op.c_bool(True))
            self.lambda_cleanAk4Jets = returnLambdaCleaningWithRespectToLeadingLeptons(0.4)
        

        self.ak4Jets            = op.select(self.ak4JetsPreSel,self.lambda_cleanAk4Jets) # Pt ordered
        self.ak4JetsByBtagScore = op.sort(self.ak4Jets, lambda j : -j.btagDeepFlavB) # Btag score ordered
    
        ############     Btagging     #############
        # The pfDeepFlavour (DeepJet) algorithm is used
        if era == "2016": 
            self.lambda_ak4BtagLoose =   lambda jet    : jet.btagDeepFlavB > 0.0614
            self.lambda_ak4Btag =   lambda jet    : jet.btagDeepFlavB > 0.3093
            self.lambda_ak4NoBtag = lambda jet    : jet.btagDeepFlavB <= 0.3093
        elif era =="2017":
            self.lambda_ak4BtagLoose =   lambda jet    : jet.btagDeepFlavB > 0.0521
            self.lambda_ak4Btag =   lambda jet    : jet.btagDeepFlavB > 0.3033
            self.lambda_ak4NoBtag = lambda jet    : jet.btagDeepFlavB <= 0.3033
        elif era == "2018":
            self.lambda_ak4BtagLoose =   lambda jet    : jet.btagDeepFlavB > 0.0494
            self.lambda_ak4Btag =   lambda jet    : jet.btagDeepFlavB > 0.2770
            self.lambda_ak4NoBtag = lambda jet    : jet.btagDeepFlavB <= 0.2770

        self.ak4BJets = op.select(self.ak4Jets, self.lambda_ak4Btag)
        self.ak4BJetsLoose = op.select(self.ak4Jets, self.lambda_ak4BtagLoose)
        self.ak4LightJetsByPt = op.select(self.ak4Jets, self.lambda_ak4NoBtag)
        self.ak4LightJetsByBtagScore = op.sort(self.ak4LightJetsByPt, lambda jet : -jet.btagDeepFlavB)
        # Sorted by btag score because for 0 and 1 btag categories, 

        # Doesn't contain the leading bTag scored Light Jet
        # --------------- not used -------------------- #
        self.remainingJets = op.select(self.ak4LightJetsByPt, lambda jet : jet.idx != self.ak4LightJetsByBtagScore[0].idx)
        self.makeJetPairs  = lambda jets : op.combine(jets, N=2, pred=lambda j1, j2 : j1.pt > j2.pt, samePred=lambda j1,j2 : j1.idx != j2.idx)
        # --------------------------------------------- #
        self.bJetsByScore        = self.ak4JetsByBtagScore[:op.min(op.rng_len(self.ak4JetsByBtagScore),op.static_cast("std::size_t",op.c_int(2)))]
        self.probableWJets       = op.select(self.ak4Jets, lambda jet : op.NOT(op.rng_any(self.bJetsByScore, lambda bjet : jet.idx == bjet.idx)))
        self.wJetsByPt           = self.probableWJets[:op.min(op.rng_len(self.probableWJets),op.static_cast("std::size_t",op.c_int(2)))] # used as real wjets, not used for VBF

        #wMassWindow              = lambda dijet : op.abs(op.invariant_mass(dijet[0].p4,dijet[1].p4)-80.4)
        #probableWJetPairs        = op.combine(self.probableWJets, N=2)
        #wJetsByPtPair            = op.combine(self.wJetsByPt, N=2)
        #self.wJetsPairs          = op.sort(op.select(probableWJetPairs, lambda dijet : wMassWindow(dijet) < op.c_float(15.0)), 
        #                                   lambda dijet : wMassWindow(dijet)) # used for VBF selection

        self.lambda_passWMassCut = lambda wjets : op.switch(op.rng_len(wjets) == 2, op.abs(op.invariant_mass(wjets[0].p4, wjets[1].p4)-80.4) < op.c_float(15.0), op.c_bool(False))

        #############################################################################
        #                                AK8 Jets                                   #
        #############################################################################
        self.ak8JetsByPt = op.sort(t.FatJet, lambda jet : -jet.pt)
        self.ak8JetsByDeepB = op.sort(t.FatJet, lambda jet : -jet.btagDeepB)
        # Preselection #
        self.lambda_ak8JetsPreSel = lambda j : op.AND(j.jetId & 1 if era == "2016" else j.jetId & 2, # Jet ID flags bit1 is loose, bit2 is tight, bit3 is tightLepVeto
                                                      j.pt >= 200.,
                                                      op.abs(j.p4.Eta())<= 2.4,
                                                      op.AND(j.subJet1.isValid, j.subJet1.pt >= 20. , op.abs(j.subJet1.eta)<=2.4,
                                                             j.subJet2.isValid, j.subJet2.pt >= 20. , op.abs(j.subJet2.eta)<=2.4),
                                                             # Fatjet subjets must exist before checking Pt and eta 
                                                      op.AND(j.msoftdrop >= 30, j.msoftdrop <= 210),
                                                      j.tau2/j.tau1 <= 0.75)

        if channel == 'SL':
            self.ak8JetsPreSel = op.select(self.ak8JetsByDeepB, self.lambda_ak8JetsPreSel)
        if channel == 'DL':
            self.ak8JetsPreSel = op.select(self.ak8JetsByPt, self.lambda_ak8JetsPreSel)

        # Cleaning #
        if channel == 'SL':
            self.lambda_cleanAk8Jets = returnLambdaCleaningWithRespectToLeadingLepton(0.8)
        if channel == 'DL':
            self.lambda_cleanAk8Jets = returnLambdaCleaningWithRespectToLeadingLeptons(0.8)
        # remove jets within cone of DR<0.8 of preselected electrons and muons
        self.ak8Jets = op.select(self.ak8JetsPreSel,self.lambda_cleanAk8Jets)

        ############     Btagging     #############
        # The DeepCSV b-tagging algorithm is used on subjets 
        if era == "2016": 
            self.lambda_subjetBtag = lambda subjet: subjet.btagDeepB > 0.6321
        elif era =="2017":
            self.lambda_subjetBtag = lambda subjet: subjet.btagDeepB > 0.4941
        elif era == "2018":
            self.lambda_subjetBtag = lambda subjet: subjet.btagDeepB > 0.4184

        self.lambda_ak8Btag = lambda fatjet : op.OR(op.AND(fatjet.subJet1.pt >= 30, self.lambda_subjetBtag(fatjet.subJet1)),
                                                    op.AND(fatjet.subJet2.pt >= 30, self.lambda_subjetBtag(fatjet.subJet2)))
        self.lambda_ak8noBtag = lambda fatjet : op.NOT(op.OR(op.AND(fatjet.subJet1.pt >= 30, self.lambda_subjetBtag(fatjet.subJet1)),
                                                       op.AND(fatjet.subJet2.pt >= 30, self.lambda_subjetBtag(fatjet.subJet2))))
        self.lambda_ak8Btag_bothSubJets = lambda fatjet : op.AND(op.AND(fatjet.subJet1.pt >= 30, self.lambda_subjetBtag(fatjet.subJet1)),
                                                                 op.AND(fatjet.subJet2.pt >= 30, self.lambda_subjetBtag(fatjet.subJet2)))

        self.ak8BJets = op.select(self.ak8Jets, self.lambda_ak8Btag)
        self.ak8nonBJets = op.select(self.ak8Jets, self.lambda_ak8noBtag)
        # Ak4 Jet Collection cleaned from Ak8b #
        self.lambda_cleanAk4FromAk8b = lambda ak4j : op.AND(op.rng_len(self.ak8BJets) > 0, op.deltaR(ak4j.p4,self.ak8BJets[0].p4) > 1.2)
        self.ak4JetsCleanedFromAk8b  = op.select(self.ak4Jets, self.lambda_cleanAk4FromAk8b)
        
        # used as a BDT input for SemiBoosted category
        self.lambda_btaggedSubJets = lambda fjet : op.switch(self.lambda_ak8Btag_bothSubJets(fjet), op.c_float(2.0), op.c_float(1.0))
        self.nMediumBTaggedSubJets = op.rng_sum(self.ak8BJets, self.lambda_btaggedSubJets)

        #############################################################################
        #                                VBF Jets                                   #
        #############################################################################
        self.lambda_VBFJets = lambda j : op.AND(j.jetId & 1 if era == "2016" else j.jetId & 2, # Jet ID flags bit1 is loose, bit2 is tight, bit3 is tightLepVeto
                                                j.pt >= 30.,
                                                op.abs(j.eta) <= 4.7,
                                                op.OR(j.pt >= 60.,
                                                      op.abs(j.eta) < 2.7, 
                                                      op.abs(j.eta) > 3.0))
        self.VBFJetsPreSel        = op.select(self.ak4JetsByPt, lambda j : op.AND(self.lambda_VBFJets(j),self.lambda_jetPUID(j)))


        self.VBFJets            = op.select(self.VBFJetsPreSel, self.lambda_cleanAk4Jets)

        self.lambda_VBFPair     = lambda j1,j2 : op.AND(op.invariant_mass(j1.p4,j2.p4) > 500.,
                                                        op.abs(j1.eta - j2.eta) > 3.)
        
        if channel == "DL":
            self.lambda_cleanVBFAk4 = lambda j : op.multiSwitch((op.rng_len(self.ak4JetsByBtagScore)>1,op.AND(op.deltaR(j.p4, self.ak4JetsByBtagScore[0].p4)>0.8,
                                                                                                              op.deltaR(j.p4, self.ak4JetsByBtagScore[1].p4)>0.8)),
                                                                (op.rng_len(self.ak4JetsByBtagScore)==1,op.deltaR(j.p4, self.ak4JetsByBtagScore[0].p4)>0.8),
                                                                op.c_bool(True))
            self.lambda_cleanVBFAk8 = lambda j : op.multiSwitch((op.rng_len(self.ak8BJets)>0,op.deltaR(j.p4, self.ak8BJets[0].p4) > 1.2),
                                                                (op.rng_len(self.ak8Jets)>0,op.deltaR(j.p4, self.ak8Jets[0].p4) > 1.2),
                                                                op.c_bool(True))


            self.VBFJetsResolved = op.select(self.VBFJets, self.lambda_cleanVBFAk4)
            self.VBFJetsBoosted  = op.select(self.VBFJets, self.lambda_cleanVBFAk8)

            self.VBFJetPairs         = op.sort(op.combine(self.VBFJets, N=2, pred=self.lambda_VBFPair), lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))
            self.VBFJetPairsResolved = op.sort(op.combine(self.VBFJetsResolved, N=2, pred=self.lambda_VBFPair), lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))
            self.VBFJetPairsBoosted  = op.sort(op.combine(self.VBFJetsBoosted,  N=2, pred=self.lambda_VBFPair), lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))

        if channel == "SL":
            #lambda_isOverlappedWithWjets  = lambda j : op.rng_any(self.wJetsPairs, lambda wjp : op.OR(op.deltaR(wjp[0].p4, j.p4) <= 0.8, op.deltaR(wjp[1].p4, j.p4) <= 0.8))
            lambda_isOverlappedWithWjets  = lambda j : op.switch(self.lambda_passWMassCut(self.wJetsByPt), 
                                                                 op.OR(op.deltaR(self.wJetsByPt[0].p4, j.p4) <= 0.8, op.deltaR(self.wJetsByPt[1].p4, j.p4) <= 0.8),
                                                                 op.c_bool(False))
            lambda_isOverlappedWithBjets  = lambda j : op.rng_any(self.bJetsByScore, lambda bj : op.deltaR(bj.p4, j.p4) <= 0.8)
            self.lambda_cleanVBFAk4       = lambda j : op.NOT(op.OR(lambda_isOverlappedWithBjets(j), lambda_isOverlappedWithWjets(j)))
            self.lambda_cleanVBFAk8       = lambda j : op.NOT(op.deltaR(j.p4, self.ak8BJets[0].p4) <= 1.2)
            
            self.VBFJetsResolved          = op.select(self.VBFJets, self.lambda_cleanVBFAk4)
            self.VBFJetsBoosted           = op.select(self.VBFJets, self.lambda_cleanVBFAk8)

            self.VBFJetPairsResolved      = op.sort(op.combine(self.VBFJetsResolved, N=2, pred=self.lambda_VBFPair), lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))
            self.VBFJetPairsBoosted       = op.sort(op.combine(self.VBFJetsBoosted,  N=2, pred=self.lambda_VBFPair), lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))


        #############################################################################
        #                                 PU Jets                                   #
        #############################################################################
        # Correction to be computed later on all the ak4 jets (regular + VBF) on which the cut is applied 
        # To avoid double counting : OR(regular jet, VBF jet)
        self.ak4ForPUID = op.select(self.ak4JetsByPt, lambda j : op.AND(j.pt<=50.,
                                                                        op.OR(self.lambda_VBFJets(j),
                                                                              self.lambda_ak4JetsPreSel(j)),
                                                                        self.lambda_cleanAk4Jets(j)))

        #############################################################################
        #                             Scalefactors                                  #
        #############################################################################
        self.SF = ScaleFactorsbbWW()    
        if self.is_MC:
            #---- PU ID SF ----#
            # Efficiency and mistags do not have uncertainties, the systematics are in the SF 
            self.jetpuid_mc_eff = self.SF.get_scalefactor("lepton", ('jet_puid_eff','eff_{}_L'.format(era)),combine="weight", defineOnFirstUse=(not forSkimmer))
            self.jetpuid_mc_mis = self.SF.get_scalefactor("lepton", ('jet_puid_eff','mistag_{}_L'.format(era)),combine="weight", defineOnFirstUse=(not forSkimmer))
                # Eff and mistag do not have systematics, only the SF do
            self.jetpuid_sf_eff = self.SF.get_scalefactor("lepton", ('jet_puid_sf','eff_{}_L'.format(era)),combine="weight", systName="jetpuid_eff", defineOnFirstUse=(not forSkimmer))
            self.jetpuid_sf_mis = self.SF.get_scalefactor("lepton", ('jet_puid_sf','mistag_{}_L'.format(era)),combine="weight", systName="jetpuid_mistag", defineOnFirstUse=(not forSkimmer))

            #---- Object SF -----# (Will take as argument the era)
            ####  Muons ####
            self.muLooseId = self.SF.get_scalefactor("lepton", 'muon_loose_{}'.format(era), combine="weight", systName="mu_loose", defineOnFirstUse=(not forSkimmer)) 
            self.lambda_MuonLooseSF = lambda mu : [op.switch(self.lambda_is_matched(mu), self.muLooseId(mu), op.c_float(1.))]
            self.muTightMVA = self.SF.get_scalefactor("lepton", 'muon_tightMVA_{}'.format(era), systName="mu_tight", combine="weight", defineOnFirstUse=(not forSkimmer))
                # -> Old too tight ID
            self.muRelaxedTightMVA = self.SF.get_scalefactor("lepton", 'muon_tightMVArelaxed_{}'.format(era), combine="weight", systName="mu_tth", defineOnFirstUse=(not forSkimmer))
                # Systematic to take into account the relaxation of the ttH ID : |SF_final - SF| / 2 added as SF_final ± SF_err_x
                # -> New relaxed ID

            self.lambda_MuonTightSF = lambda mu : [op.switch(self.lambda_muonTightSel(mu),      # check if actual tight lepton (needed because in Fake CR one of the dilepton is not)
                                                             op.switch(self.lambda_is_matched(mu), # check if gen matched
                                                                       self.muTightMVA(mu) * self.muRelaxedTightMVA(mu),
                                                                       op.c_float(1.)),
                                                             op.c_float(1.))]

            self.muPOGTightID = self.SF.get_scalefactor("lepton", 'muon_POGSF_ID_{}'.format(era), combine="weight", systName="mu_pogtightid", defineOnFirstUse=(not forSkimmer))
            self.muPOGTightISO = self.SF.get_scalefactor("lepton", 'muon_POGSF_ISO_{}'.format(era), combine="weight", systName="mu_pogtightiso", defineOnFirstUse=(not forSkimmer))
            self.elPOGTight = self.SF.get_scalefactor("lepton", 'electron_POGSF_{}'.format(era), combine="weight", systName="el_pogtight", defineOnFirstUse=(not forSkimmer))

            ####  Electrons ####
            self.elLooseId = self.SF.get_scalefactor("lepton", 'electron_looseid_{}'.format(era) , combine="weight", systName="el_loose", defineOnFirstUse=(not forSkimmer))
            self.elLooseEff = self.SF.get_scalefactor("lepton", 'electron_looseeff_{}'.format(era) , combine="weight", systName="el_loose", defineOnFirstUse=(not forSkimmer))
            if era == "2016" or era == "2017": # Electron reco eff depend on Pt for 2016 and 2017
                self.elLooseRecoPtGt20 = self.SF.get_scalefactor("lepton", ('electron_loosereco_{}'.format(era) , 'electron_loosereco_ptgt20'),
                                                                 combine="weight", systName="el_loose", defineOnFirstUse=(not forSkimmer))
                self.elLooseRecoPtLt20 = self.SF.get_scalefactor("lepton", ('electron_loosereco_{}'.format(era) , 'electron_loosereco_ptlt20'),
                                                                 combine="weight", systName="el_loose", defineOnFirstUse=(not forSkimmer))
                self.elLooseReco = lambda el: op.switch(el.pt>20,self.elLooseRecoPtGt20(el),self.elLooseRecoPtGt20(el))
            elif era == "2018": # Does not depend on pt for 2018
                self.elLooseReco = self.SF.get_scalefactor("lepton", 'electron_loosereco_{}'.format(era), combine="weight", systName="el_loose", defineOnFirstUse=(not forSkimmer))

            # reco-to-loose single systematic = reco x loose id x loose eff
            # Up/down = product up/down -> need to have the same name 

            self.lambda_ElectronLooseSF = lambda el : [op.switch(self.lambda_is_matched(el),self.elLooseId(el),op.c_float(1.)),
                                                       op.switch(self.lambda_is_matched(el),self.elLooseEff(el),op.c_float(1.)),
                                                       op.switch(self.lambda_is_matched(el),self.elLooseReco(el),op.c_float(1.))]

            self.elTightMVA = self.SF.get_scalefactor("lepton", 'electron_tightMVA_{}'.format(era) , combine="weight", systName="el_tight", defineOnFirstUse=(not forSkimmer))
                # -> Old too tight ID
            self.elRelaxedTightMVA = self.SF.get_scalefactor("lepton", 'electron_tightMVArelaxed_{}'.format(era) , combine="weight", systName="el_tth", defineOnFirstUse=(not forSkimmer))
                # Systematic to take into account the relaxation of the ttH ID : |SF_final - SF| / 2 added as SF_final ± SF_err_x
                # -> New relaxed ID

            self.lambda_ElectronTightSF = lambda el : [op.switch(self.lambda_electronTightSel(el),      # check if actual tight lepton (needed because in Fake CR one of the dilepton is not)
                                                                 op.switch(self.lambda_is_matched(el), # check if gen matched
                                                                           self.elTightMVA(el) * self.elRelaxedTightMVA(el),
                                                                           op.c_float(1.)),
                                                                 op.c_float(1.))]

            #### Ak8 btag efficiencies ####
            self.Ak8Eff_bjets = self.SF.get_scalefactor("lepton",('ak8btag_eff','eff_bjets_{}'.format(era)),combine="weight", defineOnFirstUse=(not forSkimmer),
                                                        additionalVariables={'Pt':lambda x : x.pt,'Eta':lambda x : x.eta, 'BTagDiscri':lambda x : x.btagDeepB})
            self.Ak8Eff_cjets = self.SF.get_scalefactor("lepton",('ak8btag_eff','eff_cjets_{}'.format(era)),combine="weight", defineOnFirstUse=(not forSkimmer),
                                                        additionalVariables={'Pt':lambda x : x.pt,'Eta':lambda x : x.eta, 'BTagDiscri':lambda x : x.btagDeepB})
            self.Ak8Eff_lightjets = self.SF.get_scalefactor("lepton",('ak8btag_eff','eff_lightjets_{}'.format(era)),combine="weight", defineOnFirstUse=(not forSkimmer),
                                                        additionalVariables={'Pt':lambda x : x.pt,'Eta':lambda x : x.eta, 'BTagDiscri':lambda x : x.btagDeepB})

            #----- Triggers -----# 
            #### Single lepton triggers ####
            self.ttH_singleElectron_trigSF = self.SF.get_scalefactor("lepton", 'singleTrigger_electron_{}'.format(era) , combine="weight", systName="ttH_singleElectron_trigSF", defineOnFirstUse=(not forSkimmer))
            self.ttH_singleMuon_trigSF = self.SF.get_scalefactor("lepton", 'singleTrigger_muon_{}'.format(era) , combine="weight", systName="ttH_singleMuon_trigSF", defineOnFirstUse=(not forSkimmer))
            
            #### Double lepton triggers #### (Need to split according to era) 
                # https://gitlab.cern.ch/ttH_leptons/doc/-/blob/master/Legacy/data_to_mc_corrections.md#trigger-efficiency-scale-factors -> deprecated
                # New ref (more up to date) : https://cernbox.cern.ch/index.php/s/lW2BiTli5tJR0MN 
                # Same but clearer  : http://cms.cern.ch/iCMS/jsp/iCMS.jsp?mode=single&part=publications (AN-2019/111) at Table 30, page 43
                # -> Based on subleading conept lepton (except mumu 2018), use lambda 
                # Lambdas return list to be easily concatenated with other SF in lists as well (the weight argument of op.refine requires a list)
            if era == "2016":
                # Double Electron #
                self.lambda_ttH_doubleElectron_trigSF = lambda dilep : op.multiSwitch(
                    (self.electron_conept[dilep[1].idx]<25 , op.systematic(op.c_float(0.98), name="ttH_doubleElectron_trigSF", up=op.c_float(0.98*1.02), down=op.c_float(0.98*0.98))), 
                    op.systematic(op.c_float(1.), name="ttH_doubleElectron_trigSF", up=op.c_float(1.02), down=op.c_float(0.98)))
                # Double Muon #
                self.lambda_ttH_doubleMuon_trigSF = lambda dilep : op.systematic(op.c_float(0.99), name="ttH_doubleMuon_trigSF", up=op.c_float(0.99*1.01), down=op.c_float(0.99*0.99))
                # Electron Muon #
                self.lambda_ttH_electronMuon_trigSF = lambda dilep : op.systematic(op.c_float(1.00), name="ttH_electronMuon_trigSF", up=op.c_float(1.01), down=op.c_float(0.99))
            elif era == "2017":
                # Double Electron #
                self.lambda_ttH_doubleElectron_trigSF = lambda dilep : op.multiSwitch(
                    (self.electron_conept[dilep[1].idx]<40 , op.systematic(op.c_float(0.98), name="ttH_doubleElectron_trigSF", up=op.c_float(0.98*1.01), down=op.c_float(0.98*0.99))), 
                    op.systematic(op.c_float(1.), name="ttH_doubleElectron_trigSF", up=op.c_float(1.01), down=op.c_float(0.99)))
                # Double Muon #
                self.lambda_ttH_doubleMuon_trigSF = lambda dilep : op.multiSwitch(
                    (self.muon_conept[dilep[1].idx]<40 , op.systematic(op.c_float(0.97), name="ttH_doubleMuon_trigSF", up=op.c_float(0.97*1.02), down=op.c_float(0.97*0.98))),
                    (self.muon_conept[dilep[1].idx]<55 , op.systematic(op.c_float(0.995), name="ttH_doubleMuon_trigSF", up=op.c_float(0.995*1.02), down=op.c_float(0.995*0.98))),
                    (self.muon_conept[dilep[1].idx]<70 , op.systematic(op.c_float(0.96), name="ttH_doubleMuon_trigSF", up=op.c_float(0.96*1.02), down=op.c_float(0.96*0.98))),
                    op.systematic(op.c_float(0.94), name="ttH_doubleMuon_trigSF", up=op.c_float(0.94*1.02), down=op.c_float(0.94*0.98)))
                # Electron Muon : /!\ While ElEl and MuMu is conept ordered, ElMu is not #
                self.lambda_ttH_electronMuon_trigSF = lambda dilep : op.multiSwitch(
                    (op.min(self.electron_conept[dilep[0].idx],self.muon_conept[dilep[1].idx])<40, op.systematic(op.c_float(0.98), name="ttH_electronMuon_trigSF", up=op.c_float(0.98*1.01), down=op.c_float(0.98*0.99))),
                    op.systematic(op.c_float(0.99), name="ttH_electronMuon_trigSF", up=op.c_float(0.99*1.01), down=op.c_float(0.99*0.99)))
            elif era == "2018":
                # Double Electron #
                self.lambda_ttH_doubleElectron_trigSF = lambda dilep : op.multiSwitch(
                    (self.electron_conept[dilep[1].idx]<25 , op.systematic(op.c_float(0.98), name="ttH_doubleElectron_trigSF", up=op.c_float(0.98*1.01), down=op.c_float(0.98*0.99))), 
                    op.systematic(op.c_float(1.), name="ttH_doubleElectron_trigSF", up=op.c_float(1.01), down=op.c_float(0.99)))
                # Double Muon /!\ only this one using leading conept lepton #
                self.lambda_ttH_doubleMuon_trigSF = lambda dilep : op.multiSwitch(
                    (self.muon_conept[dilep[0].idx]<40 , op.systematic(op.c_float(1.01), name="ttH_doubleMuon_trigSF", up=op.c_float(1.01*1.01), down=op.c_float(1.01*0.99))),
                    (self.muon_conept[dilep[0].idx]<70 , op.systematic(op.c_float(0.995), name="ttH_doubleMuon_trigSF", up=op.c_float(0.995*1.01), down=op.c_float(0.995*0.99))),
                    op.systematic(op.c_float(0.98), name="ttH_doubleMuon_trigSF", up=op.c_float(0.98*1.01), down=op.c_float(0.98*0.99)))
                # Electron Muon : /!\ While ElEl and MuMu is conept ordered, ElMu is not #
                self.lambda_ttH_electronMuon_trigSF = lambda dilep : op.multiSwitch(
                    (op.min(self.electron_conept[dilep[0].idx],self.muon_conept[dilep[1].idx])<25, op.systematic(op.c_float(0.98), name="ttH_electronMuon_trigSF", up=op.c_float(0.98*1.01), down=op.c_float(0.98*0.99))),
                    op.systematic(op.c_float(1.00), name="ttH_electronMuon_trigSF", up=op.c_float(1.01), down=op.c_float(0.99)))
        #----- Fake rates -----#
        FRSysts = [f'Loose_{channel}_pt_syst',f'Loose_{channel}_barrel_syst',f'Loose_{channel}_norm_syst']
        self.electronFRList = [self.SF.get_scalefactor("lepton", ('electron_fakerates_'+era, syst), combine="weight", systName="el_FR_"+syst, 
                                             additionalVariables={'Pt' : lambda obj : self.electron_conept[obj.idx]}) for syst in FRSysts]
        self.muonFRList = [self.SF.get_scalefactor("lepton", ('muon_fakerates_'+era, syst), combine="weight", systName="mu_FR_"+syst, 
                                         additionalVariables={'Pt' : lambda obj : self.muon_conept[obj.idx]}) for syst in FRSysts ] 
        if channel == 'SL': # Not needed for DL
            self.electronFRNC = self.SF.get_scalefactor("lepton", ("fakerates_nonclosure_{}".format(self.era),'Loose_Electron_SL_{}'.format(self.era)), combine="weight",
                                                        defineOnFirstUse=False,
                                                        additionalVariables={'AbsEta': lambda x : op.abs(x.eta),'Pt': lambda x : x.pt})
            self.muonFRNC     = self.SF.get_scalefactor("lepton", ("fakerates_nonclosure_{}".format(self.era),'Loose_Muon_SL_{}'.format(self.era)), combine="weight",
                                                        defineOnFirstUse=False,
                                                        additionalVariables={'AbsEta': lambda x : op.abs(x.eta),'Pt': lambda x : x.pt})

        def returnFFSF(obj,list_SF,systName):
            """ Helper when several systematics are present  """
            args = [ a(obj) for a in list_SF[0]._args ] ## get the things the SF depends on
            systArgs = {'nominal':list_SF[0].sfOp.get(*(args+[op.extVar("int", "Nominal")])),'name':systName}
            systArgs.update({SF._systName+'up':SF.sfOp.get(*(args+[op.extVar("int", "Up")])) for SF in list_SF})
            systArgs.update({SF._systName+'down':SF.sfOp.get(*(args+[op.extVar("int", "Down")])) for SF in list_SF})
            return op.systematic(**systArgs)
             
        #----- DY reweighting -----#
        mode = 'data'
        if hasattr(self,'datadrivenContributions') and "PseudoData" in self.datadrivenContributions:
            mode = 'mc'
        # Resolved 
        self.ResolvedDYReweighting1b = self.SF.get_scalefactor("lepton", ('DY_resolved_{}'.format(era),'SF_HT_{}_1b'.format(mode)), combine="weight", 
                                                              systName="dy_resolved_1b", 
                                                              defineOnFirstUse=(not forSkimmer),
                                                              additionalVariables={'Eta': lambda x : op.c_float(0.),'Pt': lambda jets: op.rng_sum(jets, lambda j : j.pt)})
        self.ResolvedDYReweighting2b = self.SF.get_scalefactor("lepton", ('DY_resolved_{}'.format(era),'SF_HT_{}_2b'.format(mode)), combine="weight", 
                                                              systName="dy_resolved_2b", 
                                                              defineOnFirstUse=(not forSkimmer),
                                                              additionalVariables={'Eta': lambda x : op.c_float(0.),'Pt': lambda jets: op.rng_sum(jets, lambda j : j.pt)})
        # Boosted
        self.BoostedDYReweighting1b =  self.SF.get_scalefactor("lepton", ('DY_boosted_{}'.format(era),'SF_fatjetsoftDropmass_{}_1b'.format(mode)), combine="weight", 
                                                              systName="dy_boosted_1b", 
                                                              defineOnFirstUse=(not forSkimmer),
                                                              additionalVariables={'Eta': lambda x : op.c_float(0.),'Pt': lambda x: x.msoftdrop})


        #############################################################################
        #                             High level lambdas                            #
        #############################################################################
        self.HLL = highlevelLambdas(self)
        
        #############################################################################
        #                             Fake Factors                                  #
        #############################################################################
        # Non closure between ttbar and QCD correction #
        # https://gitlab.cern.ch/cms-hh-bbww/cms-hh-to-bbww/-/blob/master/Legacy/backgroundEstimation.md
        if era == '2016':
            self.electronCorrFR = op.systematic(op.c_float(1.376), name="electronCorrFR",up=op.c_float(1.376*1.376),down=op.c_float(1.))
            self.muonCorrFR     = op.systematic(op.c_float(1.050), name="muonCorrFR",up=op.c_float(1.050*1.050),down=op.c_float(1.))
        elif era == '2017':
            self.electronCorrFR = op.systematic(op.c_float(1.252), name="electronCorrFR",up=op.c_float(1.252*1.252),down=op.c_float(1.))
            self.muonCorrFR     = op.systematic(op.c_float(1.157), name="muonCorrFR",up=op.c_float(1.157*1.157),down=op.c_float(1.))
        elif era == '2018':
            self.electronCorrFR = op.systematic(op.c_float(1.325), name="electronCorrFR",up=op.c_float(1.325*1.325),down=op.c_float(1.))
            self.muonCorrFR     = op.systematic(op.c_float(1.067), name="muonCorrFR",up=op.c_float(1.067*1.067),down=op.c_float(1.))

        self.lambda_FR_el = lambda el : returnFFSF(el,self.electronFRList,"el_FR")
        self.lambda_FR_mu = lambda mu : returnFFSF(mu,self.muonFRList,"mu_FR")


        if channel == "SL":
            self.ElFakeFactor = lambda el : self.lambda_FR_el(el)/(1-self.lambda_FR_el(el))
            self.MuFakeFactor = lambda mu : self.lambda_FR_mu(mu)/(1-self.lambda_FR_mu(mu))
            self.ElFakeFactorNonClosure = lambda el : self.electronFRNC(el)/(1-self.electronFRNC(el))
            self.MuFakeFactorNonClosure = lambda mu : self.muonFRNC(mu)/(1-self.muonFRNC(mu))
        if channel == "DL":
            # Correction factor only for DL #
            self.lambda_FRcorr_el   = lambda el : self.lambda_FR_el(el)*self.electronCorrFR
            self.lambda_FRcorr_mu   = lambda mu : self.lambda_FR_mu(mu)*self.muonCorrFR

            self.lambda_FF_el   = lambda el : self.lambda_FRcorr_el(el)/(1-self.lambda_FRcorr_el(el))
            self.lambda_FF_mu   = lambda mu : self.lambda_FRcorr_mu(mu)/(1-self.lambda_FRcorr_mu(mu))
            self.ElElFakeFactor = lambda dilep : op.multiSwitch((op.AND(op.NOT(self.lambda_electronTightSel(dilep[0])),op.NOT(self.lambda_electronTightSel(dilep[1]))),
                                                                 # Both electrons fail tight -> -F1*F2
                                                                 -self.lambda_FF_el(dilep[0])*self.lambda_FF_el(dilep[1])),
                                                                (op.AND(op.NOT(self.lambda_electronTightSel(dilep[0])),self.lambda_electronTightSel(dilep[1])),
                                                                 # Only leading electron fails tight -> F1
                                                                 self.lambda_FF_el(dilep[0])),
                                                                (op.AND(self.lambda_electronTightSel(dilep[0]),op.NOT(self.lambda_electronTightSel(dilep[1]))),
                                                                 # Only subleading electron fails tight -> F2
                                                                 self.lambda_FF_el(dilep[1])),
                                                                 op.c_float(1.)) # Should not happen
            self.MuMuFakeFactor = lambda dilep : op.multiSwitch((op.AND(op.NOT(self.lambda_muonTightSel(dilep[0])),op.NOT(self.lambda_muonTightSel(dilep[1]))),
                                                                 # Both muons fail tight -> -F1*F2
                                                                 -self.lambda_FF_mu(dilep[0])*self.lambda_FF_mu(dilep[1])),
                                                                (op.AND(op.NOT(self.lambda_muonTightSel(dilep[0])),self.lambda_muonTightSel(dilep[1])),
                                                                 # Only leading muon fails tight -> F1
                                                                 self.lambda_FF_mu(dilep[0])),
                                                                (op.AND(self.lambda_muonTightSel(dilep[0]),op.NOT(self.lambda_muonTightSel(dilep[1]))),
                                                                 # Only subleading muon fails tight -> F2
                                                                 self.lambda_FF_mu(dilep[1])),
                                                                 op.c_float(1.)) # Should not happen
            self.ElMuFakeFactor = lambda dilep : op.multiSwitch((op.AND(op.NOT(self.lambda_electronTightSel(dilep[0])),op.NOT(self.lambda_muonTightSel(dilep[1]))),
                                                                 # Both electron and muon fail tight -> -F1*F2
                                                                 -self.lambda_FF_el(dilep[0])*self.lambda_FF_mu(dilep[1])),
                                                                (op.AND(op.NOT(self.lambda_electronTightSel(dilep[0])),self.lambda_muonTightSel(dilep[1])),
                                                                 # Only electron fails tight -> F1
                                                                 self.lambda_FF_el(dilep[0])),
                                                                (op.AND(self.lambda_electronTightSel(dilep[0]),op.NOT(self.lambda_muonTightSel(dilep[1]))),
                                                                 # Only subleading electron fails tight -> F2
                                                                 self.lambda_FF_mu(dilep[1])),
                                                                 op.c_float(1.)) # Should not happen

        ###########################################################################
        #                    b-tagging efficiency scale factors                   #
        ###########################################################################
        if self.is_MC:
            if self.era == '2016' and (('db' in sampleCfg.keys() and not 'TuneCP5' in sampleCfg['db']) or ('files' in sampleCfg.keys() and not any(['TuneCP5' in f for f in sampleCfg['files']]))):
                csvFileNameAk4 = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data", "ScaleFactors_POG" , "deepjet_2016_oldtune.csv")
            else:
                csvFileNameAk4 = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data", "ScaleFactors_POG" , "deepjet_{}.csv".format(self.era))
            print ('Btag Ak4 CSV file',csvFileNameAk4)

            #----- Ak4 SF -----# 
            self.systTypes = ["hf","lf","hfstats1", "hfstats2","lfstats1", "lfstats2", "cferr1", "cferr2"]
            if era == "2016":
                self.systTypes += ['jesBBEC1','jesFlavorQCD','jesHF_2016','jesRelativeSample_2016','jesEC2','jesAbsolute_2016','jesAbsolute','jesRelativeBal','jesEC2_2016','jesBBEC1_2016','jesHF']
            if era == "2017":
                self.systTypes += ['jesBBEC1','jesHF_2017','jesFlavorQCD','jesEC2_2017','jesBBEC1_2017','jesEC2','jesAbsolute','jesRelativeBal','jesAbsolute_2017','jesRelativeSample_2017','jesHF']
            if era == "2018":
                self.systTypes += ['jesRelativeSample_2018','jesHEMIssue','jesBBEC1','jesAbsolute_2018','jesAbsolute','jesFlavorQCD','jesRelativeBal','jesEC2_2018','jesBBEC1_2018','jesEC2','jesHF_2018','jesHF'] 
                
            self.systTypes = ['up_'+s for s in self.systTypes]+['down_'+s for s in self.systTypes]

            def translate_btagSFVarToJECVar_invert(btagVarName, prefix="btagSF_"):
                if btagVarName.startswith("up_jes") or btagVarName.startswith("down_jes"):
                    if btagVarName.endswith("_jes"):
                        return "jesTotal{0}".format(btagVarName.split("_jes")[0])
                    elif "RelativeBal" in btagVarName or "RelativeSample_2016" in btagVarName:
                        logger.warning("B-tag SF systematic variation {0} has been inverted".format(btagVarName))
                        if btagVarName.startswith("up_jes"):
                            btagVarName = btagVarName.replace("up_jes","down_jes")
                        elif btagVarName.startswith("down_jes"):
                            btagVarName = btagVarName.replace("down_jes","up_jes")
                    return "jes{1}{0}".format(*btagVarName.split("_jes"))
                elif btagVarName.startswith("up_") or btagVarName.startswith("down_"):
                    tk = btagVarName.split("_")
                    return "".join((prefix, "_".join(tk[1:]), tk[0]))
                elif btagVarName in ("up", "down"):
                    return "".join((prefix, btagVarName))
                else:
                    return btagVarName
                

                # From https://twiki.cern.ch/twiki/bin/view/CMS/BTagShapeCalibration#Systematic_uncertainties
            self.DeepJetDiscReshapingSF = BtagSF(taggerName       = "deepjet", 
                                                 csvFileName      = csvFileNameAk4,
                                                 wp               = "reshaping",  # "loose", "medium", "tight" or "reshaping"
                                                 sysType          = "central", 
                                                 otherSysTypes    = self.systTypes,
                                                 measurementType  = "iterativefit", 
                                                 sel              = noSel, 
                                                 getters          = {'Discri':lambda j : j.btagDeepFlavB},
                                                 jesTranslate     = translate_btagSFVarToJECVar_invert,
                                                 uName            = self.sample)

        if era == '2016':
            csvFileNameAk8= os.path.join(os.path.abspath(os.path.dirname(__file__)), "data", "ScaleFactors_POG" , "subjet_DeepCSV_2016LegacySF_V1.csv")
        if era == '2017':
            csvFileNameAk8 = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data", "ScaleFactors_POG" , "subjet_DeepCSV_94XSF_V4_B_F_v2.csv")
        if era == '2018':
            csvFileNameAk8 = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data", "ScaleFactors_POG" , "subjet_DeepCSV_102XSF_V1.csv")
        
        #----- Ak8 SF -----# 
        if not os.path.exists(csvFileNameAk8):
            raise RuntimeError('Could not find Ak8 csv file %s'%csvFileNameAk8)
        print ('Btag Ak8 CSV file',csvFileNameAk8)
        self.DeepCsvSubjetMediumSF = BtagSF(taggerName       = "deepcsvSubjet", 
                                            csvFileName      = csvFileNameAk8,
                                            wp               = "medium",
                                            sysType          = "central", 
                                            otherSysTypes    = ['up','down'],
                                            measurementType  = {"B": "lt", "C": "lt", "UDSG": "incl"},
                                            getters          = {'Discri':lambda subjet : subjet.btagDeepB,
                                                                'JetFlavour': lambda subjet : op.static_cast("BTagEntry::JetFlavor",
                                                                                                    op.multiSwitch((subjet.nBHadrons>0,op.c_int(0)), # B -> flav = 5 -> BTV = 0
                                                                                                                   (subjet.nCHadrons>0,op.c_int(1)), # C -> flav = 4 -> BTV = 1
                                                                                                                   op.c_int(2)))},                  # UDSG -> flav = 0 -> BTV = 2
                                            sel              = noSel, 
                                            uName            = self.sample)


        ###########################################################################
        #                                 RETURN                                  # 
        ###########################################################################
        return noSel



    def beforeJetselection(self,sel,name=''):
        if isinstance(sel,SelectionWithDataDriven) or isinstance(sel,Selection):
            is_refine = True
        elif isinstance(sel,SelectionObject):
            is_refine = False
            selObj = sel
        else:
            raise RuntimeError(f"Could not understand type of selection {type(sel)}")

        ##############################################################################
        #                             Jets forceDefines                              #
        ##############################################################################
        from bamboo.analysisutils import forceDefine
        # Forcedefine : calculate once per event (for every event) #
        if not self.args.Synchronization:
            if is_refine:
                #forceDefine(t._Muon.calcProd, sel) # Muons for Rochester corrections
                forceDefine(self.tree._Jet.calcProd, sel)  # Jets for configureJets
                forceDefine(self.tree._FatJet.calcProd, sel)  # FatJets for configureJets
                forceDefine(getattr(self.tree, "_{0}".format("MET" if self.era != "2017" else "METFixEE2017")).calcProd,sel) # MET for configureMET
            else:
                #forceDefine(t._Muon.calcProd, sel) # Muons for Rochester corrections
                forceDefine(self.tree._Jet.calcProd, selObj.sel)  # Jets for configureJets
                forceDefine(self.tree._FatJet.calcProd, selObj.sel)  # FatJets for configureJets
                forceDefine(getattr(self.tree, "_{0}".format("MET" if self.era != "2017" else "METFixEE2017")).calcProd,selObj.sel) # MET for configureMET
        else:
            print ("No jet corrections applied")

        #############################################################################
        #                           Jet PU ID reweighting                           #
        #############################################################################
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID
        if self.is_MC and not self.args.BtagReweightingOff and not self.args.BtagReweightingOn:
            wFail = op.extMethod("scalefactorWeightForFailingObject", returnType="double")
            # method 1.a : https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#1a_Event_reweighting_using_scale
            lambda_puid_weight = lambda j : op.switch(j.genJet.isValid, 
                                                      # Is gen jet
                                                      op.switch(((j.puId >> 2) & 1), # passes jet pu id cut
                                                                self.jetpuid_sf_eff(j), # SF_eff
                                                                wFail(self.jetpuid_sf_eff(j), self.jetpuid_mc_eff(j))), # (1-SF_eff*eff) / (1-eff)
                                                      # Is PU jet
                                                      op.switch(((j.puId >> 2) & 1), # passed jet pu id cut
                                                                self.jetpuid_sf_mis(j), # SF_mis
                                                                wFail(self.jetpuid_sf_mis(j), self.jetpuid_mc_mis(j)))) # (1-SF_mis*mis) / (1-mis)
            # Eff and mistag kept separate here for the sync ntuple #
            lambda_puid_efficiency = lambda j : op.switch(j.genJet.isValid,
                                                          op.switch(((j.puId >> 2) & 1),
                                                                    self.jetpuid_sf_eff(j),
                                                                    wFail(self.jetpuid_sf_eff(j), self.jetpuid_mc_eff(j))),
                                                          op.c_float(1.))
            lambda_puid_mistag     = lambda j : op.switch(j.genJet.isValid,
                                                          op.c_float(1.),
                                                          op.switch(((j.puId >> 2) & 1),
                                                                    self.jetpuid_sf_mis(j), 
                                                                    wFail(self.jetpuid_sf_mis(j), self.jetpuid_mc_mis(j))))

            self.puid_reweighting = op.rng_product(self.ak4ForPUID, lambda j : lambda_puid_weight(j))
            # Variables used in the skims for sync ntuples purposes
            self.puid_reweighting_efficiency = op.rng_product(self.ak4ForPUID, lambda j : lambda_puid_efficiency(j))
            self.puid_reweighting_mistag = op.rng_product(self.ak4ForPUID, lambda j : lambda_puid_mistag(j))
            if is_refine: 
                sel = sel.refine("jetPUIDReweighting"+name,weight=self.puid_reweighting)
            else:
                selObj.selName += "jetPUIDReweighting"
                selObj.yieldTitle += " + Jet PU ID reweighting "
                selObj.refine(weight=self.puid_reweighting)
        else:
            if not is_refine:
                selObj.selName += "jetPUIDReweighting"

        ###########################################################################
        #                    b-tagging efficiency scale factors                   #
        ###########################################################################
        if self.is_MC:
            #----- AK4 jets -> using Method 1.d -----#
            # See https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagShapeCalibration
            # W_btag = Π_i(all jets) SD(jet_i)  which must be multiplied by r = Σ w(before)/Σ w(after) (before/after using the btag weight, no btag selection for both)
                
            self.btagAk4SF = op.rng_product(self.ak4Jets , lambda j : op.switch(j.hadronFlavour == 4,
                                                                                self.DeepJetDiscReshapingSF(j, systVars = [systType for systType in self.systTypes if "cferr" in systType]),
                                                                                self.DeepJetDiscReshapingSF(j, systVars = [systType for systType in self.systTypes if not "cferr" in systType])))

            if self.args.BtagReweightingOn and self.args.BtagReweightingOff: 
                raise RuntimeError("Reweighting cannot be both on and off") 

            if self.args.BtagReweightingOn:
                sel = sel.refine("BtagSF" , weight = self.btagAk4SF)
            elif self.args.BtagReweightingOff:
                pass # Do not apply any SF
            else:
                if 'related-sample' in self.sampleCfg.keys():
                    ReweightingFileName = os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','ScaleFactors_Btag',f'BtagReweightingRatio_jetN_{self.sampleCfg["related-sample"]}_{self.era}.json')
                    nameHint = f'bamboo_nJetsWeight_{self.sampleCfg["related-sample"]}'.replace('-','_')
                else:
                    ReweightingFileName = os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','ScaleFactors_Btag',f'BtagReweightingRatio_jetN_{self.sample}_{self.era}.json')
                    nameHint = f'bamboo_nJetsWeight_{self.sample}'.replace('-','_')
                if not os.path.exists(ReweightingFileName):
                    raise RuntimeError("Could not find reweighting file %s"%ReweightingFileName)
                print ('Reweighting file',ReweightingFileName)

                self.BtagRatioWeight = makeBtagRatioReweighting(jsonFile = ReweightingFileName,
                                                                numJets  = op.rng_len(self.ak4Jets),
                                                                nameHint = nameHint)

                if is_refine: 
                    sel = sel.refine("BtagAk4SF"+name , weight = [self.btagAk4SF,self.BtagRatioWeight])
                else:
                    selObj.selName += "BtagAk4SF"
                    selObj.yieldTitle += " + Btag Ak4 SF "
                    selObj.refine(weight = [self.btagAk4SF,self.BtagRatioWeight])

            if not self.args.BtagReweightingOff and not self.args.BtagReweightingOn:
                #----- AK8 jets -> using Method 1.a -----#
                # Reweighting #
                wFail = op.extMethod("scalefactorWeightForFailingObject", returnType="double") # 
                # double scalefactorWeightForFailingObject(double sf, double eff) {
                #  return (1.-sf*eff)/(1.-eff);
                #  }

                # Method 1.a
                # P(MC) = Π_{i tagged} eff_i Π_{j not tagged} (1-eff_j)
                # P(data) = Π_{i tagged} SF_i x eff_i Π_{j not tagged} (1-SF_jxeff_j)
                # w = P(data) / P(MC) = Π_{i tagged} SF_i Π_{j not tagged} (1-SF_jxeff_j)/(1-eff_j)
                # NB : for  SF_i, self.DeepCsvSubjetMediumSF will find the correct SF based on the true flavour
                #      however for eff_i, this comes from json SF and differs by flavour -> need multiSwitch
                lambda_subjetWeight = lambda subjet : op.multiSwitch((subjet.nBHadrons>0,      # True bjets 
                                                                      op.switch(self.lambda_subjetBtag(subjet),                                         # check if tagged
                                                                                self.DeepCsvSubjetMediumSF(subjet),                                     # Tag : return SF_i
                                                                                wFail(self.DeepCsvSubjetMediumSF(subjet),self.Ak8Eff_bjets(subjet)))),  # Not tagged : return (1-SF_jxeff_j)/(1-eff_j)
                                                                     (subjet.nCHadrons>0,      # True cjets 
                                                                      op.switch(self.lambda_subjetBtag(subjet),                                         # check if tagged
                                                                                self.DeepCsvSubjetMediumSF(subjet),                                     # Tag : return SF_i
                                                                                wFail(self.DeepCsvSubjetMediumSF(subjet),self.Ak8Eff_cjets(subjet)))),  # Not tagged : return (1-SF_jxeff_j)/(1-eff_j)
                                                                      # Else : true lightjets 
                                                                      op.switch(self.lambda_subjetBtag(subjet),                                               # check if tagged
                                                                                self.DeepCsvSubjetMediumSF(subjet),                                           # Tag : return SF_i
                                                                                wFail(self.DeepCsvSubjetMediumSF(subjet),self.Ak8Eff_lightjets(subjet))))     # Not tagged : return (1-SF_jxeff_j)/(1-eff_j)
                self.ak8BtagReweighting = op.rng_product(self.ak8Jets, lambda j : lambda_subjetWeight(j.subJet1)*lambda_subjetWeight(j.subJet2))
                if is_refine: 
                    sel = sel.refine("BtagAk8SF"+name , weight = [self.ak8BtagReweighting])
                else:
                    selObj.selName += "BtagAk8SF"
                    selObj.yieldTitle += " + Btag Ak8 SF "
                    selObj.refine(weight = [self.ak8BtagReweighting])
        else:
            if not is_refine:
                selObj.selName += "BtagAk4SFBtagAk8SF"
            
        
        if is_refine:
            return sel
       



    ###########################################################################
    #                            Lepton triggers                              #
    ###########################################################################
    def returnTriggers(self,keys):
        triggerRanges = returnTriggerRanges(self.era)
        return op.OR(*[trig for k in keys for trig in self.triggersPerPrimaryDataset[k]])

    ###########################################################################
    #                          Add Signal Reweighting                         #
    ###########################################################################
    def addSignalReweighting(self, selection):
        """
        Signal reweighting

        For the Plotter we want the many-to-many worfklow.
        This means create as many files (with the create) 
        as there are NLO weight files (+ original LO file)
        """
        baseName = selection.name
        for benchmark, reweightLO in self.signalReweightBenchmarks.items():
            selection = SelectionWithDataDriven.create(
                    parent = selection,
                    name = f"{baseName}{benchmark}",
                    ddSuffix = benchmark,
                    ddWeight = reweightLO,
                    enable = True
                    )
        return selection

    ###########################################################################
    #                             postProcess                                 #
    ###########################################################################
    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None, forSkimmer=True):
        from bamboo.root import gbl

        if not forSkimmer:
            if not self.plotList:
                self.plotList = self.getPlotList(resultsdir=resultsdir, config=config)
            ### Disable plotting for some signal samples ###
            samplesToHide = [sample for sample,sampleCfg in config['samples'].items() if 'hide' in sampleCfg.keys() and sampleCfg['hide']]
            for sample in samplesToHide:
                config['samples'][sample]['line-color'] = "#01000000" # transparent 
                config['samples'][sample]['legend'] = ""

            ### HH reweighting ###
            if self.args.HHReweighting is not None:
                print ('Starting HH reweighting postprocessing')
                LOFiles = []
                channels = ['GluGluToHHTo2B2VTo2L2Nu','GluGluToHHTo2B2WToLNu2J','GluGluToHHTo2B2VLNu2J','GluGluToHHTo2B2Tau']
                    # bbww1l signal samples have different naming depending on LO/NLO
                ddScenarios = set(self.datadrivenScenarios)
                attrToKeep = ['legend','line-color','line-type','order','type']
                for contribName in self.analysisConfig['benchmarks']['targets']:
                    print (contribName)
                    if contribName not in self.args.HHReweighting:
                        continue
                    contrib = self.datadrivenContributions[contribName]
                    contribSamples = {sample:sampleCfg for sample,sampleCfg in config['samples'].items() if contrib.usesSample(sample,sampleCfg)}
                    LOFiles.extend([sample for sample in contribSamples.keys() if sample not in LOFiles])
                    _, eras = self.args.eras
                    if eras is None:
                        eras = list(config["eras"].keys())
                    for era in eras:
                        for channel in channels:
                            contribPerChannel = {sample:sampleCfg for sample,sampleCfg in contribSamples.items() if channel in sample and sampleCfg['era']==era}
                            if len(contribPerChannel) == 0:
                                continue
                            crossSection = [sampleCfg['cross-section'] for sampleCfg in contribPerChannel.values()]
                            if len(set(crossSection)) != 1:
                                print(f"Not all samples for channel {channel} have the same cross-section, will use 1. as default")
                                crossSection = 1.
                            else:
                                crossSection = crossSection[0]
                            if all(['branching-ratio' in sampleCfg for sampleCfg in contribPerChannel.values()]):
                                branchingRatio = [sampleCfg['branching-ratio'] for sampleCfg in contribPerChannel.values()]
                                if len(set(branchingRatio)) != 1:
                                    raise RuntimeError(f"Not all samples for channel {channel} have the same branching-ratio")
                                else:
                                    branchingRatio = branchingRatio[0]
                            else:
                                branchingRatio = None 
                            attrLO = {attr:contribPerChannel[list(contribPerChannel.keys())[0]][attr]  for attr in attrToKeep 
                                        if attr in contribPerChannel[list(contribPerChannel.keys())[0]].keys()}  
                            files = {sample:gbl.TFile.Open(os.path.join(resultsdir,f"{sample}{contribName}.root")) for sample in contribPerChannel.keys()}
                            generatedEvents = {sample:self.readCounters(files[sample])[contribPerChannel[sample]['generated-events']] for sample in contribPerChannel.keys()}
                            outPath = os.path.join(resultsdir,f"{channel}_NLO{contribName}.root")
                            if os.path.exists(outPath):
                                print (f'Aggregated file {outPath} already exists, will not recompute')
                            else:
                                outFile = gbl.TFile.Open(outPath, "RECREATE")
                                hNames = [key.GetName() for key in files[list(files.keys())[0]].GetListOfKeys() if 'TH' in key.GetClassName()]
                                for hName in hNames:
                                    hist = None
                                    for sample,f in files.items():
                                        h = f.Get(hName)
                                        h.Scale(1./(generatedEvents[sample]))
                                        if hist is None:
                                            hist = copy.deepcopy(h)
                                        else:
                                            hist.Add(h)
                                    hist.Write()
                                outFile.Close()
                                print (f'Aggregated file {outPath} produced')
                                for f in files.values():
                                    f.Close()
                            
                            config['samples'][f"{channel}_NLO{contribName}"] = {'cross-section'     : crossSection,
                                                                                'era'               : era,
                                                                                'generated-events'  : len(generatedEvents)}
                            config['samples'][f"{channel}_NLO{contribName}"].update(attrLO)
                            if branchingRatio is not None:
                                config['samples'][f"{channel}_NLO{contribName}"]['branching-ratio'] = branchingRatio

                        for sc in list(ddScenarios):
                            if contribName in sc:
                                newSc = tuple(cb for cb in sc if cb != contribName)
                                ddScenarios.remove(sc)
                                ddScenarios.add(newSc)

                self.datadrivenScenarios = list(ddScenarios)

                for LOFile in LOFiles:
                    del config['samples'][LOFile]

        super(BaseNanoHHtobbWW,self).postProcess(taskList=taskList, config=config, workdir=workdir, resultsdir=resultsdir)

    ###########################################################################
    #                                  HME                                    #
    ###########################################################################
    def computeResolvedHMEAfterLeptonSelections(self, sel, l1, l2, bjets, met):
        if self.args.analysis != 'res':
            raise RuntimeError('HME was only implemented for resonant analysis')

        hme_pair = self.hmeEval.runHME(l1.p4, l2.p4, bjets[0].p4, bjets[1].p4, met.p4, self.tree.event, op.c_bool(False))
        hme_pair = op.switch(op.rng_len(bjets) >= 2, 
                             hme_pair, 
                             op.construct(op.typeOf(hme_pair), [op.c_float(0.), op.c_float(0.)]))

        hme_pair = op.forSystematicVariation(hme_pair, "nominal")   # no variations, always nominal
        forceDefine(hme_pair, sel)
        hme = hme_pair.first
        eff = hme_pair.second
        forceDefine(hme, sel)
        forceDefine(eff, sel)

        return hme,eff

    def computeBoostedHMEAfterLeptonSelections(self, sel, l1, l2, fatjets, met):
        if self.args.analysis != 'res':
            raise RuntimeError('HME was only implemented for resonant analysis')

        empty_p4 = op.construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >",([op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.)]))

        hme_pair = self.hmeEval.runHME(l1.p4, l2.p4, fatjets[0].subJet1.p4, fatjets[0].subJet2.p4, met.p4, self.tree.event, op.c_bool(True))
        hme_pair = op.switch(op.rng_len(fatjets) >= 1, 
                             hme_pair, 
                             op.construct(op.typeOf(hme_pair), [op.c_float(0.), op.c_float(0.)]))

        hme_pair = op.forSystematicVariation(hme_pair, "nominal")   # no variations, always nominal
        forceDefine(hme_pair, sel)
        hme = hme_pair.first
        eff = hme_pair.second
        forceDefine(hme, sel)
        forceDefine(eff, sel)

        return hme,eff
