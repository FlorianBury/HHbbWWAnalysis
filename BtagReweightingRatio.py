import os
import sys
from copy import copy

from itertools import chain

import logging
logger = logging.getLogger(__name__) 

import bamboo
from bamboo.analysismodules import HistogramsModule, DataDrivenBackgroundHistogramsModule

from bamboo import treefunctions as op
from bamboo.plots import CutFlowReport, Plot, EquidistantBinning, SummedPlot

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)))) # Add scripts in this directory
from BaseHHtobbWW import BaseNanoHHtobbWW
from plotDef import *
from selectionDef import *

#===============================================================================================#
#                                   BtagReweightingRatioNano                                    #
#===============================================================================================#
class BtagReweightingRatioNano(BaseNanoHHtobbWW,DataDrivenBackgroundHistogramsModule):
    def __init__(self, args):
        super(BtagReweightingRatioNano, self).__init__(args)

    def initialize(self):
        super(BtagReweightingRatioNano, self).initialize(forSkimmer=True)

    def definePlots(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(BtagReweightingRatioNano,self).prepareObjects(t, noSel, sample, sampleCfg, 'DL') 
            # Using DL tag for SL plots should not matter, as there is no jet selection
            # Cannot use SL because some lepton quantities (eg dileptons) are only defined there #

        plots = []

        era = sampleCfg['era']

        noSel = self.beforeJetselection(noSel,'noSel')

        ElSelObj,MuSelObj                = makeSingleLeptonSelection(self,noSel,use_dd=False)
        ElElSelObj,MuMuSelObj,ElMuSelObj = makeDoubleLeptonSelection(self,noSel,use_dd=False)

        ElSelObj.sel    = self.beforeJetselection(ElSelObj.sel,'El')
        MuSelObj.sel    = self.beforeJetselection(MuSelObj.sel,'Mu')
        ElElSelObj.sel  = self.beforeJetselection(ElElSelObj.sel,'ElEl')
        MuMuSelObj.sel  = self.beforeJetselection(MuMuSelObj.sel,'MuMu')
        ElMuSelObj.sel  = self.beforeJetselection(ElMuSelObj.sel,'ElMu')

        #----- Ratio reweighting variables (before lepton and jet selection) -----#
        if self.args.BtagReweightingOff or self.args.BtagReweightingOn:
            plots.append(objectsNumberPlot(channel="NoChannel",suffix='NoSelection',sel=noSel,objCont=self.ak4Jets,objName='Ak4Jets',Nmax=15,xTitle='N(Ak4 jets)'))
            plots.append(objectsNumberPlot(channel="El",suffix='El',sel=ElSelObj.sel,objCont=self.ak4Jets,objName='Ak4Jets',Nmax=15,xTitle='N(Ak4 jets)'))
            plots.append(objectsNumberPlot(channel="Mu",suffix='Mu',sel=MuSelObj.sel,objCont=self.ak4Jets,objName='Ak4Jets',Nmax=15,xTitle='N(Ak4 jets)'))
            plots.append(objectsNumberPlot(channel="ElEl",suffix='ElEl',sel=ElElSelObj.sel,objCont=self.ak4Jets,objName='Ak4Jets',Nmax=15,xTitle='N(Ak4 jets)'))
            plots.append(objectsNumberPlot(channel="MuMu",suffix='MuMu',sel=MuMuSelObj.sel,objCont=self.ak4Jets,objName='Ak4Jets',Nmax=15,xTitle='N(Ak4 jets)'))
            plots.append(objectsNumberPlot(channel="ElMu",suffix='ElMu',sel=ElMuSelObj.sel,objCont=self.ak4Jets,objName='Ak4Jets',Nmax=15,xTitle='N(Ak4 jets)'))
            plots.append(CutFlowReport("BtagReweightingCutFlowReport",noSel,printInLog=True,recursive=True))
            plots.append(CutFlowReport("BtagReweightingCutFlowReport",ElSelObj.sel,printInLog=True,recursive=True))
            plots.append(CutFlowReport("BtagReweightingCutFlowReport",MuSelObj.sel,printInLog=True,recursive=True))
            plots.append(CutFlowReport("BtagReweightingCutFlowReport",ElElSelObj.sel,printInLog=True,recursive=True))
            plots.append(CutFlowReport("BtagReweightingCutFlowReport",MuMuSelObj.sel,printInLog=True,recursive=True))
            plots.append(CutFlowReport("BtagReweightingCutFlowReport",ElMuSelObj.sel,printInLog=True,recursive=True))
        else:
            raise RuntimeError("Either --BtagReweightingOff or --BtagReweightingOn must be used")

        return plots
