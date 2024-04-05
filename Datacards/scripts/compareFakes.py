import os
import sys
import glob
import copy
import json 
import math
import argparse
import ctypes
from array import array
import numpy as np
import ROOT

from IPython import embed

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser(description='Comparing DY')
parser.add_argument('--mc', action='store', required=True, type=str,
                    help='Path to Fake MC shapes')
parser.add_argument('--cl', action='store', required=True, type=str,
                    help='Path to Fake Closure shapes')
parser.add_argument('--rebin', action='store', required=False, default=None, type=int,
                    help='Use the rebin number')
parser.add_argument('--name', action='store', required=True, type=str,
                    help='Name of the output files')
args = parser.parse_args()

factors = {}

plotDict_MC = {}
plotDict_CL = {}

path_MC  = args.mc
path_CL  = args.cl

plotNames = [os.path.basename(name) for name in glob.glob(os.path.join(path_MC,'*root'))]

def processHist(h,x):
    for i in range(1,h.GetNbinsX()+1):
        if h.GetBinContent(i) < 1e-3:
            h.SetBinContent(i,0.)
            h.SetBinError(i,0.)
    if isinstance(xn,int):
        return h.Rebin(xn)
    else:
        return h.Rebin(len(x)-1,'rebin',array('d',x))
     

for plotName in sorted(plotNames):
    f_MC = ROOT.TFile(os.path.join(path_MC,plotName))
    f_CL = ROOT.TFile(os.path.join(path_CL,plotName))

    plotName = plotName.replace('.root','')

    h_MC = f_MC.Get("TT")
    h_CL = f_CL.Get("TT")

    if args.rebin is not None:
        xn = args.rebin
    else:
        xn = np.linspace(0,1,21)

    h_MC = processHist(h_MC,xn)
    h_CL = processHist(h_CL,xn)

#    for i in range(1,h_MC.GetNbinsX()+1):
#        if h_MC.GetBinContent(i)<1e-3:
#            h_MC.SetBinContent(i,0.)
#    for i in range(1,h_CL.GetNbinsX()+1):
#        if h_CL.GetBinContent(i)<1e-3:
#            h_CL.SetBinContent(i,0.)
#        if h_MC.GetBinContent(i) == 0.:
#            h_CL.SetBinContent(i,0.)

    if plotName not in plotDict_MC.keys():
        plotDict_MC[plotName] = copy.deepcopy(h_MC)
    else:
        plotDict_MC[plotName].Add(h_MC)

    if plotName not in plotDict_CL.keys():
        plotDict_CL[plotName] = copy.deepcopy(h_CL)
    else:
        plotDict_CL[plotName].Add(h_CL)

    f_MC.Close()
    f_CL.Close()

def getY(h):
    return np.array([h.GetBinContent(i) for i in range(1,h.GetNbinsX()+1)])
def getX(h):
    return np.array([h.GetXaxis().GetBinCenter(i) for i in range(1,h.GetNbinsX()+1)])
def getEdges(h):
    return np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1,h.GetNbinsX()+2)])

categories = sorted(list(plotDict_MC.keys()))

C = ROOT.TCanvas("C","C",700,800)
pdfName = f"CompareFakes/ComparisonFakes_{args.name}.pdf"
C.Print(pdfName+"[")
data = {}
for cat in categories:
    C.Clear()

    h_MC = plotDict_MC[cat]
    h_CL = plotDict_CL[cat]


    Nerr_MC = ctypes.c_double(0.)
    Nerr_CL = ctypes.c_double(0.)
    N_MC = h_MC.IntegralAndError(0,h_MC.GetNbinsX()+1,Nerr_MC)
    N_CL = h_CL.IntegralAndError(0,h_CL.GetNbinsX()+1,Nerr_CL)
    Nerr_MC = Nerr_MC.value 
    Nerr_CL = Nerr_CL.value

    
    era = cat.split('_')[-1]
    title = ' '.join(cat.replace('HH_','').split('_')[:-1])
    h_MC.SetTitle(f"{title} category ({era});;Events / {h_MC.GetXaxis().GetBinWidth(1):0.2f}")
    hmax = max(h_MC.GetMaximum(),h_CL.GetMaximum())
    h_MC.GetYaxis().SetRangeUser(0,hmax*1.1)
    
    h_CL.SetLineColor(1)
    h_MC.SetLineColor(632)
    
    h_CL.SetMarkerColor(1)
    h_MC.SetMarkerColor(632)
    
    h_MC.SetMarkerStyle(25)
    h_CL.SetMarkerStyle(4)

    x_values = getX(h_MC)
    y_MC = getY(h_MC)
    y_CL = getY(h_CL)

    h_MC.GetXaxis().SetTitleSize(0.08)
    h_MC.GetYaxis().SetTitleSize(0.08)

    if N_CL == 0 or N_MC == 0:
        # Default params #
        nom   = 2
        slope = 1
        cog   = np.mean(x_values)

        # Plot only upper part (no ratio) #
        pad1 = C.cd(1)
        #pad1.SetLogy()
        pad1.SetTopMargin(0.1)
        pad1.SetBottomMargin(0.15)
        pad1.SetLeftMargin(0.15)
        pad1.SetRightMargin(0.05)
        pad1.Draw()

        
        h_MC.Draw("e1p")
        h_CL.Draw("e1p same")
        
        leg_up = ROOT.TLegend(0.55,0.55,0.9,0.9)
        leg_up.AddEntry(h_MC,f"#splitline{{Fakes (MC)}}{{Integral = {N_MC:0.3f} #pm {Nerr_MC:0.3f}}}")
        leg_up.AddEntry(h_CL,f"#splitline{{Fakes (Closure)}}{{Integral = {N_CL:0.3f} #pm {Nerr_CL:0.3f}}}")
        leg_up.SetBorderSize(0)
        leg_up.SetFillStyle(0)
        leg_up.SetTextSize(0.05)
        leg_up.Draw()


    else:
        cog = h_MC.GetMean() # Center of gravity
        nom = max(0,min(2,N_MC/N_CL))
        #h_MC.Scale(1./N_MC)
        #h_CL.Scale(1./N_CL)
        h_CL.Scale(N_MC/N_CL)

        
        ratio = h_MC.Clone("ratio")
        ratio.Divide(h_CL)
        ratio.SetTitle("")
        
        C.Divide(1,2)

        pad1 = C.cd(1)
        #pad1.SetLogy()
        pad1.SetTopMargin(0.1)
        pad1.SetBottomMargin(0.05)
        pad1.SetLeftMargin(0.15)
        pad1.SetRightMargin(0.05)
        pad1.Draw()
        
        h_MC.Draw("e1p")
        h_CL.Draw("e1p same")
        
        leg_up = ROOT.TLegend(0.55,0.55,0.9,0.9)
        leg_up.AddEntry(h_MC,f"#splitline{{Fakes (MC)}}{{Integral = {N_MC:0.3f} #pm {Nerr_MC:0.3f}}}")
        leg_up.AddEntry(h_CL,f"#splitline{{Fakes (Closure)}}{{Integral = {N_CL:0.3f} #pm {Nerr_CL:0.3f}}}")
        leg_up.SetBorderSize(0)
        leg_up.SetFillStyle(0)
        leg_up.SetTextSize(0.05)
        leg_up.Draw()
        
        pad2 = C.cd(2)
        pad2.SetTopMargin(0.0)
        pad2.SetBottomMargin(0.16)
        pad2.SetLeftMargin(0.15)
        pad2.SetRightMargin(0.05)
        pad2.Draw()

        # compute quantiles for the fit #
        prob = np.array([0.01,0.99])
        quantiles = np.array([0.,0.])
        h_MC.GetQuantiles(2,quantiles,prob)
        quantiles = quantiles-cog

        y_values = getY(ratio)
        edges_cog = getEdges(ratio)-cog
        ratio_cog = ROOT.TH1F("ratio_cog","ratio_cog",edges_cog.shape[0]-1,edges_cog)
        for i in range(1,ratio_cog.GetNbinsX()+1):
            ratio_cog.SetBinContent(i,ratio.GetBinContent(i))
            ratio_cog.SetBinError(i,ratio.GetBinError(i))
        ratio_cog.SetMarkerColor(1)
        ratio_cog.SetMarkerStyle(4)
        ratio_cog.SetLineColor(1)
        if getY(ratio_cog).sum() > 0:
            ratio_cog.SetTitle("")
            ratio_cog.GetYaxis().SetRangeUser(-1,3)
            ratio_cog.GetXaxis().SetTitle("DNN score (x-#bar{x})")
            ratio_cog.GetYaxis().SetTitle("Ratio")
            ratio_cog.GetXaxis().SetTitleSize(0.08)
            ratio_cog.GetYaxis().SetTitleSize(0.08)
            ratio_cog.Draw("e1p")

            ratio_cog.Fit("pol1","QRS","",quantiles[0],quantiles[1])
            fit = ratio_cog.GetFunction('pol1')
            fit.SetLineColor(1)
            slope = fit.GetParameter(1)
            k = fit.GetParameter(0)
        else:
            k = 0
            nom   = 2
            slope = 1
            cog   = np.mean(x_values)
        func = lambda x : slope*x+k

        print (f'Category {cat} : cog = {cog:5.3f}, nom = {nom:5.3f}, slope = {slope:5.3f}')

        leg_down = ROOT.TLegend(0.6,0.80,0.9,1.0)
        leg_down.AddEntry(ratio_cog,"\splitline{{#frac{{Fakes (MC)}}{{Fakes (Closure)}}}}{{Fit = {:.3f} x + {:.3f}}}".format(slope,k))
        leg_down.SetTextSize(0.04)
        leg_down.SetBorderSize(0)
        leg_down.SetFillStyle(0)
        leg_down.SetTextSize(0.05)
        leg_down.Draw()
        line = ROOT.TLine(min(edges_cog),1.,max(edges_cog),1.)
        line.Draw()
   
    # Print canvas and save data #
    C.Print(pdfName,f"Title:{cat}")
    if 'GGF' in cat or 'VBF' in cat:
        data[cat] = {'nom':nom}
    else:
        data[cat] = {'cog':cog,'nom':nom,'slope':slope}
#
    
C.Print(pdfName+"]")

jsonPath = f'CompareFakes/factors_{args.name}.json'
with open(jsonPath,'w') as handle:
    json.dump(data,handle,indent=4)
print (f"Saved {jsonPath}")
