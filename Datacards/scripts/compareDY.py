import os
import sys
import glob
import copy
import ctypes
import json
import math
import argparse
from array import array
import numpy as np
import ROOT
from IPython import embed

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser(description='Comparing DY')
parser.add_argument('--mc', action='store', required=True, type=str,
                    help='Path to DY MC shapes')
parser.add_argument('--dd', action='store', required=True, type=str,
                    help='Path to DY datadriven shapes')
parser.add_argument('--cl', action='store', required=True, type=str,
                    help='Path to DY closure shapes')
parser.add_argument('--name', action='store', required=True, type=str,
                    help='Name of the output files')
parser.add_argument('--rebin', action='store', required=False, default=None, type=int,
                    help='Use the rebin number')
parser.add_argument('--fit', action='store_true', required=False, default=False,
                    help='Use fit method for all categories')
parser.add_argument('--split_eras', action='store_true', required=False, default=False,
                    help='For splitting per eras')
args = parser.parse_args()

if '{era}' not in args.mc:
    raise RuntimeError(f'Missing {{era}} in {args.mc}')
if '{era}' not in args.dd:
    raise RuntimeError(f'Missing {{era}} in {args.dd}')
if '{era}' not in args.cl:
    raise RuntimeError(f'Missing {{era}} in {args.cl}')

factors = {}

eras = ['2016','2017','2018']

plotDict_MC = {era:{} for era in eras}
plotDict_DD = {era:{} for era in eras}
plotDict_CL = {era:{} for era in eras}

def processHist(h,x):
    for i in range(1,h.GetNbinsX()+1):
        if h.GetBinContent(i) < 1e-3:
            h.SetBinContent(i,0.)
            h.SetBinError(i,0.)
    if isinstance(x,int):
        return h.Rebin(xn)
    else:
        return h.Rebin(len(x)-1,'rebin',array('d',x))


for era in eras:
    path_DY_MC  = args.mc.format(era=era)
    path_DY_DD  = args.dd.format(era=era)
    path_DY_CL  = args.cl.format(era=era)

    plotNames = [os.path.basename(name) for name in glob.glob(path_DY_MC)]

    for plotName in sorted(plotNames):
        print (plotName,era)
        f_MC = ROOT.TFile(os.path.join(os.path.dirname(path_DY_MC),plotName))
        f_DD = ROOT.TFile(os.path.join(os.path.dirname(path_DY_DD),plotName))
        f_CL = ROOT.TFile(os.path.join(os.path.dirname(path_DY_CL),plotName))

        plotName = plotName.replace('.root','').replace(f'_{era}','')

        h_MC = f_MC.Get("DY")
        h_DD = f_DD.Get("DY")
        h_CL = f_CL.Get("DY")
        if args.rebin is not None:
            xn = args.rebin
        elif 'boosted' in plotName:
            #xn = [0.,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.]
            xn = [0.,0.20,0.40,0.60,0.80,1.00]
        else:
            xn = [0.,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.9,0.95,1.]

        h_MC = processHist(h_MC,xn)
        h_DD = processHist(h_DD,xn)
        h_CL = processHist(h_CL,xn)


        if plotName not in plotDict_MC[era].keys():
            plotDict_MC[era][plotName] = copy.deepcopy(h_MC)
        else:
            plotDict_MC[era][plotName].Add(h_MC)

        if plotName not in plotDict_DD[era].keys():
            plotDict_DD[era][plotName] = copy.deepcopy(h_DD)
        else:
            plotDict_DD[era][plotName].Add(h_DD)

        if plotName not in plotDict_CL[era].keys():
            plotDict_CL[era][plotName] = copy.deepcopy(h_CL)
        else:
            plotDict_CL[era][plotName].Add(h_CL)

        f_MC.Close()
        f_DD.Close()
        f_CL.Close()

def computeFactor(a,b,sa,sb):
    """
        factor = a/b
        a +/- sa
        b +/- sb
        -> sfactor by error propagation
    """
    factor = a/b
    sfactor = math.sqrt((1./b**2)*sa**2 + (a**2/b**4)*sb**2)
    return factor,sfactor

def decorrelatedBands(coefficients,covariance,x_min,x_max,N=101):
    coefficients = np.array(coefficients)
    covariance   = np.array(covariance)

    eigenvalues , eigenvectors = np.linalg.eigh(covariance)

    order = len(coefficients)

    xs = np.linspace(x_min,x_max,N)
    graphs = [ROOT.TGraphErrors(N) for _ in range(order)]

    for i in range(len(xs)):
        x = [xs[i]**o for o in range(order)]
        y = np.multiply(x,coefficients).sum()
        for o in range(order):
            var = np.abs(np.multiply(np.sqrt(eigenvalues[o]) * eigenvectors[:,o], x).sum())
            graphs[o].SetPoint(i,xs[i],y)
            graphs[o].SetPointError(i,0,var)
    return graphs






categories = set(plotDict_MC['2016'])
for plotDict in [plotDict_MC,plotDict_DD,plotDict_CL]:
    for era,histDict in plotDict.items():
        categories = categories.intersection(set(histDict.keys()))

def plotCategory(C,h_MC,h_DD,h_CL,title='',method=''):
    C.Clear()

    Nerr_MC = ctypes.c_double(0.)
    Nerr_DD = ctypes.c_double(0.)
    Nerr_CL = ctypes.c_double(0.)
    N_MC = h_MC.IntegralAndError(0,h_MC.GetNbinsX()+1,Nerr_MC)
    N_DD = h_DD.IntegralAndError(0,h_DD.GetNbinsX()+1,Nerr_DD)
    N_CL = h_CL.IntegralAndError(0,h_CL.GetNbinsX()+1,Nerr_CL)

    Nerr_MC = Nerr_MC.value
    Nerr_DD = Nerr_DD.value
    Nerr_CL = Nerr_CL.value

    factor_DD = computeFactor(N_MC,N_DD,Nerr_MC,Nerr_DD)
    factor_CL = computeFactor(N_MC,N_CL,Nerr_MC,Nerr_CL)

    era = title.split('_')[-1]
    title = ' '.join(title.replace('HH_','').split('_')[:-1])

    h_MC.SetTitle(f"{title} category ({era});;Events / {h_MC.GetXaxis().GetBinWidth(1):0.2f}")
    hmax = max([h.GetMaximum() for h in [h_MC,h_DD,h_CL]])
    h_MC.SetMaximum(hmax*1.8)

    h_MC.SetLineColor(632)
    h_DD.SetLineColor(600)
    h_CL.SetLineColor(1)

    h_MC.SetMarkerColor(632)
    h_DD.SetMarkerColor(600)
    h_CL.SetMarkerColor(1)

    h_MC.SetMarkerStyle(25)
    h_DD.SetMarkerStyle(8)
    h_CL.SetMarkerStyle(4)

    h_MC.SetTitleSize(0.08)
    h_MC.GetXaxis().SetTitleSize(0.08)
    h_MC.GetYaxis().SetTitleSize(0.08)

    # compute quantiles for the fit #
    prob = np.array([0.,1.])
    #prob = np.array([0.,1.])
    quantiles = np.array([0.,0.])
    h_MC.GetQuantiles(2,quantiles,prob)
    #quantiles[1] = 0.6
    #quantiles[0] = 0.6

    ratio_DY_DD = h_MC.Clone("ratio_DY_DD")
    ratio_DY_CL = h_MC.Clone("ratio_DY_CL")

    #ratio_DY_DD.Add(h_DD,-1)
    #ratio_DY_CL.Add(h_CL,-1)

    ratio_DY_DD.Divide(h_DD)
    ratio_DY_CL.Divide(h_CL)

    ratio_DY_DD.SetLineColor(600)
    ratio_DY_CL.SetLineColor(1)

    ratio_DY_DD.SetMarkerColor(600)
    ratio_DY_CL.SetMarkerColor(1)

    ratio_DY_DD.SetMarkerStyle(8)
    ratio_DY_CL.SetMarkerStyle(4)

    ratio_DY_DD.SetTitle(";DNN score;Ratio")
    ratio_DY_DD.GetXaxis().SetTitleSize(0.08)
    ratio_DY_DD.GetYaxis().SetTitleSize(0.08)

    C.Divide(1,2)

    pad1 = C.cd(1)
    #pad1.SetLogy()
    pad1.SetTopMargin(0.1)
    pad1.SetBottomMargin(0.)
    pad1.SetLeftMargin(0.15)
    pad1.SetRightMargin(0.05)
    pad1.Draw()

    h_MC.Draw("e1p")
    h_DD.Draw("e1p same")
    h_CL.Draw("e1p same")

    leg_up = ROOT.TLegend(0.45,0.55,0.9,0.9)
    leg_up.AddEntry(h_MC,"#splitline{{DY (MC)}}{{Integral = {:0.3f} #pm {:0.3f}}}".format(N_MC,Nerr_MC))
    leg_up.AddEntry(h_DD,"#splitline{{DY (data)}}{{Integral = {:0.3f} #pm {:0.3f}}}".format(N_DD,Nerr_DD))
    leg_up.AddEntry(h_CL,"#splitline{{DY (closure)}}{{Integral = {:0.3f} #pm {:0.3f}}}".format(N_CL,Nerr_CL))
    leg_up.SetBorderSize(0)
    leg_up.SetFillStyle(0)
    leg_up.SetTextSize(0.05)
    leg_up.Draw()

    pad2 = C.cd(2)
    pad2.SetTopMargin(0.0)
    pad2.SetBottomMargin(0.18)
    pad2.SetLeftMargin(0.15)
    pad2.SetRightMargin(0.05)
    pad2.Draw()

    ratio_DY_DD.Draw("e1p")
    ratio_DY_CL.Draw("e1p same")
    ratio_DY_DD.GetYaxis().SetRangeUser(0.,2.)
    pol = "pol1"
    bands_DD = None
    bands_CL = None
    if method == 'fit':
        xmins = []
        xmaxs = []
        for i in range(1,h_CL.GetNbinsX()+1):
            if h_CL.GetBinContent(i) > 1.0 and h_DD.GetBinContent(i) > 1.0 and h_MC.GetBinContent(i) > 1.0:
                xmins.append(ratio_DY_CL.GetXaxis().GetBinLowEdge(i))
                xmaxs.append(ratio_DY_CL.GetXaxis().GetBinUpEdge(i))
        ratio_DY_DD.Fit(pol,"QR","",quantiles[0],quantiles[1])
        fitResult = ratio_DY_CL.Fit(pol,"QRS","",quantiles[0],quantiles[1])
        fit_DD = ratio_DY_DD.GetFunction(pol)
        fit_CL = ratio_DY_CL.GetFunction(pol)
        fit_DD.SetLineColor(600)
        fit_CL.SetLineColor(1)
        covFit = fitResult.GetCovarianceMatrix()
        if pol == "pol1":
            params_DD = (fit_DD.GetParameter(0),fit_DD.GetParameter(1))
            params_CL = (fit_CL.GetParameter(0),fit_CL.GetParameter(1))
            formula_DD = f"Fit : {params_DD[1]:0.2f} x + {params_DD[0]:0.2f}"
            formula_CL = f"Fit : {params_CL[1]:0.2f} x + {params_CL[0]:0.2f}"
            covMatrix = [[covFit[0][0],covFit[0][1]],[covFit[1][0],covFit[1][1]]]
            if math.isnan(covFit[0][0]):
                covMatrix = [[1,0],[0,1]]
        elif pol == "pol2":
            params_DD = (fit_DD.GetParameter(0),fit_DD.GetParameter(1),fit_DD.GetParameter(2))
            params_CL = (fit_CL.GetParameter(0),fit_CL.GetParameter(1),fit_CL.GetParameter(2))
            formula_DD = f"Fit : {params_DD[2]:0.2f} x^2 + {params_DD[1]:0.2f} x + {params_DD[0]:0.2f}"
            formula_CL = f"Fit : {params_CL[2]:0.2f} x^2 + {params_CL[1]:0.2f} x + {params_CL[0]:0.2f}"
            covMatrix = [[covFit[0][0],covFit[0][1],covFit[0][2]],
                         [covFit[1][0],covFit[1][1],covFit[1][2]],
                         [covFit[2][0],covFit[2][1],covFit[2][2]]]
            if math.isnan(covFit[0][0]):
                covMatrix = [[1,0,0],[0,1,0],[0,0,1]]
        else:
            raise NotImplementedError

        bands_DD = decorrelatedBands(params_DD,covMatrix,quantiles[0],quantiles[1])
        bands_CL = decorrelatedBands(params_CL,covMatrix,quantiles[0],quantiles[1])

    leg_down = ROOT.TLegend(0.5,0.65,0.9,1.0)
    if method == 'factor':
        leg_down.AddEntry(ratio_DY_DD,"\splitline{{#frac{{DY (MC)}}{{DY (data)}}}}{{Factor = {:0.5f} #pm {:0.5f}}}".format(*factor_DD))
        leg_down.AddEntry(ratio_DY_CL,"\splitline{{#frac{{DY (MC)}}{{DY (closure)}}}}{{Factor = {:0.5f} #pm {:0.5f}}}".format(*factor_CL))
    elif method == 'fit':
        leg_down.AddEntry(ratio_DY_DD,"\splitline{{#frac{{DY (MC)}}{{DY (data)}}}}{{{:s}}}".format(formula_DD))
        leg_down.AddEntry(ratio_DY_CL,"\splitline{{#frac{{DY (MC)}}{{DY (closure)}}}}{{{:s}}}".format(formula_CL))
    else:
        leg_down.AddEntry(ratio_DY_DD,"#frac{{DY (MC)}}{{DY (data)}}")
        leg_down.AddEntry(ratio_DY_CL,"#frac{{DY (MC)}}{{DY (closure)}}}}{{Factor = {:0.5f} #pm {:0.5f}")
    leg_down.SetTextSize(0.04)
    leg_down.SetBorderSize(0)
    leg_down.SetFillStyle(0)
    leg_down.SetTextSize(0.05)
    leg_down.Draw()
    line = ROOT.TLine(0,1.,1,1.)
    line.Draw()

    if bands_DD is not None:
        for band in bands_DD:
            band.Draw("e3")
            band.SetFillColorAlpha(ratio_DY_DD.GetLineColor(),0.3)
    if bands_CL is not None:
        for band in bands_CL:
            band.Draw("e3")
            band.SetFillColorAlpha(ratio_DY_CL.GetLineColor(),0.3)

    C.Print(f"CompareDY/ComparisonDY_{args.name}.pdf","Title:"+title)

    if method == 'factor':
        return factor_CL[0]
    elif method == 'fit':
        return {'coefficients':params_CL,'covariance':covMatrix}
    else:
        return None


C = ROOT.TCanvas("C","C",700,800)
C.Print(f"CompareDY/ComparisonDY_{args.name}.pdf[")
numerics = {}
factors = {}
fit_results = {}
for cat in sorted(categories):
    if '_GGF' in cat or '_VBF' in cat or '_H' in cat:
        method = 'factor'
    else:
        method= 'fit'
    if args.fit:
        method='fit'
    title = ''
    if 'resolved' in cat or 'inclusive' in cat or args.split_eras:
        for era in eras:
            h_MC = plotDict_MC[era][cat]
            h_DD = plotDict_DD[era][cat]
            h_CL = plotDict_CL[era][cat]
            title = f'{cat}_{era}'
            numeric = plotCategory(C,h_MC,h_DD,h_CL,title=title,method=method)
            if method == 'factor':
                factors[title] = numeric
            if method == 'fit':
                fit_results[title] = numeric

    elif 'boosted' in cat:
        h_MC = plotDict_MC['2016'][cat]
        h_MC.Add(plotDict_MC['2017'][cat])
        h_MC.Add(plotDict_MC['2018'][cat])

        h_DD = plotDict_DD['2016'][cat]
        h_DD.Add(plotDict_DD['2017'][cat])
        h_DD.Add(plotDict_DD['2018'][cat])

        h_CL = plotDict_CL['2016'][cat]
        h_CL.Add(plotDict_CL['2017'][cat])
        h_CL.Add(plotDict_CL['2018'][cat])

        title = f'{cat}'
        numeric = plotCategory(C,h_MC,h_DD,h_CL,title=title,method=method)
        if method == 'factor':
            factors[title] = numeric
        if method == 'fit':
            fit_results[title] = numeric

    else:
        raise RuntimeError('Not understood category')

C.Print(f"CompareDY/ComparisonDY_{args.name}.pdf]")

with open(f'CompareDY/factors_{args.name}.json','w') as handle:
    json.dump(factors,handle,indent=4)
with open(f'CompareDY/fit_results_{args.name}.json','w') as handle:
    json.dump(fit_results,handle,indent=4)
print (f'Saved numerical factors to CompareDY/factors_{args.name}.json')
print (f'Saved fit results to CompareDY/fit_results_{args.name}.json')
