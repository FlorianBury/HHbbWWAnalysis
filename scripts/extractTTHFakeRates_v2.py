#! /bin/env python

import math
import argparse
import ROOT
import json

parser = argparse.ArgumentParser()
parser.add_argument('--file', required=True, type=str, help='ROOT file containing fake rates')
parser.add_argument('--era', required=True, type=str, help='Era')
parser.add_argument('--wp', required=True, type=str,choices=['Loose','Tight'], help='WP for the lepton selection : Loose | Tight')
parser.add_argument('--type', required=False, default='data', type=str, choices=['data','mc'], help='Type : data (actual datadriven) [default] | mc (derived from QCD, used to MC closure -> only for SL)')
parser.add_argument('--channel', required=True, type=str, choices=['DL','SL'], help='Lepton channel : DL | SL')

args = parser.parse_args()

assert args.wp == "Loose" or args.wp == "Tight"
assert args.channel == "DL" or args.channel == "SL"

f = ROOT.TFile.Open(args.file)

if args.wp == "Tight": 
    el_nom = f.Get('FR_mva080_el_data_comb_NC')
    mu_nom = f.Get('FR_mva085_mu_data_comb')
    if args.type == 'mc':
        raise RuntimeError('Tight WP with mc closure not implemented')
if args.wp == "Loose": 
    if args.type == 'data':
        el_nom = f.Get('FR_mva030_el_data_comb')
        mu_nom = f.Get('FR_mva050_mu_data_comb')
    if args.type == 'mc':
        el_nom = f.Get('FR_mva030_el_data_comb_QCD_fakes')
        mu_nom = f.Get('FR_mva050_mu_data_comb_QCD_fakes')
            
        

def getValues(h,lepton):
    xaxis = h.GetXaxis()
    yaxis = h.GetYaxis()
    idx = [_ for _ in range(1,xaxis.GetNbins()+1)]
    idy = [_ for _ in range(1,yaxis.GetNbins()+1)]
    xbinning = [xaxis.GetBinLowEdge(i) for i in range(1,xaxis.GetNbins()+2)]
    ybinning = [yaxis.GetBinLowEdge(i) for i in range(1,yaxis.GetNbins()+2)]

    data_dict = {}
    for iy in range(1,yaxis.GetNbins()+1):
        yMin = yaxis.GetBinLowEdge(iy)
        yMax = yaxis.GetBinUpEdge(iy)
        data_dict[(yMin,yMax)] = {}
        for ix in range(1,xaxis.GetNbins()+1):
            xMin = xaxis.GetBinLowEdge(ix)
            xMax = xaxis.GetBinUpEdge(ix)
            content = h.GetBinContent(ix,iy)
            error   = h.GetBinError(ix,iy)
            data_dict[(yMin,yMax)][(xMin,xMax)] = (content,error)

    for yedges in zip(ybinning[:-1],ybinning[1:]):
        for xedges in zip(xbinning[:-1],xbinning[1:]):
            data = []
            for ykey in data_dict.keys():
                ydict = {'bin':ykey,'values':[]}
                for xkey in data_dict[ykey].keys():
                    if xedges == xkey and yedges == ykey:
                        content,error = data_dict[ykey][xkey]
                        if args.channel == "SL":
                            content = min(content,0.99) # clipping
                    else:
                        content,error = 1.,0.
                    ydict['values'].append({'bin':xkey,
                                            'value':content,
                                            'error_low': error,
                                            'error_high': error})

                data.append(ydict)
            json_content = {'dimension': 2, 'variables': ['AbsEta', 'Pt'], 'binning': {'x': ybinning, 'y': xbinning}, 'data': data, 'error_type': 'absolute'}
            syst_suffix = f"abseta_{yedges[0]}_{yedges[1]}_pt_{xedges[0]}_{xedges[1]}".replace('.','p')
            filename = f'TTHFakeRates_{args.wp}MVA_{args.channel}_{lepton}_{args.era}_{syst_suffix}.json'
            with open(filename, 'w') as j:
                json.dump(json_content, j, indent=2)
            print ("Created file %s"%filename)

getValues(el_nom,"Electron")
getValues(mu_nom,"Muon")
