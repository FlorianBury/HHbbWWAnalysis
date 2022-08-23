import os
import sys
import json
import yaml

assert len(sys.argv) == 3
input_json = sys.argv[1]
output_yaml = sys.argv[2]

with open(input_json,"r") as handle:
    factors = json.load(handle)
d = {}
for cat,val in factors.items():
    u = min(0.99,max(0,abs(1-val)))
    print (cat,val,u)
    era = None
    if '2016' in cat:
        era = '2016'
    if '2017' in cat:
        era = '2017'
    if '2018' in cat:
        era = '2018'
    entry = {'cat':cat,'group':'DY','val':[round(1-u,3),round(1+u,3)]}
    if era is not None:
        entry['era'] = era
        entry['cat'] = entry['cat'].replace('_'+era,'')
    d[cat.replace('HH_','CMS_bbwwdl_DY_ncc_').replace('DL_','')] = entry
with open(output_yaml,'w') as handle:
    yaml.dump(d,handle)

from pprint import pprint
pprint (d)
