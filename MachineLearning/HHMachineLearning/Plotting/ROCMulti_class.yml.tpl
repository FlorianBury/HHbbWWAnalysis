ROC_multi:
  tree: tree
  classes:
    - GGF
    - H
    - Rare
    - ST
    - TT
    - VBF
    - WJets
  prob_branches:
    - output_GGF
    - output_H
    - output_Rare
    - output_ST
    - output_TT
    - output_VBF
    - output_WJets
  labels:
    - P(GGF)
    - P(H)
    - P(Rare)
    - P(ST)
    - P(TT)
    - P(VBF)
    - P(WJets)
  colors:
    - "#288a24"
    - "#06b894"
    - "#610596"
    - "#99053d"
    - "#cc7a16"
    - "#8f0a1e"
    - "#d95564"

  weight : '1'
  title : Multiclass
  cut : ''
  selector :
    'GGF'   : 'GGF'
    'H'     : 'H'
    'Rare'  : 'Rare'
    'ST'    : 'ST'
    'TT'    : 'TT'
    'VBF'   : 'VBF'
    'WJets' : 'WJets'
 

