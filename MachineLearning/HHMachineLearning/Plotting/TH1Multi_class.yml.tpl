multiclass:
  filename: 
  tree: tree
  weight: 'event_weight'
  name: multiclass 
  bins: 50
  xmin: 0
  xmax: 1
  title: ''
  xlabel: DNN output
  ylabel: events
  list_variable:
    - output_Ewk
    - output_GGF
    - output_H
    - output_Top
    - output_VBF
    - output_WJets
  list_color:
<<<<<<< HEAD
    - "#288a24"
    - "#06b894"
    - "#610596"
    - "#99053d"
    - "#8f0a1e"
    - "#d95564"
=======
    - '#610596'
    - '#288a24'
    - '#06b894'
    - '#cc7a16'
    - '#8f0a1e'
    - '#d95564'
>>>>>>> eed57b2aa0f195c36370ef8af9adb2772972dcc4
  list_legend:
    - node Ewk
    - node GGF
    - node H
    - node Top
    - node VBF
    - node WJets
  list_cut : '1'
  legend_pos:
    - 0.4
    - 0.60
    - 0.75
    - 0.92
  norm : True
