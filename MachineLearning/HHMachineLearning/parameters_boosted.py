# Defines the parameters that users might need to change
# Must be included manually in each script

# /!\ Two types of files need to be filled with users parameters 
#           - parameters.py (mort important one)
#           - sampleList.py (on what samples to run)
#           (optionnaly NeuralNet.py for early_stopping etc)
import multiprocessing
import os
from tensorflow.keras.losses import binary_crossentropy, mean_squared_error, logcosh, categorical_crossentropy
from tensorflow.keras.optimizers import RMSprop, Adam, Nadam, SGD            
from tensorflow.keras.activations import relu, elu, selu, softmax, tanh, sigmoid
from tensorflow.keras.regularizers import l1,l2 

from grouped_entropy import GroupedXEnt

##################################  Path variables ####################################

main_path = os.path.abspath(os.path.dirname(__file__))
path_out = '/nfs/scratch/fynu/gsaha/HHMachineLearning_Boosted_output/'
path_model = os.path.join(main_path,'model')

scaler_path = '/home/ucl/cp3/gsaha/bamboodev/HHbbWWAnalysis/MachineLearning/HHMachineLearning/scaler_boosted2016_2017_2018.pkl'
##############################  Datasets proportion   #################################
crossvalidation = True

# Classic training #
# -> For crossvalidation == False
training_ratio = 0.7    # Training set sent to keras
evaluation_ratio = 0.1  # Evaluation set sent to autom8
output_ratio = 0.2      # Output set for plotting later
assert training_ratio + evaluation_ratio + output_ratio == 1 
# Will only be taken into account for the masks generation, ignored after

# Cross-validation #
# -> For crossvalidation == True
N_models = 5  # Number of models to train
N_train  = 3  # Number of slices on which to train
N_eval   = 1  # Number of slices on which to evaluate the model
N_apply  = 1  # Number of slices on which the model will be applied for uses in analysis
N_slices = N_train+N_eval+N_apply
splitbranch = 'event' # Will split the training based on "event % N_slices" 
# /!\ Don't forget to add "splitbranch" in other_variables

if N_slices % N_models != 0: # will not work otherwise
    raise RuntimeError("N_slices [{}] % N_models [{}] should be == 0".format(N_slices,N_models))
if N_apply != N_slices/N_models: # Otherwise the same slice can be applied on several networks
    raise RuntimeError("You have asked {} models that should be applied {} times, the number of slices should be {} but is {}+{}+{}".format(N_models,N_apply,N_models*N_apply,N_train,N_eval,N_apply))

############################### Slurm parameters ######################################
partition = 'gpu'  # Def, cp3 or cp3-gpu
QOS = 'normal' # cp3 or normal
time = '0-06:00:00' # days-hh:mm:ss
mem = '50000' # ram in MB
tasks = 1 # Number of threads(as a string) (not parallel training for classic mode)
cpus = 1
gpus = 1
workers = 7

##################################  Naming ######################################
# Physics Config #
config = os.path.join(os.path.abspath(os.path.dirname(__file__)),'sampleListSL_Boosted.yml')
lumidict = {'2016':35922,'2017':41529.152060112,'2018':59740.565201546}
#eras = ['2016']
eras = ['2016','2017','2018'] # To enable or disable eras, add or remove from this list

categories = ['boosted2b2Wj','boosted2b1Wj', 'boosted2b0Wj']
channels = ['El','Mu']

# Better put them in alphabetical order
nodes = ['GGF','H','Rare','ST','TT','VBF','WJets']
group_ids = [
        (1.0, [0,5]),           # signals
        (1.0, [1,2,3,4,6]),     # backgrounds
        (1.0, [0]),             # GGF
        (1.0, [4]),             # TT
        (1.0, [5]),             # VBF
        (1.0, [6]),             # WJets
        (1.0, [1,2,3]),         # rest : H + Rare + ST       
            ]
        
#    def loss(self): # -> DL
#        return GroupedXEnt(
#            group_ids=[
#                (1.0, [0, 1]),  # Signal
#                (1.0, [2, 3, 4, 5, 6, 7, 8]),  # Bkg
#                (1.0, [0]),  # GGF
#                (1.0, [1]),  # VBF
#                (1.0, [3]),  # DY
#                (1.0, [5]),  # ttbar
#                (1.0, [2, 4, 6, 7, 8]),  # Other
#            ]
#        )

# Input plots options #
node_colors = {
            'GGF'   : '#288a24',
            'H'     : '#06b894',
            'Rare'  : '#610596',
            'ST'    : '#99053d',
            'TT'    : '#cc7a16',
            'VBF'   : '#8f0a1e',
            'WJets' : '#d95564',
             }

# Tree name #
tree_name = 'Events'

# scaler and mask names #
suffix = 'boosted' 
# scaler_name -> 'scaler_{suffix}.pkl'  If does not exist will be created 
# mask_name -> 'mask_{suffix}_{sample}.npy'  If does not exist will be created 
scaler_name = 'scaler_'+suffix+'_'.join([str(era) for era in eras])+'.pkl'
scaler_path = os.path.join(main_path,scaler_name)

# Data cache #
train_cache = os.path.join(path_out,'train_cache_'+'_'.join([str(era) for era in eras])+'.pkl')
test_cache = os.path.join(path_out,'test_cache_'+'_'.join([str(era) for era in eras])+'.pkl')

# Meta config info #
xsec_json = os.path.join(main_path,'background_{era}_xsec.json')
event_weight_sum_json = os.path.join(main_path,'background_{era}_event_weight_sum.json')

# Training resume #
resume_model = ''

# Output #
output_batch_size = 50000
split_name = 'tag' # 'sample' or 'tag' : criterion for output file splitting

##############################  Evaluation criterion   ################################

eval_criterion = "eval_error" # either val_loss or eval_error or val_acc

##############################  Model callbacks ################################
# Early stopping to stop learning after some criterion 
early_stopping_params = {'monitor'   : 'val_loss',  # Value to monitor
                         'min_delta' : 0.0001,          # Minimum delta to declare an improvement
                         'patience'  : 10,          # How much time to wait for an improvement
                         'verbose'   : 1,           # Verbosity level
                         'restore_best_weights':True,
                         'mode'      : 'min'}       # Mode : 'auto', 'min', 'max'

# Reduce LR on plateau : if no improvement for some time, will reduce lr by a certain factor
reduceLR_params = {'monitor'    : 'val_loss',   # Value to monitor
                   'factor'     : 0.1,          # Multiplicative factor by which to multiply LR
                   'min_lr'     : 1e-6,         # Minimum value for LR
                   'min_delta'  : 0.0001,       # Minimum delta to declare an improvement
                   'patience'   : 5,            # How much time to wait for an improvement
                   'cooldown'   : 1,            # How many epochs before starting again to monitor
                   'verbose'    : 1,            # Verbosity level
                   'mode'      : 'min'}         # Mode : 'auto', 'min', 'max'
    
#################################  Scan dictionary   ##################################
# /!\ Lists must always contain something (even if 0), otherwise 0 hyperparameters #
# Classification #
grouped_loss = GroupedXEnt(group_ids)
#p = { 
#    'lr' : [0.01,0.001], 
#    'first_neuron' : [64,128,256,512,1024],
#    'activation' : [relu],
#    'dropout' : [0.,0.1,0.2],
#    'hidden_layers' : [2,3,4,5,6], # does not take into account the first layer
#    'output_activation' : [softmax],
#    'l2' : [0,0.01,0.1],
#    'optimizer' : [Adam],  
#    'epochs' : [50],   
#    'batch_size' : [10000], 
#    'n_particles' : [16],
#    'loss_function' : [grouped_loss] , #  [categorical_crossentropy]
#}
p = { 
    'lr' : [0.001], 
    'first_neuron' : [128],
    'activation' : [relu],
    'dropout' : [0.],
    'hidden_layers' : [3], # does not take into account the first layer
    'output_activation' : [softmax],
    'l2' : [1e-5],
    'optimizer' : [Adam],  
    'epochs' : [1],   
    'batch_size' : [100000], 
    'n_particles' : [10],
    'loss_function' : [grouped_loss],
}


repetition = 1 # How many times each hyperparameter has to be used 

###################################  Variables   ######################################

cut = 'MC_weight > 0 && abs(total_weight)<1e5'
#cut = 'MC_weight > 0'

weight = 'total_weight'
#weight = None

# Input branches (combinations possible just as in ROOT #
#/!\ onehot variables need to be at the beginning of the list (checked later)
inputs = [
            # Onehot #
            #'$era@op_era',
            'lep_pdgId@op_pdgid',
            'lep_charge@op_charge',
            'JPAcat@op_boosted_JPAcat',
            # JPA values #
            'L2_Hbb2Wj',
            'L2_Hbb1Wj',
            'L2_Hbb0Wj',
            # LL variables #
            'METpt',               
            'lep_pt',
            'fatbj_pt',
            'fatbj_softdropMass',
            'fatbj_btagDeepB',
            'cosThetaS_Hbb',
            'mT_top_3particle',
            'wj1_pt',
            'wj1_bTagDeepFlavB',
            'wj2_pt',
            'wj2_bTagDeepFlavB ',
            'nAk8Jets',
            'nAk4JetsCleandWithAk8b',           
            'VBF_tag',
            # HL variables #
            'lepmet_DPhi',
            'lepmet_pt',
            'lep_MT',
            'MET_LD',
            'fatbj_lepDR',
            'fatbj_Wj1DR',
            'fatbj_wj2DR',
            'wj1wj2_pt',
            'wj1wj2DR',
            'wj1wj2invM',
            'WWplaneAngle_withMET', 
            'HWplaneAngle',
            'HWW_Simple_Mass',
            'HWW_dR',
            'cosThetaS_Wjj_simple',
            'cosThetaS_WW_simple_met ',
            'cosThetaS_HH_simple_met',
            'angleBetWWPlane',
            'angleBetHWPlane',
            'bj1bj2_DR',
            'bj1bj2_DPhi',
            'bj2wj1_DR',
            'bj2wj1_DPhi',
            'wj1wj2_DR',
            'wj1wj2_DPhi',
            'bj1wj2_DR',
            'bj1wj2_DPhi',
            'bj1wj1_DR',
            'bj1wj1_DPhi',
#            'VBFj1pt',
#            'VBFj2pt',
#            'VBFj1eta',
#            'VBFj2eta',
#            'VBFj1j2dEta',
#            'VBFj1j2dPhi',
#            'VBFj1j2invM',
            'zeppenfeldVar',
            'minJetDR',
            'minLepJetDR',
            'HT2_lepJetMet',
#            'HT2R_lepJetMet',
    ]

operations = [inp.split('@')[1] if '@' in inp else None  for inp  in  inputs]
check_op = [(o is not None)*1 for o in operations]
if check_op != sorted(check_op,reverse=True):
    raise RuntimeError('Onehot inputs need to be at the beginning of the inputs list')
mask_op = [len(inp.split('@'))==2 for inp  in  inputs]
inputs = [inp.split('@')[0] for inp  in  inputs]

LBN_inputs = [
              'lep_E','lep_Px','lep_Py','lep_Pz',
              'fatbj_E','fatbj_Px','fatbj_Py','fatbj_Pz',
              'wj1_E','wj1_Px','wj1_Py','wj1_Pz',
              'wj2_E','wj2_Px','wj2_Py','wj2_Pz',
             ]
# /!\ format = [E,Px,Py,Pz] per particle

assert len(LBN_inputs)%4 == 0


# Output branches #

outputs = [
            '$GGF',
            '$H',
            '$Rare',
            '$ST',
            '$TT',
            '$VBF',
            '$WJets',
          ] 

# Other variables you might want to save in the tree #
other_variables = [
    'event',
    'ls',
    'run',
    'MC_weight',
]


################################  dtype operation ##############################

def make_dtype(list_names): # root_numpy does not like . and ()
    list_dtype = [(name.replace('.','_').replace('(','').replace(')','').replace('-','_minus_').replace('*','_times_')) for name in list_names]
    return list_dtype
