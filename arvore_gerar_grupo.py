import warnings
warnings.filterwarnings("ignore")
import phylodeep
from phylodeep import encoding
import dendropy
from dendropy.simulate import treesim
import ete3
from ete3 import Tree
import sys
import argparse
import io
import numpy as np
import pandas as pd
import phylodeep.sumstats as sumstats
from collections import Counter
from ete3 import Tree
from io import StringIO
import matplotlib.pyplot as plt
from datetime import datetime

DISTANCE_TO_ROOT = "dist_to_root"
DEPTH = "depth"
HEIGHT = "height"
LADDER = "ladder"
VISITED = "visited"
TARGET_AVG_BL = 1
STOP_REASON = 'stop_reason'
STOP_UNKNOWN = 0
STOP_TRANSMISSION = 1
STOP_REMOVAL_WOS = 2
STOP_SAMPLING = 3
STOP_TIME = 4
HISTORY = 'history'
SAMPLING = 'sampling'
TRANSMISSION = 'transmission'
DIST_TO_START = 'DIST_TO_START'
PROCESSED = 'processed'

import phy
import inspect

parser = argparse.ArgumentParser(description='Generates a tree')
parser.add_argument('inputFile', help='an input file with parameters for the tree')
parser.add_argument('maxTime', type=float, help='an input float for the maximal simulation time')

args = parser.parse_args()

# read experiment design
with open(args.inputFile, 'r') as des:
    des_data = des.read()
des.close()

design = pd.read_table(io.StringIO(des_data), index_col='index')

#The .loc indexer in pandas is used to access rows and columns by label (not by numeric position).
#The colon : means "select all rows"
design = design.loc[:,
         ['R_nought', 'transmission_rate', 'removal_rate', 'sampling_proba', 'infectious_time', 'tree_size'
          ]]

nb_samples = len(design)

maxTime = args.maxTime

tree_arrays = []
parameter_arrays = []
statistics_arrays = []

# PREPARE EXPORT
# stock all trees in a list
forest = []

# col names of export statistics
col = ['tree']
forest_export = pd.DataFrame(index=design.index, columns=col)

col2 = ['total_leaves', 'removed_leaves', 'sampled_leaves', 'time_of_simulation', 'nb_trials']
#esses stats são apenas informações sobre a simulação, não estão relacionados com as SS
stats_export = pd.DataFrame(index=design.index, columns=col2)

def debug_print(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

for experiment_id in range(nb_samples):
    
    params = design.iloc[experiment_id, ] #pega a linha inteira (todas as colunas) na posição experiment_id

    # simulation
    tr, vector_counter = phy.simulate_bd_tree_gillespie(transmission_r=params[1], removal_r=params[2], sampling_p=params[3],
                                                    max_s=params[5], max_t=maxTime)
    # for display purposes
    i = 0
    for node in tr.traverse("levelorder"):
        node.name = "n" + str(i)
        i += 1

    # STOCK the tree
    # remove unsampled tips
    tr = phy.remove_certain_leaves(tr, to_remove=lambda node: getattr(node, STOP_REASON) != STOP_SAMPLING)

    #tentativa de debug:
    if tr is not None and len(tr.children) > 0:
        tuple_tree, _ = phy.encode_into_most_recent(tr, 1) #esse 1 é o sampling probability
            

    debug_print(f"Experiment {experiment_id} - tuple_tree shape: {tuple_tree.shape}")

    if tuple_tree is not None:
        # Convert the entire tuple_tree to numpy array and flatten
        tree_array = tuple_tree.values.flatten()  # Shape (402,)
        
        # Store arrays for deep learning
        tree_arrays.append(tree_array)
        
        # For the original export, convert numpy array to string representation
        # Convert array to tab-separated string
        array_string = '\t'.join(map(str, tree_array))
        forest_export.iloc[experiment_id][0] = array_string
        
        # Alternative: Keep the original tuple_tree CSV format
        # forest_export.iloc[experiment_id][0] = tuple_tree.to_csv(sep='\t', index=False)
    else:
        #caso a arvore seja nula 
        placeholder_tree = np.zeros(402)  
        tree_arrays.append(placeholder_tree)
        parameter_arrays.append(params.values)
        statistics_arrays.append(vector_counter)
        forest_export.iloc[experiment_id][0] = "NA"

    stats_export.iloc[experiment_id] = vector_counter

# Convert lists to numpy matrices for deep learning
if tree_arrays: 
    tree_matrix = np.array(tree_arrays)  # Shape: (nb_samples, 402)


    save_data = {
        'tree_data': tree_matrix,
    }

    # Save as compressed numpy format
    np.savez_compressed('deep_learning_dataset2.npz', **save_data)




# subpopulations to export as csv
#tats_export.to_csv(path_or_buf="subpopulations.txt", sep='\t', index=True, header=True)

# for the pipe : export to stdout
#sys.stdout.write(forest_export.to_csv(sep='\t', index=True, header=True))




#python treegeneratortest.py input_classificacao.txt 10000 > test4.txt

#rodem do jeito abaixo, substituindo o arquivo de entrada como desejado:
#python treegeneratortest.py input_regressao.txt 10000 


# no codigo do augusto, conferir oq acontece depois da função gillespie
