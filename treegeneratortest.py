# código do augusto para gerar as arvores

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

    #print('--------------------------------------')
    #print(vector_counter)
    #print('--------------------------------------')
    # tr é o grapho da arvore
    # vector counter é statistics on the number of branches, removed tips, sampled tips, time of simulation and number of sim trials

    #tentativa de debug:
    if tr is not None and len(tr.children) > 0:
        # a função encode_into_most_recent retorna a árvore CBLV e o fator de reescalonamento
        tuple_tree, _ = phy.encode_into_most_recent(tr, 1) #esse 1 é o sampling probability

    #print('------------------------------------')
    #print(tuple_tree)
    #print('-'*30)
    #tuple_tree é a arvore codificada em CBLV, mas parece que não está sendo salvada 

    #Processa a primeira metade da árvore (colunas 0-200)
    first_half = tuple_tree.iloc[:, :201].values.flatten()
    first_half = first_half[first_half != 0]  # Remove zeros
    if len(first_half) > 0:
        first_half = first_half[:-1]  #Remove último elemento (que é o reescaling factor)
        first_half = np.append(first_half, 0)  #Adiciona zero no final
            
    #Processa a segunda metade da árvore (colunas 201 em diante)
    second_half = tuple_tree.iloc[:, 201:].values.flatten()
    second_half = second_half[second_half != 0]  #Remove zeros
    if len(second_half) > 0:
        second_half = second_half[:-1]  #Remove último elemento

    debug_print(f"Experiment {experiment_id} - tuple_tree shape: {tuple_tree.shape}")

    if tuple_tree is not None:
        # Convert the entire tuple_tree to numpy array and flatten
        tree_array = tuple_tree.values.flatten()  # Shape (402,)
        
        # Store arrays for deep learning
        tree_arrays.append(tree_array)
        parameter_arrays.append(params.values)  # Store parameters
        statistics_arrays.append(vector_counter)  # Store statistics
        
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
    parameter_matrix = np.array(parameter_arrays)  # Shape: (nb_samples, n_params)
    statistics_matrix = np.array(statistics_arrays)  # Shape: (nb_samples, n_stats)

    # Debug information
    debug_print(f"\n=== DEEP LEARNING DATASET SUMMARY ===")
    debug_print(f"Tree matrix shape: {tree_matrix.shape}")
    debug_print(f"Parameter matrix shape: {parameter_matrix.shape}")
    debug_print(f"Statistics matrix shape: {statistics_matrix.shape}")
    debug_print(f"Tree matrix dtype: {tree_matrix.dtype}")
    debug_print(f"Parameter matrix dtype: {parameter_matrix.dtype}")

    save_data = {
        'tree_data': tree_matrix,
        'parameters': parameter_matrix,
        'statistics': statistics_matrix,
        'experiment_ids': design.index.values,
        'parameter_names': list(design.columns),
        'statistics_names': col2,
        'timestamp': datetime.now().isoformat(),
        'maxTime': maxTime,
        'nb_samples': nb_samples
    }

    # Get column names if tuple_tree exists
    if 'tuple_tree' in locals() and tuple_tree is not None:
        save_data['column_names'] = list(tuple_tree.columns)
    else:
        save_data['column_names'] = [f'col_{i}' for i in range(tree_matrix.shape[1])]

    # Save as compressed numpy format
    np.savez_compressed('deep_learning_dataset2.npz', **save_data)
    debug_print("Deep learning dataset saved as: deep_learning_dataset2.npz")
    #aqui as informacoes adicionais sao salvas
    np.save('tree_data.npy', tree_matrix)
    np.save('parameters.npy', parameter_matrix)
    np.save('statistics.npy', statistics_matrix)

    debug_print("Individual matrices saved as .npy files")
else:
    debug_print("erro")


# subpopulations to export as csv
#tats_export.to_csv(path_or_buf="subpopulations.txt", sep='\t', index=True, header=True)

# for the pipe : export to stdout
sys.stdout.write(forest_export.to_csv(sep='\t', index=True, header=True))


#python treegeneratortest.py input_classificacao.txt 10000 > test4.txt

#rodem do jeito abaixo, substituindo o arquivo de entrada como desejado:
#python treegeneratortest.py input_regressao.txt 10000 


# no codigo do augusto, conferir oq acontece depois da função gillespie
