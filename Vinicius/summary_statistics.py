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

parser = argparse.ArgumentParser(description='Generates a tree')
parser.add_argument('inputFile', help='an input file with parameters for the tree')
parser.add_argument('maxTime', type=float, help='an input float for the maximal simulation time')

args = parser.parse_args()

# read experiment design
with open(args.inputFile, 'r') as des:
    des_data = des.read()
des.close()

design = pd.read_table(io.StringIO(des_data), index_col='index')

design = design.loc[:,
         ['R_nought', 'transmission_rate', 'removal_rate', 'sampling_proba', 'infectious_time', 'tree_size'
          ]]

nb_samples = len(design)
maxTime = args.maxTime

tree_arrays = []
parameter_arrays = []

# PREPARE EXPORT
col = ['tree']
forest_export = pd.DataFrame(index=design.index, columns=col)

col2 = ['total_leaves', 'removed_leaves', 'sampled_leaves', 'time_of_simulation', 'nb_trials']
stats_export = pd.DataFrame(index=design.index, columns=col2)

def debug_print(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

successful_trees = 0
failed_trees = 0

def is_valid_for_summary_stats(tree, min_leaves=3):
    """
    Check whether a tree is valid for computing summary statistics.
    Summary statistics assume a non-degenerate tree.
    """
    if tree is None:
        return False

    # Number of leaves after pruning
    if len(tree) < min_leaves:
        return False

    # Must have at least one internal node with >= 2 children
    internal_nodes = [n for n in tree.traverse() if not n.is_leaf()]
    if len(internal_nodes) == 0:
        return False

    for n in internal_nodes:
        if len(n.children) >= 2:
            return True

    return False



for experiment_id in range(nb_samples):
    
    params = design.iloc[experiment_id, ]

    tr, vector_counter = phy.simulate_bd_tree_gillespie(
        transmission_r=params["transmission_rate"],
        removal_r=params["removal_rate"],
        sampling_p=params["sampling_proba"],
        max_s=int(params["tree_size"]),
        max_t=maxTime
    )

    # Rename nodes
    i = 0
    for node in tr.traverse("levelorder"):
        node.name = "n" + str(i)
        i += 1

    # Remove unsampled tips
    tr = phy.remove_certain_leaves(
        tr,
        to_remove=lambda node: getattr(node, STOP_REASON) != STOP_SAMPLING
    )

    SS_tree = None

    if is_valid_for_summary_stats(tr):
        try:
            SS_tree, _ = phy.encode_into_summary_statistics(tr, 1)
        except Exception as e:
            debug_print(
                f"Experiment {experiment_id} - Encoding failed: {type(e).__name__}: {e}"
            )
            SS_tree = None

    if SS_tree is not None:
        debug_print(f"Experiment {experiment_id} - SS_tree shape: {SS_tree.shape}")
        tree_arrays.append(SS_tree)
        parameter_arrays.append(params.values)
        successful_trees += 1
    else:
        debug_print(
            f"Experiment {experiment_id} - Skipped: Tree invalid for summary statistics"
        )
        failed_trees += 1


if tree_arrays: 
    tree_matrix = np.array(tree_arrays)  
    parameter_matrix = np.array(parameter_arrays) 

    save_data = {
        'tree_data': tree_matrix,
        'parameters': parameter_matrix,
        'experiment_ids': design.index.values,
        'parameter_names': list(design.columns),
        'timestamp': datetime.now().isoformat(),
        'maxTime': maxTime,
        'nb_samples': nb_samples,
        'successful_trees': successful_trees,
        'failed_trees': failed_trees
    }

    np.savez_compressed('arvores_nao_voznica2.npz', **save_data)
    debug_print(f"Deep learning dataset saved as: deep_learning_dataset2.npz")
    debug_print(f"Successfully processed: {successful_trees} trees")
    debug_print(f"Failed trees: {failed_trees}")
    
    np.save('tree_data.npy', tree_matrix)
    np.save('parameters.npy', parameter_matrix)
    debug_print("Individual matrices saved as .npy files")
else:
    debug_print("Error: No valid trees were generated")

stats_export.to_csv(path_or_buf="subpopulations.txt", sep='\t', index=True, header=True)

sys.stdout.write(forest_export.to_csv(sep='\t', index=True, header=True))


#python summary_statistics.py input_classificacao_vini.txt 1000 
#python summary_statistics2.py teste_parametros_3k.txt 3000
#python summary_statistics.py test_cod_antigo_10k.txt 10000
#python summary_statistics.py test_params_invertidos.txt 10000
#python summary_statistics.py parametros_bd_controlados_10k.txt 10000