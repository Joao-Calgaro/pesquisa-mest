# Tentativa minha de simular uma arvore unica geradora e codifica-la

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

import phy
import inspect


with open('input_classificacao_unica.txt', 'r') as des:
    des_data = des.read()
des.close()

design = pd.read_table(io.StringIO(des_data), index_col='index')


design = design.loc[:,
         ['R_nought', 'transmission_rate', 'removal_rate', 'sampling_proba', 'infectious_time', 'tree_size'
          ]]
params = design.iloc[0, ] 


tr, vector_counter = phy.simulate_bd_tree_gillespie(transmission_r=params[1], removal_r=params[2], sampling_p=params[3],
                                                    max_s=params[5], max_t=10000)

#print(tr)

i = 0
for node in tr.traverse("levelorder"):
    node.name = "n" + str(i)
    i += 1

# STOCK the tree
# remove unsampled tips
tr = phy.remove_certain_leaves(tr, to_remove=lambda node: getattr(node, 'stop_reason') != 3)

#tentativa de debug:
if tr is not None and len(tr.children) > 0:
    tuple_tree, _ = phy.encode_into_most_recent(tr, 1) #esse 1 é o sampling probability

print(tuple_tree)

tree_array = tuple_tree.values.flatten()  # Shape (402,)

print(tree_array)

save_data = {
        'tree_data': tree_array,
    }


np.savez_compressed('arvore_unica_vetorizada.npz', **save_data)
