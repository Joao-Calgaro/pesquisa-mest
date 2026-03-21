#!/usr/bin/env python
import warnings
warnings.filterwarnings("ignore")

import os
import sys
import argparse
import io
import numpy as np
import pandas as pd
from datetime import datetime

import phy
from ete3 import Tree

# -----------------------------
# ARGUMENTOS
# -----------------------------
parser = argparse.ArgumentParser(description='Generates BD trees + summary statistics')
parser.add_argument('inputFile', help='design file with parameters')
parser.add_argument('maxTime', type=float, help='maximal simulation time')
args = parser.parse_args()

# -----------------------------
# LER DESIGN
# -----------------------------
design = pd.read_table(args.inputFile, index_col='index')

design = design.loc[:, [
    'R_nought',
    'transmission_rate',
    'removal_rate',
    'sampling_proba',
    'infectious_period',
    'tree_size'
]]

nb_samples = len(design)
maxTime = args.maxTime

# -----------------------------
# PASTA PARA ÁRVORES
# -----------------------------
TREE_DIR = "simulated_trees"
os.makedirs(TREE_DIR, exist_ok=True)

# -----------------------------
# ARRAYS DE SAÍDA
# -----------------------------
tree_arrays = []
parameter_arrays = []

successful_trees = 0
failed_trees = 0

# -----------------------------
# FUNÇÕES AUXILIARES
# -----------------------------
STOP_REASON = 'stop_reason'
STOP_SAMPLING = 3

def is_valid_for_summary_stats(tree, min_leaves=3):
    if tree is None:
        return False
    if len(tree) < min_leaves:
        return False
    internal_nodes = [n for n in tree.traverse() if not n.is_leaf()]
    for n in internal_nodes:
        if len(n.children) >= 2:
            return True
    return False

def debug_print(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# -----------------------------
# LOOP PRINCIPAL
# -----------------------------
for experiment_id in range(nb_samples):

    params = design.iloc[experiment_id]

    try:
        # 1️⃣ Simula árvore completa
        tr, vector_counter = phy.simulate_bd_tree_gillespie(
            transmission_r=params["transmission_rate"],
            removal_r=params["removal_rate"],
            sampling_p=params["sampling_proba"],
            max_s=int(params["tree_size"]),
            max_t=maxTime
        )

        # Nomear nós
        for i, node in enumerate(tr.traverse("levelorder")):
            node.name = f"n{i}"

        # 2️⃣ Podar folhas não amostradas
        tr_pruned = phy.remove_certain_leaves(
            tr.copy(),
            to_remove=lambda node: getattr(node, STOP_REASON) != STOP_SAMPLING
        )

        # 3️⃣ Verificar validade
        if not is_valid_for_summary_stats(tr_pruned):
            failed_trees += 1
            continue

        # 4️⃣ Codificar estatísticas
        SS_tree, _ = phy.encode_into_summary_statistics(tr_pruned, 1)

        # 5️⃣ Se chegou até aqui, é válida → salvar tudo
        tree_path = os.path.join(TREE_DIR, f"tree_{successful_trees}.nwk")
        tr.write(outfile=tree_path, format=1)

        tree_arrays.append(SS_tree)
        parameter_arrays.append(params.values)
        successful_trees += 1

    except Exception as e:
        failed_trees += 1
        debug_print(f"Experiment {experiment_id} failed: {e}")


# -----------------------------
# SALVAR DATASET
# -----------------------------
if tree_arrays:
    tree_matrix = np.array(tree_arrays)
    parameter_matrix = np.array(parameter_arrays)

    np.savez_compressed(
        "arvores_target.npz",
        tree_data=tree_matrix,
        parameters=parameter_matrix,
        experiment_ids=design.index.values,
        parameter_names=list(design.columns),
        timestamp=datetime.now().isoformat(),
        maxTime=maxTime,
        nb_samples=nb_samples,
        successful_trees=successful_trees,
        failed_trees=failed_trees
    )

    np.save("tree_data.npy", tree_matrix)
    np.save("parameters.npy", parameter_matrix)

    debug_print(f"Dataset salvo com sucesso")
    debug_print(f"Árvores válidas: {successful_trees}")
    debug_print(f"Árvores falhas: {failed_trees}")
else:
    debug_print("Nenhuma árvore válida para summary statistics")

#python summary_statistics2.py 1k.txt 10000