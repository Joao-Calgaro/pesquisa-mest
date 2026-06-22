import numpy as np

cblv = np.load("testando_com_uma_arvore\\tree_data.npy")

import numpy as np
from ete3 import Tree


def cblv_to_ete(x):
    x = np.asarray(x).reshape(-1)

    internals = x[1:499]
    leaves = x[500:1000]

    root = Tree()
    root.dist = 0.0

    current = root

    n = min(len(internals), len(leaves))

    for i in range(n):

        # create internal node (always extend spine)
        internal = current.add_child()
        internal.dist = float(internals[i])

        # attach leaf
        leaf = internal.add_child(name=f"L{i}")
        leaf.dist = float(leaves[i])

        # move forward along backbone
        current = internal

    return root

print(cblv_to_ete(cblv))