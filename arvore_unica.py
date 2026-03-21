import pandas as pd
import numpy as np

#R_nought = [1]

R_nought = [1.5, 2, 2.5, 3, 3.5, 4, 4.5, 1.5, 2, 2.5, 3, 3.5, 4, 4.5]
removal_rate = [1/2, 1/2, 1/2, 1/3, 1/3, 1/4, 1/4, 1/5, 1/5, 1/6, 1/7, 1/8, 1/8, 1/8]
transmission_rate = np.multiply(R_nought, removal_rate)
sampling_proba = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]
infectious_time = [2, 2, 2, 3, 3, 4, 4, 5, 5, 6, 7, 8, 8, 8]  
tree_size = [150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150]



df = pd.DataFrame({
    "index": range(14),
    "R_nought": R_nought,
    "transmission_rate": transmission_rate,
    "removal_rate": removal_rate,
    "sampling_proba": sampling_proba,
    "infectious_time": infectious_time,
    "tree_size": tree_size
})

df.to_csv("input_classificacao_grupo.txt", sep="\t", index=False)


