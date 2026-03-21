import pandas as pd
import numpy as np

#R_nought = [1]

R_nought = [1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5] 
removal_rate = [1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2]
transmission_rate = np.multiply(R_nought, removal_rate)
sampling_proba = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]
infectious_time = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2]  
tree_size = [150, 150, 150, 150, 150, 150, 150, 150, 150, 150]



df = pd.DataFrame({
    "index": range(10),
    "R_nought": R_nought,
    "transmission_rate": transmission_rate,
    "removal_rate": removal_rate,
    "sampling_proba": sampling_proba,
    "infectious_time": infectious_time,
    "tree_size": tree_size
})

df.to_csv("input_classificacao_10x.txt", sep="\t", index=False)


