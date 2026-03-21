# Geração da tabela de parâmetros

import pandas as pd
import numpy as np

n = 20
rnought_list = []
r_list = []
t_list = []
s_list = []
tree_list = []
tra_list = []

for i in range(n):

    R_nought = np.random.choice([1, 5])
    rnought_list.append(R_nought)
    
    removal_rate = np.round(np.random.uniform(0.02, 0.5), 3)
    r_list.append(removal_rate)
    
    sampling_proba = np.round(np.random.uniform(0.1, 0.6), 3)
    s_list.append(sampling_proba)
    
    infectious_time = np.round(1 / removal_rate, 3)
    t_list.append(infectious_time)
    
    tree_size = np.random.randint(100, 200)
    tree_list.append(tree_size)
    
    transmission_rate = np.round(R_nought * removal_rate, 3)
    tra_list.append(transmission_rate)

df = pd.DataFrame({
    "index": np.arange(n),
    "R_nought": rnought_list,
    "transmission_rate": tra_list,
    "removal_rate": r_list,
    "sampling_proba": s_list,
    "infectious_time": t_list,
    "tree_size": tree_list
})

df.to_csv("input_classificacao.txt", sep="\t", index=False)

#r_nought = transmission_rate / removal_rate
#transmission_rate = r_nought * removal_rate

#removal rate mais tempos maiores e menos tempos menores

#https://transportgeography.org/contents/applications/transportation-pandemics/basic-reproduction-number-r0-of-major-infectious-diseases/
import matplotlib.pyplot as plt  # usar pyplot

# Supondo que df seja seu DataFrame
#plt.hist(df['removal_rate'], bins=20, color='skyblue', edgecolor='black')
#plt.xlabel('removal_rate')
#plt.ylabel('Frequency')
#plt.title('Histogram of removal_rate')
#plt.show()

#gerar dados
#gerar estatisticas descritivas


# classificacaoIF gera os parâmetros
# para gerar as arvores, utilizar o comando no terminal:
# python treegeneratortest.py input_classificacao.txt 10000