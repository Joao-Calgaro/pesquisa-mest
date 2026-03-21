# Geração da tabela de parâmetros

import pandas as pd
import numpy as np


#def simulate_bdss_tree_gillespie(tr_r11, tr_r12, tr_r21, tr_r22, removal_r, sampling_p, max_s, max_t, fraction_1):
#vou precisar simular R0, gama, a razão do tr22/tr11, tr12, 



n = 1
rnought1_list = []
rnought2_list = []
r_list = []
t_list = []
s_list = []
tree_list = []
proportion_list = []
tr_r11_list = []
tr_r22_list = []
tr_r21_list = []
tr_r12_list = []
fraction_list = []

for i in range(n):

    R_nought1 =  np.round(np.random.uniform(1, 5), 3)
    rnought1_list.append(R_nought1)
    
    # Razão entre tr22 e tr11 (o "quão mais transmissível" é a mutação)
    proportion_tr22_tr11 = np.round(np.random.uniform(1.5, 3), 3)
    proportion_list.append(proportion_tr22_tr11)

    R_nought2 = np.round(proportion_tr22_tr11 * R_nought1, 3)
    rnought2_list.append(R_nought2)

    infectious_time = np.round(np.random.uniform(1, 10), 3)
    t_list.append(infectious_time)

    removal_rate = np.round(1/infectious_time, 3)
    r_list.append(removal_rate)
    
    #sampling_proba = np.round(np.random.uniform(0.01, 1), 3)
    sampling_proba = 1
    s_list.append(sampling_proba)
    
    
    
    tree_size = np.random.randint(50, 200)
    tree_size = 450
    tree_list.append(tree_size)
    

    # Possíveis valores: [0.1; 25]
    tr_r11 = np.round(R_nought1 * removal_rate, 3)
    tr_r11_list.append(tr_r11)

   
    tr_r22 = np.round((R_nought2 * removal_rate), 3)
    tr_r22_list.append(tr_r22)

    tr_r21 = 0
    tr_r21_list.append(tr_r21)

    # a taxa de surgimento da mutação. É necessário ajustar os valores
    # aparentemente, por volta de 0.001 que fica mais realista, ocorre mutação em algumas arvores, outras não
    #tr_r12 = np.round(np.random.uniform(0.001, 0.003), 3) 
    #tr_r12 = 0.002
    tr_r12 = 0.005
    tr_r12_list.append(tr_r12)

    fraction = 1 
    fraction_list.append(fraction)



df = pd.DataFrame({
    "index": np.arange(n),
    "R_nought_1": rnought1_list,
    "R_nought_2": rnought2_list,
    "tr_11":tr_r11_list,
    "tr_22":tr_r22_list,      
    "tr_12":tr_r12_list,
    "tr_21":tr_r21_list,
    "removal_rate": r_list,
    "sampling_proba": s_list,
    "infectious_time": t_list,
    "tree_size": tree_list,
    "proportion_t22_t11": proportion_list,

    "fraction_1": fraction_list
})



#design = design.loc[:, ['R_nought', 'tr_rate_1_1', 'tr_rate_2_2', 'tr_rate_1_2', 'tr_rate_2_1', 'removal_rate',
#                        'sampling_proba', 'R_nought_1', 'R_nought_2', 'R_nought_verif', 'tree_size', 'x_transmission',
#                       'fraction_1', 'infectious_period']]

df.to_csv("tabela_BDSS_mutation.txt", sep="\t", index=False)
