import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

#['tree_data', 'rescale_factor', 'internal_nodes_order', 'ancestor', 'mutation_positions', 'number_of_mutantes', 'time_of_surgimento_mutacao', 'total_time_of_simulation']

data = np.load(f'C:\\Users\\JPC\\Documents\\MESTRADO\\Projeto\\gillespie-2\\histograms\\10000_test_plot\\tr_12(0.05)--infectuos(1-10).npz', allow_pickle=True)


num_of_mut = data['number_of_mutantes'][data['number_of_mutantes'] > 5].ravel()

temp_surg = data['time_of_surgimento_mutacao'][data['number_of_mutantes'] > 5].ravel()
temp_total = data['total_time_of_simulation'][data['number_of_mutantes'] > 5].ravel()
tempo_relativo = [x/y for x,y in zip(temp_surg, temp_total) if x != None]

proportions = [sum(data['number_of_mutantes'] <= 5), sum(data['number_of_mutantes'] > 5)]


fig, axes = plt.subplots(1, 2, figsize=(15, 5), sharey=False)
fig.suptitle('Transmission Rate 1_2 = 0.05', fontsize=16)


sns.barplot(x=["≤ 5 Type 2 Nodes", "> 5 Type 2 Nodes"], y=proportions, ax=axes[0])
axes[0].grid(alpha=0.3)
axes[0].spines["top"].set_visible(False)
axes[0].spines["right"].set_visible(False)
axes[0].set_xlabel('Count of Mutated Trees', fontsize=10)
#-----------------
num_of_mut = data['number_of_mutantes'][data['number_of_mutantes'] > 5].ravel()
sns.histplot(num_of_mut, bins=30, edgecolor="black", linewidth=0.8, kde=True, ax=axes[1])
axes[1].grid(alpha=0.3)
axes[1].spines["top"].set_visible(False)
axes[1].spines["right"].set_visible(False)
axes[1].set_xlabel('Number of Mutated Nodes (> 5)', fontsize=10)

#-----------------
temp_surg = data['time_of_surgimento_mutacao'][data['number_of_mutantes'] > 5].ravel()
temp_total = data['total_time_of_simulation'][data['number_of_mutantes'] > 5].ravel()
tempo_relativo = [x/y for x,y in zip(temp_surg, temp_total) if x != None]

#sns.histplot(tempo_relativo, bins=30, edgecolor="black", linewidth=0.8, kde=True,ax=axes[2])
#axes[2].grid(alpha=0.3)
##axes[2].spines["top"].set_visible(False)
#axes[2].spines["right"].set_visible(False)
#axes[2].set_xlabel('Relative Time of Mutation (> 5)', fontsize=10)
#----------------



for i, value in enumerate(proportions):
    axes[0].text(i, value + max(proportions) * 0.01, f"{value:,}", ha="center", va="bottom")

plt.tight_layout()
#plt.savefig("histograms\\10000_test_plot\\tr_12(0.05)_1x2.png", dpi=300, bbox_inches="tight")
#plt.savefig("histograms\\10000_test_plot\\tr_12(0.05)_1x2.pdf", bbox_inches="tight")
#plt.show()

tamanho = np.load('C:\\Users\\JPC\\Documents\\MESTRADO\\Projeto\\gillespie-2\\histograms\\grid\\data_grid\\tr12 0.003 proportion 1.5.npz', allow_pickle=True)
print(tamanho['number_of_mutantes'].shape[0])