import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#['tree_data', 'rescale_factor', 'internal_nodes_order', 'ancestor', 'mutation_positions', 'number_of_mutantes', 'time_of_surgimento_mutacao', 'total_time_of_simulation']
data = np.load('C:\\Users\\JPC\Documents\\MESTRADO\\Projeto\\gillespie-2\\mutation_deep_learning.npz',allow_pickle=True)

data1 = data['mutation_positions']
data2 = data['number_of_mutantes']
data3 = data['time_of_surgimento_mutacao'] # Já está normalizado
data4 = data['total_time_of_simulation'] # Já está normalizado

tempo_mutacao = data3
tempo_total = data4
result = [x/y for x,y in zip(tempo_mutacao, tempo_total) if x != None]

data5 = [sum(data['time_of_surgimento_mutacao'] != None), sum(data['time_of_surgimento_mutacao'] == None)]


datasets = [data1, data2, result, data5]
titles = ['Índice do Ancestral', 'Nº Mutantes', 'Tempo Relativo da Mutação', 'Proporção de Mutantes']




fig, axes = plt.subplots(2, 2, figsize=(10, 8))
axes = axes.flatten()


for i, ax in enumerate(axes):
    if i == 3:
        sns.barplot(x=['Ocorreu Mutação', 'Não Ocorreu Mutação'], y=data5, ax=ax)
        ax.bar_label(ax.containers[0])
        ax.set_title(titles[i])
        ax.set_xlabel('')
        ax.set_ylabel('')
    else:
        sns.histplot(datasets[i], bins=30, kde=True, ax=ax)
        if i == 3:
            ax.set_xlim(0, 50)
        ax.set_title(titles[i])
        ax.set_xlabel('')
        ax.set_ylabel('')

caminho = input("Nome para salvar: ")
if caminho != '':
    fig.suptitle(f"{caminho}", fontsize=16)
    plt.tight_layout()
    plt.savefig(f'histograms\\{caminho}.png', dpi=300)
    plt.savefig(f'histograms\\{caminho}.pdf', dpi=300)
    

plt.show()


print(sum(data2==1))



#tr_12 = 0.005 teve 145 árvores que 1 (um) nó mutante apenas
#tr_12 = 0.05 teve 124 árvores que 1 (um) nó mutante apenas