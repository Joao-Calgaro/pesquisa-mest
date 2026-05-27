import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

#['tree_data', 'rescale_factor', 'internal_nodes_order', 'ancestor', 'mutation_positions', 'number_of_mutantes', 'time_of_surgimento_mutacao', 'total_time_of_simulation']
arquivos = ['tr12 0.003 proportion 1.5.npz', 'tr12 0.003 proportion 2.2.npz', 'tr12 0.003 proportion 3.0.npz',
             'tr12 0.005 proportion 1.5.npz', 'tr12 0.005 proportion 2.2.npz', 'tr12 0.005 proportion 3.0.npz',
             'tr12 0.007 proportion 1.5.npz', 'tr12 0.007 proportion 2.2.npz', 'tr12 0.007 proportion 3.0.npz']
titulo_colunas = ['Proportion 1.5', 'Proportion 2.2', 'Proportion 3.0']
titulo_linhas = ['tr12 = 0.003\nFrequência', 'tr12 = 0.005\nFrequência', 'tr12 = 0.007\nFrequência']


fig, axes = plt.subplots(3, 3, figsize=(10, 10), sharey=True)
fig.suptitle('Proporção de Mutações', fontsize=16)
for idx, file in enumerate(arquivos):
    row = idx // 3
    col = idx % 3
    ax = axes[row, col]

    data = np.load(f'C:\\Users\\JPC\\Documents\\MESTRADO\\Projeto\\gillespie-2\\histograms\\grid\\{file}', allow_pickle=True)

    
    #arr = data['number_of_mutantes'][(data['number_of_mutantes'] >=10)].ravel()
    #tempo_mutacao = data['time_of_surgimento_mutacao'][(data['number_of_mutantes'] >=10)].ravel() 
    #tempo_total = data['total_time_of_simulation'][(data['number_of_mutantes'] >=10)].ravel()
    #result = [x/y for x,y in zip(tempo_mutacao, tempo_total) if x != None]
    arr = [sum(data['number_of_mutantes'] == 0),sum((data['number_of_mutantes'] > 0) & (data['number_of_mutantes'] < 10)), sum(data['number_of_mutantes'] >= 10)]

    #sns.histplot(result, bins=30, kde=False, ax=ax)
    sns.barplot(x=['0 Mutações', '(1,10)', '>=10'], y=arr, ax=ax, palette=['green','darkorange', 'steelblue'])
    ax.tick_params(labelsize=7)

    if row == 0:
        ax.set_title(titulo_colunas[col], fontsize=10)
    #if row == 2:  
    #    ax.set_xlabel('Tempo Relativo da Mutação', fontsize=10)
    if col == 0: 
        ax.set_ylabel(titulo_linhas[row], fontsize=10)

salvar = input('Deseja salvar? (y): ')
if salvar == 'y':
    #caminho = input("Nome para salvar: ")
    caminho = 'proporcao_mutacao_grid'
    plt.tight_layout()
    plt.savefig(f'histograms\\grid\\{caminho}.png', dpi=300)
    plt.savefig(f'histograms\\grid\\{caminho}.pdf', dpi=300)
plt.show()
plt.tight_layout()
