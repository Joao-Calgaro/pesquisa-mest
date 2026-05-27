import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

data = np.load(f'C:\\Users\\JPC\\Documents\\MESTRADO\\Projeto\\gillespie-2\\dados_gerados\\mutation_deep_learning_mp_10000.npz', allow_pickle=True)

tempo_mutacao = data['time_of_surgimento_mutacao'][(data['number_of_mutantes'] >=10)].ravel() 
tempo_total = data['total_time_of_simulation'][(data['number_of_mutantes'] >=10)].ravel()
result = [x/y for x,y in zip(tempo_mutacao, tempo_total) if x != None]

sns.histplot(result, bins=30, kde=False)

plt.xlabel("Tempo Relativo da Mutação")
plt.ylabel("Frequência")    
plt.title("Tempo Relativo da Mutação (>=10 Mutantes)")
plt.show()