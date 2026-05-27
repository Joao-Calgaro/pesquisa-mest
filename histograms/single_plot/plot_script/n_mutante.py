import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

data = np.load(f'C:\\Users\\JPC\\Documents\\MESTRADO\\Projeto\\gillespie-2\\dados_gerados\\mutation_deep_learning_mp_10000.npz', allow_pickle=True)

arr = data['number_of_mutantes'][(data['number_of_mutantes'] >=10)].ravel()


sns.histplot(arr, bins=30, kde=True)


plt.xlabel("Número de Mutantes (>=10)")
plt.ylabel("Frequência")
plt.title("")
plt.show()


