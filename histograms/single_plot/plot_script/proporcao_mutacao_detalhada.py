import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

data = np.load(f'C:\\Users\\JPC\\Documents\\MESTRADO\\Projeto\\gillespie-2\\dados_gerados\\mutation_deep_learning_mp_10000.npz', allow_pickle=True)
arr = [sum(data['number_of_mutantes'] == 0),sum((data['number_of_mutantes'] > 0) & (data['number_of_mutantes'] < 10)), sum(data['number_of_mutantes'] >= 10)]
sns.barplot(x=['0 Mutações', '(1,10)', '>=10'], y=arr, palette=['green','darkorange', 'steelblue'])
plt.xlabel("Número de Mutantes")
plt.ylabel("Frequência")
plt.title("Proporção de Mutações")
plt.show()