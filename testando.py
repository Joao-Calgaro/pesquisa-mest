from ete3 import Tree
import sys
import pandas as pd
import numpy as np
import phy

#t = Tree( "((H:0.3,I:0.1):0.5, A:1, (B:0.4,(C:0.5,(J:1.3, (F:1.2, D:0.1):0.5):0.5):0.5):0.5);" )
#print(t)
#ancestor = t.get_common_ancestor("C", "J", "B")
#print(ancestor)

def inorder_traversal(node):
    children = node.children

    # leaf node
    if len(children) == 0:
        yield node
        return

    # binary node
    if len(children) == 2:
        right, left = children
        yield from inorder_traversal(left)
        yield node
        yield from inorder_traversal(right)
        return

    # single child (rare but possible)
    if len(children) == 1:
        yield from inorder_traversal(children[0])
        yield node
        return

    # if more than two children → not binary
    raise ValueError("Tree is not binary. Resolve polytomies first.")

tr = Tree("testando_formato_after.nw", format=1)




#print(tr.search_nodes(stop_reason="5")[0].name)


#for node in tr.search_nodes(i_t="2"):
  #print(node.name)

#ancestor = tr.get_common_ancestor(tr.search_nodes(i_t="2"))
#print(ancestor.name)

#chamando o nome do no
#print(tr.search_nodes(name="n503")[0].name)

data = np.load('C:\\Users\\JPC\Documents\\MESTRADO\\Projeto\\gillespie-2\\mutation_deep_learning.npz',allow_pickle=True)

#print(data['tree_data'][0])

'''
t.populate(15)
print t
print t.children
print t.get_children()
print t.up
print t.name
print t.dist
print t.is_leaf()
print t.get_tree_root()
print t.children[0].get_tree_root()
print t.children[0].children[0].get_tree_root()
# You can also iterate over tree leaves using a simple syntax
for leaf in t:
  print leaf.name
'''
#for node in tr.traverse("postorder"):
#  print(node.name)

#tr.write(outfile="testando_formato_after.nw", format=3)




# entender o postorder traveral                                                     OK  
# entender como se relaciona com a lista do CBLV                                    OK
# Confirmar q de fato é o método CBLV                                               OK
# Encontrar o nome do nó da transmissao a partir do número do nó na lista do CBLV


# get_dist_to_root de cada nó interno, e comparar com o cblv para conferir como está sendo a ordem no CBLV
# tentar passar pela arvore pelo metodo inorder
# dentro do encode() no ponto phy, tentar anotar por quais nós ele está passando. Essa lista é a ordem que será salva no CBLV, e com isso,
# sabemos a posição que será salva  

#print(data['tree_data'][0])

#print(str(tr))

'''
print('Essas são as folhas da árvore:')
for leaf in tr.traverse("postorder"):
    if leaf.is_leaf():
        print(leaf.name)

print('-'*30)
print('Esses são os nós internos da árvore:')
#for node in tr.traverse("inorder"):
#    if not node.is_leaf():
#        print(node.name)
'''







#for node in inorder_traversal(tr):
#    if not node.is_leaf():
#      print(node.name)
#      print(node.get_distance(tr))
#      #print(phy.get_dist_to_root(node))
      


node1 = tr&"n24"
node2 = tr&"n14"




#
# arvore e os dist_to_root e conferir como está sendo feito as posições dos nós
# a principio, os nos internos estao na primeira metade do vetor(?)
# salva uma arvore, e aplica a função para essa unica arvore.
# descobrir a ordem relacionando os nós com as posições do vetor 




#18/3
# Testar o encode_into_most_recent() e para árvore laderizada e não laderizada, e ver se a propria função ja laderaliza
# Esperamos q a função laderalize, e q tenha o msm output
# mudar o parametro do ladderize direction = 1 para 0
# procurar onde tem a função ladderize() no codigo deles
# talvez o encode_esteja puxando a árvore nao podada, porque nao foi copiado hard 
# trabalhar com a função de encode dentro desse arquivo, q pode 
# pq eu estou passando a arvore sem os nomes do nós? já que a função atribui para nós nao nomeados


tr.ladderize(direction=0)
tuple_tree, _ = phy.encode_into_most_recent(tr, 1)

for node in tr.traverse():
    print(node.features)
