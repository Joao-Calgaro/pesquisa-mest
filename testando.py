from ete3 import Tree
import sys
import ete3
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
        left, right = children   
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

tr = Tree("tree_after_removal.nw", format=1)

print(list(inorder_traversal(tr)))


#print(tr.search_nodes(stop_reason="5")[0].name)


#for node in tr.search_nodes(i_t="2"):
  #print(node.name)

#ancestor = tr.get_common_ancestor(tr.search_nodes(i_t="2"))
#print(ancestor.name)

#chamando o nome do no
#print(tr.search_nodes(name="n503")[0].name)

#data = np.load('C:\\Users\\JPC\Documents\\MESTRADO\\Projeto\\gillespie-2\\mutation_deep_learning.npz',allow_pickle=True)

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
##    if not node.is_leaf():
 #     print(node.name)
#      print(node.get_distance(tr))
#      print(node.features)
#      #print(phy.get_dist_to_root(node))
#      print('-'*30)


#node1 = tr&"n24"
#node2 = tr&"n14"




#
# arvore e os dist_to_root e conferir como está sendo feito as posições dos nós
# a principio, os nos internos estao na primeira metade do vetor(!)
# salva uma arvore, e aplica a função para essa unica arvore.
# descobrir a ordem relacionando os nós com as posições do vetor 




#18/3
# Testar o encode_into_most_recent() e para árvore laderizada e não laderizada, e ver se a propria função ja laderaliza - OK, a função laderaliza, e nao importa a direção do input
# Esperamos q a função laderalize, e q tenha o msm output - OK
# mudar o parametro do ladderize direction = 1 para 0 - OK
# procurar onde tem a função ladderize() no codigo deles - ok?
# talvez o encode_esteja puxando a árvore nao podada, porque nao foi copiado hard - O problema é que ele estava nomeando uma árvore ja nomeada pela função name_tree
# trabalhar com a função de encode dentro desse arquivo, q pode 
# pq eu estou passando a arvore sem os nomes do nós? já que a função atribui para nós nao nomeados
# agora ver se está passando corretamenta. realizar o msm procedimento q fiz semana passada
'''

print('-'*30)

#for node in tr.traverse():
#    print(node.name)

print(tr)
tuple_tree1, _ = phy_autoral.encode_into_most_recent(tr, 1)


tr = Tree("testando_formato_before.nw", format=1)
print('-'*30)
tr.ladderize(direction=0)
print(tr)
tuple_tree0, _ = phy_autoral.encode_into_most_recent(tr, 1)


print(tuple_tree1.equals(tuple_tree0))'''

#for node in inorder_traversal(tr):
#      if not node.is_leaf():
##        print(node.name)
#        print(node.get_distance(tr))
#        print(node.features)
#        #print(phy.get_dist_to_root(node))
#        print('-'*30)


# dentro do encode(anc), printar as nodes q ele está passando de tal modo q ele tenha a msm ordenação do vetor
# extrair o ancestral comum da mutação, e relacionar a posição dele no vetor
# produzir então arvores para estudar o parametro tr_r12

#for node in tr.traverse("postorder"):
#    if node.is_leaf():
#        if hasattr(node, "i_t"):
#            print(node.name, node.i_t)


#ordem_leaf_encode = ['n97', 'n68', 'n56', 'n45', 'n34', 'n89', 'n27', 'n19', 'n21', 'n10', 'n7', 'n8', 'n12', 'n11', 'n5']

#sub_window_slide_2 = [[ordem_leaf_encode[i], ordem_leaf_encode[i+1]] for i in range(len(ordem_leaf_encode) - 1)]


#for n in sub_window_slide_2:
#    node0 = tr&n[0]
#    node1 = tr&n[1]
#    ancestor = tr.get_common_ancestor(node0, node1)
#    print(ancestor.name)

#for node in tr.traverse('postorder'):
#    if not node.is_leaf():
#        if hasattr(node, "i_t"):
#            print(node.name, node.i_t, node.name)
'''
nodes_it_2 = [node.name for node in tr.traverse() 
              if hasattr(node, "i_t") and node.i_t == '2']
#print(nodes_it_2)
if len(nodes_it_2) > 1:
    #caso o nó mutante for unicamente a folha, se passarmos pela a função, o ancestral comum de uma unico nó é a raiz
    try:
        ancestor = tr.get_common_ancestor(nodes_it_2)
        #print(f'Ancestral comum: {type(ancestor.name)}')
    except ete3.coretype.tree.TreeError:
        print('entrou no exception --- não houve mutação')
        ancestor = None
else:
    print('não houve mutação')

#print(ancestor.name)

data = np.load('C:\\Users\\JPC\Documents\\MESTRADO\\Projeto\\gillespie-2\\mutation_deep_learning.npz',allow_pickle=True)
#print(data['tree_data'][0][501] * data['rescale_factor'][0])
'''


# 02-04
# adicionar a lista que contem a POSIÇÃO de mutação de cada árvore OK ATENÇÂO o valor n no .npz será a posição n+1 no vetor CBLV, já que a primeira posição do vetor é 0
# Conferir se os filhos do ancestral da mutação são todos mutantes, e que nenhum outro mutante exista fora OK
# conferir se bate a posição do ancestral com o vetor OK
# calibrar os parametros de surgimento  
# ver na literatura a razão entre os R_0, e ir ajustanto o tr_12 e proportion_11_22, fazer umas 100 arvores com cada ajuste
#       \__ A razão entre os R_0 fica em torno de 1.5 a 3
#       \__
# e registrar número de mutantes, tempo de surgimento da mutação, tempo total da árvore, número de árvores que ocorreram a mutação 
# multiprocessing
# embedding

#print(data['time_of_surgimento_mutacao'] == None)


#08/04
#Conferir o tempo relativo da mutação no caso tr_12 = 0.05, e pq n tem valores proximos de 1? pq parou no 0.8? foi problema de gráfico?
#entender as distribuições das arvores com poucos nós mutantes 
#(EM OUTRO MOMENTO)remover as arvores q tiveram <=5 mutantes, complicado em inferir qualquer coisa
# estabelecer o intervalo (0.003, 0.007)
# Começar a dissertação: método de simulação, codificação CBLV, 
# identificar o ancestral comum relacioando com a poda, talvez não seja salvo o nó mutante propriamente, mas os seus descendentes
# recuperar os parametros de rreferencias dos r0, infectuos time

'''print(tr)
print(tr.children)
total_nodes = len(tr.get_descendants()) 
print(total_nodes)'''

'''for i in range(10):
    try:
        if i == 5:
            tuple_tree, _, leaf_list = phy.encode_into_most_recent(tr, 1)
        t = 5
    except:
        print('saiu erro')
        continue
    
    print('-'*30)
    print(t)
    print(i)
'''
#data = np.load('C:\\Users\\JPC\\Documents\\MESTRADO\\Projeto\\gillespie-2\\dados_gerados\\001_MP_mutation_deep_learning_mp_10000.npz', allow_pickle=True)
#testando = data.keys()
#for i in list(testando):
#    print(i)
#print(data['simulacoes_com_erro'])

#15/04
# implementar o multiprocessing
# fazer tr12 FIXADO {0.003, 0.005, e 0.007}
# fixar a proportion {1.5, 3}
# o que provavelmente está acontecendo é q msm com o tr12 alto, ainda ha mts arvores sem mutações
# uma das hipoteses é q ocorre a mutação e q é retirada na poda. Vamos conferir a quantidade de mutantes pré-poda
#
#para tr_12 0.05, a prof acredita que dentre as 169 arvores que nao ocorreu mutação, na vdd ocorreu mutação, olhar antes da poda para ver 
# 
#                    tr: 0.003   0.005   0.007
#   propotion   1.5     OK        OK        ok
#               2.2     OK      ok      ok
#               3.0  OK     ok


# 22-04 
# o y do grid com a proporção tem q ser normalziado 
# inverter a ordem do grid na proporção, de 0 pra [1,10) para < 10
# confirmar se o grid do tempo relativo, quando a tr12 é baixa, resultara em menos arvores mutantes, 
#assim, a área dos histogramas são diferentes, ver a versão normalizada, onde dividimos peloa quantidade total de arvores mutantes de cada um
################
#1) ESCREVER !!!
#2) Testar a rede neural, o parâmetro q ele estimar está normalziado, temos que multiplicar pelo reascale factor o infectuos time
#3) dai treinar a nossa
#4) por enquanto, achar os parametros R_0_1, infectuos time and a taxa relativa (novas árvores 100.000)
#5) em um segundo momento!, encontrar o momento da mutação