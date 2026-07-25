from ete3 import Tree
import numpy as np

def right_postorder(node):
    for child in reversed(node.children):
        yield from right_postorder(child)
    yield node

    
data = np.load('C:\\Users\\JPC\\Documents\\MESTRADO\\Projeto\\gillespie-2\\testando_1_tree_CBLV\\tree_data.npy', allow_pickle=True)
CBLV = data[0]
int_nodes = CBLV[:501]
ext_nodes = CBLV[501:]



t0 = Tree()
t0.dist=ext_nodes[0]

t1 = Tree()
t1.dist=int_nodes[1]
b = t1.add_child(dist=ext_nodes[1])

t1.add_child(t0)

t_prev = t1


for i in range(2, 19): # len(int_nodes)):
    print('-'*30)

    t_new = Tree()
    t_new.dist = int_nodes[i]
    t_new.add_child(dist=ext_nodes[i])    
    print('arvore nova:')
    print(t_new.get_ascii(show_internal=True, attributes=["dist"]))
    print('arvore velha')
    print(t_prev.get_ascii(show_internal=True, attributes=["dist"]))

    
    if t_new.dist < t_prev.dist:
        t_new.add_child(t_prev)
        print(t_new.get_ascii(show_internal=True, attributes=["dist"]))
        t_prev = t_new
        
    else:
        for node in t_prev.traverse("postorder"):
            print('nó a ser removido')
            print(node.dist)
            first_internal = node.up
            print('pai do no:')
            print(first_internal.dist)
            break
       
        for node in first_internal.traverse("postorder"): 
            print('arvore nova')
            print(t_new.get_ascii(show_internal=True, attributes=["dist"]))
            print('arvore velhaa')
            print(t_prev.get_ascii(show_internal=True, attributes=["dist"]))
            
            print('add o nó na new que vamos remover da principal')
            t_new.add_child(dist=node.dist)
            print(t_new.get_ascii(show_internal=True, attributes=["dist"]))
            

            idx = first_internal.children.index(node)
            print("agora, é para ter removido o 3.65")
            first_internal.remove_child(node)
            print(t_prev.get_ascii(show_internal=True, attributes=["dist"]))

            print('add a t_new no t_prev')
            first_internal.children.insert(idx, t_new)
            print(t_prev.get_ascii(show_internal=True, attributes=["dist"]))
            
            break
        

        

        

    
    
        



        



#print(t_prev.get_ascii(show_internal=True, attributes=["dist"]))

