from ete3 import Tree

leaf_list = []
CBLV = [6, 2, 1, 1, 4, 3, 1, 0, 3]

t0 = Tree()
#t0.add_child(name='a', dist=CBLV[0])
t0.name='a'

t1 = Tree()
A = t1.add_child(name='A', dist=CBLV[1])
b = A.add_child(name='b', dist=CBLV[2])

t2 = Tree()
B = t2.add_child(name='B', dist=CBLV[3])
c = B.add_child(name='c', dist=CBLV[4])

t3 = Tree()
C = t3.add_child(name='C', dist=CBLV[3])
d = C.add_child(name='d', dist=CBLV[4])

t4 = Tree()
D = t4.add_child(name='D', dist=CBLV[3])
e = D.add_child(name='e', dist=CBLV[4])



A.add_child(t0)
B.add_child(A)
C.add_child(c)
B.remove_child(c)
B.add_child(C)
D.add_child(B)



print(t4.get_ascii(show_internal=True))

substituir = Tree('((D,E), (F,G));')


D.add_child(substituir)

print(t4.get_ascii(show_internal=True))