import lapbuild
import numpy
import matplotlib.pyplot as mpl

n = int(raw_input("Define how many point to use: "))
h = float(raw_input("Define h: "))


L = lapbuild.laplacian(n,h)

#print(L)

E = numpy.linalg.eigh(L)  

Y = E[0]
Z = E[1] 

#print("These are the eigenvalues:")
#print(Y)
print("These are the eigenvectors:")
print(Z) 

Q = Z.T	
V = abs(numpy.dot(Q,Z)-numpy.eye(n)) 
 
#print("maximum orthogonality error:")
#print(V.max()) 

mpl.plot(Z[:,0:5],'-') 


mpl.show() 



