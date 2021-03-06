import lapbuild
import numpy as np
import matplotlib.pyplot as mpl
import scipy.sparse.linalg
import scipy as sp 


L = int(raw_input("Define the domain of this: "))
h = float(raw_input("Define h: "))


n = int(L/h)+1 
n = 2*n 

L = lapbuild.laplacian(n,h)


eig_values, eig_vectors = scipy.sparse.linalg.eigsh(L) 

idx = eig_values.argsort()   
eigenValues = eig_values[idx]
eigenVectors = eig_vectors[:,idx]


Z =  eigenVectors[:,:5]


mpl.plot(Z) 

mpl.show()




#Q = Z.T
#V = abs(numpy.dot(Q,Z)-numpy.eye(n)) 
 
#print("maximum orthogonality error:")
#print(V.max()) 

#mpl.plot(Z[:,0:5],'.') 


#mpl.show() 
