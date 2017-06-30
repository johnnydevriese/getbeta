import lapbuild
import numpy as np
import matplotlib.pyplot as mpl
import scipy.sparse.linalg
import scipy.linalg as sp 
import numpy.linalg as linalg 
import math





L = int(raw_input("Define the domain of this: ")) #you want this to be a large number (for the love of god) something like 200. 
h = float(raw_input("Define h: ")) #you wan this to be a small number. 


n = math.ceil(2*L/h)-1  
h_new = float(2.0*L/(n+1)) 

A = lapbuild.laplacian(n,h_new)


eigenValues,eigenVectors = scipy.sparse.linalg.eigsh(A,k=5,sigma=0.0,return_eigenvectors = True)

idx = eigenValues.argsort()   
eigenValues = eigenValues[idx]
eigenVectors = eigenVectors[:,idx]

#print(eigenVectors)

T = np.zeros((n+2,5))
T[1:(n+1),:] = eigenVectors


#print(eigenVectors) 

#Z = eigenVectors[:,4]

X = np.linspace(-L,L,n+2)

mpl.plot(X,T,'.') 

mpl.show()







