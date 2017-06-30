# this function computes the eigenvalues and eigenvectors of the 1d discrete Laplacian
# on the interval -L to L with step size h and plots the 5 with smallest eigenvalues

# load packages
import numpy as np
import scipy as sp
from scipy import sparse
import scipy.sparse.linalg
import fdiffsubroutines as fd
import pylab
import math


# input L and h
L = float(raw_input("Enter L:"))
h = float(raw_input("Enter h:"))

# compute N to nearest integer value
N = math.ceil(2*L/h)-1
hn = 2.0*L/(N+1)

#compute eigenvalues and eigenvectors
D,V = sparse.linalg.eigsh(fd.lap1d0(N,hn),k=5,sigma=0.0,return_eigenvectors=True)

# sort eigenvalues
minm = 10000000000000.0
for ii in range(0,4):
	pos = ii
	for jj in range(ii,5):
		if(D[jj] < minm):
			minm = D[jj]
			pos = jj

	if(pos != ii):
		teig = D[ii]
		temp = V[:,ii]

		D[ii] = D[pos]
		V[:,ii] = V[:,pos]

		D[pos] = teig
		V[:,pos] = temp



# plot k smallest eigenfunctions
k = 3
X = np.linspace(-L,L,N+2)
Y = np.zeros([N+2,k])
Y[1:(N+1),:] = V[:,0:k]**2

pylab.plot(X,Y,'.')

pylab.show()
