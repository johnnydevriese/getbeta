import numpy
import scipy
import math

#constansts h is reduced plancks constant and m is mass of electron 
m = 9.10938291*10**-31
h = 6.5821*10**-16  


def laplacian(n,h):

	L = numpy.zeros([n,n])
	L = L + numpy.diag(-2.0*numpy.ones(n),0)
	L = L + numpy.diag(numpy.ones(n-1),-1)
	L = L + numpy.diag(numpy.ones(n-1),1)
	L = L/h**2 

	return L

#build the potential function
def Potential(n,h):

	V = numpy.eye(n)
	
	V = V + numpy.diag((h)**2,0)

	return V

#add the two matrices(operators) together
H = L + V 

#multiply by constant 
H = -((h**2)/(2*m)) * H 

return H 

	
