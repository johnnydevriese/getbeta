import numpy
import scipy
import math


def laplacian(n,h):

	L = numpy.zeros([n,n])
	L = L + numpy.diag(-2.0*numpy.ones(n),0)
	L = L + numpy.diag(numpy.ones(n-1),-1)
	L = L + numpy.diag(numpy.ones(n-1),1)
	L = L/h**2 

	return L



