from pylab import * 

r = int(raw_input("How many rows?"))

c = int(raw_input("how many columns?"))

A = randn(r,c) 

#!print A  

T = zeros((r,c)) 
	
#print T 

for i in range(0, c): 
	T[:,i] = A[:,i]

	for j in range(0,i): 
	
		T[:,i] = T[:,i]- dot(A[:,i],T[:,j])*T[:,j] 

	T[:,i] = T[:,i]/norm(T[:,i]) 

print T 

I = dot(transpose(T),T) 

#!print I 

print norm(I- eye(c)) 



