def getbetaHO():
    

    import numpy as np 
    import scipy 
    import math 
    import scipy.sparse.linalg
    from matplotlib import pyplot as mpl
    from scipy import integrate 


    w = 1.  #float(raw_input("Define the frequency:")) 
    G = 10. #float(raw_input("Define the domain:")) 
    s = 500. #float(raw_input("Define the number of steps:"))
    #k = 5 #int(raw_input("Define the first eigenvector to compare(must be between 1 & 5):"))
    #j = 5 #int(raw_input("Define the second eigenvector to compare(must be between 1 & 5):"))
    
    #constansts h is reduced plancks constant(h-bar) and m is mass of electron 
    #m = 9.10938291e-31
    #h = 6.5821e-16  
    m = 1.
    h = 1.
    e = 1.0
    # /end constants 
    
    #H = C(d/dx)^2 + Kx^2
    C = -(h**2)/(2*m)  #calculated constant for derivative operator 
    
    K = (0.5)*m*(w**2)  #calculated constant for potential operator 
    
    #calculating the interval of the space(domain) 
    x = np.linspace(0,G,s,endpoint=False)
    
    #building derivative operator 
    L = ( -2.0 * np.eye(s) + np.eye(s,k=-1)+np.eye(s,k=1) ) / (G/s)**2
    
    L = np.multiply(L,C)  
    
    
    #building potential operator 
    V = np.diag(x**2)
    
    V = V * K 
    
    #addition of the two operators
    H = V + L
    
    #computing the first ____ eigenvalues and eigenvectors (k= number to compute) 
    E, psi = np.linalg.eigh(H) 
    #E, psi = scipy.sparse.linalg.eigsh(H,k=50,sigma=0.0,return_eigenvectors=True)
    
    #sorting the first ___ eigenvalues and eigenvectors 
    idx = E.argsort()   
    E = E[idx]
    psi = psi[:,idx]
    
    #normalizing the eigenvectors
    #psi = psi/ math.sqrt(1000) 
    
    #creating probability 
    #prob = abs(psi)**2
    
    #normalizing eigenvectors 
    for i in range(20):
        psi[:,i] = psi[:,i] / math.sqrt(np.trapz( psi[:,i]*psi[:,i], x))
    #print np.trapz( psi[:,0]*psi[:,0], x)
    
    
    #plotting eigenfunctions
    mpl.figure() 
    for i in range(10): 
   	mpl.plot(x, psi[:,i])  
    #mpl.figure() 
    #mpl.plot(E[0:10], '.') 
    mpl.show()
    
    #plotting probability 
    #mpl.figure 
    #for i in range(5): 
    #	mpl.plot(x, prob[:,i]) 
    #mpl.show()
    
    
    ####### Integration of eigenvectors to find transition moments ########
    
    
    ###using trapezoidal to integrate (trapz) (just testing integration of two vectors.) 
    #y = psi[:,2] * x * psi[:,2] 
    
    #x = x #same as x originally defined;just reminding myself 
        
    #y_int = integrate.trapz(y, x)
    
    #print(y_int) 
    
    #### /end testing of integration (note: this is to check because psi[:,1]*psi[:,1] should be one.
    
    
    
    ###############creating matrix of all transition moments################## 
    
        #using for loop method
    #for i in range(2)
    #    for j in range(i)):
    #        matrix[i][j] = integrate.trapz(psi[i,i]) * x * integrate.trapz(psi[j,j]) 
    #print(matrix)
    
    
    ###this section of code works(using list comprehensions) #########
    X = [ np.trapz(psi[:,k] * x * psi[:,j], x) for k in range(20) for j in range(20)]
    X = np.reshape(X, (20,20))
    #print(X)
    ###### /end section that works######
    
    #creating Omega(note that range must be previous range!!)
    Omega = E[:20] - E[0] 
    
    #Omega = range(20)
    
    #print(Omega)
    
    ###Computing Beta###
    #since X is symmetric the vectors of X[1:,0] are the same. (Generally one of them would have to be complex conjugated)
    #print(X[0,1:]) <- this works for slicing 0th(first) column starting at 1st(second) value.) 
    
    Ket = X[1:,0] / Omega[1:]
    
    Bra = X[1:,0] / Omega[1:]
    
    mu = (X[1:,1:]) - (np.eye(19)*(X[0,0]))
        
    Beta = np.dot(Bra,np.dot(mu,Ket))

    #multiplying in the constants  
    Beta = 3 * (e**3) * Beta 
    
    #print(Beta)  
        
    ###compute Beta max###
    
    B_max = (3.0**0.25)*(((e*h)/(m**0.5))**3) * 1./ (Omega[1])**(3.5)
    
    
    #print(B_max) 
    
    #calculate off resonance 
    
    Beta_int = Beta / B_max  
    
    return Beta_int
    
        
 
    