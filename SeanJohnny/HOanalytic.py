def HOAnalytic(x):
    """ Generates the energy spectrum and transition moment matrices for each 
    cartesian direction for the one D harmonic oscillatorn as a function of
    input natural frequency to characterize the harmonic potential.
    Input:
        x = domain
    Output:
        psi, E, mu 
    """
    
    import math
    import matplotlib.pyplot as plt
    import numpy as np
    from numpy.polynomial.hermite import hermval
    
    #Solution information
    NumStates = 100
    
    #Constants
    w = 1.
    m = 1.
    hbar = 1.
    e = 1.
    
    def OneDSol(n, w, x):
        b = np.sqrt(m * w / hbar)
            #So, hermval calls the Hermite polynomial corresponding to a set of
            #coefficients c. If you only want on term of a Hermite polynomial, 
            #you make an array of all 0 with the last term being 1 to grab the
            #nth term
        c = np.zeros(n+1)
        c[-1] = 1
        psi =( (b**2 / math.pi)**0.25 * 
            1 / math.sqrt(2**(n) * math.factorial(n)) * 
            np.exp(-0.5 * b**2 * x**2) * 
            hermval(b * x, c, tensor=False) )
        return psi
        
    psi = [OneDSol(i, w, x) for i in range(NumStates)]
    
    #Plot some solutions
    #for i in range(2):
    #    plt.plot(x, psi[i])
    
    #Checking normalization of the first 5 states
    norm = [np.trapz(np.square(OneDSol(i, w, x)),x) for i in range(5)]
    print norm
    
    #List comprehension method of creating the transition moment matrix
    mu = [e * np.trapz(x * OneDSol(l, w, x) * OneDSol(p, w, x), x )
            for l in range(NumStates) for p in range(NumStates)]
    mu = np.reshape(mu,(NumStates,NumStates))
                   
    E = hbar * w * ( np.arange(NumStates) + 0.5 )
    
    E = E - E[0]    
    
    return psi, E, mu