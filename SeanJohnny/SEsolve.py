def HOFinDif(x):
    import numpy as np
    from matplotlib import pyplot as plt
    
    NumStates = 40
    
    hbar = 1.
    m = 1.
    w = 1.
    e = 1.
    
    N = len(x)
    dx = x[1] - x[0]
    
    #The potential!
    V = 0.5 * m * w**2 * np.diag(x**2)
    
    T = (-2 * np.eye(N) + np.eye(N, k=1) + np.eye(N, k=-1) ) / (dx)**2
    H = -hbar**2 / ( 2 * m ) * T + V
    
    E, psi = np.linalg.eigh(H)
    
    
    for i in range(NumStates):
        psi[:,i] = psi[:,i] / np.sqrt(np.trapz( psi[:,i]*psi[:,i], x))
    
    #plt.figure()
    #for i in range(2):
    #    plt.plot(x, psi[:,i])
    #plt.figure()
    #plt.plot(E[0:10], '.')
    #plt.show()
    
    mu = [e * np.trapz(x * psi[:,l] * psi[:,p], x )
            for l in range(NumStates) for p in range(NumStates)]
    mu = np.reshape(mu,(NumStates,NumStates))
    
    E = E[1:NumStates] - E[0]    
    psi = psi.transpose()
    
    return psi, E, mu