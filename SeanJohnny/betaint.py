def BetaInt(mu, E):
    """Calculates the off-resonant (zero frequency limit) intrinsic first 
    hyperpolarizability in one dimension.
    Input: 
        mu = transition moment matrix where mu[i,j] = e<i|x|j>
        E = energy eigenvalues biased by the ground state where E[i] = E_i - E_0
    """
    import numpy as np
    
    #constants
    hbar = 1.
    m = 1.
    e = 1.
    
    beta = 3 * np.dot( mu[0,1:], 
                np.dot( mu[1:,1:] - np.eye(len(mu)-1)*mu[0,0], 
                        mu[1:,0]/E[1:] ) / E[1:] )
    
    betamax = 3**0.25 * (e * hbar / np.sqrt(m))**3 / E[1]**3.5
    
    return beta/betamax