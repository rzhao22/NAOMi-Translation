import numpy as np
import warnings
# nearest_small_prime(N,P)
# return the nearest small prime
def nearest_small_prime(N,P = 7):

    # N = nearest_small_prime(N,P)
    # 
    # Find the closest, larger number to N with no factors greater than P.
    # 
    # 2018 - Adam Charles (from code by Chris Turnes)

    ###########################################################################
    ## Check inputs

    if not P:
        P = 7

    ###########################################################################
    ## Find next largest value with a good factorization

    if np.size(N) > 1:
        for kk in range(np.size(N)):
            N[kk] = nearest_small_prime(N[kk],P)                               # If a vector of values is given, find the 
        
    else:
        if abs(N-round(N))>1e-3:                                               # Test that N is an integer (or close to it)
            raise ValueError('N is not even close to an integer!')             # If not, throw an error
        else:                   
            N = round(N)                                                       # Otherwise make sure that N is technically an integer
        
        if N > 0:
            while max(factor(N)) > P:
                N = N + 1                                                      # Iteratively increase N until the maximum factor is P
            
        else:
            warnings.warn('N is not a positive number')
        
    

    return N

    ###########################################################################
    ###########################################################################
def factor(n):
    factors = []
    for i in range(2,n + 1):
        if n % i == 0:
            count = 1
            for j in range(2,(i//2 + 1)):
                if(i % j == 0):
                    count = 0
                    break
            if(count == 1):
                factors.append(i)
    return set(factors)