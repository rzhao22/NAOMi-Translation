import numpy as np

def expression_varaiation(trace_mat, p_off, min_mod):

##function trace_mat = expression_variation(trace_mat, p_off, min_mod)

# function trace_mat = expression_variation(trace_mat, p_off, min_red)
# 
# Function to take a matrix of traces (N neurons by T time-points) and to
# modulate each neuron's time trace by a value that is 0 with probability
# p_off (that cell has no expression and is invisible to the optics) or by
# a value that is randomly chosen at uniform between min_mod and 1. 
# 
# 2018 - Adam Charles

###########################################################################
## Input parsing

    if (min_mod[0] < 0) or (min_mod[0] > 1):
        ValueError('minimum modulation level must be between zero and one.\n')

    ###########################################################################
    ## Modulate the per-trace levels

    if np.size(trace_mat) > 1:
        N = np.size(trace_mat,0)                                                 # Extract number of cells
    elif np.size(trace_mat) == 1:
        N = trace_mat                                                         # Extract number of cells
    
    if(len(min_mod)==1):
        x = min_mod + (1-min_mod)*np.random.rand(N,1)                                     # Initialize modulation variables as uniform random between min_red and 1
    else:
        x = np.random.gamma(min_mod[1],min_mod[0],N,1)
        #   x = min_mod[0] + abs(min_mod[1]*randn(N,1)+1-min_mod[0])                # Initialize modulation variables as uniform random between min_red and 1  
    

    x = x*(p_off < np.ranodm.rand(N,1))                                                # Each cells to zero fluorescence with probability p_off
    if np.size(trace_mat) > 1:
        trace_mat = trace_mat * x                              # Modulate the activity of each cell (the rows of trace_mat)
    elif np.size(trace_mat) == 1:
        trace_mat = x                                                         # Otherwise just output the multiplicative values
    
    return trace_mat

###########################################################################
###########################################################################