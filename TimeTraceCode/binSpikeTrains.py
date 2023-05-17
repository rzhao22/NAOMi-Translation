import numpy as np

def binSpikeTrains(evt, evm, N_node, dt, T):

##function S = binSpikeTrains(evt, evm, N_node, dt, T)

# function S = binSpikeTrains(evt, evm, N_node)
# 
# Function that takes a list of marked events (times in evt and marks that
# are integers in [1, N_node] in evm) and makes a spike count matrix S
# where each element is the number of spikes for each neuron (the rows) and
# a time-bin of size dt for each of the T times (the columns).
# Inputs:
#  - evt    - List of event times
#  - evm    - List of event marks (neuron numbers)
#  - N_node - Total number of possible marks
#  - dt     - Bin size (time-span) of each time-pooint
#  - T      - Total number of time-points
#
# Outputs:
#  - S      - N_node-by-T matrix of binned spike counts
#
# 2017 - Adam Charles

###########################################################################
## Input parsing

    if np.size(evt)!=np.size(evm):
        ValueError('Number of events must be the same as the number of marks!')

    if max(evm) > N_node:
        ValueError('Total number of marks must be greater than the highest seen mark!')
    

    if np.ceil(max(evt)/dt) > T:
        ValueError('Total time-span must be longer than the latest seen event!')

    ###########################################################################
    ## Get bin counts

    S = np.zeros((N_node, T))                                                    # Initialize count matrix
    for kk in range(np.size(evt)):                                                      # Iterate over events
        S[evm(kk), np.ceil(evt(kk)/dt)] = S(evm(kk), np.ceil(evt(kk)/dt)) + 1       # Add a spike to the location/mark of the event
    
    return S

###########################################################################
###########################################################################