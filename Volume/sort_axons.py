import numpy as np
import check_axon_params as check_axon_params
import check_vol_params as check_vol_params

def sort_axons(vol_params, axon_params, gp_bgvals, cell_pos):

# sort_axons(vol_params, axon_params, gp_bgvals, cell_pos)
# return bg_proc
# 
# This function ranndomly sorts length(gp_bgvals) axons into N_proc bins
#
# The inputs to this function are:
#   - vol_params - Class instance containing parameters for the volume generation
#       .vol_sz   - 3-element vector with the size (in um) of the volume to
#                   generate (default = 100x100x30um)
#       .min_dist - Minimum distance between neurons (default = 15um)
#       .N_neur   - Number of neurons to generate (default = 50)
#       .vres     - resolution to simulate volume at (default = 2
#                   samples/um)
#       .N_den    - Width of apical dendrites (default = 10)
#       .N_bg     - Number of background/neuropil components to simulate
#                   (default = 50)
#       .vol_depth- Depth of the volume under the brain surface (default =
#                   100um)
#       .verbose  - Level of verbosity in the output during the volume
#                   generation. Can be 0,1,2. 0 = no text updates, 1 = some
#                   text outputs. 2 = detailed text outputs (default = 1)
#   - axon_params   - Class instance containing parameters for background generation
#       .distsc      - Parameter to determine how directed the random walk 
#                      to generate background processes is Higher values
#                      are more directed/less random (default = 0.5)
#       .fillweight  - Maximum length for a single process branch (default
#                      = 100 um)
#       .maxlength   - Maximum length for background processes (default =
#                      200 um) 
#       .minlength   - Minimum length for background processes (default =
#                      10 um) 
#       .maxdist     - Maximum distance to the end of a process (default =
#                      100 um) 
#       .maxel       - Max number of axons per voxel (default = 8)
#       .numbranches - Number of allowable branches for a single process
#                      (default = 20) 
#       .varbranches - Standard deviation of the number of branches per
#                      process (default = 5) 
#       .maxfill     - Voxel maximum occupation: fraction of volume that is
#                      to be filled by background processes (default = 0.7)
#       .N_proc      - Number of background components (default = 10)
#       .l           - Gaussian process length scale for correllation
#                      background processes over the image (default = 25) 
#       .rho         - Gaussian process variance parameter for correllation
#                      background processes over the image (default = 0.1) 
#   - gp_bgvals     - Nested list with number of rows equal to the number of
#                     background processes, which contains the locations and
#                     values of the processes at those locations
#   - cell_pos         - Position of cells within the volume
#
# 
# The ouptut to this function is:
#     bg_proc - The volume Class instance provided by simulate_neural_vol,
#               modified to have the background components.
# 
# 2017 - Adam Charles and Alex Song
#
###########################################################################
## Parse inputs

    vol_params  = check_vol_params(vol_params)                                # Check volume parameters
    axon_params   = check_axon_params(axon_params)                            # Check background parameters

    ###########################################################################
    ## Calculate background

    N_proc    = axon_params.N_proc                                            # Extract the number of correlated background components to create
    vol_sz    = vol_params.vol_sz*vol_params.vres
    if (vol_params.verbose > 0):                                                  # Optional verbose output
        print('Sorting axons...')
    

    bg_proc  = [[None for j in range(2)] for i in range(N_proc)]                                                 # Initialize the cell array to store the background processes
    if(len(bg_proc)>vol_params.N_neur+vol_params.N_den):
        N_comps = vol_params.N_neur+vol_params.N_den
        gp_bgpos = np.zeros((len(gp_bgvals),3))
        for kk in range(len(gp_bgvals)):
            del TMP_pos
            if gp_bgvals[kk][0]:
                [TMP_pos[:,0],TMP_pos[:,1],TMP_pos[:,2]] = np.unravel_index(gp_bgvals[kk][0], vol_sz)
                gp_bgpos[kk,:] = np.mean(TMP_pos)
            
        
    
        cell_pos2 = cell_pos[:N_comps,:]
        dist_mat = np.sqrt(((cell_pos2[:,0] - gp_bgpos[:,0].H)**2) +
            ((cell_pos2[:,1] - gp_bgpos[:,1].H)**2) + ((cell_pos2[:,2] - gp_bgpos[:,2].H)**2))
        
        idxlist = np.zeros((N_comps,1))
        for ii in range(N_comps):
            [_,idx] = min(dist_mat[ii,:])
            dist_mat[:,idx] = np.Inf
            bg_proc[ii][0] = gp_bgvals[idx][0]
            bg_proc[ii][1] = gp_bgvals[idx][1]
            idxlist[ii] = idx
        

        for kk in range(len(gp_bgvals)):
            if(kk not in idxlist):
                index = N_comps+np.ceil((N_proc-N_comps)*np.random.rand)
                bg_proc[index][0] = np.concatenate(bg_proc[index][0],gp_bgvals[kk][0], axis=0)
                bg_proc[index][1] = np.concatenate(bg_proc[index][1],gp_bgvals[kk][1], axis=0)
            
        
    else:
        for kk in range(len(gp_bgvals)):
            index = np.ceil(N_proc*np.random.rand)
            bg_proc[index][0] = np.concatenate(bg_proc[index][0],gp_bgvals[kk][0], axis=0)
            bg_proc[index][1] = np.concatenate(bg_proc[index][1],gp_bgvals[kk][1], axis=0)
        
    
    if (vol_params.verbose > 0):    
        print('done.\n')
    

    
    return bg_proc

###########################################################################
###########################################################################