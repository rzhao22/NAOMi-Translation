import time

def generate_axons(vol_params, axon_params, neur_vol, neur_num, gp_vals, gp_nuc, neur_vol_flag=1):

# [bg_pix,neur_vol,gp_bgvals, axon_params, vol_params] = ...
#                      generate_axons(vol_params, axon_params, neur_vol, ...
#                                               neur_num, gp_vals, gp_nuc)
# This function simulates background components. The inputs to this
# function are:
#   - vol_params  - Struct with parameters for the volume generation
#       .vol_sz       - 3-element vector with the size (in um) of the 
#                       volume to generate (default = 100x100x30um)
#       .min_dist     - Minimum distance between neurons (default = 15um)
#       .N_neur       - Number of neurons to generate (default = 50)
#       .vres         - resolution to simulate volume at (default = 2
#                       samples/um)
#       .N_den        - Width of apical dendrites (default = 10)
#       .N_bg         - Number of background/neuropil components to 
#                       simulate (default = 50)
#       .vol_depth    - Depth of the volume under the brain surface
#                       (default = 100um)
#       .dendrite_tau - Dendrite decay strength exponential distance
#                       (default = 5)
#       .verbose      - Level of verbosity in the output during the volume
#                       generation. Can be 0,1,2. 0 = no text updates, 1 =
#                       some text outputs. 2 = detailed text outputs
#                       (default = 1)
#   - axon_params   - Struct containing parameters for background generation
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
#   - neur_vol    - An 3D array with the overall base fluorescence at each
#                   3D location (includes cells/vasculature)
#   - neur_num    - An array where the k^th neuron's locations in the
#                   volume (both the soma and dendrites) are deliniated by
#                   the value 'k'
#   - gp_vals     - Cell array that contains the locations and fluorescence
#                   values 
#   - gp_nuc      - Locations of the nucleus voxels in the volume
#
# The ouptut to this function is:
#   - neur_vol   - Updated 3D array with the overall base fluorescence at
#                  each 3D location (includes cells/vasculature and now also
#                  background/neuropil)
#   - gp_bgvals  - Cell array with number of rows equal to the number of
#                  background processes, which contains the locations and
#                  values of the processes at those locations
#   - axon_params  - Potentially updated background parameter struct
#   - vol_params - Potentially updated volume parameter struct
#
# 2016 - Adam Charles and Alex Song
#
###########################################################################
  ## Input parsing

  vol_params  = check_vol_params(vol_params)                               # Check volume parameters
  axon_params   = check_axon_params(axon_params)                                 # Check background parameters

  ###########################################################################
  ## Calculate background

  if vol_params.verbose == 1:
    print('Generating background fluorescence.')
  elif vol_params.verbose >1:
    print('Generating background fluorescence...\n')

  bg_pix = (neur_num == 0)                                                # Get the locations where the background can be
  for kk in range(0,gp_nuc.shape[0]):
    bg_pix[gp_nuc[kk,0]] = 0 # CHANGED CELL BRAKCETS {kk,1} TO NORMAL MATRIX HERE         # Remove any zeros that may be in the nuclii
    
  fillnum   = np.rint((axon_params.maxfill)*(axon_params.maxvoxel)*sum(bg_pix.T.flatten()))   # Set the fill-number for filling the background with processes
  volsize   = vol_params.vol_sz*vol_params.vres                              # Extract the size of the neural volume
  N_bg      = vol_params.N_bg                                               # Get the number of neurons in the volume
  gp_bgvals = np.array(N_bg*2).reshape((N_bg,2))                                                  # Initialize an array of background components

  if vol_params.verbose >1:                                                   # Optional verbose output
    print('Initializing volume')
  if(neur_vol_flag):
    neur_vol  = np.zeros(neur_vol.shape)                               # Initialize baseline neural volume
    for kk in range(gp_vals.shape[0]):
      neur_vol(gp_vals[kk,1]) = gp_vals[kk,2]                                 # Initialize the new full neural volume to the set neural fluorescences. Do this for somas...
      if kk <= gp_nuc.shape[0]:
        neur_vol[gp_nuc[kk,1]]  = gp_nuc[kk,2]                              # ... and nuclei (if needed)
      if vol_params.verbose >=1:                                                # Optional verbose output
        print('.')
        
  if vol_params.verbose >1:                                                   # Optional verbose output
    print('\n')

  padsize = axon_params.padsize
  volpad = volsize+2*padsize

  M = np.random.rand(volpad)                         # 
  M(padarray(bg_pix==0,padsize*[1 1 1],false,'both')) = realmax('single');#IDK WHAT TO DO WITH THIS 

  if vol_params.verbose >1:                                                  # Optional verbose output
    tic = time.time()

  j      = 1                                                               # Initialize background process count
  numit2 = 0
  nummax = 10000
  while((fillnum>0)and(j<=N_bg)and(numit2<nummax)):
      bgpts  = []                                                          # Initialize the list of points in the background as an empty vector
      numit2 = 0                                                            # Set up a counter to test for stuck processes with nowhere to grow
      while((len(bgpts)<axon_params.minlength)and(numit2<nummax)):
          numit2 = numit2+1                                                 # Increment counters of number of trials where there is nowhere to grow
          root   = np.ceil((volpad-2)*np.random.rand(1,3)+1)
          while(M[root[2][root[0],root[1]]>(axon_params.fillweight*...
                                                           axon_params.maxvoxel)):
              root = np.ceil((volpad-2)*np.random.rand(1,3)+1)
          ends   = np.ceil(root + 2*axon_params.maxdist*vol_params.vres*...
                                                           (rand(1,3)-0.5)) # 
          if(ends[0]>volpad[0]):
                ends[0] = volpad[0]
          if(ends[1]>volpad[1]):
                ends[1] = volpad[1]
          if(ends[2]>volpad[2]):
                ends[2] = volpad[2]
          if(ends[0]<1):
                ends[0] = 1
          if(ends[1]<1):
                ends[1] = 1
          if(ends[2]<1):
                ends[2] = 1
          bgpts = dendrite_randomwalk2(M,root,ends,axon_params.distsc,...
                           axon_params.maxlen,axon_params.fillweight,...
                                     axon_params.maxvoxel,axon_params.minlength)   # 
      if (bgpts.any()):
        nbranches  = np.max(0,round(axon_params.numbranches + ...
                                             axon_params.varbranches*np.random.randn()))
        for i in range(nbranches):
          bgpts2   = []
          numit = 0
          while(len(bgpts2)<axon_params.minlength and numit<100):
            numit = numit+1
            root  = bgpts[np.ceil(np.random.rand()*len(bgpts))]
            while(root[0] == 1 || root[0] == volpad[0] || root[1] == 1 ||...
                                  root[1] == volpad[1] || root[2] == 1 ||...
                                                       root[2] == volpad[2]):
              root = bgpts[np.ceil(np.random.rand()*len(bgpts))]
              ends=np.ceil(root + np.round(2*axon_params.maxdist*...
                                          vol_params.vres*(np.random.rand(1,3)-0.5)))
            if(ends[0]>volpad[0]):
              ends[0] = volpad[0]
            if(ends[1]>volpad[1]):
              ends[1] = volpad[1]
            if(ends[2]>volpad[2]):
              ends[2] = volpad[2]
            if(ends[0]<1):
                ends[0] = 1
            if(ends[1]<1):
                ends[1] = 1
            if(ends[2]<1):
                ends[2] = 1

            bgpts2 = dendrite_randomwalk2(M,root,ends,axon_params.distsc,...
                              axon_params.maxlength,axon_params.fillweight,...
                                       axon_params.maxvoxel,axon_params.minlength)
        bgpts = np.array(bgpts,bgpts2)
        bgpts = bgpts-padsize
        TMPidxs = (bgpts[:,0]<=0)|(bgpts[:,0]>volsize[0]) ...
            |(bgpts[:,1]<=0)|(bgpts[:,1]>volsize[1]) ...
            |(bgpts[:,2]<=0)|(bgpts[:,2]>volsize[2])
        
        bgpts[TMPidxs,:] = [] #THIS DOES NOT WORK, NEED TO FIND SOMETHING ELSE TO FILL WITH      
        if(bgpts.any()):
           gp_bgvals[j,1] = bgpts[:,0]+(bgpts[:,1]-1)*volsize[0]+...
                                           (bgpts[:,2]-1)*volsize[0]*volsize[1]
          gp_bgvals[j,2] = (1/axon_params.maxel)*np.ones(bgpts.shape[0],1)* ... 
                            np.max(0,1+axon_params.varfill*np.random.randn())
          fillnum = fillnum-bgpts.shape[0]
          if(neur_vol_flag):
            neur_vol[gp_bgvals[j,1]] = neur_vol[gp_bgvals[j,1]]+gp_bgvals[j,2]
          j = j+1;                                                          # Increment background process count
          if vol_params.verbose >1:                                             # Optional verbose output
          if(j%1000)==0):
            toc = tim.time()
            Tdone = tic - toc
            print('#d (#f seconds).\n',j,Tdone)
  if(j>N_bg):
    j = N_bg
    
  vol_params.N_bg = j                                                      # Save the number of background pieces generated
  gp_bgvals  = gp_bgvals[1:j,:]                                       # Only output the number of components actually generated (initialized to a much more optomistic number)

  if not neur_vol_flag:
    neur_vol = []
  if vol_params.verbose >= 1:
    print('done.\n')


  return neur_vol, gp_bgvals, axon_params, vol_params
