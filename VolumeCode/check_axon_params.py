import numpy as np
from setParams import setParams
def check_axon_params(axon_params):

##function axon_params = check_axon_params(axon_params)

# axon_params = check_axon_params(axon_params)
#  
# This function checks the elements of the struct axon_params to ensure
# that all fields are set. Fields not set are set to a list of default
# parameters. The struct checked is:
# 
#   - axon_params   - Struct containing parameters for axon generation
#       .flag        - Flag for generation of background dendrites
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
#       .maxvoxel    - Max number of allowable axon elements per voxel
#                      (default = 6)
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
#
# 2017 - Adam Charles and Alex Song
#
###########################################################################
## Run the checks

    if not axon_params:                                                    # Make sure that axon_params is a struct
        del axon_params
        class Axon_params:

            def __init__(self,flag,distsc,fillweight,maxlength,minlength,maxdist,maxel,varfill,maxvoxel,padsize,numbranches,varbranches,maxfill,N_proc,l,rho):
                self.flag        = flag
                self.distsc      = distsc
                self.fillweight  = fillweight
                self.maxlength   = maxlength
                self.minlength   = minlength
                self.maxdist     = maxdist
                self.maxel       = maxel
                self.varfill     = varfill                                               # variation in filling weight (std around 1)
                self.maxvoxel    = maxvoxel                                                   # Maximum number of elements per voxel
                self.padsize     = padsize                                                  # Background padding size (for smoothness in background)
                self.numbranches = numbranches
                self.varbranches = varbranches
                self.maxfill     = maxfill
                self.N_proc      = N_proc
                self.l           = l
                self.rho         = rho
        axon_params = Axon_params(1,0.5,100,200,10,100,8,0.3,6,20,20,5,0.5,10,25,0.1)
    

    # dParams.flag        = 1
    # dParams.distsc      = 0.5
    # dParams.fillweight  = 100
    # dParams.maxlength   = 200
    # dParams.minlength   = 10
    # dParams.maxdist     = 100
    # dParams.maxel       = 8
    # dParams.varfill     = 0.3                                                 # variation in filling weight (std around 1)
    # dParams.maxvoxel    = 6                                                   # Maximum number of elements per voxel
    # dParams.padsize     = 20                                                  # Background padding size (for smoothness in background)
    # dParams.numbranches = 20
    # dParams.varbranches = 5
    # dParams.maxfill     = 0.5
    # dParams.N_proc      = 10
    # dParams.l           = 25
    # dParams.rho         = 0.1
    dParams = Axon_params(1,0.5,100,200,10,100,8,0.3,6,20,20,5,0.5,10,25,0.1)
    axon_params = setParams(dParams,axon_params)
        
    
    return axon_params

###########################################################################
###########################################################################
