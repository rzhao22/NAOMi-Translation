import numpy as np
from setParams import setParams
def check_vasc_params(vasc_params):

##function vasc_params = check_vasc_params(vasc_params)

# vasc_params = check_vasc_params(vasc_params)
#  
# This function checks the elements of the struct vasc_params to ensure
# that all fields are set. Fields not set are set to a list of default
# parameters. The struct checked is:
# 
#   - vasc_params - Struct containing parameters for vasculature simulation
#       .flag            - On/off flag for  vasculature simulation
#       .ves_shift       - 3-vector of the amount of wobble allowed for
#                          blood vessels (default = [5 15 5] um)
#       .depth_vasc      - Depth into tissue for which vasculature is
#                          simulated (default = 200um) 
#       .depth_surf      - Depth into tissue of surface vasculature
#                          (default = 15um) 
#       .distWeightScale - Scaling factor for weight of node distance
#                          scaling (default = 2) 
#       .randWeightScale - scaling factor for weight of nodes (additional
#                          variability) (default = 0.1) 
#       .cappAmpScale    - scaling factor for weights (capillaries,
#                          lateral) (default = 0.5) 
#       .cappAmpZscale   - scaling factor for weights (capillaries, axial)
#                          (default = 0.5) 
#       .vesSize         - vessel radius (surface, axial, capillaries) in
#                          um (default = [15 10 3] um) 
#       .vesFreq         - blood vessel freqency in microns (default = 
#                          [125 200 50] um) 
#       .vesNumScale     - blood vessel number random scaling factor
#                          (default = 0.2) 
#       .sourceFreq      - Rate of generation of source nodes (default =
#                          1000 um/node) 
#       .sepweight       - Set weight that guides how far, on average, the
#                          vasculature nodes are placed (value is between 0
#                          and 1 default = 0.75)
#       .distsc          - How strongly local capillary connections are.
#                          Higher numbers indicate more local connections
#                          (default = 4)
#       .node_params     - Set of parameters to guide node placement:
#          .maxit        - Maximum iteration to place nodes (default = 25)
#          .lensc        - Average distance between vasculature branch
#                          points (default = 50 um) 
#          .varsc        - Standard deviation of distances between
#                          vascualure branch points (default = 15 um) 
#          .mindist      - Minimum inter-node distance (default = 10 um)
#          .varpos       - Standard deviation of vasculature placement
#                          (default = 5 um) 
#          .dirvar       - The maximum branching angle (default = pi/8)
#          .branchp      - Probability of branching surface vasculature
#                          (default = 0.02) 
#          .vesrad       - Radius of surface vasculature (default = 25 um)
# 
# 2017 - Adam Charles and Alex Song
#
###########################################################################
## Run the checks

    if not vasc_params:                                                    # Make sure that vasc_params is a struct
        del vasc_params
        class Vasc_params:
            def __init__(self,flag,ves_shift,depth_vasc,depth_surf,distWeightScale,randWeightScale
                ,cappAmpScale,cappAmpZscale,vesSize,vesFreq,sourceFreq,vesNumScale,sepweight,distsc):
                self.flag = flag                                                          # amount of wobble allowed for blood vessels
                self.ves_shift       = ves_shift                                        # amount of wobble allowed for blood vessels
                self.depth_vasc      = depth_vasc                                             # depth into tissue for which vasculature is simulated
                self.depth_surf      = depth_surf                                              # depth into tissue of surface vasculature
                self.distWeightScale = distWeightScale                                               # scaling factor for weight of node distance scaling
                self.randWeightScale = randWeightScale                                             # scaling factor for weight of nodes (additional variability)
                self.cappAmpScale    = cappAmpScale                                             # scaling factor for weights (capillaries, lateral)
                self.cappAmpZscale   = cappAmpZscale                                             # scaling factor for weights (capillaries, axial)
                self.vesSize         = vesSize                                        # vessel radius (surface, axial, capillaries) in um
                self.vesFreq         = vesFreq                                    # vessel freqency in microns
                self.sourceFreq      = sourceFreq                                            # 
                self.vesNumScale     = vesNumScale                                             # vessel number random scaling factor
                self.sepweight       = sepweight                                            #
                self.distsc          = distsc
        vasc_params = Vasc_params(1,np.array([5,15,5]),200,15,2,0.1,0.5,0.5,np.array([15,9,2]),np.array([125,200,50]),1000,0.2,0.75,4)
                 

    # dParams.flag            = 1                                               # amount of wobble allowed for blood vessels
    # dParams.ves_shift       = [5, 15, 5]                                        # amount of wobble allowed for blood vessels
    # dParams.depth_vasc      = 200                                             # depth into tissue for which vasculature is simulated
    # dParams.depth_surf      = 15                                              # depth into tissue of surface vasculature
    # dParams.distWeightScale = 2                                               # scaling factor for weight of node distance scaling
    # dParams.randWeightScale = 0.1                                             # scaling factor for weight of nodes (additional variability)
    # dParams.cappAmpScale    = 0.5                                             # scaling factor for weights (capillaries, lateral)
    # dParams.cappAmpZscale   = 0.5                                             # scaling factor for weights (capillaries, axial)
    # dParams.vesSize         = [15, 9, 2]                                        # vessel radius (surface, axial, capillaries) in um
    # dParams.vesFreq         = [125, 200, 50]                                    # vessel freqency in microns
    # dParams.sourceFreq      = 1000                                            # 
    # dParams.vesNumScale     = 0.2                                             # vessel number random scaling factor
    # dParams.sepweight       = 0.75                                            #
    # dParams.distsc          = 4                                               #
    dParams = Vasc_params(1,np.array([5,15,5]),200,15,2,0.1,0.5,0.5,np.array([15,9,2]),np.array([125,200,50]),1000,0.2,0.75,4)
    vasc_params = setParams(dParams,vasc_params)

    if (not hasattr(vasc_params,'node_params')) or not vasc_params.node_params: # Check if the node parameter sub-struct exists. If not - make it
        class Node_params:
            def __init__(self,maxit,lensc,varsc,mindist,varpos,dirvar,branchp,vesrad):
                self.maxit   = maxit
                self.lensc   = lensc
                self.varsc   = varsc
                self.mindist = mindist
                self.varpos  = varpos
                self.dirvar  = dirvar
                self.branchp = branchp
                self.vesrad  = vesrad
        vasc_params.node_params = Node_params(25,50,15,10,5,np.pi/8,0.02,25)
        # vasc_params.node_params.maxit   = 25                                                # Maximum iteration to place nodes
        # vasc_params.node_params.lensc   = 50                                                # 
        # vasc_params.node_params.varsc   = 15                                                # 
        # vasc_params.node_params.mindist = 10                                                # Minimum inter-node distance
        # vasc_params.node_params.varpos  = 5                                                 # 
        # vasc_params.node_params.dirvar  = np.pi/8                                              # 
        # vasc_params.node_params.branchp = 0.02                                              # 
        # vasc_params.node_params.vesrad  = 25                                                # 
    else:                                                                       # Check if the node parameter sub-struct exists
        if (not hasattr(vasc_params.node_params,'maxit')) or not vasc_params.node_params.maxit:
            vasc_params.node_params.maxit = 25                                        
        
        if (not hasattr(vasc_params.node_params,'lensc')) or not vasc_params.node_params.lensc:
            vasc_params.node_params.lensc = 50                                        
        
        if (not hasattr(vasc_params.node_params,'varsc')) or not vasc_params.node_params.varsc:
            vasc_params.node_params.varsc = 15                                        
        
        if (not hasattr(vasc_params.node_params,'mindist')) or not vasc_params.node_params.mindist:
            vasc_params.node_params.mindist = 10                                        
        
        if (not hasattr(vasc_params.node_params,'varpos')) or not vasc_params.node_params.varpos:
            vasc_params.node_params.varpos = 5                                        
        
        if (not hasattr(vasc_params.node_params,'dirvar')) or not vasc_params.node_params.dirvar:
            vasc_params.node_params.dirvar = np.pi/8                                        
        
        if (not hasattr(vasc_params.node_params,'branchp')) or not vasc_params.node_params.branchp:
            vasc_params.node_params.branchp = 0.02                                        
        
        if (not hasattr(vasc_params.node_params,'vesrad')) or not vasc_params.node_params.vesrad:
            vasc_params.node_params.vesrad = 25                                        
    
    return vasc_params

###########################################################################
###########################################################################
