import time

def simulatebloodvessels(vol_params,vasc_params): 

# [neur_ves,vasc_params] = simulatebloodvessels(vol_params,vasc_params)
#
# Function to simulate vasculature in and above a neural volume. The inputs
# to this function are:
#   - vol_params      - Struct with parameters for the volume generation
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
#   - vasc_params - Struct containing parameters for vasculature simulation
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
#                          um (default = [15 6 2] um)
#       .vesFreq         - blood vessel freqency in microns (default = 
#                          [125 200 50] um) 
#       .vesNumScale     - blood vessel number random scaling factor
#                          (default = 0.2) 
#       .sourceFreq      - Rate of generation of source nodes (default =
#                          1000 um/node) 
#       .sepweight       - Set weight that guides how far, on average, the
#                          vasculature nodes are placed (value is between 0
#                          and 1; default = 0.75)
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
# The outputs of this function are
#   - neur_ves     - Array of the blood vessel locations
#   - vasc_params  - Updated scruct of parameters (with defaults filled in
#                    if needed)
#
# 2017 - Alex Song and Adam Charles
#
###########################################################################
## Input parsing

    vol_params  = check_vol_params(vol_params)                                # Check volume parameters
    vasc_params = check_vasc_params(vasc_params)                              # Check vasculature parameters

###########################################################################
## Simulating Blood Vessels

    if vol_params.verbose == 1:
      print('Generating in-volume blood vessels...')
    elif vol_params.verbose >1:
      print('Generating in-volume blood vessels...\n')
      tic = time.time()

    # Setup parameters for simulation
    vres = vol_params.vres                                              # Pull out volume resolution

    vp = vasc_params
    vp.depth_surf = vasc_params.depth_surf*vres                                 # Calculate depth beneath the surface
    vp.mindists = vasc_params.vesFreq*vres/2                                  # Set minimum distance between nodes
    vp.maxcappdist = 2*vasc_params.vesFreq[2]*vres                               # Maximum capilliary distance
    vp.vesSize = vasc_params.vesSize*vres                                        # Calculate the blood-vessel size

    node_params = vasc_params.node_params
    nop = node_params
    nop.lensc   = node_params.lensc*vres
    nop.varsc   = node_params.varsc*vres
    nop.mindist = node_params.mindist*vres
    nop.varpos  = node_params.varpos*vres
    nop.vesrad  = node_params.vesrad*vres

    if (not hasattr(vol_params,'vasc_sz')) or not bool(vol_params.vasc_sz)  #  or nargout<3
      nv.vol_sz  = vol_params.vol_sz+np.array([0 0 1])*vol_params.vol_depth             # Extract volume size
    else
      nv.vol_sz  = vol_params.vasc_sz

    nv.size = nv.vol_sz*vres                                                   #
    nv.szum = nv.vol_sz                                                   #
    nv.nsource = np.maximum(round((2*(nv.vol_sz[0]+nv.vol_sz[1])/(vp.sourceFreq))*...
                                     abs(1+vp.vesNumScale*np.random.randn(1))),0) # 
    nv.nvert = np.maximum(round((nv.vol_sz[0]*nv.vol_sz[1]/(vp.vesFreq[1]^2))*...
                                     abs(1+vp.vesNumScale*np.random.randn(1))),0) #
    nv.nsurf = np.maximum(round((nv.vol_sz[0]*nv.vol_sz[1]/(vp.vesFreq[0]^2))*...
                                     abs(1+vp.vesNumScale*np.random.randn(1))),0) # 
    nv.ncapp = np.maximum(round((prod(nv.vol_sz)/(vp.vesFreq[2]^3))*...
                                     abs(1+vp.vesNumScale*np.random.randn(1))),0) #
                                   
    ## Initialize a few points for vertical vessels, Initialize some points in surface for surface vessels
    [nodes,nv] = growMajorVessels(nv,nop,vp)

    ## Convert node structure to a connection structure
    conn = nodesToConn(nodes)
    nv.nconn = len(conn)

    # Shift surface vessel location to adjust for vessel diameter
    for i in range(len(conn)):
        if(sum([nodes[conn[i].start].type in {'edge', 'surf', 'sfvt'} for i in range(len(conn))])):
            nodes[conn[i].start].pos[2] = np.minimum(nodes[conn[i].start].pos[2]+ ...
            np.ceil(conn[i].weight/len(nodes[conn[i].start].conn)),nv.size[2])
        if(sum([nodes[conn[i].ends].type in {'edge', 'surf', 'sfvt'} for i in range(len(conn))])):
            nodes(conn[i].ends).pos[2] = np.minimum(nodes[conn[i].ends].pos[2]+ ...
            np.ceil(conn[i].weight/len(nodes[conn[i].start].conn)),nv.size[2])

    ## Create initial volume with major vessels
    [neur_ves,conn] = connToVol(nodes,conn,nv)

    ## Initialize and connect capillaries
    [nodes,conn,nv] = growCapillaries(nodes,conn,neur_ves,nv,vp,vres)

    ## Add capillaries to rest of volume
    cappidxs = np.where([len(c.locs) == 0 for c in conn])[0]
    neur_ves = connToVol(nodes,conn,nv,cappidxs,neur_ves)[0]

    if hasattr(vol_params,'vasc_sz')and bool(vol_params.vasc_sz):  # &&nargout==3
        neur_ves_all = neur_ves
        sz = [vol_params.vol_sz[0:1], vol_params.vol_depth+vol_params.vol_sz[2]]*vres
        sz_diff  = np.ceil((vol_params.vasc_sz*vres-sz)/2)
        neur_ves = neur_ves[(0:sz[2])+sz_diff[2],(0:sz[0])+sz_diff[0],(0:sz[1])+sz_diff[1]];
    '''
    elif (nargout==3)
      neur_ves_all = neur_ves;
    '''

    if vol_params.verbose == 1:
        print('done.\n');
    elif vol_params.verbose >1:
        toc = time.time()
        Tdone = toc-tic
        print('done (#f seconds).\n',Tdone)

    return [neur_ves,vasc_params,neur_ves_all]

###########################################################################
###########################################################################
