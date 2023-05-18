import numpy as np
import math
import scipy
from scipy.sparse.linalg import eigs
from check_neur_params import check_neur_params
from teardrop_poj import teardrop_poj

def generateNeuralBody(neur_params, eo, Vs):

# generateNeuralBody(neur_params)
# return Vcell,Vnuc,Tri,a
#
# This function samples a neural shape using an Isotropic Gaussian process
# on a sphere. The function basically draws a set of points whos
# covariance matrix depends on the geodescic distance along the sphere.
# Also provided is a sphere embedded in the 
# 
# Inputs to this function are
#   - neur_params - Struct containing parameters for neuron generation
#       .n_samps     - Number of sphere samples to use in creating the mesh
#                      for generating soma and nucleus shapes (default =
#                      200) 
#       .l_scale     - length-scale for the isotropic GP of the soma
#                      shapes. This controls the shape `bumpiness' (default
#                      = 90) 
#       .p_scale     - Overall variance of the isotropic GP of the soma 
#                      shape. (default = 90) 
#       .avg_rad     - Average radius of each neuron in um (default =
#                      6.8 um) 
#       .nuc_fluorsc - Potential fluorescence in the nucleus (default =
#                      0)
#       .min_thic    - Minimum cytoplasmic thickness (default = 0.8)
#       .eccen       - Maximum eccentricity of neuron (default = 0.35)
#       .exts        - Parameters dictating the max/min of the soma radii
#                      (Default = [0.75,1.7])
#       .nexts       - Parameters dictating the extent to shrink and smooth
#                      the nucleus (Default = [60,20])
#       .neur_type   - Option for neuron type (Default 'pyr')
# 
# The outputs of this function are
#   - Vcell:       Nx3 matrix indicating points on the neural surface
#   - Vnuc         Nx3 matrix indicating points on the nucleus surface
#   - Tri          Triangulation matrix indicating the faces of the neural
#                  and nucleus surfaces (same triangulation for both)
#   - a            Vector with the rotation angle of the cell (Rx,Ry,Rz)
# 
# 2016 - Adam Charles
#
##########################################################################
## Input Parsing

    if not neur_params:                                                    # Make sure the neur_params struct exists at all
        neur_params.TMP = []
    neur_params = check_neur_params(neur_params)                              # Check neuron parameters
    pwr       = 1                                                             # GP power: DON'T CHANGE - VERY SENSITIVE
    nucoff    = 3                                                             # Nucleus offset

    ###########################################################################
    ## Get sampling on a sphere

    if (((not hasattr(neur_params,'S_samp')) or not neur_params.S_samp) 
          and ((not hasattr(neur_params,'Tri')) or not neur_params.Tri)):              # Check if a sampling is already provided
        V,Tri=SpiralSampleSphere(neur_params.n_samps)                 # This function uses a spiral method to sample uniformly on a sphere
    else:
        V   = neur_params.S_samp                                              # If provided, use the provided sampling and triangulation
        Tri = neur_params.Tri

    ###########################################################################
    ## Calculate covariance based on distances

    if (neur_params.neur_type == 'pyr'):                                     # In the case of pyramidal neurons, add a tear-dropped mean to bias the GP
        Vtear = teardrop_poj(V,1)                                          # Change samples to samples on a teardrop
    elif (neur_params.neur_type == 'peanut'):                              # Take me out to the ball-game
        Vtear = teardrop_poj(V,2)                                          # Change samples to samples on a teardrop
    else:
        # Do nothing
        pass

    if ((not hasattr(neur_params,'dists')) or not neur_params.dists):           # If 'dists' is not provided, calculate arc-lengths between points
        dists = np.reshape(V, [np.size(V,0),1,np.size(V,1)]) - np.reshape(V, 
            [1,np.size(V,0),np.size(V,1)]) # Calculate element-wise distance between points
        dists = np.sqrt(sum(np.square(dists),3))                                         # Calculate full distance between points
        dists = 2*math.asin(dists/2)                                            # The geodesic distance is the arc-length
        dists = neur_params.p_scale*np.exp(-(dists/(neur_params.l_scale))**pwr) # The actual covariance elements are e^{-dist/l}
    else:
        dists = neur_params.dists                                             # If provided, use the provided distances
        dists = neur_params.p_scale*np.exp(-(dists/(neur_params.l_scale))**pwr) # The actual covariance elements are e^{-dist/l}

    if ((not hasattr(neur_params,'Rtear')) or not neur_params.Rtear):
        if (neur_params.neur_type == 'pyr'):                                 # In the case of pyramidal neurons, add a tear-dropped mean to bias the GP
            Rtear = 1*np.sqrt(sum(Vtear**2,2))                                   # Get radii of points on the teardrop
        elif (neur_params.neur_type == 'peanut'):                          # Take me out to the ball-game
            Vtear = Vs                                                        # Initialize points to the points on a sphere
            Rtear = 1*np.sqrt(sum(Vtear**2,2))                                   # Get radii of points on the teardrop
        else: 
            Rtear = 1
    else:
        Rtear = neur_params.Rtear                                             # If provided, use the provided means
    # eo.issym  = 1                                                             # Use symmetric flag for eigenvalue finding 
    eo.isreal = 1                                                             # Use real flag for eigenvalue finding
    min_eig   = 1.03*eigs(dists,1,'SR',eo) #### check here                                   # To combat ill conditioning from numerical errors, find the minimum eigenvalue
    if min_eig < 0:
        dists = dists + abs(min_eig)*np.eye(dists.shape)                         # Find the value that makes sure the covariance is PSD
    x_bnds = neur_params.exts*neur_params.avg_rad                             # Set bounds for the maximum and minimum radii of the soma shape
    x_base = abs(np.random.multivariate_normal(0*Rtear,dists).T)                                     # Sample from a Gaussian with the correct covariance matrix
    x      = x_base - np.mean(x_base) + neur_params.avg_rad*Rtear                # Re-center the samples 
    xmin   = min(min(x),x_bnds[0])                                            # Calculate the minimum radius of the current shape
    x      = (x_bnds[1] - x_bnds[0])*(x - xmin)/(
        max(max(x),x_bnds[1]) - xmin) + x_bnds[0]  # Renormalize the radii so that the cell isn't ever too big or too small in any direction
    if (neur_params.neur_type == 'pyr'): 
        x2   = x_base - np.mean(x_base) + neur_params.avg_rad                    # Re-center the samples (now for the nucleus)
        xmin = min(min(x2),x_bnds[0])                                         # Calculate the minimum radius of the current shape (now for the nucleus)
        x2   = (x_bnds[1] - x_bnds[0])*(x2 - xmin)/(
                max(max(x2),x_bnds[1]) - xmin) + x_bnds[0] # Renormalize the radii so that the nucleus isn't ever too big or too small in any direction (now for the nucleus)
    else:
        x2 = x

    ###########################################################################
    ## Create Elliptical Nucleus

    # eccens = ones(1,3)+neur_params.eccen                      # Get ellipse eccentricities
    eccens = np.ones((1,3))+neur_params.eccen * (np.random.rand(1,3)-np.array([0.5,0.5,0]))                      # Get ellipse eccentricities
    # eccens = ones(1,3)+neur_params.eccen*(rand(1,3)-0.5)                      # Get ellipse eccentricities
    # eccens = eccens./prod(eccens)
    eccens = eccens/(np.prod(eccens)**(1/3))

    if (neur_params.neur_type == 'pyr'):  
        Vetear = Vtear * eccens                                  # Create an elliptical teardrop
        Vetear = Vetear/np.sqrt(np.mean(sum(Vetear**2,2)))                          # Normalize the elliptical teardrop
    else:
        Vetear = V * eccens                                      # Create an elliptical shape
        Vetear = Vetear/np.sqrt(np.mean(sum(V**2,2)))                               # Normalize the elliptical shape

    Vcell   = Vetear* x.flatten()                              # Multiply through the sampled points to get the samples on the neural surface
    Vcell   = Vcell + np.array([0, 0, -nucoff])                           # Shift the soma down (for easier modulation of the nucleus)
    Vnorms  = np.sqrt(sum(Vcell**2,2))                                           # Get soma surface point extents
    # Vnuc    = bsxfun(@times, bsxfun(@times, V, [1,1,-1].*eccens), x2(:))      # Initialize the nucleus to an off-set cell shape
    # Vnorms2 = sqrt(sum(Vnuc.^2,2))                                            # Get soma surface point extents
    # Vnorms2 = neur_params.nexts(2)*(neur_params.nexts(1)*(Vnorms2 ...
    #                 - min(Vnorms2)) + (1-neur_params.nexts(1))*max(Vnorms2))  # Shrink & smooth the nucleus
    # Vnorms2 = Vnorms2 + min(Vnorms - Vnorms2) - neur_params.min_thic(1)          # Make sure the nucleus fits inside the soma with the minimum thickness
    # Vnuc    = bsxfun(@times,Vnuc,Vnorms2./ sqrt(sum(Vnuc.^2,2)))              # Apply shrinkage values to soma surface points

    Vnuc    = (V * np.array([1,1,-1])) * x2.flatten()      # Initialize the nucleus to an off-set cell shape
    Vnorms2 = np.sqrt(sum(Vnuc**2,2))                                            # Get soma surface point extents
    Vnorms2 = neur_params.nexts[1]*(neur_params.nexts[0]*(Vnorms2
                    - min(Vnorms2)) + (1-neur_params.nexts[0])*max(Vnorms2))  # Shrink & smooth the nucleus
    Vnorms2 = Vnorms2 + min(Vnorms - Vnorms2) - neur_params.min_thic[0]          # Make sure the nucleus fits inside the soma with the minimum thickness
    Vnuc    = (Vnuc * eccens) * Vnorms2 / np.sqrt(sum(Vnuc**2,1))              # Apply shrinkage values to soma surface points

    lat_ang = np.random.rand(1)*2*np.pi
    lat_shft = (1-abs(np.random.rand(1)-np.random.rand(1)))*neur_params.min_thic[1]*[math.sin(lat_ang), math.cos(lat_ang)]
    # lat_shft = rand(1)*neur_params.min_thic(2)*[sin(lat_ang), cos(lat_ang)]

    Vcell   = Vcell + np.array([0, 0, nucoff])                            # Shift the soma up
    Vnuc    = Vnuc + np.array([lat_shft[0], lat_shft[1], nucoff])                             # Shift the nucleus up

    # Vcell   = bsxfun(@plus, Vcell, [0, 0, nucoff])                            # Shift the soma up
    # Vnuc    = bsxfun(@plus, Vnuc, [0, 0, nucoff])                             # Shift the nucleus up

    # Optional scaling of nucleus size to adjust to normalized radius
    if (hasattr(neur_params,'nuc_rad') and not (neur_params.nuc_rad)):
        hull = scipy.spatial.ConvexHull(Vnuc[:,0],Vnuc[:,1],Vnuc[:,2])
        VnucSz = hull.volume
        nucsz = (4/3)*np.pi*(neur_params.nuc_rad[0]**3)
        if (max(neur_params.nuc_rad.shape)>1):
            Vnuc = Vnuc*(((nucsz/VnucSz)**(1/3))**(1/neur_params.nuc_rad[1]))
        else:
            Vnuc = Vnuc*((nucsz/VnucSz)**(1/3))


    if ((not hasattr(neur_params,'max_ang')) or not (neur_params.max_ang)):
        max_ang = 20
    else:
        max_ang = neur_params.max_ang
    a     = -abs(max_ang)+2*abs(max_ang)*np.random.rand(1,3)                            # Choose a random rotation angle
    Rx    = np.array([[1, 0, 0], [0, cosd(a[0]), -sind(a[0])], [0, sind(a[0]), cosd(a[0])]])        # ... Rotation around x
    Ry    = np.array([[cosd(a[1]), 0, sind(a[1])], [0, 1, 0], [-sind(a[1]), 0, cosd(a[1])]])        # ... Rotation around y
    Rz    = np.array([[cosd(a[2]), -sind(a[2]), 0], [sind(a[2]), cosd(a[2]), 0], [0, 0, 1]])        # ... Rotation around z
    Vnuc  = Vnuc*Rx*Ry*Rz                                                     # Apply the three rotation matrices to the nucleus shape 
    Vcell = Vcell*Rx*Ry*Rz                                                    # Apply the three rotation matrices to the soma shape
                                
    return Vcell, Vnuc, Tri, a

    ##########################################################################
    ##########################################################################

def sind(x):
    I = x/180.
    y = np.sin(I * np.pi)
    mask = (I == np.trunc(I)) & np.isfinite(I)
    y[mask] = 0
    return y

def cosd(x):
    I = x/180.
    y = np.cos(I * np.pi)
    mask = (I == np.trunc(I)) & np.isfinite(I)
    y[mask] = 0
    return y

def SpiralSampleSphere(n, flag=True):
    inc = np.pi * (3 - np.sqrt(5))
    off = 2 / float(n)
    phi = 0
    theta = 0
    points = np.zeros((n, 3))
    for i in range(n):
        y = i * off - 1 + (off / 2)
        r = np.sqrt(1 - y*y)
        phi = phi + inc
        if phi > 2*np.pi:
            phi = phi - (2*np.pi)
            theta = theta + inc
        x = np.cos(phi) * r
        z = np.sin(phi) * r
        points[i] = [x, y, z]
    tri = np.zeros(((n-2)*3, 3), dtype=int)
    k = 0
    for i in range(1, n-1):
        tri[k] = [i, i+1, 0]
        k = k + 1
    tri[k] = [n-1, 1, 0]
    k = k + 1
    for i in range(1, n-2):
        tri[k] = [i, i+1, i+2]
        k = k + 1
    if flag:
        return points, tri
    else:
        return points