import numpy as np
import math
from gaussianBeamSize import gaussianBeamSize
from generateGaussianProfile import generateGaussianProfile
from generateZernike import generateZernike
from applyZernike import applyZernike
def generateBA(vol_params, psf_params):

##function Uout2 = generateBA(vol_params,psf_params)

# Uout2 = generateBA(vol_params,psf_params)
#
# This function generates a Gaussian back aperture intensity profile with 
# the specified aberrations applied. The inputs are
# 
# - vol_params          - Struct with parameters for the volume generation
#          .vol_sz      - 3-element vector with the size (in um) of the 
#                         volume to generate (default = 100x100x30um)
#          .vres        - resolution to simulate volume at (default = 2
#                         samples/um)
# - psf_params          - Struct contaning the parameters for the PSF
#          .NA          - Numerical aperture of Gaussian beam
#          .objNA       - Numerical aperture of objective lens
#          .n           - Refractive index of propagation material
#          .lambda      - Two-photon excitation wavelength (um)
#          .obj_fl      - Objective focal length (mm)
#          .ss          - Subsampling factor for fresnel propagation
#          .sampling    - Spatial sampling for tissue occlusion mask
#          .zernikeDst  - Microscope aberration weights (Zernike basis) as
#                         a function of position
# The output is
# - Uout2               - Output scalar field
#
# 2017 - Alex Song

###########################################################################

    if ((not hasattr(vol_params,'vasc_sz')) or not vol_params.vasc_sz):
        vol_params.vasc_sz = gaussianBeamSize(psf_params,vol_params.vol_depth+
        vol_params.vol_sz[2]/2)+vol_params.vol_sz+[0, 0, 1]*vol_params.vol_depth# Get beam size
    

    fl     = (psf_params.obj_fl/1000).astype(np.float32)                                   # focal length [m]
    D2     = (1e-6*(1/vol_params.vres).astype(np.float32) /psf_params.ss)                   # observation grid spacing [m]
    N      = (1e-6*(vol_params.vasc_sz[:1]-vol_params.vol_sz[:1])/D2).astype(np.float32)  # gridsize
    D1     = (max(gaussianBeamSize(psf_params,fl*1e6)/1e6)/min(N)).astype(np.float32)       # source grid spacing [m]
    nre    = (psf_params.n).astype(np.float32)                                              # immersion numerical index
    rad    = (math.tan(math.asin(psf_params.NA/nre))*fl).astype(np.float32)                           # source radius [m]
    objrad = (math.tan(math.asin(psf_params.objNA/nre))*fl).astype(np.float32)                        # source radius [m]
    k      = 2*nre*np.pi/(psf_params.Lambda*1e-6).astype(np.float32)                           # optical wavenumber [rad/m]
    X,Y  = np.meshgrid(np.arange(-N[0]/2,N[0]/2)*D1,np.arange(-N[1]/2,N[1]/2)*D1)             
    Uout   = generateGaussianProfile(X,Y,rad,objrad,k,fl)                     # generate gaussian wavefront with apodization 1

    imax = round(vol_params.vol_sz[0]/psf_params.sampling)+1
    jmax = round(vol_params.vol_sz[1]/psf_params.sampling)+1

    if (imax*jmax==1 or (not hasattr(psf_params,'zernikeDst')) or (not psf_params.zernikeDst)):
        abb = generateZernike(psf_params)
        Uout2 = applyZernike(Uout,X/objrad,Y/objrad,k,abb)
    else:
        Uout2 = [[None for j in range(jmax)] for i in range(imax)]
        for i in range(imax):
            for j in range(jmax):
                abb = generateZernike(psf_params)
                Uout2[i][j] = applyZernike(Uout,X,Y,k,abb)
    
    return Uout2