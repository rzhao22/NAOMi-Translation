
def generateZernike(psf_params,offset=[0,0]):

# lambda called Lambda here to avoid built in python lambda

# abb = generateZernike(psf_params,offset)
#  
# Generates the Zernike aberration weights given a input distribution and
# offset. The inputs are
# 
# - psf_params        - Struct contaning the parameters for the PSF
#        .Lambda      - Two-photon excitation wavelength (um)
#        .zernikeWt   - Microscope aberration weights (Zernike basis)
#        .zernikeDst  - Microscope aberration weights (Zernike basis) as
#                       a function of position
# - offset            - offset for calculating aberrations at a distance
#                       away from the center of field
#
# The output is
# - abb               - Zernike abberation weights
# 
# 2017 - Alex Song

###########################################################################

    if((not hasattr(psf_params,'zernikeWt')) or (len(psf_params.zernikeWt)==0)):
        abb = 0
    else:
        if((not hasattr(psf_params,'zernikeDst')) or (len(psf_params.zernikeDst)==0)):
            abb = psf_params.zernikeWt
        else:
            abb = psf_params.zernikeDst[offset]*psf_params.zernikeWt
            
    abb = abb*psf_params.Lambda*1e-6

    return abb
