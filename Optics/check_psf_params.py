import numpy as np
import sys
import inspect

# adding misc folder to the system path
sys.path.insert(0, './misc')
from setParams import setParams
def check_psf_params(psf_params):

# check_psf_params(psf_params)
# return psf_params
#  
# This function checks the elements of the Class instance psf_params to ensure
# that all fields are set. Fields not set are set to a list of default
# parameters. The Class instance checked is:
# 
#   - psf_params        - Class instance contaning the parameters for the PSF
#          .NA          - Numerical aperture of Gaussian beam
#          .n           - Refractive index of propagation material
#          .n_diff      - Shift in refractive index from vessels to tissue
#          .lambda      - Two-photon excitation wavelength (um)
#          .obj_fl      - Objective focal length (mm)
#          .ss          - Subsampling factor for fresnel propagation
#          .sampling    - Spatial sampling for tissue occlusion mask
#          .psf_sz      - Default two-photon PSF size simulated (um)
#          .prop_sz     - Fresnel propagation length outside of volume (um)
#          .blur        - PSF lateral blurring (um)
#          .scatter_sz  - Scattering object sizes (um), column vector
#          .scatter_wt  - Scattering object weights, column vector
#          .zernikeWt   - Microscope aberration weights (Zernike basis)
#
# 2017 - Adam Charles and Alex Song

###########################################################################
## Run the checks
    class Psf_params:
        def __init__(self,NA,objNA,n,n_diff,Lambda,obj_fl,ss,sampling,
                    psf_sz,prop_sz,blur,scatter_sz,scatter_wt,zernikeWt,
                    taillength,type,scaling,hemoabs,propcrop,fastmask):
            self.NA         = NA                                                  # Default excitation numerical aperture
            self.objNA      = objNA                                                  # Default objective numerical aperture
            self.n          = n                                                 # Default index of refraction in tissue
            self.n_diff     = n_diff                                                 # Default shift in index of refraction from vessels to tissue
            self.Lambda     = Lambda                                                 # Default two-photon excitation wavelength (microns)
            self.obj_fl     = obj_fl                                                  # Default objective focal length (mm)
            self.ss         = ss                                                    # Default subsampling factor for fresnel propagation (from volume voxel size)
            self.sampling   = sampling                                                   # Default spatial sampling for tissue occlusion mask
            self.psf_sz     = psf_sz                                           # Default two-photon PSF size simulated (microns)
            self.prop_sz    = prop_sz                                                   # Default fresnel propagation length outside of volume (microns)
            self.blur       = blur                                                    # Default psf lateral blurring (microns)
            self.scatter_sz = scatter_sz                              # Default scattering object sizes (microns), column vector
            self.scatter_wt = scatter_wt                               # Default scattering object weights, column vector
            self.zernikeWt  = zernikeWt                         # Default microscope aberration weights (Zernike basis). In units of wavelength, a small amount of spherical aberration and astigmatism added as uncorrected "system" aberrations
            self.taillength = taillength                                                   # Distance from edge of PSF_sz to estimate tailweight (um)    
            self.type       = type                                           # Default PSF type ('gaussian','vtwins','bessel')          
            self.scaling    = scaling                                         # Default PSF scaling type ('two-photon','three-photon','temporal-focusing')
            self.hemoabs    = hemoabs                                      # Hemoglobin absorbance scaling factor. Default assumes 150mg/ml Hb, 64500 g/mol Hb, 2.9 (abs/um)/(mol/L) in units of abs/um. Absorbance calculated from Scott Prahl's Hb curve and eGFP emission spectrum
            self.propcrop   = propcrop                                                # Flag to crop scanned beam during optical propagation (default true) 
            self.fastmask   = fastmask
    if not inspect.isclass(psf_params):
        del psf_params
        
        psf_params = Psf_params(0.6,0.8,1.35,0.02,0.92,4.5,2,50,[20,20,50],10,3,
                                np.matrix([0.51, 1.56, 4.52, 14.78]).H,np.matrix([0.57, 0.29, 0.19, 0.15]).H,
                                [0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 0.12],50,'gaussian',
                                'two-photon',0.00674*np.log(10),True,True)
    dParams = Psf_params(0.6,0.8,1.35,0.02,0.92,4.5,2,50,[20,20,50],10,3,
                                np.matrix([0.51, 1.56, 4.52, 14.78]).H,np.matrix([0.57, 0.29, 0.19, 0.15]).H,
                                [0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 0.12],50,'gaussian',
                                'two-photon',0.00674*np.log(10),True,True)
    psf_params = setParams(dParams, psf_params)

    class FM:
        pass

    if psf_params.fastmask:
        psf_params.FM = FM()
        psf_params.FM.sampling = 10
        psf_params.FM.fineSamp = 2
        psf_params.FM.ss       = 1
    

    
    return psf_params

###########################################################################
###########################################################################
# testing
if __name__ == "__main__":
    print(check_psf_params(np.array([1])).FM.ss)
