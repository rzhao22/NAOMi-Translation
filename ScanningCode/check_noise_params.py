import numpy as np
import inspect
def check_noise_params(noise_params):

# check_noise_params(noise_params)
# return noise_params 
#  
# This function checks the elements of the Class instance noise_params to ensure
# that all fields are set. Fields not set are set to a list of default
# parameters. The Class instance checked is:
# 
#   - noise_params - Class instance contaning the parameters for the noise model
#          .mu       - Mean measurement increase per photon (default = 100)
#          .mu0      - Electronics offset (default = 0)
#          .sigma    - Variance increase per photon (default = 2300)
#          .sigma0   - Electronics base noise variance (default = 2.7)
#          .sigscale - (default = 2.0000e-07)
#          .bleedp   - (default = 0.3000)
#          .bleedw   - (default = 0.4000)
#    
# 2017 - Adam Charles and Alex Song

###########################################################################
## Run the checks
    class Noise_params:
        pass
    if not inspect.isclass(noise_params):
        del noise_params
        
        noise_params = Noise_params()
    

    if (not hasattr(noise_params,'mu')) or not noise_params.mu:                 # Default mean measurement increase per photon is 100
        noise_params.mu = 100
    
    if (not hasattr(noise_params,'mu0')) or not noise_params.mu0:               # Default electronics offset is 0
        noise_params.mu0 = 0
    
    if (not hasattr(noise_params,'sigma')) or not noise_params.sigma:           # Default variance increase per photon is 2300
        noise_params.sigma = 2300
    
    if (not hasattr(noise_params,'sigma0')) or not noise_params.sigma0:         # Default electronics base noise variance is 2.7
        noise_params.sigma0 = 2.7
    
    if (not hasattr(noise_params,'darkcount')) or not noise_params.darkcount:   # Default PMT dark count rate (and excess photon rate from light leakage and autofluoresence)
        noise_params.darkcount = 0.05
    
    if (not hasattr(noise_params,'sigscale')) or not noise_params.sigscale:     # Default signal magnitude scale is 2e-7
        noise_params.sigscale = 2e-7
    
    if (not hasattr(noise_params,'bleedp')) or not noise_params.bleedp:         # Default electronics pixel bleed through probability is 0.3
        noise_params.bleedp = 0.3
    
    if (not hasattr(noise_params,'bleedw')) or not noise_params.bleedw:         # Default electronics average pixel bleed-through (if it occurs) is 0.4
        noise_params.bleedw = 0.4
    
    if (not hasattr(noise_params,'bleedp')) or not noise_params.bleedp:      
        noise_params.bleedp = 0.2
    
    if (not hasattr(noise_params,'bleedw')) or not noise_params.bleedw:      
        noise_params.bleedw = 0.5
    

    
    return noise_params

###########################################################################
###########################################################################

# testing
if __name__ == "__main__":
    print(check_noise_params(np.array([1])).bleedw)