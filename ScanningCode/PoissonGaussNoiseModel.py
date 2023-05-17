import numpy as np
from check_noise_params import check_noise_params

def PoissonGaussNoiseModel(clean_in, noise_params):
# PoissonGaussNoiseModel(clean_in, noise_params)
# return noisy_out
#
# Draw noisy outputs given clean fluorescence rates. The noise model first
# draws random i.i.d. Poisson photon counts, x, and then simulates the
# electronic measurement noise. The inputs to this function are
#   - clean_in     - Raw fluorescence value at each location in the
#                    field-of-view
#   - noise_params - Struct contaning the parameters for the noise model
#          .mu     - Mean measurement increase per photon (default = 150)
#          .mu0    - Electronics offset (default = 0)
#          .sigma  - Variance increase per photon (default = 5000)
#          .sigma0 - Electronics base noise variance (default = 3.5)
#
# The output is
#   - noisy_out - the currupted version of clean_in
#
# 2017 - Adam Charles and Alex Song

###########################################################################
    ## Input Parsing

    noise_params = check_noise_params(noise_params)                           # Check the noise parameter struct for missing elements

    ###########################################################################
    ## Draw noisy outputs

    cnt_tmp = np.random.poisson(clean_in+noise_params.darkcount)                       # Draw Poisson counts from the distribution plus the dark count rate

    m      = cnt_tmp*noise_params.mu                                          # Mean count is scaled Poisson count
    v      = cnt_tmp*noise_params.sigma                                       # Standard deviation is scaled Poisson count
    mu2    = np.log((m**2)/np.sqrt(v+m**2))
    sigma2 = np.sqrt(np.log(v/(m**2)+1))

    noisy_out = np.random.lognormal(mu2,sigma2)                                           # Sample counts from a log-normal random variable with Poisson-based mean and variance
    noisy_out[np.isnan(noisy_out)] = 0                                           # Get rid of pesky NaN values (those shoule be zero)
    noisy_out = round(noisy_out+np.random.normal(noise_params.mu0,noise_params.sigma0,
    np.size(noisy_out,0),np.size(noisy_out,1)))                                   # Add Gaussian noise that stems from the electronics

    return noisy_out
###########################################################################
###########################################################################
