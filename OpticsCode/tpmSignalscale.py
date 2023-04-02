import numpy as np

def tpmSignalscale(tpm_params, psf_params):
##function Ftavg = tpmSignalscale(tpm_params,psf_params)

# function Ftavg = tpmSignalscale(tpm_params,psf_params)
# 
# Function to estimate the 
# <F(t)> = 1/2*phi*eta*C*delta*gp/f/tau*8*nidx*<P(t)>^2/pi/lambda
# gives the average photons collected at a particualr concentration.
# Adapted from Chris Xu, Watt Webb (1996) JOSA B
#
#   - tpm_params - Struct containing parameters for estimating photon flux
#     .nidx = 1.33      - refractive index [], water
#     .nac = 0.8        - objective NA []
#     .phi = 0.8*SA*0.4 - 80 percent transmission and 0.8NA solid angle
#                         with average 40 perfect PMT QE []
#     .eta = 0.6        - eGFP quantum yield, estimate for 2P QY []
#     .conc = 10        - fluorophore concentration, average of literature 
#                         (Huber et al 2012, Zariwala et al 2012), [uM]
#     .delta = 35       - two-photon abs. cross section, estimate from
#                          Harris lab (Janelia), saturated GCaMP [GM]
#     .delta = 2        - estimate of 2pCS at resting calcium levels [GM]
#     .gp = 0.588       - pulse-shape temporal coherence [], sech 
#     .f = 80           - Ti:S laser rep rate [MHz]
#     .tau = 150        - Ti:S laser pulse-width [fs]
#     .pavg = 15        - laser power [mW]
#     .lambda = 0.92    - wavelength used for GCaMP excitation
#   - psf_params - Struct containing PSF parameters for calculating flux
#
#   - Ftavg - average number of photons collected per second
# 
# 2017 - Alex Song and Adam Charles
#
###########################################################################

    if(psf_params):
        nidx   = float(psf_params.n)
        nac    = float(psf_params.objNA)
        Lambda = float(psf_params.Lambda)
    else:
        nidx   = float(tpm_params.nidx)                                        # refractive index []
        nac    = float(tpm_params.nac)                                         # collection NA []
        Lambda = float(tpm_params.Lambda)                                      # wavelength [um]
    
    phi   = float(tpm_params.phi)                                            # collection efficiency []
    eta   = float(tpm_params.eta)                                            # fluorophore QE []
    conc  = float(tpm_params.conc)                                           # fluorophore concentration [uM]
    delta = float(tpm_params.delta)                                          # two-photon cross-section [GM]
    gp    = float(tpm_params.gp)                                             # pulse-shape temporal coherence []
    f     = float(tpm_params.f)                                              # laser repetition rate [MHz]
    tau   = float(tpm_params.tau)                                            # laser pulse-width [fs]
    pavg  = float(tpm_params.pavg)                                           # average power [mW]

    conc   = conc*1e-6*6.02e23*1e3                                            # molar concentration to number of molecules
    delta  = delta*1e-58
    f      = f*1e6
    tau    = tau*1e-15
    Lambda = Lambda*1e-6
    pavg   = 1e-3*pavg/(6.626e-34*3e8/Lambda) # times hc/Lambda

    Ftavg  = phi*eta*conc*delta*gp*8*nidx*(pavg**2)/(2*f*tau*np.pi*Lambda) # about 1 photon per pulse (close to 80M photons/s) - or average 10 photons/pixel
    return Ftavg