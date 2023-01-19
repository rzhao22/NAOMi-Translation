# spike_opts = check_spike_opts(spike_opts) 
#  
# This function checks the elements of the struct bg_params to ensure
# that all fields are set. Fields not set are set to a list of default
# parameters. The struct checked is:
# 
#   spike_opts - Struct containing the main options used for generating
#                time traces. Includes
#     .nt         - Number of time-steps to simulate (default = 1000)
#     .K          - Number of neurons to generate time traces for (default =
#                   30)
#     .rate       - Spiking rate for each neuron (default = 0.001 t/spike)
#     .dt         - Time difference between time-points (inverse of
#                   frame-rate to generate samples at) (default = 1/30 s) 
#     .mu         - (default = 0)
#     .sig        - (default = 1)
#     .dyn_type   - Type of dynamics to simulate. Can be 'AR1'
#                   (autoregressive model with 1 time lag), 'AR2' (AR model
#                   with 2 time lags), 'single' (single-compartment Ca2+ 
#                   <-> protein dynamics), 'Ca_DE' (single-compartment Ca2+
#                    <-> double-expo protein dynamics) (default = 'Ca_DE')
#     .N_bg       - Number of background/neuropil time-traces to generate
#                   (default = 0)
#     .prot       - Protein selection for Ca2+ dynamics and fluorescence
#                   simulation. Options are 'GCaMP3' and 'GCaMP6' (default
#                   = 'GCaMP6')
#     .alpha      - Gamma distribution parameter that modulates the rate
#                   distribution. Defaults to an to an Exponential
#                   distribution (alpha=1) 
#     .burst_mean - Mean of the Poisson distribution for number of spikes
#                   per burst to be (default = 4)
# 
# 2017 - Adam Charles and Alex Song

###########################################################################
## Run the checks

import numpy as np

def check_spike_opts(spike_opts):
    
    if (!isfield(spike_opts,'K'))|isempty(spike_opts.K):
        spike_opts.K = 30;                                                     # Default to 30 neurons
    
    if (!isfield(spike_opts,'mu'))|isempty(spike_opts.mu)
        spike_opts.mu = 0;                                                     # Default the mean of the normal r.v. used in the log-normal to zero
    
    if (!isfield(spike_opts,'sig'))|isempty(spike_opts.sig)
        spike_opts.sig = 1;                                                    # Default the standard-deviation of the normal r.v. used in the log-normal to one
    
    if (!isfield(spike_opts,'dyn_type'))|isempty(spike_opts.dyn_type)         # Current options are 'AR1', 'AR2', 'single' or 'Ca_AR2'
        spike_opts.dyn_type = 'Ca_DE';                                         # Default to a single-compartment Ca2+ model with double-exponential smoothing
    
    if (!isfield(spike_opts,'rate_dist'))|isempty(spike_opts.rate_dist)       # Current options are 'gamma' or 'uniform'
        spike_opts.rate_dist = 'gamma';                                        # Default to a gamma distribution of firing rates
    
    if (!isfield(spike_opts,'dt'))|isempty(spike_opts.dt)         
        spike_opts.dt = 1/30;                                                  # Default sampling rate is 30 Hz
    
    if (!isfield(spike_opts,'nt'))|isempty(spike_opts.nt)         
        spike_opts.nt = 1000;                                                  # Default number of time-steps to simulate is 1000
    
    if (!isfield(spike_opts,'rate'))|isempty(spike_opts.rate)                 
        spike_opts.rate = 1e-3;                                                # Default inverse average of 1s between bursts (1000 1/100 s units) WAS: spike_opts.rate = 0.16;
    
    if (!isfield(spike_opts,'N_bg'))|isempty(spike_opts.N_bg)         
        spike_opts.N_bg = 0;                                                   # Default to only one background component
    
    if (!isfield(spike_opts,'prot'))|isempty(spike_opts.prot)                 # Current options are 'GCaMP6', 'GCaMP3'
        spike_opts.prot = 'GCaMP6';                                            # Default to 'GCaMP6' 
    
    if !isfield(spike_opts,'alpha')
        spike_opts.alpha = 1;                                                  # Default the Gamma distribution to an Exponential distribution (alpha=1)
    
    if !isfield(spike_opts,'burst_mean')
        spike_opts.burst_mean = 10;                                            # Default the mean of the Poisson distribution for number of spikes per burst to be 4
    
    if !isfield(spike_opts,'smod_flag')
        spike_opts.smod_flag = 'hawkes';                                       # Default simulation model is the Hawke's model
    
    if !isfield(spike_opts,'p_off')
        spike_opts.p_off = 0.2;                                           
    
    if !isfield(spike_opts,'selfact')
        spike_opts.selfact = 1.2;                                           
    
    if !isfield(spike_opts,'min_mod')
        spike_opts.min_mod = [0.4 2.53];                                       # shape parameter estimate for values estimated from Thy1 paper, visual cortex gp 5.17 WAS: spike_opts.min_mod = 0.3;  
    
    if !isfield(spike_opts,'spikeflag')|isempty(spike_opts.spikeflag)  
        spike_opts.spikeflag = 1;                                              # Flag to save spike data when simulating activity traces
    
    if !isfield(spike_opts,'dflag')|isempty(spike_opts.dflag)  
        spike_opts.dflag = 1;                                               # Flag to simulate drite traces per component
    
    if !isfield(spike_opts,'axonflag')|isempty(spike_opts.axonflag)  
        spike_opts.axonflag = 1;                                               # Flag to simulate axon traces per component
    
    
    return spike_opts

###########################################################################
###########################################################################