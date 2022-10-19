import scipy.io

def reSampleMovie(fullPath, varargin):

    if nargin > 1:
        boostLvl = varargin{1}
    else:
        boostLvl = []


    print('Loading movie parts...')
    result = fileparts(fullPath)
    save_path = result[0]
    save_filename = result[1]
    m = scipy.io.loadmat(fullPath)                                     # Access the mat file with writable permission

    if isempty(boostLvl):
        print('No boost level provided. Auto-normalizing to ~1 photon/pixel.')
        vol_params = m.vol_params
        boostLvl   = prod(vol_params.vol_sz(1:2))/(250)**2  
        print('done.\n')

    if not boostLvl == 1:
        noise_params = m.noise_params                                             # Load the parameter struct that needs to be adjusted
        Fsim_clean   = m.Fsim_clean                                               # Load the clean video

    print('done.\n')


    if not boostLvl == 1:
        print('Boosting clean movie and updating signal scale parameter...')
        noise_params.sigscale = boostLvl*noise_params.sigscale                    # Boost the signal scale parameter
        Fsim_clean            = boostLvl*Fsim_clean                               # Boost the actual signal
        print('done.\n')
        print('Resampling from the noise model with new signal scale...')
        Fsim                  = applyNoiseModel(Fsim_clean, noise_params)         # Rerun the noise model with the boosted singal
        print('done.\n')

        print('Saving new parameters and boosted movies...')
        m.Fsim_clean          = Fsim_clean                                        # Write back to the mat file the adjusted clean video 
        m.Fsim                = Fsim                                              # Write back to the mat file the adjusted sampled video 
        m.noise_params        = noise_params                                      # Write back the parameter struct that needs to be adjusted
        print('done.\n')


    tifFiles = strcat(save_path,'/',allFilesOfType(save_path, '.tif'))        # Find all the TIF files to delete for replacement
    aviFiles = strcat(save_path,'/',allFilesOfType(save_path, '.avi'))        # Find all the AVI files to delete for replacement
    if not isempty(tifFiles): 
        delete(tifFiles{:})                            # Delete all the TIF files
    if not isempty(aviFiles):
        delete(aviFiles{:})                            # Delete all the AVI files


    print('Creating movies...')
    if not boostLvl == 1:
        make_avi(Fsim,       [save_path,'/',save_filename, '.avi'],      0.2) # Make an avi of the noisy video
        make_avi(Fsim_clean, [save_path,'/',save_filename, 'clean.avi'], 0.2) # Make an avi of the clean video
    else
        make_avi(m.Fsim,       [save_path,'/',save_filename, '.avi'],     0.2)# Make an avi of the noisy video
        make_avi(m.Fsim_clean, [save_path,'/',save_filename, 'clean.avi'],0.2)# Make an avi of the clean video
    print('done.\n')

    print('Saving as TIF stacks...')
    if not boostLvl == 1
        write_TPM_movie(Fsim,       [save_path,'/',save_filename, '_mov.tif'])
        write_TPM_movie(Fsim_clean, [save_path,'/',save_filename, '_movClean.tif'])
        clear Fsim Fsim_clean
    else
        write_TPM_movie(m.Fsim,       [save_path,'/',save_filename, '_mov.tif'])
        write_TPM_movie(m.Fsim_clean, [save_path,'/',save_filename, '_movClean.tif'])
    print('done.\n')


    print('Creating separate model parts...')
    saveSimulationParts([save_path,'/',save_filename,'.mat'])
    print('done.\n')
    print('Finished!\n')

