import numpy as np
from check_scan_params import check_scan_params
from psf_fft import psf_fft
import inspect

def setup_scan_volume_frame(neur_vol,PSF_struct,scan_params):

# function setup_scan_volume_frame(neur_vol,PSF_struct,scan_params)
# return scan_vol 

    scan_params  = check_scan_params(scan_params)                             # Check the scanning parameter class instance for missing elements

    scan_avg  = scan_params.scan_avg                                          # Get scanning stepping amount (how many sub-resolution steps between each sample point in the FOV)
    sfrac     = scan_params.sfrac                                             # Subsampling factor

    N1     = scan_params.vol_sz[0]                                          # Get length dimension
    N2     = scan_params.vol_sz[1]                                          # Get width dimension
    N3     = scan_params.vol_sz[2]                                          # Get depth dimension

    if not hasattr(PSF_struct, 'psf'):
        ValueError('Must provide PSF to scan!')
    else:
        PSF    = PSF_struct.psf                                               # Extract point spread function from PSF class instance

    if not inspect.isclass(PSF):                                                          # Get size of the PSF
        Np1,Np2,Np3 = PSF.shape                                              # warning! this line is not working, as the "shape" won't work for class instance!
    else:
        TypeError('Unknown input configuration!')

    if (N1 < Np1) or (N2 < Np2):
        ValueError('PSF extent is bigger than the volume!')                         # Check that the PSF fits inside the volume (transversally)
    if (N3 < Np3):
        ValueError('PSF depth is larger than the volume depth!')                    # Check that the PSF fits inside the volume (axially)


    if (hasattr(PSF_struct,'psfT')) and (PSF_struct.psfT):
        psfT = PSF_struct.psfT
        psfB = PSF_struct.psfB
        psfT.mask = psfT.mask/np.mean(psfT.mask[:])
        psfB.mask = psfB.mask/np.mean(psfB.mask[:])
        psfT.convmask = psfT.convmask/sum(psfT.convmask[:])
        psfB.convmask = psfB.convmask/sum(psfB.convmask[:])
        psfT.freq_psf = psf_fft([N1, N2, N3], psfT.convmask)  
        psfB.freq_psf = psf_fft([N1, N2, N3], psfB.convmask)  
    else:
        psfT = []
        psfB = []

    if not hasattr(PSF_struct, 'mask'):
        t_mask = []                                                           # No masking if mask is not supplied
    else:
        t_mask = PSF_struct.mask                                              # Extract mask from PSF struct
        t_thresh = 1e-5
        t_mask[t_mask<t_thresh] = t_thresh
    
    if hasattr(PSF_struct, 'colmask'):
        t_mask = t_mask*PSF_struct.colmask                                   # No masking if mask is not supplied

    if(t_mask):
        top_mask = 1/t_mask
        bot_mask = 1/t_mask
    else:
        top_mask = np.ones((N1/sfrac,N2/sfrac),'single')
        bot_mask = np.ones((N1/sfrac,N2/sfrac),'single')
    
    if(hasattr(psfT,'mask')):
        top_mask = top_mask*psfT.mask
        bot_mask = bot_mask*psfB.mask
    
    # create an empty class for scan_vol
    class Scan_vol:
        pass
    scan_vol = Scan_vol()
    scan_vol.top_mask = top_mask
    scan_vol.bot_mask = bot_mask

    del top_mask, bot_mask

    scan_vol.psfT = psfT
    scan_vol.psfB = psfB

    del psfT, psfB

    scan_vol.freq_psf = psf_fft([N1, N2, N3], PSF, scan_avg)                           # Pre-calculate the FFT of the PSF for faster scanning


    if not hasattr(PSF_struct, 'g_blur'):
        scan_vol.g_blur = []                                                            # No additional blurring if transversal blur function not supplied
    else:
        scan_vol.g_blur = PSF_struct.blur                                              # Extract point transversal blur function from PSF struct

    somaVol = [[None, None] for _ in range(len(neur_vol.gp_vals))]
    dendVol = [[None, None] for _ in range(len(neur_vol.gp_vals))]

    for i in range(len(neur_vol.gp_vals)):
        somaVol[i][0] = neur_vol.gp_vals[i][0][neur_vol.gp_vals[i][2]]      #Pre-alocate separate soma/dendrite indexing for faster activity modulation in the volume
        dendVol[i][0] = neur_vol.gp_vals[i][0][~neur_vol.gp_vals[i][2]]
        
        if t_mask is None:
            somaVol[i][1] = neur_vol.gp_vals[i][1][neur_vol.gp_vals[i][2]]
            dendVol[i][1] = neur_vol.gp_vals[i][1][~neur_vol.gp_vals[i][2]]
        else:
            somaVol[i][1] = neur_vol.gp_vals[i][1][neur_vol.gp_vals[i][2]] * t_mask[max(1, (somaVol[i][0] % (N1 * N2)))]
            dendVol[i][1] = neur_vol.gp_vals[i][1][~neur_vol.gp_vals[i][2]] * t_mask[max(1, (dendVol[i][0] % (N1 * N2)))]


    scan_vol.somaVol = somaVol
    del somaVol
    scan_vol.dendVol = dendVol
    del dendVol
 
    if neur_vol.bg_proc:  # t_mask multiplied by gp_vals{i,2}
        if t_mask is None:
            axonVol = neur_vol.bg_proc
        else:
            axonVol = [[None, None] for _ in range(len(neur_vol.bg_proc))]
            for i in range(len(neur_vol.bg_proc)):
                axonVol[i][0] = neur_vol.bg_proc[i][0]
                axonVol[i][1] = neur_vol.bg_proc[i][1] * t_mask[max(1, (axonVol[i][0] % (N1 * N2)))]
    else:
        axonVol = []


    scan_vol.axonVol = axonVol
    del axonVol

    if(hasattr(scan_params,'nuc_label')) and (scan_params.nuc_label>=1):
        if t_mask is None:
            nucVol = neur_vol.gp_nuc
        else:
            nucVol = [[None, None] for _ in range(len(neur_vol.gp_nuc))]
            for i in range(len(neur_vol.gp_nuc)):
                nucVol[i][0] = neur_vol.gp_nuc[i][0]
                nucVol[i][1] = t_mask[max(1, (nucVol[i][0] % (N1 * N2)))]
        scan_vol.nucVol = nucVol
    del nucVol
    
    return scan_vol
