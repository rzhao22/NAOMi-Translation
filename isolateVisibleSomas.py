
import numpy as np

def isolateVisibleSomas(vol, psf, vol_params, neur_params, *args):

    if len(args)!=0:
        thresh = args[0]
    else:
        thresh = 0
        
    # Get midpoint of volume (plane of imaging)
    vMiid = vol_params['col_sz'][2]/2

    # Get number of neural somas in the volume
    nNeur = len(vol['gp_nuc'])

    # Get the z (axial) locations of the neurons
    zLocs = vol['locs'][0:nNeur,-1]

    # Calculate the half-width of the point-spread function
    psfHalfWidth = widthestimate3D(psf['psf'])

    # The threshold is calculated based on the psf half-width + neural radius
    thresh = thresh + psfHalfWidtch[-1]/2 + neur_params['avg_rad'][0]

    # Find all the somas no more than 'thresh' away from the scan plane
    idList = np.where(abs(zLocs - vMid)<thresh)

    
    return idList
