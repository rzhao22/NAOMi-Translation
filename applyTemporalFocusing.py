def applyTemporalFocusing(psf, length, dz, offset=0):
    
# psf = applyTemporalFocusing(psf,length,dz,offset)
#
# This function applies temporal focusing restriction on the PSF by
# assuming the constraint on the axial profiles follows 1/sqrt(1+(dz/zR)^2)
# relationship, where dz is the offset axially and zR is the constrained
# width, and is related to the FWHM by a factor of 2*sqrt(3)
# This approximates the theoretical valeus from Durst 2006 (Optics Express)
# 
# Inputs:
#     psf     - 3D psf matrix (with axial scale dz)
#     length  - length of FWHM axially from temporal focusing (um)
#     dz      - step size in z per voxel
#     offset  - offset from center of psf matrix in z (in um)
# 
# Outputs:
#     psf     - temporally focused psf
#
# 2017 - Alex Song
#
###########################################################################


    zR = length/(dz*2*np.sqrt(3))
    z0 = np.ceil(len(psf)/2)

    zprofile = 1/np.sqrt(1+((z0-np.arange(1,len(psf)+1)+(offset*dz))/zR)**2)
    psf = np.multiply(np.reshape(zprofile, (1, 1, -1), psf)

    return psf
