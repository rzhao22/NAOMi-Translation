import numpy as np
def generateGaussianProfile(X, Y, rad, aper, k, fl, offset = [0,0]):
# generateGaussianProfile(X,Y,rad,aper,k,fl,offset)
# return Uout
#
# This function generates a Gaussian back aperture intensity profile with a
# fixed aperture size with a focused phase fl. The inputs are
#
# - X           - lateral position X [m]
# - Y           - lateral position X [m]
# - rad         - radius [m]
# - aper        - aperture distance [m]
# - k           - optical wavenumber [rad/m]
# - fl          - focal lenth of lens [m]
# - offset      - offset of Gaussian position [m]
#
# The output is
# - Uout        - scalar field for a gaussian beam with apodization 1
#
# 2017 - Alex Song

###########################################################################
    
    rho2 = X**2+Y**2
    Uout = np.exp(-((X-offset[0])**2+(Y-offset[1])**2)/(rad**2))                  # Gaussian intensity profile
    Uout = Uout*(rho2<(aper**2))                                              # Fixed aperture
    Uout = Uout*np.exp(-1j*k/(2*fl)*rho2)                                       # Apply ideal phase
    return Uout
