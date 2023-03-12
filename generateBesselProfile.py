import warnings

def generateBesselProfile(X,Y,rad,width,k,fl,type_string = 'gaussian',aper=float('inf'),offset=[0,0]):

# Uout = generateBesselProfile(X,Y,rad,aper,k,fl,offset)
#
# This function generates a Gaussian back aperture intensity profile with a
# fixed aperture size with a focused phase fl. The inputs are
#
# - X           - lateral position X [m]
# - Y           - lateral position X [m]
# - rad         - radius [m]
# - width       - width of annulus (FWHM) on back aperture [m]
# - k           - optical wavenumber [rad/m]
# - fl          - focal lenth of lens [m]
# - type        - one of 'gaussian', 'uniform', 'posax', 'negax' for shape
#                 of profile. Gaussian is symmetric, uniform is flat across
#                 the width, posax and negax are x*gauss(x) and are the
#                 profile formed by a positive or negative strength axicon
# - aper        - aperture distance [m]
# - offset      - offset of Gaussian position [m]
#
# The output is
# - Uout        - scalar field for a gaussian beam with apodization 1
#
# 2017 - Alex Song

###########################################################################


    X2 = X-offset[0]
    Y2 = Y-offset[1]
    rho = np.sqrt(X2**2+Y2**2)
    
    if (type_string=='gaussian'):
        width = width/(2*np.sqrt(np.log(2)))
        Uout = np.exp(-((rho-rad)/width)**2)  
    elif(type_string=='uniform'):
        Uout = ((rho-rad)<(width/2))*((rho-rad)>-(width/2))   
    elif((type_string=='posax') or (type_string=='negax')):
        width = width/1.133
        Uout = (rho-rad+width/np.sqrt(2))*np.exp(-((rho-rad+width/np.sqrt(2))/width)**2)*(rho-rad+width/np.sqrt(2)>0) 
    else:
        warnings.warn('Specified type of Bessel beam unavailable, defaulting to Gaussian')
      
    rho2 = X**2+Y**2
    Uout = Uout*(rho2<(aper**2)) # Fixed aperture
    Uout = Uout*np.exp(-1j*k/(2*fl)*rho2)                                      # Apply ideal phase

    return Uout
