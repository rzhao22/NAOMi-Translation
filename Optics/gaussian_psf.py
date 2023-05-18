import numpy as np

def gaussian_psf(psflen = 10, Lambda = 0.92, sampling = [0.1,0.1,0.1], matSize = [101,101,1001], theta = 0):


##function gaussian_psf(psflen,lambda,sampling,matSize,theta)
# return psf,x,y,z
# length is the FWHM (1/e2 over 1.7) of the PSF, lambda is the wavelength,
# sampling is the spacing between pixels, and matSize is the output psf
# size. psf is the output, x y and z are the coordinate positions
# for FWHM axial length, we describe this as where the signal for a
# particular plane is 1/2 of the central plane
    

    if len(sampling)<3:
        sampling = [sampling[0], sampling[0], sampling[0]]
    

    n = 1
    # zr = 1/(sqrt(sqrt[1]-1))*psflen/2
    zr = psflen/2

    x = (np.arange(1,matSize[0]+1)-round(matSize[0]/2))*sampling[0]
    y = (np.arange(1,matSize[1]+1)-round(matSize[1]/2))*sampling[1]
    z = (np.arange(1,matSize[2]+1)-round(matSize[2]/2))*sampling[2]
    xg,yg,zg = np.meshgrid(x,y,z)

    R = np.array([[cosd(theta), -sind(theta)], [sind(theta), cosd(theta)]])
    xg2 = R[0,0]*xg+R[0,1]*zg
    zg2 = R[1,0]*xg+R[1,1]*zg
    xg = xg2
    zg = zg2

    psf = (np.exp(-2*np.pi*n*(xg**2+yg**2)/(zr*Lambda*(1+(zg/zr)**2)))/(1+(zg/zr)**2))**2
    return psf, x, y, z

def sind(x):
    I = x/180.
    y = np.sin(I * np.pi)
    mask = (I == np.trunc(I)) & np.isfinite(I)
    y[mask] = 0
    return y

def cosd(x):
    I = x/180.
    y = np.cos(I * np.pi)
    mask = (I == np.trunc(I)) & np.isfinite(I)
    y[mask] = 0
    return y