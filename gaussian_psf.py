# length is the FWHM (1/e2 over 1.7) of the PSF, lambda is the wavelength,
# sampling is the spacing between pixels, and matSize is the output psf
# size. psf is the output, x y and z are the coordinate positions
# for FWHM axial length, we describe this as where the signal for a
# particular plane is 1/2 of the central plane
import numpy as np

def gaussian_psf(psflen = None, lambdaWave = None,sampling = None, matSize = None,theta = None):

  nargin = 0

  if (psflen != None):
    nargin += 1
  if (lambdaWave != None):
    nargin += 1
  if (sampling != None):
    nargin += 1
  if (matSize != None):
    nargin += 1
  if (theta != None):
    nargin += 1

  if nargin<5:
    theta = 0

  if nargin<4:
    matSize = [101, 101, 1001]

  if nargin<3:
    sampling = [0.1, 0.1, 0.1]

  if len(sampling)<3:
    sampling = [sampling(1), sampling(1), sampling(1)]

  if nargin<2:
    lambdaWave = 0.92

  if nargin<1:
    psflen = 10

  n = 1
  zr = psflen/2

  x = (range(1,matSize(1))-round(matSize(1)/2))*sampling(1)
  y = (range(1,matSize(2))-round(matSize(2)/2))*sampling(2)
  z = (range(1,matSize(3))-round(matSize(3)/2))*sampling(3)
  [xg,yg,zg] = np.meshgrid(x,y,z)

  R = [np.cosd(theta), -1 * np.sind(theta), np.sind(theta), np.cosd(theta)]
  xg2 = R(1,1)*xg+R(1,2)*zg
  zg2 = R(2,1)*xg+R(2,2)*zg
  xg = xg2
  zg = zg2
  psf = np.exp(-2 * np.pi * n * (np.power(xg,xg) + np.power(yg,yg)) / (zr *lambdaWave * (1+np.power((zg / zr))))) / (1 + np.power((zg / zr),(zg / zr)))
  psf = np.power(psf,psf)
  return [psf,x,y,z]
