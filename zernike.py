# function Z = zernike(i, x, y)
#
# Creates the Zernike polynomial with mode index i. The inputs are
#
# - i     - zernike index
# - x     - x position of evaluated points
# - y     - y position of evaluated points
#
# The output is
# - Z     - zernike polynomial profile
#
# Adapted from Numerical Simulation of Optical Wave Propagation by Jason D. Schmidt (2010)
#
# 2017 - Alex Song

###########################################################################
import numpy as np
from scipy.special import gamma

def cart2pol(x, y):
  rho = np.sqrt(x ** 2 + y ** 2)
  phi = np.arctan2(y, x)
  return (rho, phi)

def pol2cart(rho, phi):
  x = rho * np.cos(phi)
  y = rho * np.sin(phi)
  return (x, y)

def zrf(n, m, r):
  R = 0
  for s in range(0 , (n-m)/2):
    R = R + (-1) ^ s * gamma(n-s+1) / (gamma(s+1) * gamma((n+m)/2-s+1) * gamma((n-m)/2-s+1)) * r ^ (n-2*s)
  return R

# Zernike index function
def zidx(i,j):
  idx = np.ceil(np.sqrt(0.25 + 2 * i) - 1.5)
  if j==1:
    out = idx
  elif j==2:
    out = np.mod(idx+1,2) *2.* np.floor((i-(idx+1) * idx/2)/2) + np.mod(idx,2) * (2 * np.ceil((i-(idx+1) * idx/2)/2)-1)
  return out

def zernike(i,x,y):

  [theta, r] = cart2pol(x, y)
  n = zidx(i,1)
  m = zidx(i,2)
  if m == 0:
    Z = np.sqrt(n + 1) * zrf(n, 0, r)
  else:
    if np.mod(i,2) == 0:
      Z = np.sqrt(2*(n+1))*zrf(n,m,r) * np.cos(m*theta)
    else:
      Z = np.sqrt(2*(n+1)) * zrf(n,m,r) * np.sin(m*theta)
  return Z


