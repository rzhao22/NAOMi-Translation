# [pos,pdf] = pseudoRandSample2D(sz, nsamps, width, weight, pdf, maxit)
# 
# Function to sample pseudo randomly within a 2D matrix, with partial
# exclusion of previously sampled locations. The inputs are:
#
#  - sz         - Size of input matrix
#  - nsamps     - Number of locations to sample
#  - width      - Width of Gaussian exclusionary zone
#  - weight     - Weight of Gaussian (0-1, 1 fully excludes previously
#                 sampled locations)
#  - pdf        - Probability density function of allowable locations to
#                 sample
#  - maxit      - Maximum number of iterations to look for new positions
#
# The outputs are:
#  - pos        - Sampled positions
#  - pdf        - Resulting probability density function of future
#                 positions that may be sampled
# 
# 2017 - Alex Song

###########################################################################
## Input parsing
import numpy as np

def pseudoRandSample2D(sz=None, nsamps=None, width=None, weight=None, pdf=None, maxit=None):

  nargin = 0

  if (sz != None):
    nargin += 1
  if (nsamps != None):
    nargin += 1
  if (width != None):
    nargin += 1
  if (weight != None):
    nargin += 1
  if (pdf != None):
    nargin += 1
  if (maxit != None):
    nargin += 1

  if(nargin<6):
    maxit  = 1000
  if(nargin<5):
    pdf    = np.ones(sz,'single')
  if(nargin<4):
    weight = 1
  if(nargin<3):
    width  = 2

  ###########################################################################
  ## Main function

  [X,Y] = np.meshgrid(range(-1 * np.ceil(2*width), np.ceil(2*width)))                             # Set up a mesh-grid
  gpdf  = np.single(-weight*np.exp(-(np.power(X,2) + np.power(Y,2))/(width^2)))                        # Greate a Gaussian PDF (proportional to)
  pos   = np.zeros(nsamps,2)                                                    # Initialize the position matrix
  i     = 1
  numit = 0

  while(i<=nsamps & numit<maxit):
    numit = numit+1
    rndpt = np.ceil(np.random.rand(1,2)*sz)
    if((pdf(rndpt(1),rndpt(2))-np.random.rand)>0):
      xc = int([max(0,2*width+1-rndpt(1)),min(0,sz(1)-rndpt(1)-2*width)+4*width]+1)
      yc = int([max(0,2*width+1-rndpt(2)),min(0,sz(2)-rndpt(2)-2*width)+4*width]+1)

      xi = int([max(1,rndpt(1)-2*width),min(sz(1),rndpt(1)+2*width)])
      yi = int([max(1,rndpt(2)-2*width),min(sz(2),rndpt(2)+2*width)])

      [xc,xi] = correctOffsetRoundError(xc,xi)
      [yc,yi] = correctOffsetRoundError(yc,yi)

      pdf(range(xi(1),xi(2)),range(yi(1),yi(2))) = pdf(range(xi(1),xi(2)),range(yi(1),yi(2)))+gpdf(range(xc(1),xc(2)),range(yc(1),yc(2)))
      pos(i,range(len(pos))) = rndpt
      i        = i+1                                                        # Add one to the # samples selected
      numit    = 0

  return [pos, pdf]

def correctOffsetRoundError(xc,xi):

  if np.diff(xc) > np.diff(xi):
    xc[1] = xc(1) + (np.diff(xc)-np.diff(xi))
  elif np.diff(xc) < np.diff(xi):
    xi[1] = xi(1) + (np.diff(xi)-np.diff(xc))

  return [xc,xi]


###########################################################################
###########################################################################
