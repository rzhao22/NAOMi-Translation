# Uout = applyZernike(Uin,X,Y,k,abb)
#
# Applies the Zernike aberration given an input field, position and
# aberrations. The inputs are:
#
#   - Uin            - Input scalar field
#   - X              - X field positions
#   - Y              - Y field positions
#   - k              - Wavenumber [1/m]
#   - abb            - Zernike aberrations
#
# The outpus is:
#   - Uout           - Output scalar field
#
# 2017 - Alex Song

###########################################################################

import numpy as np

def applyZernike(Uin,X,Y,k,abb):
  phase = np.zeros(np.shape(X),'single')
  idxs = np.argwhere(abb != 0)
  if (idxs.size != 0):
    for i in idxs:
      phase = phase + abb(i) * zernike(i,X,Y)
  phase = phase - phase.mean()
  Uout = Uin * np.exp(1j * k * phase)
  return Uout