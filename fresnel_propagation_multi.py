from scipy.fft import fft2, ifft2, fftshift, ifftshift

def fresnel_propagation_multi(Uin, Lambda, dx, z, phi, nidx, nargout): #add nargout as input for python

# [Uout, UoutAll] = fresnel_propagation_multi(Uin, Lambda, dx, z, phi, nidx)
#
# Simulates optical propagation through a medium with phase masks applied
# at specified points. The inputs are
#
#  - Uin          - incident field
#  - Lambda       - wavelength [m]
#  - dx           - lateral pixel width at each step of propagatio [m]
#  - z            - position of each step of propagation [m]
#  - phi          - phaseshift at each step of propagation [rad]
#  - nidx         - index of refraction
#
# The outputs are
#  - Uout         - output field at end of propgation
#  - UoutAll      - output field at all steps of propagation
#
# Adapted from Numerical Simulation of Optical Wave Propagation by Jason D.
# Schmidt (2010)
#
# 2017 - Alex Song

###########################################################################

    Lambda   = Lambda/nidx
    sz       = Uin.shape                                                           
    [nx, ny] = np.meshgrid(np.arange(-sz[1]/2,sz[1]/2),np.arange(-sz[2]/2,sz[2]/2))
    nx = nx.T
    ny = ny.T
    k   = 2*np.pi/Lambda                                                          
    n   = len(z)                                                          
    df  = 1/(sz[0]*dx)
    dz  = z[1:n] - z[0:n-1]                                                 
    sc  = dx[1:n]/dx[0:n-1]                                                
    Uin = Uin*np.exp(1j*k*((nx*dx[0])**2 + (ny*dx[0])**2)*(1-sc[0])/(2*dz[0])).*phi[:,:,1]

    if (nargout>1):
        UoutAll = np.zeros_like(Uin, shape=(sz[0],sz[1],n))
        UoutAll[:,:,0] = Uin
    tol = 1e-12
    Q = np.exp(-1j*np.pi*Lambda*dz[0]*((nx*df[0])**2+(ny*df[0])**2)/sc[0])
    for i in range(n):
        if(i>0 and not(abs(dz[i]-dz[i-1])<tol and abs(df[i]-df[i-1])<tol and all(abs(sc[i]-sc[i-1])<tol))):
            Q = np.exp(-1j*np.pi*Lambda*dz[i]*((nx*df[i])**2+(ny*df[i])**2)/sc[i])
        
      if(i==n-1 and phi.shape[2]<n):
          Uin = ifftshift(ifft2(ifftshift(Q * fftshift(fft2(fftshift(Uin/sc[i]))))))
      else:
          Uin = phi[:, :, i+1] * ifftshift(ifft2(ifftshift(Q * fftshift(fft2(fftshift(Uin/sc[i]))))))

      if (nargout>1):
          UoutAll[:,:,i+1] = Uin

    Uout = Uin*np.exp(1j*k/2*(sc[n-2]-1)/(sc[n-2]*dz[i])*((nx*dx[n])**2+(ny*dx[n])**2))

    return [Uout, UoutAll]

###########################################################################
###########################################################################
