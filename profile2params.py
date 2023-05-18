

#A has to be a np array set up like this: [[[a,b],[c,d],[e,f]],[[g,h],[i,j],[k,l]]]

def profile2params(A, z_tilt=1):

    A_params = np.zeros((len(A),4))
    paramsM = np.zeros((len(A),4))
    
    for i in range(len(A)):
        im = A[i]
        imbin = (im>0).astype(int))
        box - 

    return A_params, paramsM
