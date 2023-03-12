
def widthestimate3D(matrix,fraction=0.5):
    #calculates fwhm if fraction is 0.5, fw 1/e for 1/e, etc., uses linear
    #interpolation

    vec1 = np.sum(np.sum(matrix, axis=2), axis=1)
    vec1 = np.squeeze(vec1)
    
    vec2 = np.sum(np.sum(matrix, axis=2), axis=0)
    vec2 = np.squeeze(vec2)
    
    vec3 = np.sum(np.sum(matrix, axis=1), axis=0)
    vec3 = np.squeeze(vec3)

    width1 = widthestimate(vec1,fraction)
    width2 = widthestimate(vec2,fraction)
    width3 = widthestimate(vec3,fraction)
    
    width = [width1, width2, width3]

    return width
