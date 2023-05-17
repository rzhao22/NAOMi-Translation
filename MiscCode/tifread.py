import numpy as np
def tifread(fileTif,subImages=0):
# function [finalImage, infoImage] = tifread(fileTif,subImages)
    #Tifread - FileTif is a string specifying file location,
    #subimages is len2 vector specifying start and end images (default all)

    infoImage=imfinfo(fileTif)
    mImage=infoImage[0].Width
    nImage=infoImage[0].Height

    if subImages == 0:
        subImages = [1, len(infoImage)]
        

    tifLink = Tiff(fileTif, 'r')
    tifLink.setDirectory(subImages[0])
    im = tifLink.read()
    finalImage=np.zeros((nImage,mImage,subImages[2]-subImages[1]+1),'like',im)
    finalImage[:,:,0]=im
    for i in range(subImages[1]-subImages[0]):
        tifLink.setDirectory(i+subImages[1])
        finalImage[:,:,i+1]=tifLink.read()
    
    tifLink.close()
    return finalImage, infoImage
