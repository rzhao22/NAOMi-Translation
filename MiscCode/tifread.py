import numpy as np
import tifffile

def tifread(fileTif,subImages=0):
# tifread(fileTif,subImages)
# return finalImage, infoImage
#Tifread - FileTif is a string specifying file location,
#subimages is len2 vector specifying start and end images (default all)

    infoImage = tifffile.TiffFile(fileTif).pages
    mImage = infoImage[0].width
    nImage = infoImage[0].height

    if subImages is None:
        subImages = [1, len(infoImage)]

    tifLink = tifffile.TiffFile(fileTif)
    im = tifLink.pages[subImages[0] - 1].asarray()
    finalImage = np.zeros((nImage, mImage, subImages[1] - subImages[0] + 1), dtype=im.dtype)
    finalImage[:, :, 0] = im

    for i in range(subImages[1] - subImages[0]):
        im = tifLink.pages[i + subImages[0]].asarray()
        finalImage[:, :, i + 1] = im

    tifLink.close()

    return finalImage, infoImage
