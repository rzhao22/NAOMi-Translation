import numpy as np
import os
import math
import tifffile

def tifappend(tiflink,image,tag,blocksize,imnum,filename):

# tifappend(tiflink,image,tag,blocksize,imnum,filename)
# return tiflink
#
# Function to write append the next scanned image to disk
#
# 2019 - Alex Song

    if imnum % blocksize == 1 and imnum > 1:
        tiflink.close()
        del tiflink
        i = int(math.ceil(imnum / blocksize))
        path, name, _ = os.path.split(filename)
        filename2 = os.path.join(path, f"{name}_{i:05d}.tif")
        tiflink = tifffile.TiffWriter(filename2)

        tiflink.setTag(tag)
        tiflink.write(image)
        tiflink.writeDirectory()
    else:
        tiflink.setTag(tag)
        tiflink.write(image)
        tiflink.writeDirectory()

    return tiflink
