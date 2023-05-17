import numpy as np
def tifappend(tiflink,image,tag,blocksize,imnum,filename):
# function tiflink = tifappend(tiflink,image,tag,blocksize,imnum,filename)

    # tiflink = tifappend(tiflink,image,tag,blocksize,imnum,filename)
    #
    # Function to write append the next scanned image to disk
    #
    # 2019 - Alex Song

    if (np.mod(imnum,blocksize)==1 and imnum>1):
        tiflink.close()
        clear tiflink
        i = np.ceil(imnum/blocksize)
        [path,name,_] = fileparts(filename)
        filename2 = fullfile(path,sprintf([name '_#05d.tif'],i))
        tiflink = Tiff(filename2,'w')
        
        tiflink.setTag(tag)
        tiflink.write(image)
        tiflink.writeDirectory();  
    else:
        tiflink.setTag(tag)
        tiflink.write(image)
        tiflink.writeDirectory();  
    return tiflink
