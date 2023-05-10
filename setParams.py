import numpy as np
def setParams(dParams, paramsIn):
# function paramsOut = setParams(dParams, paramsIn)

    # paramsOut = setParams(dParams, paramsIn)
    #
    #  paramsOut = setParams(dParams, paramsIn)
    #
    #  The output params struct will include the fields of dParams and paramsIn.
    #  If the same field exists in both of them, the value will be taken from 'paramsIn'
    #
    #

    ###########################################################################

    paramsOut = dParams
    if not paramsIn:
        return

    fnames = paramsIn.__dict__.keys()
    for fname in fnames:
        if hasattr(paramsOut, fname):
            attr = getattr(paramsOut, fname)
            if isinstance(attr, dict):
                val = setParams(attr, getattr(paramsIn, fname))
            else:
                val = getattr(paramsIn, fname)
            setattr(paramsOut, fname, val)
        else:
            setattr(paramsOut, fname, getattr(paramsIn, fname))
    

    return paramsOut

    ###########################################################################
    ###########################################################################
