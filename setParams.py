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

    fnames = fieldnames(paramsIn)
    for f in range(len(fnames)):
        if isfield(paramsOut, fnames{f}):
            if isstruct(paramsOut.(fnames{f})):
                val = setParams(paramsOut.(fnames{f}), paramsIn.(fnames{f}))
            else:
                val = paramsIn.(fnames{f})
            
            paramsOut.(fnames{f}) = val
        else:
            paramsOut.(fnames{f}) = paramsIn.(fnames{f})
        
    

    return paramsOut

    ###########################################################################
    ###########################################################################
