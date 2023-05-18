def vec(x):

# vec(x)
# return y
# 
# This function just returnes the vectorized version of the input x into
# the vector y. If x is a order-N tensor with M_j elements in the j^th
# dimension, then y is of size [ \prod_{j=1}^N M_j ]x 1.
# 
# 2017 - Adam Charles

    y = x.flatten()
    return y

###########################################################################
###########################################################################