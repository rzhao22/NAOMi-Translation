import numpy as np

def teardrop_poj(V, varargin):

    if nargin > 1:                                                             # Check if the tear-drop parameter p was provided
        p = varargin{1}
    else:
        p = 1

    if nargin > 2:                                                              # Check if the tear-drop parameter p was provided :
        plot_opt = varargin{2}
    else:
        plot_opt = false

    rr = sqrt(V(:,1).^2+V(:,2).^2)                                            # Get xy-plane radii
    tt = pi - atan(rr./abs(V(:,3))) - (V(:,3)>0)*pi                           # Create azymouth angles

    Vtear(:,1) = (V(:,1)./rr).*sin(tt).*(sin(0.5*tt).**p)                      # Set x-dimension
    Vtear(:,2) = (V(:,2)./rr).*sin(tt).*(sin(0.5*tt).**p)                      # Set y-dimension
    Vtear(:,3) = -cos(tt)                                                     # Set z-dimension

    Vtear(isnan(Vtear)) = 0                                                   # Sometimes NaNs happen: but only when a value should be zero (here at least)

    if plot_opt:                                                                # Optional plotting for debugging purposes
        scatter3(Vtear(:,1), Vtear(:,2), Vtear(:,3))
        axis equal
        view(3)

    return Vtear
