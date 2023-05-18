

def check_cal_params(cal_params, prot_type):

# vol_params = check_cal_params(vol_params)
#  
# This function checks the elements of the struct vol_params to ensure
# that all fields are set. Fields not set are set to a list of default
# parameters. The struct checked is:
# 
#   - cal_params  - Struct with parameters for the calcium simulation
#     .ext_rate   - Extrusion rate for the (default = 265.73)
#     .ca_bind    - Calcium binding constant (default = 110)
#     .ca_rest    - Resting-state calcium concentration (default = 50e-9)
#     .ind_con    - Indicator concentration (default = 200e-6)
#     .ca_dis     - Calcium disassociation constant (default = 290e-9)
#     .ca_sat     - Optional calcium saturation parameter(default = 1)
#     .sat_type   - Type of dynamics to simulate (default = 'double')                                                                 
#     .dt         - Time-trace sampling rate - should be at least the video
#                   sampling rate (default = 1/30) 
#     .ca_amp     - Calcium transient amplitude (default = 130.917 for 
#                   GCaMP6; default = 0.05 for GCaMP3) 
#     .t_on       - Rising time-constant for calcium transient (default =
#                   3.1295 for GCaMP6; default = 1 for GCaMP3) 
#     .t_off      - Falling time-constant for calcium transient (default =
#                   0.020137 for GCaMP6; default = 1 for GCaMP3) 
#     .a_bind     - Binding rate for more detailed simulation (default =
#                   3.5) 
#     .a_ubind    - Unbinding rate for more detailed simulation (default =
#                   7) 
# 
# 2017 - Adam Charles and Alex Song

###########################################################################
## Run the checks

    if cal_params.is_empty(): #requires a method in class cal_params for is_empty()                                                   # Make sure that cal_params is a struct
        cal_params = {}

    # if (isfield(cal_params,'sat_type'))&&strcmp(cal_params.sat_type,'Ca_DE')
    #   cal_params.ext_rate  = 265.73;                                      # Somas needs an even higher extrusion rate
    #   cal_params.t_on      = 3.1295;
    #   cal_params.t_off     = 0.020137;
    #   cal_params.ca_amp    = 130.917;                                       # Somas needs an even higher extrusion rate
    # end

    # if (not hasattr(cal_params,'ext_rate')) or not bool(cal_params.ext_rate)         # Default extrusion rate
    #     cal_params.ext_rate = 265.73;                                          # Default for mouse layer 2/3: from Helmchen & Tank 2011 
    # end

    if (not hasattr(cal_params,'ca_bind')) or not bool(cal_params.ca_bind):           # Default calcium binding ratio
        cal_params.ca_bind = 110                                              # Default for mouse layer 2/3: from Helmchen & Tank 2011

    if (not hasattr(cal_params,'ca_rest')) or not bool(cal_params.ca_rest):           # Default calcium resting level
        cal_params.ca_rest = 50e-9                                            # ~30-100 n mol: Helmchen et al 1996, Maravall et al. 2000 via Helmchen & Tank 201

    if (not hasattr(cal_params,'ind_con')) or not bool(cal_params.ind_con):           # Default indicator concentration
        cal_params.ind_con =  200e-6

    if (not hasattr(cal_params,'ca_dis')) or not bool(cal_params.ca_dis):             # Default calcium dissociation constant
        cal_params.ca_dis = 290e-9                                           # Default for GCaMP6f: from Badura et al. 2014

    if (not hasattr(cal_params,'ca_sat')) or not bool(cal_params.ca_sat):             # Default to no saturation
        cal_params.ca_sat = 1

    if (not hasattr(cal_params,'sat_type')) or not bool(cal_params.sat_type):         # Default to a single cacium binding/unbinding equation
        cal_params.sat_type = 'double'

    if (not hasattr(cal_params,'dt')) or not bool(cal_params.dt):                     # Default to a a sampling rate of 30 Hz
        cal_params.dt = 1/30                                               # Default is 30Hz

    ###########################################################################
    ## Bind-unbind constants for full simulation

    if (not hasattr(cal_params,'a_bind')) or not bool(cal_params.a_bind):             # Default binding rate of 3.5
        cal_params.a_bind = 3.5

    if (not hasattr(cal_params,'a_ubind')) or not bool(cal_params.a_ubind)           # Default unbinding rate of 7
        cal_params.a_ubind = 7


    ###########################################################################
    ## Double exponential parameters for ca + double exp simulation

    # define dict to store vals to substitute for switch statements
    
    params_dict = {
    "gcamp6": [76.1251, 0.8535, 98.6173 292.3],
    "gcamp6f": [76.1251, 0.8535, 98.6173 292.3],
    "gcamp6s": [54.6943, 0.4526, 68.5461, 299.0833],
    "gcamp7": [230.917, 0.020137, 3.1295, 265.73],
    "gcamp3": [0.05, 1, 1, 265.73]
    }
    
    if (not hasattr(cal_params,'ca_amp')) or not bool(cal_params.ca_amp):             # Default Amplitude of 1

        cal_params.ca_amp = params_dict.get(prot_type.lower(),76.1251)[0]

        '''
        switch lower(prot_type)                                                # First calculate the dF/F curve (as in Badura et al.)
            case {'gcamp6','gcamp6f'};    cal_params.ca_amp = 76.1251;         # Fit using data from Chen et al. 2013
            case 'gcamp6s';               cal_params.ca_amp = 54.6943;         # Fit using data from Chen et al. 2013
            case 'gcamp7';                cal_params.ca_amp = 230.917;         
            case 'gcamp3';                cal_params.ca_amp = 0.05;            
            otherwise;                    cal_params.ca_amp = 76.1251;         # Fit using data from Chen et al. 2013
 '''
    
    if (not hasattr(cal_params,'t_on')) or not bool(cal_params.t_on):                 # Default on rate of 0.1

        cal_params.t_on = params_dict.get(prot_type.lower(),0.8535)[1]

        '''
        switch lower(prot_type)                                                # First calculate the dF/F curve (as in Badura et al.)
            case {'gcamp6','gcamp6f'};    cal_params.t_on = 0.8535;            # Fit using data from Chen et al. 2013
            case 'gcamp6s';               cal_params.t_on = 0.4526;            # Fit using data from Chen et al. 2013
            case 'gcamp7';                cal_params.t_on = 0.020137;
            case 'gcamp3';                cal_params.t_on = 1;
            otherwise;                    cal_params.t_on = 0.8535;            # Fit using data from Chen et al. 2013
'''
        
    if (not hasattr(cal_params,'t_off')) or not bool(cal_params.t_off):               # Default off rate of 1

        cal_params.t_off = params_dict.get(prot_type.lower(),98.6173)[2]

        '''
        switch lower(prot_type)                                                # First calculate the dF/F curve (as in Badura et al.)
            case {'gcamp6','gcamp6f'};    cal_params.t_off = 98.6173;          # Fit using data from Chen et al. 2013
            case 'gcamp6s';               cal_params.t_off = 68.5461;          # Fit using data from Chen et al. 2013
            case 'gcamp7';                cal_params.t_off = 3.1295;
            case 'gcamp3';                cal_params.t_off = 1;
            otherwise;                    cal_params.t_off = 98.6173;          # Fit using data from Chen et al. 2013
'''
        
    if (not hasattr(cal_params,'ext_rate')) or not bool(cal_params.ext_rate):         # extrusion parameter

        cal_params.ext_rate = params_dict.get(prot_type.lower(),265.73)[3]
        
        '''
        switch lower(prot_type)                                                # First calculate the dF/F curve (as in Badura et al.)
            case {'gcamp6','gcamp6f'};    cal_params.ext_rate = 292.3;         # Fit using data from Chen et al. 2013
            case 'gcamp6s';               cal_params.ext_rate = 299.0833;      # Fit using data from Chen et al. 2013
            case 'gcamp7';                cal_params.ext_rate = 265.73;        # Default for mouse layer 2/3: from Helmchen & Tank 2011 
            case 'gcamp3';                cal_params.ext_rate = 265.73;        # Default for mouse layer 2/3: from Helmchen & Tank 2011 
            otherwise;                    cal_params.ext_rate = 265.73;        # Default for mouse layer 2/3: from Helmchen & Tank 2011 
'''

    return cal_params
###########################################################################
###########################################################################
