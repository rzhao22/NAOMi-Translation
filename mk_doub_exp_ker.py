
def mk_doub_exp_ker(t_on, t_off, A, dt, *args):

    if len(args)!=0:
        dext_type = args[0]
    else:
        dext_type = 'mult'

    # invert the fisrt constant t_on to make the rest of the code simpler
    t_on = 1/t_on

    #invert the second constant t_off to make the rest of the code simpler
    t_off = 1/t_off

    
