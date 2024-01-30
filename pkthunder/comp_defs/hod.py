from halomod.hod import HOD, HODBulk


class HODVN(HODBulk):
    _defaults = {'M_min':2e12, 'M_0':4.3e10, 'alpha':0.24} # in SM/h units, z = 0 parameters
    def __init__(self, **model_parameters):
        return
    