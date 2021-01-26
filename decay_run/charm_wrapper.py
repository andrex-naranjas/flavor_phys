'''
---------------------------------------------------------------
 Python wrapper for the CharmDecayWidths C++ shared library
 Authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)
 ---------------------------------------------------------------
'''
import ctypes
from ctypes import cdll
lib = cdll.LoadLibrary('./DecayWidthsPyBind/libcharmdecay.so')


class decay(object):
    def __init__(self):
        self.obj = lib.charm_new()

    def decay_width(self, MA_val, MB_val, MC_val, SA_val,
                          LA_val, JA_val, SL_val, AL_val, AR_val,
                          baryon, excMode, prodDecay):
        
        lib.charm_execute.restype = ctypes.c_double
        MA_val = ctypes.c_double(MA_val)
        MB_val = ctypes.c_double(MB_val)
        MC_val = ctypes.c_double(MC_val)
        SA_val = ctypes.c_double(SA_val)
        LA_val = ctypes.c_double(LA_val)
        JA_val = ctypes.c_double(JA_val)
        SL_val = ctypes.c_double(SL_val)
        AL_val = ctypes.c_double(AL_val)
        AR_val = ctypes.c_double(AR_val)
        baryon = ctypes.c_int(baryon)
        excMode = ctypes.c_int(excMode)
        prodDecay = ctypes.c_int(prodDecay)
        
        decay_value = lib.charm_execute(self.obj, MA_val, MB_val, MC_val, SA_val,
                                        LA_val, JA_val, SL_val, AL_val, AR_val,
                                        baryon,  excMode, prodDecay)        
        return decay_value
