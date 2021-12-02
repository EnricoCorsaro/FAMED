from .configuring_parameters import ConfiguringParameters
from .star_object import *
from .global_object import *
from .utils import *


__all__ = ['GLOBAL','CHUNK','ECHELLE','COMPLETE']

def GLOBAL(catalog_id, star_id, teff, force=True):
    """
    Helper function to run all steps, including plotting, of GLOBAL modality.

    Parameters
    ----------
    catalog_id : str
        Catalogue ID of the star (e.g. 'KIC' for Kepler).
    star_id : str
        ID of the star as a string (e.g. '0012008916' or '7037405').
    teff : float
        Effective temperature of the star in Kelvin.
    force : bool
        Flag to force the computation of the sliding pattern fit.
    """
    famed_obj = Global(catalog_id,star_id,teff)
    famed_obj.make_islands(force=force)
    famed_obj.find_islands(force=force)
    famed_obj.make_global_plots()
    return famed_obj

def CHUNK():
    print('This function is not yet implemented')

def ECHELLE():
    print('This function is not yet implemented')

def COMPLETE():
    print('This function is not yet implemented')





