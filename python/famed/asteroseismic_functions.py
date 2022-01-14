import numpy as np
from scipy.interpolate import interp1d

from .utils import *

__all__ = ['assess_freq_asymptotic',
           'compute_scaling_dnu',
           'get_asymptotic_parameters',
           'get_linewidth',
           'get_modeid',
           'interpolate_epsilon']


def assess_freq_asymptotic(freq_list, enn, ell, dnu, epsilon, alpha, d01, numax, tolerance):
    """
    Verify that the input frequencies follow the expected asymptotic pattern.

    This function verifies that the input list of frequencies for radial and 
    dipole modes follows the expected asymptotic pattern (as set up by the 
    asymptotic parameters given as input) and discards those frequencies that 
    appear to exceed from the expected asymptotic pattern by an input tolerance
    value. 
         
    WARNING: the input frequency list is expected to include only frequencies 
    of the same angular degree (``ell``). The tolerance is provided as a 
    percentage of ``dnu``.

    Parameters
    ----------
    freq_list : array_like
        The list of `N` input frequencies to verify.
    enn : array_like
        The corresponding length `N` list of radial orders.
    ell : int
        Degree of the frequencies to be verified, either 0 or 1.
    dnu : float
        The large freqency spacing in microHz.
    epsilon : float
        The asymptotic phase term.
    alpha : float
        The curvature term.
    d01 : float
        The small fequency separation between l=0 and l=1 modes.
    numax : float
        The frequency of maximum oscillation power in microHz.
    tolerance: float
        Maximum allowed difference from expected as a percentage of ``dnu``

    Returns
    -------
    good_freq_index : list 
        List of indices of the input ``freq_list`` that follow the expected 
        asymptotic pattern
    
    """
    n_freq = len(freq_list)
    tolerance_threshold = dnu*tolerance
    good_freq_index = []
    residuals = np.zeros(n_freq)

    # When assessing the asymptotic position, also consider the curvature term
    # alpha. Assume for simplicity that alpha is the same for both l=1 and l=0.
    for i in range(0,n_freq):
        if ell == 0:
            residuals[i] = np.abs(freq_list[i] - dnu*(enn[i] + epsilon + alpha/2.*(enn[i] - numax/dnu)**2))
        else:
            residuals[i] = np.abs(freq_list[i] - dnu*(enn[i] + 0.5 + epsilon + alpha/2.*(enn[i] - numax/dnu)**2) + d01)
        if residuals[i] <= tolerance_threshold:
            good_freq_index.append(i)

    # Search for double l=0 peaks within the same radial order
    if ell == 0:
        duplicate_freq_index = []
        for i in range(0,n_freq):
            actual_enn = enn[i]
            tmp_match = np.where(enn == actual_enn)[0]
            if len(tmp_match) > 1:
                asymp_freq = dnu*(actual_enn + epsilon + alpha/2.*(actual_enn - numax/dnu)**2)
                closest_index = closest(asymp_freq,freq_list,index=True)
                duplicate_freq_index.append(tmp_match[tmp_match!=closest_index][0])

        if len(duplicate_freq_index) > 0:
            duplicate_freq_index = np.unique(np.array(duplicate_freq_index))

            # Remove these peaks from the list of good frequencies
            good_freq_index_final = []
            for j in range(0,len(good_freq_index)):
                actual_good_freq_index = good_freq_index[j]
                tmp_match = np.where(duplicate_freq_index == actual_good_freq_index)[0]
                if len(tmp_match) == 0:
                    good_freq_index_final.append(actual_good_freq_index)

            if len(good_freq_index_final) > 0:
                good_freq_index = good_freq_index_final

    return good_freq_index

def asymptotic_relation_radial(enn,numax,param):
    """
    Compute the asymptotic frequency for a radial mode.

    Parameters
    ----------
    enn : int
        The radial order of the mode.
    numax : float
        The frequency of maximum oscillation power
    param : array_like 
        Collection of the asypmtotic parameters in the order `[dnu, epsilon, 
        alpha]`.

    Returns
    -------
    freq_asymp : float
        The radial asymptotic frequency of given radial order `enn`.
    """
    freq_asymp = param[0]*(enn + param[1] + param[2]/2.*(enn - numax/param[0])**2)
    return freq_asymp
    

def asymptotic_relation_nonradial(enn,numax,ell,param):
    """
    Compute the asymptotic frequency for a non-radial mode.

    Parameters
    ----------
    enn : int
        The radial order of the mode.
    numax : float
        The frequency of maximum oscillation power in microHz.
    ell : int
        The angular degree of the mode
    param : array_like 
        Collection of the asypmtotic parameters in the order `[dnu, epsilon, 
        d0ell, alpha (dnu curvature), beta (small spacing curvature)]`.

    Returns
    -------
    freq_asymp : float
        The non-radial asymptotic frequency of given radial order `enn` and 
        degree `ell` in microHz.
    """
    freq_asymp = param[0]*(enn + param[1] + ell/2. + param[3]/2.*(enn - numax/param[0])^2 - param[4]*(enn - numax/param[0])) - param[2]
    return freq_asymp


def compute_scaling_dnu(numax, numax_threshold=300, numax_coeff_low=0.267, numax_coeff_high=0.22, numax_exponent_low=0.76, numax_exponent_high=0.797):
    """
    Compute the large frequency separation (`dnu`).

    This function computes a value of `dnu` from an input `numax` as based on 
    literature scaling relations. The relations are slightly different above 
    and below a given `numax_threshold`

    Parameters
    ----------
    numax : float
        The frequency of maximum oscillation power in microHz.
    numax_threshold : float, default: 300
        Threshold in `numax` which roughly separates red giants stars from 
        subgiant and main sequence stars. 
    numax_coeff_low : float, default: 0.267 
        Power law coefficient used for stars below *numax_threshold*.
    numax_coeff_high : float, default: 0.22 
        Power law coefficient used for stars above *numax_threshold*.
    numax_exponent_low : float, default: 0.76
        Power law exponent used for stars below *numax_threshold*.
    numax_exponent_high : float, default0.797
        Power law exponent used for stars above *numax_threshold*.  

    Returns
    -------
    dnu : float
        The calculated `dnu` value in microHz

    Notes
    -----
    The calculation takes the form:
    .. math:: \Delta\nu = a\nu_{\mathrm{max}}^{b}
    """
    
    # nuMax has to be in microHz. Following scaling relations calibrated by
    # Huber et al. 2011
    if numax < numax_threshold:
        dnu =  numax_coeff_low*numax** numax_exponent_low
    else:
        dnu =  numax_coeff_high*numax** numax_exponent_high
    return dnu


def get_asymptotic_parameters(numax, dnu, teff, d01_mass_offset=-0.073, d01_mass_slope=0.044, d01_offset=-0.063, d02_mass_offset=0.138, d02_mass_slope=-0.014, d02_offset=0.035, d03_slope=0.282, d03_offset=0.16, numax_sun=3150., dnu_sun=134.9, teff_sun=5777):
    """
    Compute relevant asymptotic parameters from `numax`, `dnu`, and `teff`.

    This function provides estimates of stellar mass and radius, as well as the 
    small frequency spacings d02, d01, d03.

    NOTE: The asymptotic parameters d02, d01, d03, are reliable only for RG 
    stars. When using this function for SG and MS, make sure that proper 
    corrections are applied.

    Parameters
    ----------
    numax : float
        The frequency of maximum oscillation power in microHz.
    dnu : float
        The large frequency spacing in microHz.
    teff : float
        Efffective temperature of the star in Kelvin.
    d01_mass_offset : float, default: -0.073
    d01_mass_slope : float, default: 0.044
    d01_offset : float, default: -0.063
    d02_mass_offset : float, default: 0.138
    d02_mass_slope : float, default: -0.014, 
    d02_offset : float, default: 0.035
    d03_slope : float, default: 0.282
    d03_offset : float, defualt: 0.16
    numax_sun : float, default: 3150.
        Solar reference value for frequency of maximum oscillations (`numax`).
    dnu_sun : float, default: 134.9 
        Solar reference value for large frequency spacing (`dnu`).
    teff_sun : float, default: 5777
        Solar reference vaule for effective temperature (`teff`).

    Returns
    -------
    asymptotic_parameters : dict
        A dictionary of the computed asymptotic values. The keys are `'mass'`, 
        `'radius'`, `'d01'`, `'d02'`, and `'d03'`. 
    """
    # Raw mass estimate from scaling relation
    mass = (numax/numax_sun)**3 * (dnu/dnu_sun)**-4 * (teff/teff_sun)**1.5

    # Raw radius estimate from scaling relation
    radius = (numax/numax_sun) * (dnu/dnu_sun)**-2 * (teff/teff_sun)**0.5

    # Scaling relation Dnu - d02 (Corsaro et al. 2012b)
    slope = d02_mass_offset + d02_mass_slope*mass
    d02 = slope*dnu + d02_offset

    # Scaling relation Dnu - d01 (Corsaro et al. 2012b)
    slope = d01_mass_offset + d01_mass_slope*mass
    d01 = slope*dnu + d01_offset

    # Scaling relation Dnu - d03 (Huber et al. 2010)     [1,20] microHz
    d03 = d03_slope*dnu + d03_offset

    asymptotic_parameters = {'d02':     d02,    
                             'd01':     d01,   
                             'd03':     d03,   
                             'mass':    mass,  
                             'radius':  radius,}

    return asymptotic_parameters


def get_linewidth(freq, teff, numax, numax_threshold=300):
    """
    Estimate the FWHM of the Lorentzian line profile.

    This function provides an empirical estimate of the linewidth of a radial 
    mode, based on the input value of the stellar Teff. It currently 
    incorporates different relations depending on which evolutionary stage the
    star belongs to. If `numax` < `numax_threshold` then the linewidth returned 
    is given by the linewidth-Teff relation from Corsaro et al. (2015a). 
    Otherwise it is given by the linewidth-numax-Teff relation from Ball et al.
    (2018).

    Parameters
    ----------
    freq : float
        The frequency at which to calculate linewidth in microHz.
    teff : float
        Effective temperature of the star in Kelvin.
    numax : float 
        Frequency of maximum oscillation power in microHz.

    Returns
    -------
    float
        Linewidth of Lorentzian profile for given stellar parameters in microHz.
    """
    if numax > numax_threshold:
        # Bilinear fits coefficients from Ball et al. 2018
        offset = np.array([-3.71*1.e0,-7.209*1.e1,-2.266*1.e-1,-2.190*1e3,-5.639*1e-1])
        teff_coeff = np.array([1.073*1.e-3,1.543*1.e-2,5.083*1e-5,4.302*1e-1,1.138*1e-4])
        numax_coeff = np.array([1.883*1.e-4,9.101*1e-4,2.715*1e-6,8.427*1e-1,1.312*1e-4])

        parameter = offset + teff_coeff*teff + numax_coeff*numax
        ln_fwhm = parameter[0]*np.log(freq/numax) + np.log(parameter[1]) + np.log(parameter[2])/(1 + (2*np.log(freq/parameter[3])/np.log(parameter[4]/numax))**2)
    else:
        # Polynomial fit from Corsaro et al. 2015b
        if teff < 5500.0:
            if teff > 4900.0:
                ln_gamma = np.log(51.6371 -0.0295244*teff + 5.61255e-06*teff**2 -3.53773e-10*teff**3)
            else:
                teff_local = 4900.
                ln_gamma = -2.12026
        else:
            # Total fit from RG to MS
            ln_gamma = 1463.49 -1.03503*teff + 0.000271565*teff**2 -3.14139e-08*teff**3 + 1.35524e-12*teff**4        
        if np.isscalar(freq):
            ln_fwhm = ln_gamma
        else:
            ln_fwhm = np.zeros(len(freq))
            ln_fwhm += ln_gamma

    return np.exp(ln_fwhm)



def get_modeid(freq, dnu, epsi_input, d01, numax, alpha):
    """
    Identifies the angular degree and radial order of a given mode frequency.

    This function obtains the order number and angular degree of an input 
    frequency value. It assumes that the input epsilon parameter is reliable 
    and incorporates the curvature term, which plays no effect if set to zero.

    Parameters
    ----------
    freq : float
        The frequency to identify.
    dnu : float
        The large frequency spacing in microHz.
    epsi_input : float
        The asymptotic phase term.
    d01 : float
        The small frequency separation between l=0 and l=1 modes.
    numax : float
        The frequency of maximum oscillation power in microHz.
    alpha : float
        The curvature term.

    Returns
    -------
    dict
        Dictionary with the keys `'order'` and `'degree'`.
    """
    if not isinstance(freq,np.ndarray): # Put single frequency into array form
        freq=np.array([freq])
       
    radial_order_array = np.tile(np.arange(60),(len(freq.T),1)).T
    epsi_radial_array = (freq - dnu*radial_order_array - dnu*alpha/2.*(radial_order_array - numax/dnu)**2)/dnu
    diff = epsi_radial_array - epsi_input
    #if no valid value, returns inf, probably breaks stuff later if it gets this
    min_vals = np.min(diff,axis=0,where=diff>=-.25,initial=np.inf)
    min_index = np.where(diff==min_vals)
    good_order = radial_order_array[min_index]
    
    epsi_local_radial = (freq - dnu*good_order - dnu*(alpha/2.*(good_order - numax/dnu)**2))/dnu
    epsi_local_dipole = (freq - dnu*(good_order+0.5) - dnu*(alpha/2.*(good_order - numax/dnu)**2) + d01)/dnu
    diff_radial = np.abs(epsi_local_radial - epsi_input)
    diff_dipole = np.abs(epsi_local_dipole - epsi_input)
    diff_radial = np.reshape(diff_radial,len(freq.T))
    diff_dipole = np.reshape(diff_dipole,len(freq.T))
    angular_degree = np.zeros(len(freq.T),dtype=int)
    angular_degree[diff_radial>diff_dipole]=1

    mode_id = {'order':      good_order,
               'degree':     angular_degree}
    
    return mode_id


def interpolate_epsilon(teff, dnu, dnu_threshold=30, epsilon_offset=0.601, epsilon_slope=0.632):
    """
    Calculate the asymptotic phase term.

    This function computes the value of epsilon for either a main-sequence, 
    subgiant star or a red giant star by exploiting empirical relations linking
    epsilon with the stellar parameters `teff` and `dnu`.

    Parameters
    ----------
    teff : float
        Effective temperature of the star in Kelvin.
    dnu : float
        The large frequency spacing in microHz.
    dnu_threshold : float, default: 30
        Threshold at which to switch from the `epsilon-teff` relation to 
        `epsilion-dnu` relation
    epsilon_offset : float, default: 0.601
        The offset of the `epsilon-dnu` relation calibrated by Corsaro et al. 
        (2012)b.
    epsilon_slope : float, default: 0.632
        The slope of the `epsilon-dnu` relation calibrated by Corsaro et al. 
        (2012)b.
    
    Returns
    -------
    epsilon : float
        The asymptotic phase term.
    dnu_array : ndarray
        The dnu values computed in the `epsilon-dnu` relation. Returns `None` if
        `dnu` > `dnu_threshold`. Used primarily for plotting. 
    epsi_array : ndarray
        The `epsilon` values from either the `epsilon-dnu` relation or the 
        data used in the `epsilon-teff` relation. Used primarily for plotting. 
    teff_array : ndarray
        The teff values of the data used in the `epsilon-teff` relation. Returns
        `None` if `dnu` < `dnu_threshold`. Used primarily for plotting. 
    fit : ndarray
        The output polynomial fit of data used in the `epsilon-teff` relation. 
        Returns `None` if `dnu` < `dnu_threshold`. Used for primarily plotting. 
    """

    if dnu > dnu_threshold:
        # Values from Lund+17 LEGACY
        epsi_array = np.array([1.114, 0.911, 1.356,0.988,1.114,1.445,1.325,1.377,1.374,
                               1.077,1.431,1.343,1.336,1.225,1.006,1.492,0.880,1.319,0.978,
                               1.392,1.054,1.358,1.112,1.368,1.117,1.504,1.075,1.455,1.547,
                               1.163,1.153,1.158,1.311,1.267,1.517,1.113,1.400,1.444,1.475,
                               1.439,1.337,1.007,0.958,1.095,1.343,1.045,1.067,1.529,1.139,
                               1.131,1.350,1.106,1.206,1.318,1.313,1.032,1.275,1.020,0.920,
                               1.516,1.200,1.061,1.437,1.461,1.281,0.928])

        teff_array = np.array([6326,6614,6045,6384,6193,5668,6107,5805,5846,6130,5853,6037,
                               6033,6313,6331,5674,6479,5832,6344,6068,6305,5775,6171,5811,
                               6248,5501,6235,5309,5488,6173,6343,6122,6067,6143,5719,6246,
                               5873,5677,5270,5852,6302,6400,6538,6278,6047,6253,6321,5457,
                               5860,6132,5949,6146,6177,5964,6045,6150,6140,6548,6642,5180,
                               6179,6276,5825,5750,5964,6580],dtype=float)

        epsi_array = epsi_array[np.argsort(teff_array)]
        teff_array = teff_array[np.argsort(teff_array)]
        result = np.polyfit(teff_array,epsi_array,4)
        fit = result[0]*teff_array**4 + result[1]*teff_array**3 + result[2]*teff_array**2 + result[3]*teff_array + result[4]

        # Need to remove duplicated temperature/fit points for interp
        teff_arr_uniq,fit_uniq=np.unique(np.array([[x,y] for x,y in zip(teff_array,fit)]),axis=0).T
        good_epsi = interp1d(teff_arr_uniq,fit_uniq,'cubic',fill_value='extrapolate')(teff)
        if good_epsi <= 0 :
            good_epsi = interp1d(teff_arr_uniq,fit_uniq,'linear',fill_value='extrapolate')(teff)
        dnu_array=None
    else:
        # epsilon-dnu relation by Corsaro et al. 2012b
        good_epsi = epsilon_offset + epsilon_slope*np.log10(dnu)
        dnu_array = np.linspace(0,1,300) * dnu_threshold + 0.15
        epsi_array = epsilon_offset + epsilon_slope*np.log10(dnu_array)
        teff_array = None
        fit = None
        
    return good_epsi, dnu_array, epsi_array, teff_array, fit
