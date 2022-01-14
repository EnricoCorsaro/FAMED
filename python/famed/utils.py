import numpy as np
import os


def gaussian_func(x, a, x0, sigma, c1, c2, c3):
    """
    Compute Gaussian function with up to quadratic background terms.

    Parameters
    ----------
    x : array
        Array of values at which to evaluate the function.
    a : float
        Height of the Gaussian function.
    x0 : float
        Location of Gaussian function.
    sigma : float
        Width of Gaussian function.
    c1 : float
        Zeroth-order background coefficient.
    c2 : float
        First-order background coefficient.
    c3 : float
        Second-order background coefficient.
    
    Returns
    -------
    array
        Gaussian function
    """
    return a*np.exp(-(x-x0)**2/(2*sigma**2)) + c1 + c2*x + c3*x**2


def smooth(x, window_len=11, window='hanning'):
    """
    Smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) at both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    Parameters
    ----------
    x : array
        The input signal 
    window_len : int
        The dimension of the smoothing window; should be an odd integer
    window : str
        The type of window from 'flat', 'hanning', 'hamming', 'bartlett', 
        'blackman'. Flat window will produce a moving average smoothing.

    Returns
    -------
    array
        The smoothed signal  
    """
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len<3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
    
    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    
    return y[window_len//2:-window_len//2+1]

def closest(val,arr,index=False):
    """
    Find either the closest value or index of closest value in an array.

    Parameters
    ----------
    arr : array
        Array to find the closest value in.
    val : float or int
        Input value to find the closest array value to.
    index : bool, default: False
        If True, return the index of the closest value instead of the value.

    Returns
    -------
    float or int
        If `index=False`, this is the closest value
    int
        If `index=True`, returns the index of the closest value isntead,
    """
    idx = np.argmin(np.abs(arr-val))
    if index:
        return idx
    else:
        return arr[idx]
