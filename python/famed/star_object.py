import numpy as np
import subprocess
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from pathlib import Path
import scipy.stats

from .configuring_parameters import ConfiguringParameters
from .utils import *

from . import diamonds_functions as diamonds
from . import asteroseismic_functions as astero

__all__ = ['FamedStar',]

class FamedStar(object):
    """
    Generic class to represent a star processed by the FAMED pipeline.

    Longer description

    Parameters
    ----------
    catalog_id : str
        Catalogue ID of the star (e.g. 'KIC' for Kepler ot 'TIC' for TESS).
    star_id : str
        ID of the star as a string (e.g. '0012008916' or '7037405').
    teff : float
        Effective temperature of the star in Kelvin.
    background_run_number : str or int
        Number of the background subfolder that contains the results of
        the background fit. 

    """

    def __init__(self, catalog_id, star_id, teff=None, background_run_number=None):
        self.catalog_id = catalog_id
        self.star_id = star_id
        if teff:
            self.teff = teff

        self.cp = ConfiguringParameters()
        peakbagging_results_dir = self.cp.diamonds_path/'PeakBagging'/'results'
        self.star_dir = peakbagging_results_dir/(self.catalog_id+self.star_id)

        if background_run_number:
            self.cp.background_run_number = background_run_number

        # Read in the fitted background parameters
        if self.cp.background_data_dir == '-99':
            self.cp.background_data_dir = self.cp.diamonds_path/'Background'/'data'
        else:
            self.cp.background_data_dir = Path(self.cp.background_data_dir)

        if self.cp.background_results_dir == '-99':
            self.cp.background_results_dir = self.cp.diamonds_path/'Background'/'results'
        else:
            self.cp.background_results_dir = Path(self.cp.background_results_dir)

        if self.cp.external_background_results_dir != '-99':
            self.cp.external_background_results_dir = Path(self.cp.external_background_results_dir)
 
        self.bgp = diamonds.get_background(self.catalog_id, self.star_id, self.cp.background_results_dir, self.cp.background_run_number, self.cp.external_background_results_dir, self.cp.external_background_filename_suffix)


    ### Methods to be used by mupltiple modalities (e.g. GLOBAL and CHUNK)
    
    def compute_acf_dnu(self, scaling_dnu, par_hist, asef_hist):
        """
        Compute `dnu` from the autocorrelation of the timeseries.

        This function computes the ACF**2 of the multi-modal sampling obtained 
        with DIAMONDS to evaluate an estimate of the large frequency separation
        `dnu`. The final value is estimated from a Gaussian fit to the ACF**2.

        Parameters
        ----------
        scaling_dnu : float
            The large frequency separation as estimated from the scaling laws.
        par_hist : array
            
        asef_hist : array
            

        Returns
        -------
        acf_dnu : float
            The peak of the Gaussian curve fit to the ACF**2 of the timeseries.
        interpolated_dnu : array
            Interpolated array of the `dnu` values. Used primarily for plotting.
        interpolated_acf : array
            Interpolated array of the ACF**2. Used primarily for plottoing.
        
        """
        # Set the limits based on scaling_dnu
        dnu_range_side = scaling_dnu * self.cp.dnu_acf_range_side
        top_dnu = scaling_dnu + dnu_range_side
        bottom_dnu = scaling_dnu - dnu_range_side

        # Autocorrelate the ASEF
        freqbin = par_hist[1] - par_hist[0]
        dnu_range_bins = round((top_dnu-bottom_dnu)/freqbin)+1
        bottom_dnu_bins = round(bottom_dnu/freqbin)
        lag = np.arange(dnu_range_bins, dtype='int') + bottom_dnu_bins
        lag = lag.astype(int)
        temp = asef_hist-np.mean(asef_hist)
        norm = np.sum(temp**2)
        result = np.correlate(temp,temp,mode='full')/norm
        result = result[len(result)//2:][lag]
        result = result-np.min(result)

        # Interpolate the results
        n_interpol = 101
        interpolated_dnu = np.arange(n_interpol)/(n_interpol-1)*(top_dnu - bottom_dnu) + bottom_dnu
        interpolated_acf = interp1d(lag*freqbin,result,'cubic',fill_value='extrapolate',bounds_error=False)(interpolated_dnu)
        j = np.argmax(interpolated_acf)
        best_acf = interpolated_acf[j]
        acf_dnu = interpolated_dnu[j]

        # Perform a Gaussian fit to the ACF^2 from the smoothed PSD
        x = interpolated_dnu
        y = interpolated_acf**2
        coeff,cov = curve_fit(gaussian_func,x,y,p0=[max(y),acf_dnu,acf_dnu/10,min(y),0,0])
        acf_dnu = coeff[1]
        return acf_dnu, interpolated_dnu, interpolated_acf


    def compute_asef(self, param, distr, n_bins=100):
        """
        Compute an Averaged Shifted Envelope Function (ASEF) 

        This function computes an Averaged Shifted Envelope Function (ASEF) 
        from an input sampling distribution. 

        Parameters
        ----------
        param : array
            The parameter over which we are computing the ASEF. In the case of 
            FAMED, this is the frequency.
        distr : array
            The distribution of the parameter to compute the ASEF for. In the 
            case of FAMED, this is the nested sampling iterations.
        n_bins : int, default: 100 
            The number of bins used to compute the Average Shifted Histogram.

        Returns
        -------
        param : array
            The rebinned parameter distribution. 
        distr : array
            The rebinned ASEF of the original distribution.
        """
        min_par = min(param)
        max_par = max(param)
        binwidth = (max_par-min_par)/n_bins

        # Now perform average shifting for ASH
        n_sim = 20    # Number of combined histograms
        simsize = binwidth/n_sim

        distr_rebin  = np.zeros(shape=(n_sim,n_bins))
        param_rebin = np.zeros(shape=(n_sim, n_bins))

        for k in range(n_sim):

            # Define bin-edges. 
            bins = [min_par + j*binwidth + k*simsize for j in range(n_bins+1)]

            rebin = scipy.stats.binned_statistic( 
                        x = param,
                        values = distr,
                        statistic='max',
                        bins = bins,
                        )
            bin_edges = rebin[1]
            bin_centres = (bin_edges[1:] + bin_edges[:-1]) / 2
            distr_rebin[k,:] = np.nan_to_num(rebin[0])
            param_rebin[k,:] = bin_centres

        param = param_rebin.mean(axis=0)
        distr = distr_rebin.mean(axis=0)

        return param, distr

    
    def evaluate_sampling_frequencies(self, par0, par_hist, freq, spsd, maximum, range_maximum):
        """
        Compute frequencies and uncertainties from the multi-modal fit.

        This function computes the frequencies and their uncertainties from 
        the multi-modal fit of DIAMONDS. The extracted frequencies do not 
        necessarily correspond to real oscillation peaks but they correspond to
        the extracted local maxima from the ASEF histogram. The function also
        provides the sampling counts, SPSD local maximum for each frequency, as
        well as an additional improvement to the definition of the frequency 
        range for each ASEF peak.

        Parameters
        ----------
        par0 : array
            The sampled frequency from DIAMONDS multi-modal fit.
        par_hist : array
            The rebinned sampled frequency as returned by the ASEF
        freq : array
            Frequency axis of the smoothed PSD.
        spsd : array
            The smoothed PSD.
        maximum : array
            The length N array of frequencies where maxima occur in the ASEF as
            returned by the hill climbing algorithm.
        range_maximum : array 
            The 2xN array with the range of frequency to consider for evaluating
            each maximum.

        Returns
        -------
        sampled_estimates : dict
            A dictionary containing the following keys: `freq1` for the  
            frequency, `freq1_sig` for the uncertainty, `sampling_counts` for 
            the total nested iterations values, and `spsd_maximum` for the 
            maxima of the smoothed PSD. 
        """
        n_maxima = len(maximum)
        freq1 = np.zeros(n_maxima)
        freq_sig1 = np.zeros(n_maxima)
        sampling_counts = np.zeros(n_maxima)
        spsd_maximum = np.zeros(n_maxima)
        nest_iter = np.arange(len(par0))
        range_maximum_old = range_maximum
        freqbin = freq[1]-freq[0]

        for i in range(0, n_maxima):
            # Consider the frequency range around the local maximum
            iterations = 0

            while iterations < self.cp.max_iterations_frequency:
                upper_bound = range_maximum[1,i]
                lower_bound = range_maximum[0,i]
                tmp_freq_range = np.where((freq >= lower_bound) & (freq <= upper_bound))[0]

                # Make sure that within each range, at least one frequency
                # point is found.
                while len(tmp_freq_range)<1:
                    upper_bound += freqbin/2
                    lower_bound -= freqbin/2
                    tmp_freq_range = np.where((freq >= lower_bound) & (freq <= upper_bound))[0]
                    range_maximum[1,i] = upper_bound
                    range_maximum[0,i] = lower_bound

                tmp_range = np.where((par0 <= upper_bound) & (par0 >= lower_bound))[0]
                par0_range = par0[tmp_range]
                spsd_range = spsd[tmp_freq_range]

                # Count the total of the nested iteration values falling in this
                # maximum bin. Note this is not the number of nested sampling
                # points, but it includes the actual nested iteration value from
                # each point (which can be up to several thousands for an
                # individual sampling point). In this way it is possible to
                # better weight the local maxima by the actual level of
                # likelihood that they have reached during the sampling.
                sampling_counts[i] = np.sum(nest_iter[tmp_range])

                # Save the maximum smoothed PSD in the region of this maximum
                spsd_maximum[i] = np.max(spsd_range)

                # Weighted by nested iteration value
                freq1[i] = np.sum(par0_range*tmp_range**2)/np.sum(tmp_range**2)
                freq_sig1[i] = np.sqrt(np.sum((par0_range-freq1[i])**2*tmp_range**2)/np.sum(tmp_range**2))
                
                if iterations == self.cp.max_iterations_frequency - 1:
                    break

                # Improve frequency ranges for computation of uncertainties and
                # evaluation of the number of sampling points. Do not exceed
                # ranges by more than estimated sigma * X times (usually 2)
                if upper_bound > freq1[i]+freq_sig1[i]*self.cp.max_sigma_range:
                    range_maximum[1,i] = freq1[i]+freq_sig1[i]*self.cp.max_sigma_range
                if lower_bound < freq1[i]-freq_sig1[i]*self.cp.max_sigma_range:
                    range_maximum[0,i] = freq1[i]-freq_sig1[i]*self.cp.max_sigma_range

                # Try to make ranges extend by at least sigma * Y times (usually 1) on each side
                left_freq = freq1[i]-freq_sig1[i]*self.cp.min_sigma_range
                right_freq = freq1[i]+freq_sig1[i]*self.cp.min_sigma_range
                
                if lower_bound > left_freq:
                    if i == 0:
                        if left_freq <= np.min(par_hist):
                            range_maximum[0,i] = np.min(par_hist)
                        if left_freq > np.min(par_hist):
                            range_maximum[0,i] = left_freq
                    else:
                        if left_freq <= range_maximum[1,i-1]:
                            range_maximum[0,i] = range_maximum[1,i-1]
                        if left_freq > range_maximum[1,i-1]:
                            range_maximum[0,i] = left_freq

                if upper_bound < right_freq:
                    if i == n_maxima-1:
                        if right_freq >= np.max(par_hist):
                            range_maximum[1,i] = np.max(par_hist)
                        if right_freq < np.max(par_hist):
                            range_maximum[1,i] = right_freq
                    else:
                        if right_freq >= range_maximum[0,i+1]:
                            range_maximum[1,i] = range_maximum[0,i+1]
                        if right_freq < range_maximum[0,i+1]:
                            range_maximum[1,i] = right_freq

                # Make sure to also include the original ASEF maximum frequency
                if range_maximum[0,i] > maximum[i]:
                    range_maximum[0,i] = maximum[i]
                if range_maximum[1,i] < maximum[i]:
                    range_maximum[1,i] = maximum[i]

                iterations +=1

        sampled_estimates = {'freq1':             freq1,             
                             'freq_sig1':         freq_sig1,   
                             'sampling_counts':   sampling_counts,
                             'spsd_maximum':      spsd_maximum}

        return sampled_estimates


    def get_minimum_freq_separation(self, dnu, global_flag=True):
        """
        Set the minimum separation of frequencies depending on modality.

        This function computes the minimum separation in frequency to define 
        the number of bins in the ASEF histogram. This is based on an estimate 
        of the minimum width required to obtain a good sampling of the actual 
        frequency peaks obtained in the multi-modal fits.

        Parameters
        ----------
        dnu : float
            The large frequency separation.
        global_flag : bool
            Flag to identify if running GLOBAL or CHUNK modality.

        """
        if global_flag:
            min_separation = dnu/self.cp.min_sep_scaling_global
        else:
            if dnu < self.cp.dnu_rg:
                min_separation = dnu/self.cp.min_sep_scaling_chunk_rg
            else:
                d02 = self.cp.d02_unique_slope*dnu + self.cp.d02_unique_offset
                min_separation = d02/self.cp.min_sep_scaling_chunk_ms
        return min_separation

    
    def get_range_divisions(self, par_hist, asef_hist, index_maximum, chunk=False):
        """
        Compute both the ranges and divisions from the ASEF histogram

        This routine obtains both ranges and divisions for each local maximum 
        peak found in the ASEF histogram. The ranges are defined as those 
        boundaries where the ASEF of the local maximum peak starts to rise 
        (left edge) and stops to decrease (right edge), while the divisions 
        are defined as the mid points between adjacent local maxima.

        Parameters
        ----------
        par_hist : array
            The rebinned parameter distribution from the ASEF. 
        asef_hist : array
            The rebinned distribution from the ASEF.
        index_maximum : array
            The length N array of indices of maxima returned from the hill 
            climbing algorithm.
        chunk : bool
            Flag to identify if running GLOBAL or CHUNK modality

        Returns
        -------
        range_maximum : array
            The 2xN array of ranges. 
        divisions_maximum
            The 2xN array of divisions.
        """
        maximum = par_hist[index_maximum]
        asef_maximum = asef_hist[index_maximum]
        n_maxima = len(maximum)
        n_bins_tot = len(par_hist)
        n_bins_max = round(n_bins_tot/self.cp.n_bins_max_fraction)
        asef_threshold = (self.cp.dp_isla['max_nested_it']+self.cp.dp_isla['n_live'])*3/4

        range_maximum = np.zeros((2,n_maxima))
        divisions_maximum = np.zeros((2,n_maxima)) 

        # At least 2 maxima should be present to define a division array.
        # Otherwise, take as first and last division only the first and last
        # frequency value of the given frequency range.
        if n_maxima >= 2:
            for i in range(0, n_maxima):
                if i == 0:
                    divisions_maximum[0,i] = np.min(par_hist)
                else:
                    divisions_maximum[0,i] = divisions_maximum[1,i-1]

                if i == n_maxima-1:
                    divisions_maximum[1,i] = np.max(par_hist)
                else:
                    divisions_maximum[1,i] = (maximum[i+1]+maximum[i])/2.0

            avg_spacing = np.mean(np.diff(maximum))
            first_division = maximum[0] - avg_spacing
            last_division = maximum[-1] + avg_spacing
            
            if first_division > divisions_maximum[0,0]:
                divisions_maximum[0,0] = first_division
            
            if last_division < divisions_maximum[1,-1]:
                divisions_maximum[1,-1] = last_division
        else:
            divisions_maximum[0,0] = np.min(par_hist)
            divisions_maximum[1,0] = np.max(par_hist)

        # Find the indices for the ranges around each local maximum. Use a X%
        # drop tolerance on the ASEF for those maxima corresponding to peaks
        # that are close to the maximum value of the ASEF. This is done in order
        # to optimize the identification of ranges for possible plateau regions.
        for i in range(0, n_maxima):
            index_max_act = index_maximum[i]

            # Left bound
            j = 0
            while (asef_hist[index_max_act-j] >= asef_hist[index_max_act-j-1]) and (index_max_act-j-1 > 0):
                j += 1

            if index_max_act-j < 0:
                range_maximum[0,i] = np.min(par_hist)
            else:
                range_maximum[0,i] = par_hist[index_max_act-j]
                if i >= 1:
                    if range_maximum[0,i] <= range_maximum[1,i-1]:
                        range_maximum[0,i] = range_maximum[1,i-1]

            # Right bound
            if chunk:
                if maximum[i] > self.cp.dnu_sg:
                    drop_tolerance = asef_hist[index_max_act]*self.cp.drop_tolerance_chunk_ms
                else:
                    drop_tolerance = asef_hist[index_max_act]*self.cp.drop_tolerance_chunk_rg
            else:
                drop_tolerance = asef_hist[index_max_act]*self.cp.drop_tolerance_global

            j = 0
            for k in range(index_max_act, len(par_hist)-1):
                if (asef_hist[k] >= (asef_hist[k+1] - drop_tolerance)):
                    j+=1
                else:
                    break

            if index_max_act+j >= len(par_hist)-1:
                range_maximum[1,i] = np.max(par_hist)
            else:
                range_maximum[1,i] = par_hist[index_max_act+j]

            # Check that actual frequency range does not exceeed the division
            if i <= n_maxima-2:
                if range_maximum[1,i] > divisions_maximum[1,i]:
                    range_maximum[1,i] = divisions_maximum[1,i]

                tmp_right_exceed = np.where((par_hist > range_maximum[1,i]) & (par_hist < maximum[i+1]))[0]

                if len(tmp_right_exceed)>0:
                    par_hist_exceed = par_hist[tmp_right_exceed]
                    index_right_maximum = closest(range_maximum[1,i],par_hist,index=True)
                    asef_right_maximum = asef_hist[index_right_maximum]
                    asef_hist_exceed = asef_hist[tmp_right_exceed]

                    min_asef_index=np.argmin(asef_hist_exceed)
                    if asef_hist_exceed[min_asef_index] < asef_right_maximum:
                        if par_hist_exceed[min_asef_index] <= divisions_maximum[1,i]:
                            range_maximum[1,i] = par_hist_exceed[min_asef_index]
                        else:
                            range_maximum[1,i] = divisions_maximum[1,i]


            # Extend right range beyond if the local maxiumu is the last one
            # of the list and there is an ASEF value smaller that that at
            # the right bound.
            if i == n_maxima-1:
                tmp_right_exceed = np.where(par_hist > range_maximum[1,i])[0]
                if len(tmp_right_exceed)>0:
                    par_hist_exceed = par_hist[tmp_right_exceed]
                    index_right_maximum = closest(range_maximum[1,i],par_hist,index=True)
                    asef_right_maximum = asef_hist[index_right_maximum]
                    asef_hist_exceed = asef_hist[tmp_right_exceed]

                    min_asef_index=np.argmin(asef_hist_exceed)
                    if asef_hist_exceed[min_asef_index] < asef_right_maximum:
                        range_maximum[1,i] = par_hist_exceed[min_asef_index]


            # In the CHUNK modality verify that if the local maxiumum is not
            # a prominent one, i.e. its ASEF is below 3/4, the ranges do not
            # extend over a given maximum number of allowed bins. This will
            # prevent from having very wide ranges if the peak is small
            # (surrounded by flat ASEF region).
            if chunk:
                tmp_bins_range = np.where((par_hist >= range_maximum[0,i]) & (par_hist <= range_maximum[1,i]))[0]
                if (len(tmp_bins_range) > n_bins_max) & (asef_maximum[i] < asef_threshold):
                    # If here then find the range that is closest to the
                    # local maximum and make the other boundary not exceed
                    # the maximum number of bins allowed starting from the
                    # selected bound.
                    index_range_closest = closest(range_maximum[:,i],maximum[i],index=True)
                    print(index_range_closest)
                    if index_range_closest == 0:
                        index_left_range = closest(range_maximum[0,i],par_hist,index=True)
                        if index_left_range + n_bins_max <= len(par_hist)-1:
                            range_maximum[1,i] = par_hist[index_left_range+n_bins_max]
                    else:
                        index_right_range = closest(range_maximum[1,i],par_hist,index=True)
                        if index_right_range - n_bins_max >= 0:
                            range_maximum[0,i] = par_hist[index_right_range-n_bins_max]

        return range_maximum,divisions_maximum

    
    def hill_climbing(self, par_hist, asef_hist, threshold, minimum_bin_separation):
        """
        Locate local maxima of input distrubition via hill climbing algorithm.

        Parameters
        ----------
        par_hist : array
            The x-axis of our distribution. For FAMED, this is the frequency.
        asef_hist : array
            The distribution. In the case of FAMED, this is the ASEF histogram.
        threshold : float
            The threshold in terms of distribution units to consider a 
            significant peak. 
        minimum_bin_separation : float
            The minimum number of adjacent bins required to consider two local
            maxima separated. 

        Returns
        -------
        index_maximum : array
            The indices of the maxima found by the hill climbing algorithm.

        """
        n_bins = len(par_hist)
        index_maximum = [0]
        new_i = 2

        for i in range(2,n_bins-1):
            if i != new_i:
                continue
            distr_start = asef_hist[i]
            new_i = i + 1
            distr_next = asef_hist[new_i]

            # Check if we are in a descending phase. If so, update the starting
            # bin at each step, until a local minimum is reached.
            while (distr_next < distr_start) & (new_i < n_bins-1):
                new_i = new_i + 1
                distr_start = distr_next
                distr_next = asef_hist[new_i]

            # Now check if we are in an ascending phase.
            # If so, update the final bin until a local maximum is reached.
            distr_act = distr_start
            while (distr_next > distr_act) & (new_i < n_bins-1):
                new_i = new_i + 1
                distr_act = distr_next      # Actual value of the distribution
                distr_next = asef_hist[new_i]   # Next value of the distribution

                # Store the local maximum only if the rising phase has an amount
                # of distribution variation more than given threshold.
                if (distr_act - distr_start) >= threshold:
                    index_maximum.append(new_i-1)

        index_maximum = index_maximum[1:]
        good_index_maximum = [0]
        index_bad = -1

        # Remove local maxima from list if within minimum_bin_separation bins
        # of one another
        for i in range(0,len(index_maximum)-1):
            if i == index_bad:
                continue
            index1 = index_maximum[i]      # where(par_hist eq maximum(i))
            index2 = index_maximum[i+1]    # where(par_hist eq maximum(i+1))
            if (index2 - index1) <= minimum_bin_separation:
                distr1 = asef_hist[index1]
                distr2 = asef_hist[index2]
                # Consider only the maximum with the larger peak
                if distr1 >= distr2:
                    index_bad = i + 1      # Right peak is bad, don't save it
                    good_index_maximum.append(index_maximum[i])
            else:
                good_index_maximum.append(index_maximum[i])

        # Add the last local maximum of the group if it is good
        if index_bad != len(index_maximum)-1:
            good_index_maximum.append(index_maximum[len(index_maximum)-1])

        index_maximum = good_index_maximum[1:]

        return index_maximum

