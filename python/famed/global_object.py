import numpy as np
import os
import shutil
import subprocess
import glob
import pickle
from statistics import median_high
import matplotlib.pyplot as plt

from .star_object import *
from .utils import *

from . import diamonds_functions as diamonds
from . import asteroseismic_functions as astero
from . import plotting_functions as famed_plots

__all__ = ['Global']

class Global(FamedStar):
    """
    FamedStar sub-class specific to the running of the GLOBAL modality.    

    Parameters
    ----------
    catalog_id : str
        Catalogue ID of the star (e.g. 'KIC' for Kepler).
    star_id : str
        ID of the star as a string (e.g. '0012008916' or '7037405').
    teff : float
        Effective temperature of the star in Kelvin.
    background_run_number : str or int
        Number of the background subfolder that contains the results of
        the background fit.
    load_islands : bool, default: False
        Flag to read the pickled object if `make_islands` or `find_islands` has 
        already been ran.
    """
    def __init__(self, catalog_id, star_id, teff, background_run_number=None, load_islands=False):
        FamedStar.__init__(self, catalog_id, star_id, teff, background_run_number)

        if self.cp.print_on_screen:
            print('-------------------------------------------------')
            print(' Performing GLOBAL modality for ' + catalog_id + star_id + '.')
            print('-------------------------------------------------')

        # Create output directories if not already present
        if not os.path.isdir(self.star_dir/self.cp.isla_subdir):
            os.makedirs(self.star_dir/self.cp.isla_subdir,exist_ok=True)
        if not os.path.isdir(self.star_dir/self.cp.pb_subdir):
            os.makedirs(self.star_dir/self.cp.pb_subdir,exist_ok=True)
        if not os.path.isdir(self.star_dir/self.cp.as_subdir):
            os.makedirs(self.star_dir/self.cp.as_subdir,exist_ok=True)
        if not os.path.isdir(self.star_dir/self.cp.figs_subdir):
            os.makedirs(self.star_dir/self.cp.figs_subdir,exist_ok=True)
        if not os.path.isdir(self.star_dir/self.cp.summary_subdir):
            os.makedirs(self.star_dir/self.cp.summary_subdir,exist_ok=True)
        if not os.path.isdir(self.star_dir/self.cp.as_subdir/'data'):
            os.makedirs(self.star_dir/self.cp.as_subdir/'data',exist_ok=True)

        if load_islands and self.cp.save_progress_pickle:
            if self.cp.print_on_screen:
                print('Loading saved information from the pickled object.')
            # Load any variables into attributes that should haven been
            # loaded or created during make_islands. Will ignore new config
            temp = pickle.load(open(self.star_dir/(self.catalog_id+self.star_id+'.pickle'),'rb'))
            self.__dict__ = temp.__dict__.copy()

        else:
            self.modality='GLOBAL'
            # Copy the configuring parameters to the summary/ directory
            target_dir = self.star_dir/self.cp.summary_subdir
            if self.cp.local_config:
                shutil.copy('famed_config.yml',target_dir/('famed_config_'+self.catalog_id + self.star_id + '_' + self.cp.isla_subdir + '_' + self.cp.global_subdir+'_GLOBAL.yml'))
            else:
                shutil.copy(self.cp.famed_path/'famed_config.yml',target_dir/('famed_config_'+self.catalog_id + self.star_id + '_' + self.cp.isla_subdir + '_' + self.cp.global_subdir+'_GLOBAL.yml'))
            shutil.copy(self.cp.configuring_parameters_file,target_dir/('famed_configuring_parameters_'+self.catalog_id + self.star_id + '_' + self.cp.isla_subdir + '_' + self.cp.global_subdir+'_GLOBAL.txt'))

            
    def make_islands(self,force=False):
        """ 
        Compute a global multi-modal fit with DIAMONDS.

        This function creates the initial multi-modal fit used in the GLOBAL
        modality of FAMED. It generates the directory structure to be used for 
        future steps as well.
        """
        diamonds.set_peakbagging(self.catalog_id, self.star_id, self.bgp, self.cp.diamonds_path, self.cp.background_data_dir, self.cp.background_results_dir, self.cp.external_background_results_dir, self.cp.dnu_cl, self.cp.dnu_tip, self.cp.n_dnu_envelope, self.cp.n_sigma_envelope, self.cp.n_sigma_envelope_cl, self.cp.n_sigma_envelope_tip, self.cp.numax_threshold, self.cp.numax_coeff_low, self.cp.numax_coeff_high, self.cp.numax_exponent_low, self.cp.numax_exponent_high,force=force)
        
        # Read input PSD and global asteroseismic parameters
        peakbagging_data_dir = self.cp.diamonds_path/'PeakBagging'/'data'
        freq, psd = np.loadtxt(peakbagging_data_dir/(self.catalog_id + self.star_id + '.txt'), unpack=True)
        gauss_par = np.loadtxt(self.star_dir/'gaussianEnvelopeParameters.txt')
        numax = gauss_par[1]
        freqbin = freq[1]-freq[0]
        dnu = astero.compute_scaling_dnu(numax, self.cp.numax_threshold, self.cp.numax_coeff_low, self.cp.numax_coeff_high, self.cp.numax_exponent_low, self.cp.numax_exponent_high)

        # Run DIAMONDS in a global multi-modal fit in order to identify the
        # chunks and estimate global frequencies and the value of DeltaNu
        avg_fwhm = np.mean(astero.get_linewidth(np.array([np.min(freq),np.max(freq)]), self.teff,numax, self.cp.numax_threshold))
        smth_bins = int(avg_fwhm/freqbin)
        spsd = smooth(psd, window_len=smth_bins, window='flat')

        # Write the prior file for the global multi-modal fit
        with open(self.star_dir/self.cp.isla_subdir/(self.cp.prior_filename + '_' + self.cp.global_subdir + '.txt'),'w') as f:
            f.write('{}    {}\n{}    {}\n'.format(np.min(freq),np.max(freq),0,np.mean(spsd)))
            
        # Evaluate linewidth at nuMax and use it for the global multi-modal fit
        if dnu <= self.cp.dnu_threshold:
            if dnu <= self.cp.dnu_tip:
                linewidth = astero.get_linewidth(numax, self.teff, numax, self.cp.numax_threshold)/self.cp.fwhm_global_scaling_tip
            else:
                linewidth = astero.get_linewidth(numax, self.teff, numax, self.cp.numax_threshold)/self.cp.fwhm_global_scaling
        else:
            linewidth = astero.get_linewidth(numax, self.teff, numax, self.cp.numax_threshold)
    
        linewdith = float(linewidth)
        
        # Perform the multi-modal fit. When performing the global fit, also
        # save the background level as an output file
        peakbagging_parameters = {'subdir':     self.cp.isla_subdir,
                                  'run':        self.cp.global_subdir,
                                  'background': self.bgp['name'], 
                                  'fwhm':       linewidth,
                                  'duplet':     False}

        flag_computation_completed = diamonds.run_peakbagging(self.catalog_id, self.star_id, peakbagging_parameters,0,0,1, self.cp.dp_isla, self.cp.diamonds_path, self.cp.n_threads, self.cp.prior_filename)

        # Set some attributes to keep around
        self.freq = freq
        self.psd = psd
        self.numax = numax
        self.scaling_dnu = dnu
        self.linewidth_numax = linewidth

        # pickle to save stuff
        if self.cp.save_progress_pickle:
            pickle.dump(self,open(self.star_dir/(self.catalog_id+self.star_id+'.pickle'),'wb'))

    def find_islands(self, force=False):
        """
        Compute initial radial and dipole mode frequencies and uncertainties.

        This routine processes the global multi-modal sampling obtained with 
        DIAMONDS in order to identify meaningful radial orders and a first 
        proxy for the location of dipole and radial mode frequencies. This is 
        also used to obtain a rather accurate estimate of `dnu` from ACF of 
        the multi-modal sampling obtained with DIAMONDS. For evolved stars, an
        additional fit to locate the central radial mode is performed, here also
        performed with DIAMONDS to gain computational speed. This is because the
        epsilon term from the asymptotic relation may change depending on the 
        evolutionary stage of the star.

        """
        
        run = self.cp.global_subdir
        peakbagging_filename_global = self.star_dir/self.cp.summary_subdir/(self.catalog_id + self.star_id + self.cp.peakbagging_filename_label + self.cp.isla_subdir + '_' + self.cp.global_subdir + '_GLOBAL.txt')

        # Read sampled frequency from DIAMONDS multi-modal fit
        par0 = np.loadtxt(self.star_dir/self.cp.isla_subdir/run/'peakbagging_parameter000.txt')

        self.nested_iters = par0
        self.run = run
        
        # Read prior height (upper limit)
        prior_down, prior_up = np.loadtxt(self.star_dir/self.cp.isla_subdir/run/'peakbagging_hyperParametersUniform.txt', unpack=True)
        left_bound = prior_down[0]
        right_bound = prior_up[0]
        upper_height = prior_up[1]

        # Read input linewidth and background name of the run
        config = np.loadtxt(self.star_dir/self.cp.isla_subdir/run/'peakbagging_computationParameters.txt',dtype='str')    
        fit_linewidth = float(config[-1])
        background_name = config[-2]

        # Read star PSD
        freq, psd = self.freq,self.psd
        freqbin = freq[1]-freq[0]

        # Read Nyquist frequency
        nyq = np.loadtxt(self.star_dir/'NyquistFrequency.txt')

        if self.cp.print_on_screen:
            print('---------------------------------------------------')
            print(' Parameter range (microHz): [',str(np.min(par0)),', ', str(np.max(par0)),']')
            print('---------------------------------------------------\n')
            print('-------------------------------------------------')
            print(' Minimum PSD frequency: ',str(np.min(freq)),' muHz')
            print(' Maximum PSD frequency: ',str(np.max(freq)),' muHz')
            print('-------------------------------------------------\n')


        # Load the nuMax information
        numax = self.numax
        scaling_dnu = self.scaling_dnu
        teff = self.teff
        

        # Set a binwidth proportional to the expected separation between
        # adjacent peaks (normally taken as d02/2) if MS or smaller if RG star.
        par_range = np.max(par0) - np.min(par0)
        min_separation = self.get_minimum_freq_separation(scaling_dnu)
        binwidth = 1.0*min_separation/self.cp.min_n_bins
        n_bins = int(np.ceil(par_range/binwidth))
        nest_iter = np.arange(len(par0))


        
        # Compute an Average Shifted Histogram (ASH) of the distribution of
        # nested iterations. Compute the ASEF with high resolution first in
        # order to determine dnu with high accuracy
        n_bins_acf = round(len(nest_iter)/self.cp.n_bins_acf_scaling)
        par_hist, asef_hist = self.compute_asef(par0,nest_iter,n_bins_acf)

        acf_dnu, interpolated_dnu, interpolated_acf = self.compute_acf_dnu(scaling_dnu,par_hist,asef_hist)

        # Save some attributes
        self.acf_dnu = acf_dnu
        self.interpolated_dnu = interpolated_dnu
        self.interpolated_acf = interpolated_acf
        
        # Compute the ASEF with input resolution for extracting the local maxima
        par_hist, asef_hist = self.compute_asef(par0,nest_iter,n_bins)

        # Save some attributes
        self.par_hist = par_hist
        self.asef_hist = asef_hist
        self.asef_bins = n_bins
        
        # Find an appropriate epsilon based on the input temperature.
        interp_epsi, dnu_array, epsi_array, teff_array, epsi_fit_array  = astero.interpolate_epsilon(teff, acf_dnu, self.cp.dnu_threshold, self.cp.epsilon_offset, self.cp.epsilon_slope)

        # Save stuff for plotting later
        self.epsi_fit = interp_epsi
        self.epsi_fit_dnu_arr = dnu_array
        self.epsi_fit_epsi_arr = epsi_array
        self.epsi_fit_teff_arr = teff_array
        self.epsi_fit_arr = epsi_fit_array
        
        # Find the local maxima using a hill-climbing algorithm on the ASEF.
        # Give as input the minimum number of adjacent bins required to consider
        # two local maxima separated.
        threshold = self.cp.threshold_asef_global*np.max(asef_hist)
        index_maximum = self.hill_climbing(par_hist,asef_hist,threshold,self.cp.min_bin_separation)
        maximum = par_hist[index_maximum]
        asef_maximum = asef_hist[index_maximum]

        # Identify the divisions among the local maxima. These will be used to
        # select chunks of parameter range that locate an optimal frequency
        # region to perform a fit around the peak. Also compute ranges for each
        # maximum such that an ASEF peak is considered up to the points of its
        # tails, i.e. the last ones in the decreasing phase of the peak. Ranges
        # are used to compute the actual frequency estimatse from the sampling.
        n_maxima = len(maximum)
        range_maximum, divisions_maximum = self.get_range_divisions(par_hist, asef_hist, index_maximum)

        self.maxima = maximum
        self.ranges = range_maximum
        self.divisions = divisions_maximum

        # Based on the identified frequency ranges, estimate the oscillation
        # frequencies using the sampling information obtained by DIAMONDS.
        n_freq = n_maxima

        # Compute a smoothed PSD by some average FWHM to evaluate its values at
        # each local maxima of the ASEF. In case of RG, adopt a finer smoothing
        # window to overcome the problem of very narrow mixed modes.
        avg_fwhm = np.mean(astero.get_linewidth(maximum,teff,numax,self.cp.numax_threshold))
        if acf_dnu <= self.cp.dnu_tip:
            avg_fwhm = fit_linewidth*self.cp.smoothing_fwhm_factor_rg

        smth_bins = int(avg_fwhm/freqbin)
        spsd = smooth(psd,window_len=smth_bins,window='flat')
        tmp_good = np.where((freq <= np.max(par_hist)) & (freq >= np.min(par_hist)))[0]
        psd = psd[tmp_good]
        freq = freq[tmp_good]
        spsd = spsd[tmp_good]
        
        # Load the background level of the star
        bg_level = np.loadtxt(self.star_dir/'backgroundLevel.txt', usecols=(1,))
        bg_level_local = bg_level[tmp_good]

        # Evaluate SNR of the dataset in the case of global modality
        snr = np.max(spsd)/np.mean(bg_level_local)

        sampled_estimates = self.evaluate_sampling_frequencies(par0, par_hist ,freq, spsd, maximum, range_maximum)
        freq1 = sampled_estimates['freq1']
        freq_sig1 = sampled_estimates['freq_sig1']
        sampling_counts = sampled_estimates['sampling_counts']
        spsd_maximum = sampled_estimates['spsd_maximum']

        # Define some weights useful for identification of proper frequency peaks during mode identification
        sampling_weights = sampling_counts/np.sum(sampling_counts)
        spsd_weights = spsd_maximum/np.sum(spsd_maximum)
        asef_weights = asef_maximum/np.sum(asef_maximum)

        # Save some attributes to our object
        self.snr = snr
        self.psd = psd
        self.freq = freq
        self.spsd = spsd
        self.bg_level = bg_level_local
        self.snr = snr
        self.freq1 = freq1
        self.freq_sig1 = freq_sig1
        
        # Save the total list of frequencies, uncertainties, ASEF maxima and
        # sampling counts for reference.
        if self.cp.save_complete_lists:
            with open(self.star_dir/self.cp.summary_subdir/(self.catalog_id + self.star_id + self.cp.peakbagging_filename_label + 'global.all.txt'),'w') as f:
                f.write('# Frequency (microHz), 1-sigma frequency (microHz), ASEF maximum (iterations), sampling counts\n')
            for i in range(0,n_freq):
                f.write('%.5f\t%.5f\t%i\t%i\n'%(freq1[i],freq_sig1[i],asef_maximum[i],sampling_counts[i]))

        if self.cp.print_on_screen:
            print(' Total number of local maxima found: ',n_maxima)

        # Sliding Pattern Analysis
        
        # Obtain epsilon from an échelle value of the l=0 ridge échelle
        # frequency position.
        self.bad_epsi = False
        fit_dnu = acf_dnu
        central_freq = numax

        run_subdir = 'sliding'
        run_names = [run_subdir+str(x) for x in range(self.cp.n_sliding_test)]

        flag_evolved_star = False
        flag_depressed_dipole = False

        ap = astero.get_asymptotic_parameters(numax, acf_dnu, teff, self.cp.d01_mass_offset, self.cp.d01_mass_slope, self.cp.d01_offset, self.cp.d02_mass_offset, self.cp.d02_mass_slope, self.cp.d02_offset, self.cp.d03_slope,self.cp.d03_offset, self.cp.numax_sun, self.cp.dnu_sun, self.cp.teff_sun)

        # Start by veryfying whether the stars has depressed dipole modes
        asef_depressed_modes_threshold = (self.cp.dp_isla['max_nested_it']+self.cp.dp_isla['n_live'])*self.cp.asef_threshold_fraction

        if fit_dnu <= self.cp.dnu_threshold:
            central_indicies = np.where((freq1 > central_freq - fit_dnu*self.cp.n_central_orders_side/2) & (freq1 < central_freq + fit_dnu*self.cp.n_central_orders_side/2))[0]
        else:
            central_indicies = np.where((freq1 > central_freq - fit_dnu*self.cp.n_central_orders_side) & (freq1 < central_freq + fit_dnu*self.cp.n_central_orders_side))[0]

        depressed_indicies = np.where(asef_maximum[central_indicies]<asef_depressed_modes_threshold)[0]
        n_central_freq = len(central_indicies)
        n_depressed_freq = len(depressed_indicies)

        if n_depressed_freq >= n_central_freq*self.cp.depressed_dipole_fraction:
            flag_depressed_dipole = True
            if self.cp.print_on_screen:
                print('\nThe star is likely to have depressed dipole modes.\n')

        # Set up priors for the sliding pattern fit
        if fit_dnu <= self.cp.dnu_threshold:
            flag_evolved_star = True
            dipole_radial_height_ratio = self.cp.dipole_radial_height_ratio_rg
            quadrupole_radial_height_ratio = self.cp.quadrupole_radial_height_ratio_rg
            dipole_radial_fwhm_ratio = self.cp.dipole_radial_fwhm_ratio_rg
            n_orders_side_prior = self.cp.n_orders_side_prior_rg
            n_orders_side_data = n_orders_side_prior

            # Set the range of the dataset and frequency prior for the sliding
            # pattern fit.
            data_freq_boundaries = [numax - n_orders_side_data*fit_dnu, numax + n_orders_side_data*fit_dnu]
            tmp_central = np.where((freq <= numax + n_orders_side_prior*fit_dnu) & (freq >= numax - n_orders_side_prior*fit_dnu))[0]
            spsd_central = spsd[tmp_central]
            freq_prior = [numax - n_orders_side_prior*fit_dnu, numax + n_orders_side_prior*fit_dnu]

            dnu_prior = [fit_dnu*self.cp.dnu_prior_lower_fraction,fit_dnu*self.cp.dnu_prior_upper_fraction]
            d02_prior = [ap['d02'],ap['d02']]

            if fit_dnu < self.cp.dnu_tip:
                d01_prior = [ap['d01'],ap['d01']]
            else:
                d01_prior = [99.0,99.0]

            d13_prior = [99.0,99.0]
            rot_split_prior = [0.0,0.0]
            cosi_prior = [0.0,0.0]

            n_orders_model = round(2*n_orders_side_prior)
            if (n_orders_model%2) == 0:
                n_orders_model += 1
        else:
            if flag_depressed_dipole:
                flag_evolved_star = True
                dipole_radial_height_ratio = 0
                dipole_radial_fwhm_ratio = 1

                dnu_prior = [fit_dnu*self.cp.dnu_prior_lower_fraction,fit_dnu*self.cp.dnu_prior_upper_fraction]
                d02_prior = [1.5,self.cp.d02_prior_upper_sg]
                d01_prior = [99.0,99.0]
                d13_prior = [99.0,99.0]
                rot_split_prior = [0.0,0.0]
                cosi_prior = [0.0,0.0]

                quadrupole_radial_height_ratio = self.cp.quadrupole_radial_height_ratio_sg
                n_orders_side_prior = self.cp.n_orders_side_prior_sg
                n_orders_side_data = n_orders_side_prior*1.33333

                #Set the range of the dataset and frequency prior forthe sliding pattern fit
                right_freq_bound = central_freq + n_orders_side_data*fit_dnu
                left_freq_bound = central_freq - n_orders_side_data*fit_dnu
                data_freq_boundaries = [left_freq_bound,right_freq_bound]
                tmp_central = np.where((freq <= right_freq_bound) & (freq >= left_freq_bound))[0]
                spsd_central = spsd[tmp_central]
                freq_prior = [left_freq_bound,right_freq_bound]

                n_orders_model = round(4*n_orders_side_prior)
                if (n_orders_model%2) == 0:
                    n_orders_model += 1
            else:
                # Evaluate the maximum spread in % of DeltaNu of each ridge
                # (both l=2,0 and l=1) with respect to their corresponding
                # median value. If such deviation for any of the two ridges is
                # larger than a given threshold, then consider that mixed modes
                # are present and classify the star as a subgiant.
                central_indices = np.where((freq1 > central_freq - fit_dnu*self.cp.n_central_orders_side) & (freq1 < central_freq + fit_dnu*self.cp.n_central_orders_side) & (asef_maximum >= asef_depressed_modes_threshold))[0]
                freq1_right = freq1[central_indices[1::2]]     # Odd frequencies
                freq1_left = freq1[central_indices[0::2]]      # Even frequencies

                # Check that extracted frequencies are approximately varying in
                # steps of dnu
                flag_double_step = False

                if len(freq1_right) >= 2:
                    freq1_right_diff = np.diff(freq1_right)
                    if median_high(freq1_right_diff) < fit_dnu*self.cp.dnu_ridge_threshold:
                        flag_double_step = True

                if len(freq1_left) >= 2:
                    freq1_left_diff = np.diff(freq1_left)
                    if median_high(freq1_left_diff) < fit_dnu*self.cp.dnu_ridge_threshold:
                        flag_double_step = True
                                                
                # If a double step is required then consider three ridges
                # (left, central, and right).
                if flag_double_step:
                    freq1_left = freq1[central_indicies[::3]]
                    freq1_right = freq1[central_indicies[2::3]]

                    if self.cp.print_on_screen:
                        print('\nUsing double step to check frequency ridge')

                freq1_left_modulo = freq1_left%fit_dnu
                freq1_right_modulo = freq1_right%fit_dnu

                tmp = np.where(freq1_left_modulo > fit_dnu/2.)[0]
                if len(tmp) > 0: 
                    freq1_left_modulo[tmp] = fit_dnu - freq1_left_modulo[tmp]

                tmp = np.where(freq1_right_modulo > fit_dnu/2.)[0]
                if len(tmp) > 0:
                    freq1_right_modulo[tmp] = fit_dnu - freq1_right_modulo[tmp]

                freq1_left_modulo_median = median_high(freq1_left_modulo)
                freq1_right_modulo_median = median_high(freq1_right_modulo)
                left_modulo_median_index = np.where(freq1_left_modulo == freq1_left_modulo_median)[0]
                right_modulo_median_index = np.where(freq1_right_modulo == freq1_right_modulo_median)[0]
                freq1_left_median = freq1_left[left_modulo_median_index]
                freq1_right_median = freq1_right[right_modulo_median_index]
                deviation_left = np.abs(freq1_left_modulo - freq1_left_modulo_median) / fit_dnu*100.
                deviation_right = np.abs(freq1_right_modulo - freq1_right_modulo_median) / fit_dnu*100.
                max_dev_left = np.max(deviation_left)
                max_dev_right = np.max(deviation_right)
                max_dev = np.max([max_dev_left,max_dev_right])

                if flag_double_step:
                    freq1_central = freq1[central_indicies[1::3]]
                    freq1_central_modulo = freq1_central%fit_dnu

                    tmp = np.where(freq1_central_modulo > fit_dnu/2)[0]
                    if len(tmp) > 0:
                        freq1_central_modulo[tmp] = fit_dnu - freq1_central_modulo[tmp]

                    freq1_central_modulo_median = median_high(freq1_central_modulo)
                    central_modulo_median_index = np.where(freq1_central_modulo == freq1_central_modulo_median)[0]
                    freq1_central_median = freq1_central[central_modulo_median_index]


                    deviation_central = np.abs(freq1_central_modulo - freq1_central_modulo_median) / fit_dnu *100
                    max_dev_central = np.max(deviation_central)
                    max_dev = max(max_dev_left,max_dev_central,max_dev_right)
                else:
                    max_dev = max(max_dev_left,max_dev_right)


                if max_dev >= self.cp.dnu_echelle_threshold:
                    flag_evolved_star = True

                if flag_evolved_star:
                    if self.cp.print_on_screen:
                        print('\n The star likely contains modes that have undergone avoided crossings, so it is classified as a subgiant.\n')

                    # Now find the ridge with the smallest maximum spread and
                    # from its modulo median value find epsilon and the value of
                    # radial mode frequency.
                    if flag_double_step:
                        freq1_median = np.r_[freq1_left_median,freq1_central_median,freq1_right_median]
                        freq1_modulo_median = np.r_[freq1_left_modulo_median,feq1_central_modulo_median,freq1_right_modulo_median]
                        ridge_index = np.argmin(np.r_[max_dev_left,max_dev_central,max_dev_right])
                        min_max_dev = min(max_dev_left,max_dev_central,max_dev_right)
                    else:
                        freq1_median = np.r_[freq1_left_median,freq1_right_median]
                        freq1_modulo_median = np.r_[freq1_left_modulo_median,freq1_right_modulo_median]
                        ridge_index = np.argmin(np.r_[max_dev_left,max_dev_right])
                        min_max_dev = min(max_dev_left,max_dev_right)

                    freq_radial_echelle = freq1_median[ridge_index]
                    dipole_radial_height_ratio = self.cp.dipole_radial_height_ratio_sg
                    quadrupole_radial_height_ratio = self.cp.quadrupole_radial_height_ratio_sg
                    dipole_radial_fwhm_ratio = self.cp.dipole_radial_fwhm_ratio_sg
                    n_orders_side_prior = self.cp.n_orders_side_prior_sg
                    n_orders_side_data = n_orders_side_prior*1.33333

                    # Set the range of the dataset and frequency prior for the
                    # sliding pattern fit
                    right_freq_bound = freq_radial_echelle + n_orders_side_data*fit_dnu
                    left_freq_bound = freq_radial_echelle - n_orders_side_data*fit_dnu
                    data_freq_boundaries = [left_freq_bound,right_freq_bound]
                    tmp_central = np.where((freq <= right_freq_bound) & (freq >= left_freq_bound))[0]
                    spsd_central = spsd[tmp_central]
                    freq_prior = [left_freq_bound, right_freq_bound]

                    dnu_prior = [fit_dnu,fit_dnu]
                    d02_prior = [self.cp.d02_prior_lower_sg,self.cp.d02_prior_upper_sg]
                    d01_prior = [99.0,99.0]
                    d13_prior = [99.0,99.0]
                    rot_split_prior = [0.0,0.0]
                    cosi_prior = [0.0,0.0]
                else:
                    if self.cp.print_on_screen:
                        print('\n The star could be a main sequence. \n')

                    dipole_radial_height_ratio = self.cp.dipole_radial_height_ratio_ms
                    quadrupole_radial_height_ratio = self.cp.quadrupole_radial_height_ratio_ms
                    dipole_radial_fwhm_ratio = self.cp.dipole_radial_fwhm_ratio_ms
                    n_orders_side_prior = self.cp.n_orders_side_prior_ms
                    n_orders_side_data = n_orders_side_prior

                    # Set the range of the dataset and frequency prior for the
                    # sliding pattern fit
                    data_freq_boundaries = [numax - n_orders_side_data*fit_dnu, numax + n_orders_side_data*fit_dnu]
                    tmp_central = np.where((freq <= numax + n_orders_side_prior*fit_dnu) & (freq >= numax - n_orders_side_prior*fit_dnu))[0]
                    spsd_central = spsd[tmp_central]
                    freq_prior = [numax - n_orders_side_prior*fit_dnu, numax + n_orders_side_prior*fit_dnu]

                    dnu_prior = [fit_dnu*self.cp.dnu_prior_lower_fraction,fit_dnu*self.cp.dnu_prior_upper_fraction]

                    d02_prior_upper_ms = self.cp.d02_prior_upper_ms
                    if d02_prior_upper_ms >= fit_dnu/4:
                        d02_prior_upper_ms = fit_dnu/4
                    d02_prior = [self.cp.d02_prior_lower_ms,d02_prior_upper_ms]

                    d01_prior_upper = self.cp.d01_prior_upper_ms
                    if d01_prior_upper >= fit_dnu/4:
                        d01_prior_upper = fit_dnu/4
                    d01_prior = [self.cp.d01_prior_lower_ms,d01_prior_upper]


                    d13_prior_upper = d01_prior_upper
                    if d13_prior_upper < (3./8.*fit_dnu):
                        d13_prior_upper = fit_dnu/4. - d01_prior_upper
                    d13_prior = [self.cp.d13_prior_lower_ms,abs(d13_prior_upper)]

                    rot_split_prior = [self.cp.rot_split_prior_lower_ms,self.cp.rot_split_prior_upper_ms]
                    cosi_prior = [self.cp.cosi_prior_lower,self.cp.cosi_prior_upper]

                n_orders_model = round(4*n_orders_side_prior)
                if (n_orders_model%2) == 0:
                    n_orders_model += 1

        flag_repeat_sliding_fit = True
        sliding_iteration = 0
        
        while flag_repeat_sliding_fit & (sliding_iteration <= 1):

            # Check if run already exists
            if not os.path.isfile(self.star_dir/self.cp.as_subdir/(run_subdir+'0')/'peakbagging_parameter000.txt') or force or (sliding_iteration > 0):
                # Make sure that the number of orders to compute the sliding
                # pattern model is an odd number. This will make sure that the
                # selected range is symmetric with respect to nuMax

                print(' Total number of radial orders in the sliding model: %i'%n_orders_model)

                asymp_param = [n_orders_model, dipole_radial_height_ratio, quadrupole_radial_height_ratio, 
                               self.cp.octupole_radial_height_ratio, dipole_radial_fwhm_ratio]

                fwhm_asymp = astero.get_linewidth(numax,teff,numax,self.cp.numax_threshold)
                asymp_filename = 'asymptoticParameters'
                filename = self.star_dir/(asymp_filename + '.txt')

                header=("""
# Asymptotic parameters to set up the asymptotic pattern model.
# Row 1: Norders (spanning orders range around nuMax for asymptotic pattern)
# Row 2: l=1/l=0 height
# Row 3: l=2/l=0 height
# Row 3: l=3/l=0 height
# Row 4: l=1/l=0 fwhm
""")
                np.savetxt(filename,np.array(asymp_param).T,fmt='%.3f',header=header)
               
                height_prior = [np.max(spsd_central)*0.1, np.max(spsd_central)*1.4]
                boundaries = [freq_prior, height_prior, dnu_prior, d02_prior, d01_prior, d13_prior, rot_split_prior, cosi_prior]

                self.bb = boundaries
                # Prepare frequency range and prior filenames for each run
                prior_filenames = []
                data_range_filenames = []

                for k in range(0, self.cp.n_sliding_test):
                    data_range_filenames.append(self.star_dir/self.cp.as_subdir/('frequencyRange_'+run_names[k] + '.txt'))
                    prior_filenames.append(self.star_dir/self.cp.as_subdir/(self.cp.prior_filename+'_'+run_names[k] + '.txt'))
                    np.savetxt(data_range_filenames[k],np.array(data_freq_boundaries),fmt='%.5f')
                    diamonds.write_uniform_prior(prior_filenames[k],np.array(boundaries))

                if self.cp.print_on_screen:
                    print('\n Performing asymptotic pattern fit with DIAMONDS. \n')

                peakbagging_parameters = {'subdir':          self.cp.as_subdir,
                                          'run':             run_names,
                                          'background':      background_name,
                                          'fwhm':            fwhm_asymp,
                                          'filename_run':    run_subdir,
                                          'duplet':          False}

                flag_computation_completed = diamonds.run_peakbagging(self.catalog_id, self.star_id, peakbagging_parameters, 0, 1, 0, self.cp.dp_slid, self.cp.diamonds_path, self.cp.n_threads, self.cp.prior_filename)
            else:
                if self.cp.print_on_screen:
                    print(' Load information from asymptotic pattern fit with DIAMONDS.')

            echelle_epsi_array = np.zeros(self.cp.n_sliding_test)
            radial_freq_reference_array = np.zeros(self.cp.n_sliding_test)
            d01_array = np.zeros(self.cp.n_sliding_test)

            for k in range(0, self.cp.n_sliding_test):
                # Read sampled frequency from multi-modal fit for nu0 central
                par_nu0 = np.loadtxt(self.star_dir/self.cp.as_subdir/(run_names[k]+'/peakbagging_parameter000.txt'))

                # Read sampled posterior distribution from multi-modal fit
                post = np.loadtxt(self.star_dir/self.cp.as_subdir/(run_names[k]+'/peakbagging_posteriorDistribution.txt'))
                post /= np.max(post)
                radial_freq_reference = np.sum(par_nu0*post)
                if self.cp.input_radial_freq_reference > 0:
                    radial_freq_reference = self.cp.input_radial_freq_reference

                radial_freq_reference_array[k] = radial_freq_reference
                modulo_reference = radial_freq_reference%fit_dnu
                echelle_epsi = modulo_reference/fit_dnu

                if (echelle_epsi < self.cp.epsilon_threshold) & (fit_dnu > self.cp.dnu_lower_threshold_epsilon):
                    echelle_epsi += 1
                echelle_epsi_array[k] = echelle_epsi

                # if star is flagged as MS, then also retrieve the value of
                # the small spacing d01
                if not flag_evolved_star:
                    par_d01 = np.loadtxt(self.star_dir/self.cp.as_subdir/(run_names[k]+'/peakbagging_parameter004.txt'))
                    d01_array[k] = np.sum(par_d01*post)/np.sum(post)

            median_echelle_epsi = median_high(echelle_epsi_array)

            if fit_dnu <= self.cp.dnu_threshold:
                epsilon_upper_limit  = self.cp.upper_epsilon_rg_slope * np.log(fit_dnu) + self.cp.upper_epsilon_rg_offset
                epsilon_lower_limit  = self.cp.lower_epsilon_rg_slope * np.log(fit_dnu) + self.cp.lower_epsilon_rg_offset

                # The sliding pattern without including the l=1 mode peak has
                # failed in providing a reliable epsilon. Therefore repeat the
                # fit by including the l=1 model peak, having the position fixed
                # to the p-mode frequency of the asymptotic relation
                if (median_echelle_epsi >= epsilon_upper_limit) or ((median_echelle_epsi <= epsilon_lower_limit) & (fit_dnu >= self.cp.dnu_cl2)):
                    if self.cp.print_on_screen:
                        if sliding_iteration>0:
                            print('Repeating the sliding-pattern fit did not solve the issue. Epsilon is likely to be wrong for this star')
                            self.bad_epsi = True
                        else:
                            if median_echelle_epsi >= epsilon_upper_limit: 
                                print('Repeating the sliding-pattern fit because epsilon exceeds the upper limit for RGs')
                            else:
                                print('Repeating the sliding-pattern fit because epsilon exceeds the lower limit for RGs')
                    d01_prior = [ap['d01'],ap['d01']]
                else:
                    # Epsilon from the sliding-pattern fit has been validated.
                    # Therefore exit the while loop.
                    flag_repeat_sliding_fit = False
            else:
                # No need to repeat the sliding pattern fit because the star is
                # not a red giant.
                flag_repeat_sliding_fit = False
            sliding_iteration +=1

        median_index = np.where(echelle_epsi_array == median_echelle_epsi)[0]
        radial_freq_reference = radial_freq_reference_array[median_index]
        
        if not flag_evolved_star:
            fit_d01 = median_high(d01_array)
        else:
            if d01_prior[0] != 99.0:
                fit_d01 = ap['d01']
            else:
                fit_d01 = 0

        # Control check on epsilon for late subgiants/early RGB
        interp_epsi_flag = False

        if self.cp.force_epsilon_dnu_value:
            if (fit_dnu <= self.cp.dnu_threshold) & (fit_dnu >= self.cp.dnu_cl):
                # If the star is an evolved subgiant/early RGB check that
                # epsilon is in agreement with the epsilon-dnu relation
                radial_freq_reference2 = closest(radial_freq_reference,freq1)
                modeid_sliding = astero.get_modeid(radial_freq_reference,fit_dnu,median_echelle_epsi,0,numax,0)
                modeid_epsi_dnu = astero.get_modeid(radial_freq_reference2,fit_dnu,interp_epsi,0,numax,0)
                diff_radial = abs(radial_freq_reference2 - radial_freq_reference)

                if (modeid_sliding['degree'] != modeid_epsi_dnu['degree']) or (diff_radial >= fit_dnu/4):
                    median_echelle_epsi = interp_epsi
                    interp_epsi_flag = True
                    self.bad_epsi = True
                    if self.cp.print_on_screen:
                        print('Sliding mode identification did not match epsilon-dnu relation.\n')
                    
                    if modeid_epsi_dnu['degree'] == 1:
                        radial_freq_reference = closest(radial_freq_reference2+fit_dnu/2,freq1)
                        if self.cp.print_on_screen:
                            print('Applying correction to epsilon from epsilon-Dnu relation.\n')
                    else:
                        radial_freq_reference = radial_freq_reference2

        # Control check on epsilon for F-type stars
        if teff > self.cp.teff_sg:
            # Due to the confusion arising from the strong blending of the modes
            # for F-type stars, the sliding pattern fit may be unreliable.
            # In this case check the obtained mode identification against the
            # one using the epsilon-Teff relation 
            radial_freq_reference2 = closest(radial_freq_reference,freq1)
            modeid_sliding = astero.get_modeid(radial_freq_reference,fit_dnu,median_echelle_epsi,fit_d01,numax,0)
            modeid_teff = astero.get_modeid(radial_freq_reference2,fit_dnu,interp_epsi,fit_d01,numax,0)
            diff_radial = abs(radial_freq_reference2 - radial_freq_reference)

            if (modeid_sliding['degree'] != modeid_teff['degree']) or (diff_radial >= fit_dnu/4):
                    median_echelle_epsi = interp_epsi
                    interp_epsi_flag = True
                    self.bad_epsi = True
                    if self.cp.print_on_screen:
                        print('Sliding mode identification did not match epsilon-Teff relation well.\n')
                    
                    if modeid_teff['degree'] == 1:
                        radial_freq_reference = closest(radial_freq_reference2+fit_dnu/2,freq1)
                        if self.cp.print_on_screen:
                            print('Applying correction to epsilon from epsilon-Teff relation well.\n')
                    else:
                        radial_freq_reference = radial_freq_reference2

        if self.cp.print_on_screen:
            print('\n Epsilon values from the sliding fit: {}'.format(echelle_epsi_array))
            print(' Final epsilon: {}'.format(median_echelle_epsi))
            print(' The reference radial mode is at: {} muHz\n'.format(radial_freq_reference))

        fit_epsi = median_echelle_epsi 
        fit_alpha = 0.        
        flag_dnu_fit = True
        iterations = 0

        # Find possible bad frequencies and remove them from the list. Then
        # compute a large frequency separation. Finally obtain a simple mode
        # identification (either l=0 or l=1) for the full set of frequencies.
        while (flag_dnu_fit) & (iterations < self.cp.max_skim_iterations_global):
            n_freq = n_maxima
            order_number = np.zeros(n_freq,dtype='int')
            angular_degree = np.zeros(n_freq,dtype='int')

            if self.cp.print_on_screen:
                print(' Iteration: %s' % iterations)

            # Perform a first mode identification based on the ACF value of dnu.
            mode_id = astero.get_modeid(freq1,fit_dnu,fit_epsi,fit_d01,numax,fit_alpha)
            enn = mode_id['order']
            ell = mode_id['degree']          
            order_number = mode_id['order']
            angular_degree = mode_id['degree']

            # Using the obtained mode identification, divide the frequency set
            # into dipole and radial modes.
            tmp_radial = np.where(angular_degree == 0)[0]
            tmp_dipole = np.where(angular_degree == 1)[0]

            # Select only l=0 modes
            if len(tmp_radial) > 0:
                n_radial = len(tmp_radial)
                freq1_radial = freq1[tmp_radial]
                freq1_radial_org = freq1[tmp_radial]
                freq_sig1_radial = freq_sig1[tmp_radial]
                order_radial = order_number[tmp_radial]
                asef_maximum_radial = asef_maximum[tmp_radial]
            else:
                n_radial = 0

            # Select only l=1 modes
            if len(tmp_dipole) > 0:
                n_dipole = len(tmp_dipole)
                freq1_dipole = freq1[tmp_dipole]
                freq_sig1_dipole = freq_sig1[tmp_dipole]
                order_dipole = order_number[tmp_dipole]
                asef_maximum_dipole = asef_maximum[tmp_dipole]
            else:
                n_dipole = 0

            # If requested, apply a small correction to the radial mode
            # frequencies that takes into account the difference between the
            # reference radial mode from the sliding pattern and the one
            # estimated from the ASEF. If the sliding-pattern fit was
            # unsuccesful, this difference will be 0.
            if self.cp.correct_radial_frequencies:
                delta_nu0 = radial_freq_reference - closest(radial_freq_reference,freq1_radial)
                ap = astero.get_asymptotic_parameters(numax,fit_dnu,teff, self.cp.d01_mass_offset, self.cp.d01_mass_slope, self.cp.d01_offset, self.cp.d02_mass_offset, self.cp.d02_mass_slope, self.cp.d02_offset, self.cp.d03_slope,self.cp.d03_offset, self.cp.numax_sun, self.cp.dnu_sun, self.cp.teff_sun)

                if delta_nu0 < ap['d02']: 
                    freq1_radial = freq1_radial + delta_nu0

            # Verify the position of each frequency by comparing it to the
            # expected asymptotic value. Discard those frequencies that deviate
            # from their asymptotic value by more than a given tolerance in dnu.
            if (fit_dnu < self.cp.dnu_rg) or ((fit_dnu < self.cp.dnu_sg) & (self.teff < self.cp.teff_sg)):
                if fit_dnu > self.cp.dnu_rg:
                    tolerance = self.cp.skim_frequency_tolerance_sg
                else:
                    if fit_dnu <= self.cp.dnu_tip:
                        tolerance = self.cp.skim_frequency_tolerance_tip
                    else:
                        tolerance = self.cp.skim_frequency_tolerance_rg
            else:
                tolerance = self.cp.skim_frequency_tolerance_ms
    
            if n_radial != 0:
                # Select only good l=0 frequencies
                good_freq_index_radial = astero.assess_freq_asymptotic(freq1_radial,order_radial,0,fit_dnu,fit_epsi,fit_alpha,fit_d01,numax,tolerance)

                if self.cp.print_on_screen:
                    if len(good_freq_index_radial) < n_radial:
                        print(' Returning fewer frequencies than input for l=0. Input: {},  Actual: {}'.format(n_radial,len(good_freq_index_radial)))

                freq1_radial = freq1_radial[good_freq_index_radial]
                freq1_radial_org = freq1_radial_org[good_freq_index_radial]
                freq_sig1_radial = freq_sig1_radial[good_freq_index_radial]
                order_radial = order_radial[good_freq_index_radial]
                asef_maximum_radial = asef_maximum_radial[good_freq_index_radial]

            # Collect all frequencies after removing the bad ones
            n_dipole = len(freq1_dipole)
            n_radial = len(freq1_radial)
            n_freq = n_dipole + n_radial

            if (n_dipole != 0) & (n_radial != 0):
                freq1_final = np.append(freq1_radial,freq1_dipole)
                freq_sig1_final =  np.append(freq_sig1_radial,freq_sig1_dipole)
                angular_degree =  np.append(np.zeros(n_radial,dtype='int'),np.ones(n_dipole,dtype='int'))
                order_number =  np.append(order_radial,order_dipole)
                asef_maximum_final =  np.append(asef_maximum_radial,asef_maximum_dipole)
            elif n_dipole == 0 and n_radial != 0:
                freq1_final = freq1_radial
                freq_sig1_final = freq_sig1_radial
                angular_degree = np.zeros(n_radial,dtype='int')
                order_number = order_radial
                asef_maximum_final = asef_maximum_radial 

            sorted_index = np.argsort(freq1_final)
            freq1_final = freq1_final[sorted_index]
            freq_sig1_final = freq_sig1_final[sorted_index]
            angular_degree = angular_degree[sorted_index]
            order_number = order_number[sorted_index]
            asef_maximum_final = asef_maximum_final[sorted_index]

            # Compute optimal value for dnu from individual frequencies only if
            # at least three modes are present. Otherwise, keep as best dnu and
            # epsilon, those from the ACF and dnu-epsilon diagram.
            flag_dnu_fit = False

            if n_radial >= 2:
                run_subdir = 'radial_global'
                filename = self.star_dir/self.cp.as_subdir/(self.cp.prior_filename + '_' + run_subdir + '.txt')
                dnu_prior = [fit_dnu*self.cp.dnu_prior_lower_fraction_as,fit_dnu*self.cp.dnu_prior_upper_fraction_as]

                # Use epsilon fixed from the sliding pattern fit, if successful
                epsi_prior = [fit_epsi,fit_epsi]
                if (n_radial == 2) or (n_radial == 3):
                    # Only DeltaNu can be estimated as a free parameter. Then fix alpha to 0 and epsilon to
                    # its former value.
                    alpha_prior = [0.0,0.0]
                else:
                    alpha_prior = [self.cp.alpha_prior_lower_as,self.cp.alpha_prior_upper_as]
                    if interp_epsi_flag:
                        epsi_prior = [fit_epsi*self.cp.epsi_prior_lower_fraction_as,fit_epsi*self.cp.epsi_prior_upper_fraction_as]

                boundaries = np.array([dnu_prior,epsi_prior,alpha_prior])
                diamonds.write_uniform_prior(filename,boundaries)
                data = np.array([order_radial,freq1_radial,freq_sig1_radial]).T
                data_dir = self.star_dir/self.cp.as_subdir/'data'
                data_filename = data_dir/(run_subdir + '.txt')
                np.savetxt(data_filename, data, fmt=['%.1f','%12.4f','%12.4f'])

                run_parameters = {'subdir':     self.cp.as_subdir,
                                  'run':        run_subdir}
              
                flag_computation_completed = diamonds.run_asymptotic(self.catalog_id, self.star_id, run_parameters, numax, 0, self.cp.dp_pb, self.cp.diamonds_path)

                parameter_filenames = np.sort(glob.glob(str(self.star_dir/self.cp.as_subdir/run_subdir/'asymptotic_parameter0*.txt')))
                n_par = len(parameter_filenames) 
                par = np.loadtxt(parameter_filenames[0])
                post = np.loadtxt(self.star_dir/self.cp.as_subdir/run_subdir/'asymptotic_posteriorDistribution.txt')
                post = post/np.max(post)
                best_dnu = np.sum(par*post)/np.sum(post)

                if n_par == 1:
                    best_epsi = epsi_prior[0]
                    best_alpha = alpha_prior[0]

                if n_par == 2:
                    par = np.loadtxt(parameter_filenames[1])
                    best_epsi = epsi_prior[0]
                    best_alpha = np.sum(par*post)/np.sum(post)
                        
                if n_par == 3:
                    par = np.loadtxt(parameter_filenames[1])
                    best_epsi = np.sum(par*post)/np.sum(post)
                    par = np.loadtxt(parameter_filenames[2])
                    best_alpha = np.sum(par*post)/np.sum(post)
                        
                # Update local values for asymptotic parameters
                fit_dnu = best_dnu
                fit_epsi = best_epsi
                fit_alpha = best_alpha
                flag_dnu_fit = True
                iterations +=1
            else:
                # Not enough radial mode frequencies were found to perform
                # an asymptotic fit
                best_dnu = acf_dnu
                best_epsi = fit_epsi
                best_alpha = self.cp.alpha_radial_universal
                if self.cp.print_on_screen:
                    print(' Not enough radial mode frequencies to perform an asymptotic fit.\n')
                    
        # Save asymptotic radial mode frequencies
        freq_radial_asymptotic = best_dnu*(best_epsi + order_radial + best_alpha/2.*(order_radial - numax/best_dnu)**2)

        # Select only l=1 modes
        tmp_dipole = np.where(angular_degree==1)[0]
        if len(tmp_dipole)>0:
            n_dipole = len(tmp_dipole)
            freq1_dipole = freq1_final[tmp_dipole]
        else:
            n_dipole = 0

        n_freq = len(freq1_final)

        # Calculate radial order positions to divide the PSD into chunks
        # Compute radial mode positions from universal pattern for RG stars
        # (Mosser et al. 2011, A&A, 525, L9). This will improve the trimming of
        # chunks with more accurate frequency positions.
        
        # Check for an additional dipole frequency after the last radial mode.
        # If so, add one more chunk to include the last dipole separately.
        n_chunks = max(order_number) - min(order_number) + 1
        unique_order_number = np.arange(n_chunks) + min(order_number)

        if max(freq1_final) > freq1_radial[n_radial-1]:
            if (max(freq1_final) + best_dnu) <= (max(par_hist) - best_dnu/2):
                n_chunks = n_chunks + 1
                unique_order_number = np.append(unique_order_number,max(unique_order_number)+1)

        if min(freq1_final) < freq1_radial[0]:
            if (min(freq1_final) - best_dnu) <= (min(par_hist) + best_dnu/2):
                n_chunks = n_chunks - 1
                unique_order_number = unique_order_number[1:]

        n_separations = n_chunks + 1
        separations = np.zeros(n_separations)
        freq_asymptotic_separation = best_dnu*(best_epsi + unique_order_number + best_alpha/2.*(unique_order_number - numax/best_dnu)**2)

        for i in range(1, n_separations-1):
            if (best_dnu < self.cp.dnu_sg) & (teff < self.cp.teff_sg):
                separations[i] = freq_asymptotic_separation[i-1] + best_dnu*self.cp.separations_dnu_tolerance_rg
                self.dnu_tol = self.cp.separations_dnu_tolerance_rg
            else:
                separations[i] = freq_asymptotic_separation[i-1] + best_dnu*self.cp.separations_dnu_tolerance_ms
                self.dnu_tol = self.cp.separations_dnu_tolerance_ms

        if separations[1] - best_dnu < min(par_hist):
            separations[0] = min(par_hist)
        else:
            separations[0] = separations[1] - best_dnu

        if separations[n_separations-2] + best_dnu > max(par_hist):
            separations[n_separations-1] = max(par_hist)
        else:
            separations[n_separations-1] = separations[n_separations-2] + best_dnu

        if self.cp.save_asymptotic_radial:
            freq1_final = np.sort(np.r_[freq1_dipole, freq_radial_asymptotic])

        
        # Save final outputs
        # Evaluate empirical linewidth for each oscillation frequency
        linewidth = astero.get_linewidth(freq1_final,teff,numax)

        # Save the value of dnu from ACF and the value of epsilon from diagram
        with open(peakbagging_filename_global, 'w') as f:
            f.write('# nuMax (microHz), DeltaNu_ACF (microHz), DeltaNu_fit (microHz), epsilon, alpha, Teff (K), N_chunks, Flag depressed dipole\n')
            f.write('{:.4f}  {:.4f}  {:.4f}  {:.4f}  {:.4f}  {:.1f}  {:d}  {:d}\n'.format(numax,acf_dnu,best_dnu,best_epsi,best_alpha,teff,n_chunks,flag_depressed_dipole))

            # Save the frequency positions of each chunk identified in GLOBAL.
            f.write('# Chunk index, start and end frequency values for each chunk (one per line), SNR\n')
            for i in range(0, n_separations-1):
                tmp_chunk = np.where((freq <= separations[i+1]) & (freq >= separations[i]))[0]
                spsd_chunk = spsd[tmp_chunk]
                bg_level_chunk = bg_level_local[tmp_chunk]
                snr_chunk = max(spsd_chunk)/np.mean(bg_level_chunk)
                f.write('{:<3d} {:>12.3f} {:>12.3f} {:>12.2f}\n'.format(i,separations[i],separations[i+1],snr_chunk))

            # Save the individual oscillation frequencies from the global fit,
            # their uncertainties, mode identification and associated FWHM from
            # Ball+18.
            f.write('# n, l, frequency (microHz), 1-sigma frequency (microHz), ASEF maximum (iterations), FWHM from predictions (microHz)\n')
            for i in range(0,n_freq):
                f.write('{:<3d} {:>2d} {:>12.5f} {:>12.5f} {:>7d} {:>12.5f}\n'.format(order_number[i],angular_degree[i],freq1_final[i],freq_sig1_final[i],int(asef_maximum_final[i]),linewidth[i]))

        # Save some quantities as object attributes.
        self.separations = separations
        self.n_freq = n_freq
        self.n_chunks = n_chunks
        self.orders = order_number
        self.degrees = angular_degree
        self.freqs = freq1_final
        self.freqs_sig = freq_sig1_final
        self.asef_maxima = asef_maximum_final
        self.linewidths = linewidth

        self.hmax_prior = upper_height
        self.alpha = best_alpha
        self.dnu = best_dnu
        self.epsilon = best_epsi
        self.flag_depressed_dipole = flag_depressed_dipole

        # Save stuff into pickle for later steps...
        if self.cp.save_progress_pickle:
            pickle.dump(self,open(self.star_dir/(self.catalog_id+self.star_id+'.pickle'),'wb'))

    def make_global_plots(self):
        """
        Produce all plots related to the GLOBAL modality and save as desired.
        """
        plt.style.use(self.cp.famed_path/self.cp.mplstyle)
        famed_plots.global_plot(self)

        if self.cp.save_png:
            plt.savefig(self.star_dir/self.cp.figs_subdir/(self.catalog_id+self.star_id+'_'+self.cp.isla_subdir+'_'+self.cp.global_subdir+'_GLOBAL.png'))
        if self.cp.save_eps:
            plt.savefig(self.star_dir/self.cp.figs_subdir/(self.catalog_id+self.star_id+'_'+self.cp.isla_subdir+'_'+self.cp.global_subdir+'_GLOBAL.eps'))

