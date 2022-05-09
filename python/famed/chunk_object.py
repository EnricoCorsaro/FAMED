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

__all__ = ['Chunk']

class Chunk(FamedStar):
    """
    FamedStar sub-class specific to the running of the CHUNK modality.    

    Parameters
    ----------
    catalog_id : str
        Catalogue ID of the star (e.g. 'KIC' for Kepler).
    star_id : str
        ID of the star as a string (e.g. '0012008916' or '7037405').
    teff : float
        Effective temperature of the star in Kelvin.
    load_islands : bool, default: False
        Flag to read the pickled object if `make_islands` or `find_islands` has 
        already been ran.
    """
    def __init__(self, catalog_id, star_id, background_run_number=None, load_islands=False):
        FamedStar.__init__(self, catalog_id, star_id, background_run_number=background_run_number)
        
        if load_islands and self.cp.save_progress_pickle:
            if self.cp.print_on_screen:
                print('Loading saved information from the pickled object.')
            # Load any variables into attributes that should haven been
            # loaded or created during make_islands. Will ignore new config
            temp = pickle.load(open(self.star_dir/(self.catalog_id+self.star_id+'_chunk.pickle'),'rb'))
            self.__dict__ = temp.__dict__.copy()
        else:
            self.modality='CHUNK'


    def make_islands(self,run,force=False,merge=False):
        """
        Compute a chunk multi-modal fit with DIAMONDS.

        This function computes a chunk multi-modal fit in order to retrieve the 
        frequencies, their uncertainties, and mode identification. This function
        cannot be executed if the GLOBAL modality was not previously performed.
        """
        peakbagging_filename_global = self.star_dir/self.cp.summary_subdir/(self.catalog_id + self.star_id + self.cp.peakbagging_filename_label + self.cp.isla_subdir + '_' + self.cp.global_subdir + '_GLOBAL.txt')
        
        # Read input PSD and global asteroseismic parameters
        peakbagging_data_dir = self.cp.diamonds_path/'PeakBagging'/'data'
        freq, psd = np.loadtxt(peakbagging_data_dir/(self.catalog_id + self.star_id + '.txt'), unpack=True)
        gauss_par = np.loadtxt(self.star_dir/'gaussianEnvelopeParameters.txt')
        numax = gauss_par[1]
        freqbin = freq[1]-freq[0]
        dnu = astero.compute_scaling_dnu(numax, self.cp.numax_threshold, self.cp.numax_coeff_low, self.cp.numax_coeff_high, self.cp.numax_exponent_low, self.cp.numax_exponent_high)
        maxpower = max(psd)
        maxpower *= 1.1

        # Run DIAMONDS to perform the multi-modal fit in the different chunks indicated by the global fit.
        if not os.path.isfile(peakbagging_filename_global):
            print('Cannot perform multi-modal fit for individual chunks. Missing global fit peakbagging summary file.')
            return False
        else:
            # Load global parameters
            acf_dnu,best_dnu,best_epsi,best_alpha,teff,n_chunks = np.loadtxt(peakbagging_filename_global,max_rows=1,usecols=(1,2,3,4,5,6))
            n_chunks = int(n_chunks)

            self.n_chunks = n_chunks
            self.acf_dnu = acf_dnu
            self.best_dnu = best_dnu
            self.best_epsi = best_epsi
            self.best_alpha = best_alpha
            self.teff = teff
            
            # Load the frequency positions for each chunk
            chunk_number,freq_left,freq_right,snr = np.loadtxt(peakbagging_filename_global,max_rows=n_chunks,skiprows=3,unpack=True)
            chunk_number = np.array(chunk_number,dtype=int)
            min_linewidths = np.zeros(n_chunks)
            bg_names = np.zeros(n_chunks,dtype='U30')
            run_labels = np.arange(n_chunks)

            # Load the background level of the star
            bg_level = np.loadtxt(self.star_dir/'backgroundLevel.txt', usecols=(1,))
            
            self.freq = [None]*n_chunks
            self.psd = [None]*n_chunks
            self.spsd = [None]*n_chunks
            self.bg_level = [None]*n_chunks

            if (run == -1) or (run >= n_chunks):
                first_it = 0
                last_it = n_chunks
            else:
                first_it = run
                last_it = run + 1

            for i in range(first_it,last_it):
                run_subdir = str(i)
   
                # Evaluate maximum PSD in the given chunk
                freq_index = np.where((freq <= freq_right[i]) & (freq >= freq_left[i]))[0]
                psd_chunk = psd[freq_index]
                freq_chunk = freq[freq_index]
                max_psd_chunk = max(psd_chunk)
                min_freq = min(freq_chunk)

                if best_dnu <= self.cp.dnu_rg:
                    # In case the star is classified as a RG, distinguish among RGB, clump and RGB-tip
                    if best_dnu <= self.cp.dnu_tip:
                        min_linewidths[i] = astero.get_linewidth(min_freq,self.teff,numax)/self.cp.fwhm_chunk_scaling_tip
                        avg_fwhm = freqbin
                    else:
                        if best_dnu <= self.cp.dnu_cl:
                            min_linewidths[i] = astero.get_linewidth(min_freq,self.teff,numax)/self.cp.fwhm_chunk_scaling_cl
                        else:
                            min_linewidths[i] = astero.get_linewidth(min_freq,self.teff,numax)/self.cp.fwhm_chunk_scaling_rg                
                        avg_fwhm = min_linewidths[i]*self.cp.smoothing_fwhm_factor_rg
                else:
                    if (best_dnu < self.cp.dnu_sg) & (self.teff < self.cp.teff_sg):
                        min_linewidths[i] = astero.get_linewidth(np.mean(freq_chunk),self.teff,numax)/self.cp.fwhm_chunk_scaling_sg
                    else:
                        min_linewidths[i] = astero.get_linewidth(np.mean(freq_chunk),teff,numax)/self.cp.fwhm_chunk_scaling_ms
                    avg_fwhm = np.mean(astero.get_linewidth([min(freq_chunk),max(freq_chunk)],self.teff,numax))

                bg_names[i] = self.bgp['name']
       
                # Evaluate smoothed PSD to obtain more accurate prior height
                smth_bins = int(avg_fwhm/freqbin)
                spsd = smooth(psd,window_len=smth_bins,window='flat')
                spsd_chunk = spsd[freq_index]
                bg_level_chunk = bg_level[freq_index]
               
                
                with open(self.star_dir/self.cp.isla_subdir/(self.cp.prior_filename + '_' + run_subdir + '.txt'),'w') as f:
                    f.write('#\n')
        
                    # Apply overlapping chunks if the star is a subgiant or redgiant.
                    # This should solve the problem of having mixed modes that happen to fall very close to the
                    # radial mode of the previous chunk.
                    if i != 0:
                        if (best_dnu < self.cp.dnu_rg) or ((best_dnu < self.cp.dnu_sg) & (teff < self.cp.teff_sg)):
                            f_left = freq_left[i] - best_dnu*self.cp.dnu_overlap_fraction_rg
                            f.write('{:.5f}  \t{:.5f}\n'.format(f_left,freq_right[i]))
                        else:
                            f_left = freq_left[i] - best_dnu*self.cp.dnu_overlap_fraction_ms
                            f.write('{:.5f}  \t{:.5f}\n'.format(f_left,freq_right[i]))
                    else:
                        f_left = freq_left[i]
                        f.write('{:.5f}  \t{:.5f}\n'.format(freq_left[i],freq_right[i]))
        
                    f.write('{:.5f}  \t{:.5f}\n'.format(0,max(spsd_chunk)-np.mean(bg_level_chunk)))

                # Saving locally the full range that is used for the sampling. Not just the limits of the separations. Useful for plotting the right range
                f_temp_index = np.where((freq <= freq_right[i]) & (freq >= f_left))[0]
                self.freq[i] = freq[f_temp_index]
                self.psd[i] = psd[f_temp_index]
                self.spsd[i] = spsd[f_temp_index]
                self.bg_level[i] = bg_level[f_temp_index]
                

            if (run >= 0) and (run < n_chunks):
                run_labels = str(run_labels[first_it])
                bg_names = bg_names[first_it]
                min_linewidths = min_linewidths[first_it]

            peakbagging_parameters = { 'subdir':     self.cp.isla_subdir,
                                        'run':        run_labels,
                                        'background': bg_names,
                                        'fwhm':       min_linewidths,
                                        'duplet':     0}

            flag_computation_completed = diamonds.run_peakbagging(self.catalog_id,self.star_id,peakbagging_parameters,0,0,0,self.cp.dp_isla,self.cp.diamonds_path, self.cp.n_threads, self.cp.prior_filename,merge=True)


        # Give as input the minimum number of adjacent bins required to consider two local maxima separated. Don't need to calculate for each chunk. Just save number now.
        if best_dnu < self.cp.dnu_rg:
            threshold_asef = self.cp.threshold_asef_chunk_rg
        else:
            if (best_dnu < self.cp.dnu_sg) & (self.teff < self.cp.teff_sg):
                threshold_asef = self.cp.threshold_asef_chunk_sg
            else:
                threshold_asef = self.cp.threshold_asef_chunk_ms

        # Save some attributes to keep around
        self.numax = numax
        self.scaling_dnu = dnu

        self.chunk_number = chunk_number
        self.bg_names = bg_names
        self.min_linewidths = min_linewidths
        self.snr = snr
        self.threshold_asef = threshold_asef

        # Possibly useful later?
        self.maxima = [None]*n_chunks
        self.ranges = [None]*n_chunks
        self.divisions = [None]*n_chunks

        self.freqs = [None]*n_chunks
        self.freqs_sig = [None]*n_chunks
        self.degrees = [None]*n_chunks
        self.orders = [None]*n_chunks
        self.freqs_global = [None]*n_chunks
        
        self.nested_iters = [None]*n_chunks
        self.par_hist = [None]*n_chunks 
        self.asef_hist = [None]*n_chunks
        self.asef_bins = [None]*n_chunks

        self.local_epsi = [None]*n_chunks
        self.local_d02 = [None]*n_chunks
        self.local_dp = [None]*n_chunks
        self.low_cut_frequency = [None]*n_chunks

        self.n_radial_chunk = [None]*n_chunks
        self.fwhm_radial_fit = [None]*n_chunks
        self.avg_fwhm = [None]*n_chunks
        self.upper_height = [None]*n_chunks
        self.n_freqs = [None]*n_chunks
        self.freqs_radial_chunk = [None]*n_chunks
        self.octupole_freq_asymp = [None]*n_chunks
        self.octupole_freq_lower = [None]*n_chunks
        self.octupole_freq_upper = [None]*n_chunks
        self.fit_linewidth = [None]*n_chunks
        
        # pickle to save stuff
        if self.cp.save_progress_pickle:
            pickle.dump(self,open(self.star_dir/(self.catalog_id+self.star_id+'_chunk.pickle'),'wb'))

        return self.snr,self.chunk_number


    def find_islands(self,run,force=False):
        run = str(run)
        run_number = int(run)
        peakbagging_filename_global = self.star_dir/self.cp.summary_subdir/(self.catalog_id + self.star_id + self.cp.peakbagging_filename_label + self.cp.isla_subdir + '_' + self.cp.global_subdir + '_GLOBAL.txt')
        peakbagging_filename_chunk = self.star_dir/self.cp.summary_subdir/(self.catalog_id + self.star_id + self.cp.peakbagging_filename_label + self.cp.isla_subdir + '_')
        
        # Copy FAMED configuring parameters used in this run. JM: Do this here? or in init?
        target_dir = self.star_dir/self.cp.summary_subdir
        if self.cp.local_config:
                shutil.copy('famed_config.yml',target_dir/('famed_config_'+self.catalog_id + self.star_id + '_' + self.cp.isla_subdir + '_' + run +'_CHUNK.yml'))
        else:
            shutil.copy(self.cp.famed_path/'famed_config.yml',target_dir/('famed_config_'+self.catalog_id + self.star_id + '_' + self.cp.isla_subdir + '_' + run +'_CHUNK.yml'))
        shutil.copy(self.cp.configuring_parameters_file,target_dir/('famed_configuring_parameters_'+self.catalog_id + self.star_id + '_' + self.cp.isla_subdir + '_' + run +'_CHUNK.txt'))

            
        # Read sampled frequency from DIAMONDS multi-modal fit
        par0 = np.loadtxt(self.star_dir/self.cp.isla_subdir/run/'peakbagging_parameter000.txt')
        self.nested_iters[run_number] = par0
        
        # Read prior height (upper limit) and frequency range
        prior_down,prior_up = np.loadtxt(self.star_dir/self.cp.isla_subdir/run/'peakbagging_hyperParametersUniform.txt',unpack=True)
        freq_left_chunk = prior_down[0]
        freq_right_chunk = prior_up[0]
        upper_height = prior_up[1]

        # Read input linewidth and background name of the run
        config = np.loadtxt(self.star_dir/self.cp.isla_subdir/run/'peakbagging_computationParameters.txt',dtype=str)
        fit_linewidth = float(config[len(config)-1])
        background_name = config[len(config)-2]

        # Read star PSD
        peakbagging_data_dir = self.cp.diamonds_path/'PeakBagging'/'data'
        freq,psd = np.loadtxt(peakbagging_data_dir/(self.catalog_id + self.star_id + '.txt'), unpack=True)
        freqbin = freq[1]-freq[0]

        # Read Nyquist frequency
        
        nyq = np.loadtxt(self.star_dir/'NyquistFrequency.txt')

        if self.cp.print_on_screen:
            print('---------------------------------------------------')
            print(' Parameter range (microHz): [{:.2f}, {:.2f}]'.format(min(par0),max(par0)))
            print('---------------------------------------------------\n')
            print('-------------------------------------------------')
            print(' PSD frequency range (microHz): [{:.2f}, {:.2f}]'.format(min(freq),max(freq)))
            print('-------------------------------------------------\n')

        # Load the information available from the global fit
        # Load asymptotic parameters
        acf_dnu,best_dnu,best_epsi,best_alpha,teff,n_chunks,flag_depressed_dipole = np.loadtxt(peakbagging_filename_global, max_rows=1, usecols=(1,2,3,4,5,6,7))
        n_chunks = int(n_chunks)
        flag_depressed_dipole = int(flag_depressed_dipole)

        # Load the SNR for the given chunk
        chunk_number,freq_left,freq_right,snr_chunks = np.loadtxt(peakbagging_filename_global, max_rows=n_chunks, skiprows=3, unpack=True)
        snr = snr_chunks[run_number]
        left_bound = freq_left[run_number]
        right_bound = freq_right[run_number]

        # Load global frequencies
        enn_global,ell_global,freq_global,freq_sig_global,fwhm_global = np.loadtxt(peakbagging_filename_global,skiprows=3+n_chunks,usecols=(0,1,2,3,5), unpack=True)
        enn_global = np.array(enn_global,dtype=int)
        ell_global = np.array(ell_global,dtype=int)

        # Select only global frequencies relevant for the chunk
        tmp_chunk = np.where((freq_global <= right_bound) & (freq_global >= left_bound))[0]
        if len(tmp_chunk) > 0:
            freq_global = freq_global[tmp_chunk]
            freq_sig_global = freq_sig_global[tmp_chunk]
            enn_global = enn_global[tmp_chunk]
            ell_global = ell_global[tmp_chunk]
            fwhm_global = fwhm_global[tmp_chunk]
            self.freqs_global[run_number] = freq_global
        else:
            print('This chunk is empty according to global fit. Skip chunk.')
            return False


        # Select the radial mode. If the global modality was ran correctly, there should be only one radial mode, but there could be none as well.
        # If more than one radial mode is present, then pick up the one closest to the right frequency bound of the chunk.

        tmp_radial = np.where(ell_global==0)[0]
        if len(tmp_radial) > 0:
            if len(tmp_radial) == 1:
                enn_radial = enn_global[tmp_radial][0]
                freq_radial = freq_global[tmp_radial][0]
                freq_sig_radial = freq_sig_global[tmp_radial][0]
                fwhm_radial = fwhm_global[tmp_radial][0]
            else:
                freq_radial = freq_global[tmp_radial]
                good_radial_index = closest(right_bound,freq_radial,index=True)
                good_radial_index = good_radial_index + min(tmp_radial)
                enn_radial = enn_global[good_radial_index]
                freq_radial = freq_global[good_radial_index]
                freq_sig_radial = freq_sig_global[good_radial_index]
                fwhm_radial = fwhm_global[good_radial_index]
            n_radial_chunk = 1
        else:
            n_radial_chunk = 0


        # Select the dipole mode(s). There might be none as well, depending if it was selected during the global fit.
        tmp_dipole = np.where(ell_global==1)[0]
        if len(tmp_dipole) > 0:
            # If more than one dipole mode frequency is selected, take the one closest to nu0 - best_dnu/2.
            if len(tmp_dipole) > 1:
                freq_dipole = freq_global[tmp_dipole]

                if n_radial_chunk==0:
                    # If the radial mode is not present, take the dipole frequency closest to the central frequency of the chunk
                    dipole_global_index = closest((left_bound + right_bound)/2.,freq_dipole,index=True)
                else:
                    dipole_global_index = closest(freq_radial - best_dnu/2.,freq_dipole,index=True)
        

                best_dipole_index = dipole_global_index + min(tmp_dipole)
                freq_dipole = freq_global[best_dipole_index]
                enn_dipole = enn_global[best_dipole_index]
                freq_sig_dipole = freq_sig_global[best_dipole_index]
                
                fwhm_dipole = fwhm_global[best_dipole_index]
            else:
                enn_dipole = enn_global[tmp_dipole][0]
                freq_dipole = freq_global[tmp_dipole][0]
                freq_sig_dipole = freq_sig_global[tmp_dipole][0]
                
                fwhm_dipole = fwhm_global[tmp_dipole][0]
            n_dipole_chunk = 1
        else:
            n_dipole_chunk = 0

        # Load the nuMax information as obtained from a previous background fit
        gauss_par = np.loadtxt(self.star_dir/'gaussianEnvelopeParameters.txt')
        numax = gauss_par[1]
        scaling_dnu = astero.compute_scaling_dnu(numax)


        # Set a binwidth proportional to the expected separation between adjacent peaks (normally taken as d02/2) if MS
        # or smaller if RG star.
        par_range = max(par0) - min(par0)
        min_separation = self.get_minimum_freq_separation(best_dnu,global_flag=False)
        binwidth = 1.0*min_separation/self.cp.min_n_bins
        n_bins = int(np.ceil(par_range/binwidth))
        nest_iter = np.arange(len(par0))

        # Compute an Average Shifted Histogram (ASH) of the distribution of nested iterations
        # Compute the ASEF with input resolution for extracting the local maxima
        par_hist, asef_hist = self.compute_asef(par0,nest_iter,n_bins)
        
        # Find the local maxima using a hill-climbing algorithm on the envelope function        
        index_maximum = self.hill_climbing(par_hist,asef_hist,self.threshold_asef,self.cp.min_bin_separation)
        maximum = par_hist[index_maximum]
        asef_maximum = asef_hist[index_maximum]

        # Identify the divisions among the local maxima. These will be used to select chunks of parameter range that locate
        # an optimal frequency region to perform a fit around the peak.
        # Also compute ranges for each maximum such that an ASEF peak is considered up to the points of its tails, i.e. 
        # the last ones in the decreasing phase off the peak. Ranges are used to compute the actual frequency estimates
        # from the sampling.
        n_maxima = len(maximum)
        range_maximum, divisions_maximum = self.get_range_divisions(par_hist,asef_hist,index_maximum,chunk=True)

        self.maxima[run_number] = maximum
        self.ranges[run_number] = range_maximum
        self.divisions[run_number] = divisions_maximum
        
        # Based on the identified frequency ranges, estimate the oscillation frequencies using the sampling information 
        # obtained by DIAMONDS.
        n_freq = n_maxima

        # Compute a smoothed PSD by some average FWHM to evaluate its values at each local maxima of the ASEF.
        # In case of RG, adopt a finer smoothing window to overcome the problem of very narrow mixed modes.

        avg_fwhm = np.mean(astero.get_linewidth(maximum,self.teff,numax,self.cp.numax_threshold))
        if best_dnu <= self.cp.dnu_rg:
            avg_fwhm = fit_linewidth*self.cp.smoothing_fwhm_factor_rg
 

        smth_bins = int(avg_fwhm/freqbin)
        spsd = smooth(psd,window_len=smth_bins)
        tmp_good = np.where((freq <= max(par_hist)) & (freq >= min(par_hist)))[0]
        psd_total = psd
        freq_total = freq
        psd = psd[tmp_good]
        freq = freq[tmp_good]
        spsd = spsd[tmp_good]

        # Load the background level of the star

        bg_level = np.loadtxt(self.star_dir/'backgroundLevel.txt',usecols=(1,))
        bg_level_local = bg_level[tmp_good]
        sampled_estimates = self.evaluate_sampling_frequencies(par0, par_hist, freq, spsd, maximum, range_maximum)
        freq1 = sampled_estimates['freq1']
        freq_sig1 = sampled_estimates['freq_sig1']
        sampling_counts = sampled_estimates['sampling_counts']
        spsd_maximum = sampled_estimates['spsd_maximum']

        # Define some weights useful for identification of proper frequency peaks during mode identification
        sampling_weights = sampling_counts/np.sum(sampling_counts)
        spsd_weights = spsd_maximum/np.sum(spsd_maximum)
        asef_weights = asef_maximum/np.sum(asef_maximum)

        if self.cp.print_on_screen:
            print(' Total number of local maxima found: ',n_maxima)

        # Define some average asymptotic frequency spacings useful for the computation
        ap = astero.get_asymptotic_parameters(numax, best_dnu, self.teff, self.cp.d01_mass_offset, self.cp.d01_mass_slope, self.cp.d01_offset, self.cp.d02_mass_offset, self.cp.d02_mass_slope, self.cp.d02_offset, self.cp.d03_slope,self.cp.d03_offset, self.cp.numax_sun, self.cp.dnu_sun, self.cp.teff_sun)

        angular_degree = np.ones(n_freq,dtype=int)
        d02 = ap['d02']
        d01 = ap['d01']
        d03 = ap['d03']

        # If possible, update the d02 spacing from existing chunk outputs
        filename_summary = np.sort(glob.glob(str(peakbagging_filename_chunk) + '*' + self.modality + '.txt'))
        
        flag_median_d02_active = 0
        if len(filename_summary) > 0:
            d02_array = np.zeros(len(filename_summary))
            for jj in range(0, len(filename_summary)):
                d02_chunk = np.loadtxt(filename_summary[jj],usecols=(1,), max_rows=1)
                d02_array[jj] = d02_chunk
    
            # Remove zeros from possible solutions (e.g. if a chunk had no l=2 mode detected)
            d02_array = d02_array[np.where(d02_array>0)[0]]

            median_d02 = np.median(d02_array)
            max_d02 = max(d02_array)

            if (max_d02 < d02) or (max_d02==0):
                max_d02 = d02
            if median_d02==0:
                median_d02 = d02
            flag_median_d02_active = 1
        else:
            median_d02 = d02
            max_d02 = d02

        # Radial mode and Quadrupole mode

        # Start by obtaining a more accurate position of the current radial mode by using the (possible) radial mode of the
        # previous (or next) chunk. The previous (or next) radial mode will correspond to a chunk with a higher SNR (also assuming that the chunks
        # are being calculated by decreasing SNR order). If no value could be found, double the step.

        previous_run_number = run_number - 1
        next_run_number = run_number + 1
        previous_run_number2 = run_number - 2
        next_run_number2 = run_number + 2

        flag_previous_radial_mode_found = 0
        flag_next_radial_mode_found = 0
        
        reference_central_freq = (freq_left_chunk + freq_right_chunk)/2.

        if reference_central_freq >= numax:
            if previous_run_number >= 0:
                if os.path.isfile(str(peakbagging_filename_chunk) + str(previous_run_number) + '_' + self.modality + '.txt'):
                    enn_previous_chunk, ell_previous_chunk, freq_previous_chunk, freq_sig_previous_chunk = np.loadtxt(str(peakbagging_filename_chunk) + str(previous_run_number) + '_' + self.modality + '.txt', usecols=(0,1,3,4),skiprows=7,unpack=True)
                    if isinstance(enn_previous_chunk,float):
                        enn_previous_chunk = [int(enn_previous_chunk)]
                        ell_previous_chunk = [ell_previous_chunk]
                        freq_previous_chunk = [freq_previous_chunk]
                        freq_sig_previous_chunk = [freq_sig_previous_chunk]
                    else:                        
                        enn_previous_chunk = np.array(enn_previous_chunk,dtype=int)
                    tmp_previous_radial = np.where(ell_previous_chunk==0)[0]
           
                    # Check whether the selected chunk contains a l=0 mode, otherwise move to the one before it.
                    if len(tmp_previous_radial) > 0:
                        if self.cp.print_on_screen:
                            print('\n Resuming l=0 frequency from previous chunk...\n')
                    
                        freq_previous_radial = freq_previous_chunk[tmp_previous_radial]
                        freq_sig_previous_radial = freq_sig_previous_chunk[tmp_previous_radial]
                        enn_previous_radial = enn_previous_chunk[tmp_previous_radial][0]
                        enn_radial = enn_previous_radial + 1
                        freq_previous_radial = freq_previous_radial[0]
                        freq_sig_previous_radial = freq_sig_previous_radial[0]
                        flag_previous_radial_mode_found = 1
               
                        # Set new global radial mode frequency
                 
                        freq_radial = freq_previous_radial + best_dnu*(1.0 + best_alpha*(enn_radial - 0.5 - numax/best_dnu))
                        fwhm_radial = astero.get_linewidth(freq_radial,self.teff,numax)
               
                        if n_radial_chunk==0:
                            if freq_radial <= max(par0):
                                n_radial_chunk = 1
                                freq_sig_radial = freq_sig_dipole
                                
                        if n_dipole_chunk==0:
                            freq_dipole = freq_radial - best_dnu/2. - d01
                            freq_sig_dipole = freq_sig_radial
                            n_dipole_chunk = 1
                
                    else:
                        if previous_run_number2 >= 0:
                            if os.path.isfile(str(peakbagging_filename_chunk) + str(previous_run_number2) + '_' + self.modality + '.txt'):
                                enn_previous_chunk, ell_previous_chunk, freq_previous_chunk, freq_sig_previous_chunk = np.loadtxt(str(peakbagging_filename_chunk) + str(previous_run_number2) + '_' + self.modality + '.txt', usecols=(0,1,3,4),skiprows=7,unpack=True)
                                if isinstance(enn_previous_chunk,float):
                                    enn_previous_chunk = [int(enn_previous_chunk)]
                                    ell_previous_chunk = [ell_previous_chunk]
                                    freq_previous_chunk = [freq_previous_chunk]
                                    freq_sig_previous_chunk = [freq_sig_previous_chunk]
                                else:
                                    enn_previous_chunk = np.array(enn_previous_chunk,dtype=int)
                                tmp_previous_radial = np.where(ell_previous_chunk==0)[0]
                        
                                if len(tmp_previous_radial) > 0:
                                    if self.cp.print_on_screen:
                                        print('\n Resuming l=0 frequency from second previous chunk...\n')
                                
                                    freq_previous_radial = freq_previous_chunk[tmp_previous_radial]
                                    freq_sig_previous_radial = freq_sig_previous_chunk[tmp_previous_radial]
                                    enn_previous_radial = enn_previous_chunk[tmp_previous_radial][0]
                                    enn_radial = enn_previous_radial + 2
                                    freq_previous_radial = freq_previous_radial[0]
                                    freq_sig_previous_radial = freq_sig_previous_radial[0]
                                    flag_previous_radial_mode_found = 2
                                
                                    # Set new global radial mode frequency
                                    
                                    freq_radial = freq_previous_radial + 2*best_dnu*(1.0 + best_alpha*((2*enn_radial-1)/2. - 0.5 - numax/best_dnu))
                                    fwhm_radial = astero.get_linewidth(freq_radial,self.teff,numax)
                           
                                    if n_radial_chunk==0:
                                        if freq_radial <= max(par0):
                                            n_radial_chunk = 1
                                            freq_sig_radial = freq_sig_dipole
                                            
                                    if n_dipole_chunk==0:
                                        freq_dipole = freq_radial - best_dnu/2. - d01
                                        freq_sig_dipole = freq_sig_radial
                                        n_dipole_chunk = 1
        else:
            if next_run_number <= n_chunks-1:
                if os.path.isfile(str(peakbagging_filename_chunk) + str(next_run_number) + '_' + self.modality + '.txt'):
                    enn_next_chunk, ell_next_chunk, freq_next_chunk, freq_sig_next_chunk = np.loadtxt(str(peakbagging_filename_chunk) + str(next_run_number) + '_' + self.modality + '.txt', usecols=(0,1,3,4),skiprows=3,unpack=True)
                    if isinstance(enn_next_chunk,float):
                        enn_next_chunk = [int(enn_next_chunk)]
                        ell_next_chunk = [ell_next_chunk]
                        freq_next_chunk = [freq_next_chunk]
                        freq_sig_next_chunk = [freq_sig_next_chunk]
                    else:                        
                        enn_next_chunk = np.array(enn_next_chunk,dtype=int)
                    tmp_next_radial = np.where(ell_next_chunk==0)[0]
           
                    # Check whether the selected chunk contains a l=0 mode, otherwise move to the one next to it.

                    if len(tmp_next_radial) > 0:
                        if self.cp.print_on_screen:
                            print('\n Resuming l=0 frequency from next chunk...\n')

                        freq_next_radial = freq_next_chunk[tmp_next_radial]
                        freq_sig_next_radial = freq_sig_next_chunk[tmp_next_radial] 
                        enn_next_radial = enn_next_chunk[tmp_next_radial][0]
                        enn_radial = enn_next_radial - 1
                        freq_next_radial = freq_next_radial[0]
                        freq_sig_next_radial = freq_sig_next_radial[0]
                        flag_next_radial_mode_found = 1
               
                        # Set new global radial mode frequency
                        freq_radial = freq_next_radial - best_dnu*(1.0 + best_alpha*(enn_next_radial - 0.5 - numax/best_dnu))
                        fwhm_radial = astero.get_linewidth(freq_radial,self.teff,numax)
                
                        if n_radial_chunk==0:
                            if freq_radial >= min(par0):
                                n_radial_chunk = 1
                                freq_sig_radial = freq_sig_dipole
                                                                
                        if n_dipole_chunk==0:
                            freq_dipole = freq_radial - best_dnu/2. - d01
                            freq_sig_dipole = freq_sig_radial
                            n_dipole_chunk = 1
                
                    else:
                        if next_run_number2 <= n_chunks-1:
                            if os.path.isfile(str(peakbagging_filename_chunk) + str(next_run_number2) + '_' + self.modality + '.txt'):
                                enn_next_chunk, ell_next_chunk, freq_next_chunk, freq_sig_next_chunk = np.loadtxt(str(peakbagging_filename_chunk) + str(next_run_number2) + '_' + self.modality + '.txt', usecols=(0,1,3,4),skiprows=7,unpack=True)
                                if isinstance(enn_next_chunk,float):
                                    enn_next_chunk = [int(enn_next_chunk)]
                                    ell_next_chunk = [ell_next_chunk]
                                    freq_next_chunk = [freq_next_chunk]
                                    freq_sig_next_chunk = [freq_sig_next_chunk]
                                else:                        
                                    enn_next_chunk = np.array(enn_next_chunk,dtype=int)
                                tmp_next_radial = np.where(ell_next_chunk==0)[0]
                        
                                if len(tmp_next_radial) > 0:
                                    if self.cp.print_on_screen:
                                        print('\n Resuming l=0 frequency from second next chunk...\n')

                                    freq_next_radial = freq_next_chunk[tmp_next_radial]
                                    freq_sig_next_radial = freq_sig_next_chunk[tmp_next_radial]
                                    enn_next_radial = enn_next_chunk[tmp_next_radial][0]
                                    enn_radial = enn_next_radial - 2
                                    freq_next_radial = freq_next_radial[0]
                                    freq_sig_next_radial = freq_sig_next_radial[0]
                                    flag_next_radial_mode_found = 2
                           
                                    # Set new global radial mode frequency
                                    freq_radial = freq_next_radial - 2*best_dnu*(1.0 + best_alpha*((2*enn_next_radial-1)/2. - 0.5 - numax/best_dnu))
                                    fwhm_radial = astero.get_linewidth(freq_radial,self.teff,numax)
                            
                                    if n_radial_chunk==0:
                                        if freq_radial >= min(par0):
                                            n_radial_chunk = 1
                                            freq_sig_radial = freq_sig_dipole
                                                                                        
                                    if n_dipole_chunk==0:
                                        freq_dipole = freq_radial - best_dnu/2. - d01
                                        freq_sig_dipole = freq_sig_radial
                                        n_dipole_chunk = 1
                            
        flag_quadrupole_found = 0

        if n_radial_chunk != 0:
            # ----------------------------------------------------------------------
            # CASE 1: One radial mode is found in the chunk from the global modality
            # ----------------------------------------------------------------------

            # Radial mode

            # There must be only one reference radial mode in the chunk for the routine to work properly. 
            # If the global modality was run successfully, this should not constitute a problem.
    
            order_number = np.zeros(n_freq,dtype=int) + enn_radial - 1 
   
            # Select the frequency corresponding to the closest and highest (in counts) peak to the global l=0. 
            # This is the most likely one to be the correct radial mode frequency because of its large linewidth as
            # compared to other dipole modes (in the case of RG) or spurious peaks. 
            # This applies only if more than one frequency is present.
            # Possible contamination by l=2 mode should be avoided.
  
            freq_diff = np.abs(freq1 - freq_radial)
            freq_weights = 1.0/freq_diff
            freq_weights /= np.sum(freq_weights)
            freq_ww = freq_weights/max(freq_weights)
            asef_ww = asef_weights/max(asef_weights)
            spsd_ww = spsd_weights/max(spsd_weights)
            sampling_ww = np.log(sampling_counts)/max(np.log(sampling_counts))

            # If we have a radial mode solution from previous (or next) chunks, then make the frequency position of the global radial mode a much stronger constraint
            weight_freq_fraction = self.cp.weight_freq_fraction
            if (flag_previous_radial_mode_found != 0) or (flag_next_radial_mode_found != 0):
                weight_freq_fraction = self.cp.weight_freq_fraction_enhanced

            total_ww = weight_freq_fraction*freq_ww + self.cp.weight_asef_fraction*asef_ww + self.cp.weight_spsd_fraction*spsd_ww + self.cp.weight_sampling_fraction*sampling_ww
            total_ww /= np.sum(total_ww)

            upper_limit_freq_radial = freq_radial + max_d02*self.cp.d02_factor_search_range

            # Make sure that the upper limit frequency for radial mode search is not larger than chunk upper limit frequency
            if upper_limit_freq_radial > right_bound:
                upper_limit_freq_radial = right_bound

            # Force input value if required 
            if self.cp.upper_limit_freq_radial != 0:
                upper_limit_freq_radial = self.cp.upper_limit_freq_radial
    
            if self.cp.print_on_screen:
                print(' Upper limit frequency for l=0: {:.3f}  muHz'.format(upper_limit_freq_radial))

            additional_freq_sig_radial = 0.
            if flag_previous_radial_mode_found != 0:
                additional_freq_sig_radial = freq_sig_previous_radial
            if flag_next_radial_mode_found != 0:
                additional_freq_sig_radial = freq_sig_next_radial     

            if n_dipole_chunk != 0:
                tmp_radial = np.where((freq1 >= (freq_radial - freq_sig_radial - additional_freq_sig_radial - d02)) & (freq1 < upper_limit_freq_radial) & (freq1 > freq_dipole))[0]
            else:
                tmp_radial = np.where((freq1 >= (freq_radial - freq_sig_radial - additional_freq_sig_radial - d02)) & (freq1 < upper_limit_freq_radial))[0]
                
            if len(tmp_radial) > 0:
                total_ww_radial = total_ww[tmp_radial]
                max_ww = max(total_ww_radial)
                index = np.where(total_ww_radial == max_ww)[0][0]
                radial_index = index + min(tmp_radial)
                min_ww = min(total_ww)
            else:
                print(' The l=0 mode could not be located based on the position from the global fit. ')
                print(' This could be caused by the low SNR of the chunk.')
                return False   
    
            # Check if the selected peak is not the adjacent l=2 by assessing the total weight of the subsequent modes.
            if radial_index < n_freq-1:
                candidate_radial_index = np.where((freq1 > freq1[radial_index]) & (freq1 < upper_limit_freq_radial))[0]
                if len(candidate_radial_index) > 0:
                    for kk in range(0,len(candidate_radial_index)):
                        local_index = candidate_radial_index[kk]
                        if total_ww[local_index] > (abs(max_ww - min_ww)*self.cp.max_ratio_search_radial + min_ww):
                            radial_index = local_index
                            
            freq_radial_chunk = freq1[radial_index]
            freq_sig_radial_chunk = freq_sig1[radial_index]
            low_cut_frequency = freq_radial_chunk - best_dnu*(1.0 + best_alpha*(enn_radial - 0.5 - numax/best_dnu)) + freq_sig_radial/2.0 
  
            # Make sure that the previous l=0 mode is not included, in case the chunk belongs to the right-end tail of the oscillation envelope.
            # First verify that if a radial mode has been found from previous chunks (to the left side), the low cut frequency is not including it.
            if flag_previous_radial_mode_found==1:
                if low_cut_frequency < freq_previous_radial + freq_sig_radial/2.:
                    low_cut_frequency = freq_previous_radial + freq_sig_radial/2.
            if flag_previous_radial_mode_found==2:
                if low_cut_frequency < freq_previous_radial + best_dnu*(1.0 + best_alpha*(enn_radial-1 - 0.5 - numax/best_dnu)) + freq_sig_radial/2.:
                    low_cut_frequency = freq_previous_radial + best_dnu*(1.0 + best_alpha*(enn_radial-1 - 0.5 - numax/best_dnu)) + freq_sig_radial/2.
    
            # However, the l=0 mode from the previous chunk may turnout to be wrong (e.g. if confused with an adjacent l=1 mixed mode). 
            # To verify this do some additional check to test whether the selected low cut frequency can be reliable. 
            if freq_radial_chunk > numax:
                first_chunk_indices = np.where(freq1 < min(par_hist) + (max(par_hist) - min(par_hist))/self.cp.previous_radial_range_fraction)[0]
        
                if len(first_chunk_indices) > 0:
                    first_chunk_sampling_counts = sampling_counts[first_chunk_indices]
                    first_chunk_asef = asef_maximum[first_chunk_indices]
            
                    # Do a first check on the sampling counts
                    previous_radial_mode_index = np.argmax(first_chunk_sampling_counts)
                    max_first_chunk_sampling_counts = first_chunk_sampling_counts[previous_radial_mode_index]
                    if sampling_counts[radial_index] < max_first_chunk_sampling_counts*self.cp.sampling_counts_fraction:
                        # Here update the lower cutting frequency of the chunk, as well as the radial mode index
                        low_cut_frequency2 = freq1[previous_radial_mode_index] + freq_sig_radial_chunk
                        if low_cut_frequency2 > low_cut_frequency:
                            low_cut_frequency = low_cut_frequency2

                        radial_index_new = closest(freq1[previous_radial_mode_index] + best_dnu*(1.0 + best_alpha*(enn_radial - 0.5 - numax/best_dnu)),freq1,index=True)
                        if sampling_counts[radial_index_new] >= sampling_counts[radial_index]:
                            radial_index = radial_index_new
                            freq_radial_chunk = freq1[radial_index]
                
                    # Do a second check on the ASEF maximum
                    previous_radial_mode_index = np.argmax(first_chunk_asef)
                    max_first_chunk_asef = first_chunk_asef[previous_radial_mode_index]
                    if max_first_chunk_asef*self.cp.asef_saturation_fraction >= asef_maximum[radial_index]:
                        # Here update the lower cutting frequency of the chunk, as well as the radial mode index, if adequate.
                        low_cut_frequency2 = freq1[previous_radial_mode_index] + freq_sig_radial_chunk
                        if low_cut_frequency2 > low_cut_frequency:
                            low_cut_frequency = low_cut_frequency2

                        radial_index_new = closest(freq1[previous_radial_mode_index] + best_dnu*(1.0 + best_alpha*(enn_radial - 0.5 - numax/best_dnu)),freq1,index=True)
                        asef_threshold = (self.cp.dp_isla['max_nested_it'] + self.cp.dp_isla['n_live'])/self.cp.asef_threshold_scaling_radial 
                        if asef_maximum[radial_index_new] > asef_threshold:
                            radial_index = radial_index_new
                            freq_radial_chunk = freq1[radial_index]
                       
            # Apply a final check on the upper bound for the low cut frequency
            if low_cut_frequency > freq_radial_chunk - best_dnu*self.cp.dnu_lower_cut_fraction:
                low_cut_frequency = freq_radial_chunk - best_dnu*self.cp.dnu_lower_cut_fraction

            angular_degree[radial_index] = 0
            order_number[radial_index] = enn_radial

            # Compute a local value for epsilon
            local_epsi = (freq_radial_chunk - best_dnu * enn_radial)/best_dnu


            # Quadrupole mode
    
            # Identify the quadrupole mode by assuming a small spacing d02.
            # Compute more reliable values for d02 from high SNR chunks, if available.
            # Mostly useful for the case of MS stars.
            freq_radial_chunk_org = freq_radial_chunk
            freq_sig_radial_chunk_org = freq_sig_radial_chunk

            flag_quadrupole_found = 1
            flag_duplet_fit = 0

            if best_dnu < self.cp.dnu_rg:
                # Evolutionary stage: RG
                quadrupole_freq_asymp = freq_radial_chunk - d02
                quadrupole_index = closest(quadrupole_freq_asymp,freq1,index=True)
                if quadrupole_index != radial_index:
                    freq_quadrupole_chunk = freq1[quadrupole_index]
                    freq_sig_quadrupole_chunk = freq_sig1[quadrupole_index]
                    range_maximum_quadrupole = range_maximum[:,quadrupole_index]
                    range_maximum_radial = range_maximum[:,radial_index]
                    angular_degree[quadrupole_index] = 2
                    order_number[quadrupole_index] = enn_radial - 1

                    # Find the l=1 mixed modes closest to l=2 on each side of the peak.
                    # If the mixed mode falls within d02/2 from the l=2 frequency position then merge it with the l=2
                    candidate_mixed = np.where((freq1 <= freq_quadrupole_chunk + d02/self.cp.d02_scaling_merge_mixed) & (freq1 >= freq_quadrupole_chunk - d02/self.cp.d02_scaling_merge_mixed) & (freq1 != freq_quadrupole_chunk) & (freq1 < freq_radial_chunk))[0]
                    tmp_mask = np.zeros(len(freq1),dtype=bool)
                    tmp_mask[candidate_mixed] = True 
                    good_freq = np.arange(len(freq1),dtype=int)[~tmp_mask]
    
                    if len(candidate_mixed) > 0:
                        # Remove the merged l=1 frequencies from the frequency set
                        sampling_counts[quadrupole_index] = sampling_counts[quadrupole_index] + np.sum(sampling_counts[candidate_mixed])
                        merged_freq = freq1[candidate_mixed]
                   
                        new_lower_bound = range_maximum[0, min(candidate_mixed)]
                        if new_lower_bound < range_maximum[0,quadrupole_index]:
                            range_maximum[0,quadrupole_index] = new_lower_bound
                
                        new_upper_bound = range_maximum[1, max(candidate_mixed)]
                        if new_upper_bound > range_maximum[1,quadrupole_index]:
                            range_maximum[1,quadrupole_index] = new_upper_bound
               
                        new_lower_bound = divisions_maximum[0, min(candidate_mixed)]
                        if new_lower_bound < divisions_maximum[0,quadrupole_index]:
                            divisions_maximum[0,quadrupole_index] = new_lower_bound
                
                        new_upper_bound = divisions_maximum[1, max(candidate_mixed)]
                        if new_upper_bound > divisions_maximum[1,quadrupole_index]:
                            divisions_maximum[1,quadrupole_index] = new_upper_bound
                
                        if len(good_freq) > 0:
                            # Deprecate one range and division between the lowest and highest merged frequencies
                            range_maximum = range_maximum[:,good_freq]
                            divisions_maximum = divisions_maximum[:,good_freq]
                            freq1 = freq1[good_freq]
                            freq_sig1 = freq_sig1[good_freq]
                            angular_degree = angular_degree[good_freq]
                            order_number = order_number[good_freq]
                            sampling_counts = sampling_counts[good_freq]
                            sampling_weights = sampling_counts/np.sum(sampling_counts)
                            asef_maximum = asef_maximum[good_freq]
                            spsd_maximum = spsd_maximum[good_freq]
                            spsd_weights = spsd_maximum/np.sum(spsd_maximum)
                            asef_weights = asef_maximum/np.sum(asef_maximum)
                            n_freq = len(freq1)
                
                        quadrupole_index = np.where(angular_degree==2)[0]
                        upper_bound = range_maximum[1,quadrupole_index]
                        lower_bound = range_maximum[0,quadrupole_index]
                        tmp_range = np.where((par0 < upper_bound[0]) & (par0 >= lower_bound[0]))[0]
                        par0_range = par0[tmp_range]
                
                        freq_quadrupole_chunk = np.sum(par0_range*tmp_range**2)/np.sum(tmp_range**2)
                        freq_sig_quadrupole_chunk = np.sqrt(np.sum((par0_range-freq_quadrupole_chunk)**2*tmp_range**2)/np.sum(tmp_range**2))

                        freq1[quadrupole_index] = freq_quadrupole_chunk
                        freq_sig1[quadrupole_index] = freq_sig_quadrupole_chunk
            

                else: 
                    # No quadrupole has been found.
                    flag_quadrupole_found = 0
                        
            else:
                # Evolutionary stage: MS, SG
                
                # Estimate the position of the l=2 and l=0 peaks. 
                # Consider as an upper limit for d02 given by the d02 from RG scaling in the case of stars in the subgiant regime. 
                # This is because here d02_RG  constitutes an upper limit for the actual d02 (see White et al. 2011). 
                # For MS stars, take as an upper limit for d02 the value Dnu/4 in order to incorporate also cases with very large d02.
                if flag_median_d02_active==1: 
                    max_local_d02 = median_d02*self.cp.d02_factor_search_range
                else:
                    if best_dnu <= self.cp.dnu_sg:
                        # Here the star is considered a SG
                        max_local_d02 = max_d02
                    else:
                        # Here the star is considered a MS
                        max_local_d02 = best_dnu/4.
            
    
                # Set up and run the multi-modal fit for the double peak, after checking whether the run already exists
                run_subdir = run + 'A'
                flag_duplet_fit = 1 
       
                if not os.path.isfile(self.star_dir/self.cp.isla_subdir/run_subdir/'peakbagging_computationParameters.txt') or force:
                    # Set up prior
                    freq_prior = [freq_radial_chunk - d02,upper_limit_freq_radial]
                    height_prior = [prior_down[1],prior_up[1]]

                    if self.cp.d02_prior_upper_duplet_fit != 0:
                        max_local_d02 = self.cp.d02_prior_upper_duplet_fit
            
                    d02_prior = [self.cp.d02_prior_lower_duplet_fit,max_local_d02]
                    boundaries = [freq_prior,height_prior,d02_prior]

                    filename = self.star_dir/self.cp.isla_subdir/(self.cp.prior_filename + '_' + run_subdir + '.txt')
                    diamonds.write_uniform_prior(filename,np.array(boundaries))

                    if self.cp.print_on_screen:
                        print(' Performing l=2,0 duplet fit with DIAMONDS...')
            
           
                    peakbagging_parameters = { 'subdir':       self.cp.isla_subdir,
                                               'run':          run_subdir,
                                               'background':   background_name,
                                               'fwhm':         astero.get_linewidth(freq_radial_chunk,self.teff,numax),
                                               'duplet':       1}

                    flag_computation_completed = diamonds.run_peakbagging(self.catalog_id,self.star_id,peakbagging_parameters,0,0,0, self.cp.dp_isla, self.cp.diamonds_path, self.cp.n_threads,self.cp.prior_filename)
                else:
                    if self.cp.print_on_screen:
                        print(' Load information from l=2,0 duplet fit with DIAMONDS.')

                # Read sampled frequency from DIAMONDS multi-modal fit for nu_0
                par_nu0 = np.loadtxt(self.star_dir/self.cp.isla_subdir/run_subdir/'peakbagging_parameter000.txt')
     
                # Read sampled frequency from DIAMONDS multi-modal fit for d02
                par_d02 = np.loadtxt(self.star_dir/self.cp.isla_subdir/run_subdir/'peakbagging_parameter002.txt')

                # Read posterior distribution from DIAMONDS multi-modal fit
                post = np.loadtxt(self.star_dir/self.cp.isla_subdir/run_subdir/'peakbagging_posteriorDistribution.txt')

                # Compute parameter estimates, using a weighted mean from the posterior
                post /= max(post) 
                nest_iter = np.arange(len(par_nu0)) 
                nu0 = np.sum(post*par_nu0)/np.sum(post)
                nu0_sig = np.sqrt(np.sum(post*(par_nu0 - nu0)**2)/np.sum(post))
                d02 = np.sum(post*par_d02)/np.sum(post)
                d02_sig = np.sqrt(np.sum(post*(par_d02 - d02)**2)/np.sum(post))

                freq_quadrupole_chunk = nu0 - d02
                freq_sig_quadrupole_chunk = np.sqrt(nu0_sig**2 + d02_sig**2)
                freq_radial_chunk = nu0
                freq_sig_radial_chunk = nu0_sig

                # Now make frequency uncertainties consistent in definition with the remainder of the extracted set of frequencies
                midpoint_02 = (freq_radial_chunk + freq_quadrupole_chunk)/2.
                upper_bound_quadrupole = midpoint_02
                best_quadrupole_index = max(np.where(range_maximum[0,:] <= freq_quadrupole_chunk)[0])
                lower_bound_quadrupole = range_maximum[0,best_quadrupole_index]
                #lower_bound_quadrupole = lower_bound_quadrupole[0]
       
                # Put safe conditions on the lower bound for l=2
                if lower_bound_quadrupole >= upper_bound_quadrupole:
                    print(' Adjusting lower frequency bound for l=2 due to swap with upper bound.')
                    lower_bound_quadrupole = freq_quadrupole_chunk - d02

                if lower_bound_quadrupole < freq_quadrupole_chunk - d02:
                    print(' Adjusting lower frequency bound for l=2 due to large range.')
                    lower_bound_quadrupole = freq_quadrupole_chunk - d02
        
                lower_bound_radial = midpoint_02
                best_radial_index = max(np.where(range_maximum[0,:] <= freq_radial_chunk)[0])
                upper_bound_radial = range_maximum[1,best_radial_index]
                #upper_bound_radial = upper_bound_radial[0]
        
                # Put safe conditions on the lower bound for l=0
                if lower_bound_radial >= upper_bound_radial:
                    print(' Adjusting upper frequency bound for l=0 due to swap with lower bound.')
                    upper_bound_radial += abs(upper_bound_radial - lower_bound_radial)*2 
        
                # Compute l=2 frequency uncertainty
                tmp_range = np.where((par0 < upper_bound_quadrupole) & (par0 >= lower_bound_quadrupole))[0]
                tmp_freq_range = np.where((freq < upper_bound_quadrupole) & (freq >= lower_bound_quadrupole))[0]
                tmp_hist_range = np.where((par_hist < upper_bound_quadrupole) & (par_hist >= lower_bound_quadrupole))[0]
                if len(tmp_range) == 0:
                    print(' Not enough sampling points in the specified frequency range for l=2. Quitting program.')
                    return False
         
                par0_range = par0[tmp_range]
                spsd_range = spsd[tmp_freq_range]
                spsd_maximum_quadrupole = max(spsd_range)
                asef_maximum_quadrupole = max(asef_hist[tmp_hist_range])
                sampling_counts_quadrupole = np.sum(nest_iter[tmp_range])
                freq_sig_quadrupole_chunk = np.sqrt(np.sum((par0_range-freq_quadrupole_chunk)**2*tmp_range**2)/np.sum(tmp_range**2))
       
                # Compute l=0 frequency uncertainty
                tmp_range = np.where((par0 < upper_bound_radial) & (par0 >= lower_bound_radial))[0]
                tmp_freq_range = np.where((freq < upper_bound_radial) & (freq >= lower_bound_radial))[0]
                tmp_hist_range = np.where((par_hist < upper_bound_radial) & (par_hist >= lower_bound_radial))[0]
                if len(tmp_range) == 0:
                    print(' Not enough sampling points in the specified frequency range for l=0. Quitting program.')
                    return False
        

                iteration = 0
                while len(tmp_hist_range) == 0:
                    best_radial_index += 1

                    if best_radial_index <= len(freq1):
                        upper_bound_radial = range_maximum[iteration%2,best_radial_index]
                        #upper_bound_radial = upper_bound_radial[0]
                        tmp_hist_range = np.where((par_hist < upper_bound_radial) & (par_hist >= lower_bound_radial))[0]
                    else:
                        print('Not enough bins in the specified ASEF range for l=0. Quitting program.')
                        return False
                    iteration +=1

                par0_range = par0[tmp_range]
                spsd_range = spsd[tmp_freq_range]
                spsd_maximum_radial = max(spsd_range)
                asef_maximum_radial = max(asef_hist[tmp_hist_range])
                sampling_counts_radial = np.sum(nest_iter[tmp_range])
                freq_sig_radial_chunk = np.sqrt(np.sum((par0_range-freq_radial_chunk)**2*tmp_range**2)/np.sum(tmp_range**2))

                # Check whether there are old frequencies (from local maxima) inside this region of the l=2,0 duplet
                tmp_inside = np.where((freq1 >= lower_bound_quadrupole) & (freq1 <= upper_bound_radial))[0]
                tmp_mask = np.zeros(len(freq1),dtype=bool)
                tmp_mask[tmp_inside] = True 
                tmp_outside = np.arange(len(freq1),dtype=int)[~tmp_mask]
                    
                tmp_radial_outside = np.where(angular_degree[tmp_outside]==0)[0]
                if len(tmp_radial_outside) > 0: 
                    angular_degree[tmp_radial_outside] = 1
                    order_number[tmp_radial_outside] -= 1
        
                if (len(tmp_inside) > 0) & (len(tmp_outside) > 0):
                    # If at least one frequency is found inside the l=2,0 range, first remove it (them) from the sample of frequencies
                    range_maximum = range_maximum[:,tmp_outside]
                    divisions_maximum = divisions_maximum[:,tmp_outside]
               
                    freq1 = freq1[tmp_outside]
                    freq_sig1 = freq_sig1[tmp_outside]
                    angular_degree = angular_degree[tmp_outside]
                    order_number = order_number[tmp_outside]
                    sampling_counts = sampling_counts[tmp_outside]
                    sampling_weights = sampling_counts/np.sum(sampling_counts)
                    asef_maximum = asef_maximum[tmp_outside]
                    spsd_maximum = spsd_maximum[tmp_outside]

                # Then add up the new l=2,0 frequencies to the list
                if len(tmp_outside) > 0:
                    # In this case there are other frequencies outside the l=2,0 range, so take them into account
                    range_maximum = np.hstack((range_maximum,[[lower_bound_quadrupole],[upper_bound_quadrupole]],[[lower_bound_radial],[upper_bound_radial]]))
                    divisions_maximum = np.hstack((divisions_maximum,[[lower_bound_quadrupole],[upper_bound_quadrupole]],[[lower_bound_radial],[upper_bound_radial]]))
                    
                    freq1 = np.append(freq1,[freq_quadrupole_chunk,freq_radial_chunk])
                    freq_sig1 = np.append(freq_sig1,[freq_sig_quadrupole_chunk,freq_sig_radial_chunk])
                    angular_degree = np.append(angular_degree,[2,0])
                    order_number = np.append(order_number,[enn_radial-1,enn_radial])
                    sampling_counts = np.append(sampling_counts,[sampling_counts_quadrupole,sampling_counts_radial])
                    sampling_weights = sampling_counts/np.sum(sampling_counts)
                    asef_maximum = np.append(asef_maximum,[asef_maximum_quadrupole,asef_maximum_radial])
                    spsd_maximum = np.append(spsd_maximum,[spsd_maximum_quadrupole,spsd_maximum_radial])

                    # Resort by increasing frequency order (useful if there are frequencies above l=0)
                    tmp_sort = np.argsort(freq1)
                    range_maximum = range_maximum[:,tmp_sort]
                    divisions_maximum = divisions_maximum[:,tmp_sort]
                    
                    freq1 = freq1[tmp_sort]
                    freq_sig1 = freq_sig1[tmp_sort]
                    angular_degree = angular_degree[tmp_sort]
                    order_number = order_number[tmp_sort]
                    sampling_counts = sampling_counts[tmp_sort]
                    sampling_weights = sampling_counts/np.sum(sampling_counts)
                    asef_maximum = asef_maximum[tmp_sort]
                    spsd_maximum = spsd_maximum[tmp_sort]
                    spsd_weights = spsd_maximum/np.sum(spsd_maximum)
                    asef_weights = asef_maximum/np.sum(asef_maximum)
                    n_freq = len(freq1)
                else:
                    # Only the candidate l=2,0 mode frequencies are available
                    range_maximum = np.array([[lower_bound_quadrupole,upper_bound_quadrupole],[lower_bound_radial,upper_bound_radial]]).T
                    divisions_maximum = np.array([[lower_bound_quadrupole,upper_bound_quadrupole],[lower_bound_radial,upper_bound_radial]]).T
                    
                    freq1 = [freq_quadrupole_chunk,freq_radial_chunk]
                    freq_sig1 = [freq_sig_quadrupole_chunk,freq_sig_radial_chunk]
                    angular_degree = [2,0]
                    order_number = [enn_radial-1,enn_radial]
                    sampling_counts = [sampling_counts_quadrupole,sampling_counts_radial]
                    sampling_weights = sampling_counts/np.sum(sampling_counts)
                    asef_maximum = [asef_maximum_quadrupole,asef_maximum_radial]
                    spsd_maximum = [spsd_maximum_quadrupole,spsd_maximum_radial]
                    spsd_weights = spsd_maximum/np.sum(spsd_maximum)
                    asef_weights = asef_maximum/np.sum(asef_maximum)
                    n_freq = len(freq1)
                    
                # Obtain the radial and quadrupole mode indices from frequency array
                quadrupole_index = np.where(angular_degree==2)[0][0]
                radial_index = np.where(angular_degree==0)[0][0]
                    
            # Compute local small spacing d02 for this chunk.
            if flag_quadrupole_found != 0: 
                local_d02 = freq_radial_chunk - freq_quadrupole_chunk
            else:
                local_d02 = d02
                
            actual_d02 = local_d02

            if flag_median_d02_active==1:
                if median_d02 > local_d02:
                    actual_d02 = median_d02
    
            # Update the lower frequency limit for this chunk
            low_cut_frequency2 = freq_radial_chunk - best_dnu*(1.0 + best_alpha*(enn_radial - 0.5 - numax/best_dnu)) + freq_sig_radial
            low_cut_frequency3 = freq_radial_chunk - best_dnu + freq_sig_radial
            if low_cut_frequency2 > low_cut_frequency:
                low_cut_frequency = low_cut_frequency2
            if low_cut_frequency3 > low_cut_frequency:
                tmp_lower_order = np.where(freq1 <= low_cut_frequency3)[0]
                if len(tmp_lower_order) > 0:
                    low_cut_frequency = low_cut_frequency3

        else:
            # ---------------------------------------------------------------------
            # CASE 2: No radial mode is found in the chunk from the global modality
            # ---------------------------------------------------------------------
            
            # Compute the estimated radial mode frequency from the asymptotic relation for radial modes.
            # Assume no quadrupole mode in this case.
            enn_radial = enn_dipole + 1
            
            freq_radial_chunk = astero.asymptotic_relation_radial(enn_radial,numax,[best_dnu,best_epsi,best_alpha])
            fwhm_radial_fit = fwhm_dipole
            order_number = np.zeros(n_freq,dtype=int) + enn_dipole

            local_epsi = 0
            local_d02 = 0
            actual_d02 = d02

            freq_quadrupole_chunk = freq_radial_chunk - d02
            low_cut_frequency = freq_radial_chunk - best_dnu + freq_sig_dipole

        # Force input value for lower limit frequency if required
        if self.cp.low_cut_frequency != 0:
            low_cut_frequency = self.cp.low_cut_frequency

        if self.cp.print_on_screen:
            print(' Lower limit for chunk frequency range: {:.3f} muHz\n'.format(low_cut_frequency))

        # Now check that the previous order l=0 frequency was not selected as a potential dipole mode by the ASEF.
        # If it is selected, then remove it from the sample of identified frequencies.
        # This could happen for RG and SG stars, where a chunk overlapping with the previous one is used during the computation of the fits.

        tmp_below_radial_double = np.where(freq1 <= low_cut_frequency)[0]
        tmp_mask = np.zeros(len(freq1),dtype=bool)
        tmp_mask[tmp_below_radial_double] = True 
        good_freq1_indices = np.arange(len(freq1),dtype=int)[~tmp_mask]
        if len(tmp_below_radial_double) > 0:
            freq1 = freq1[good_freq1_indices]
            freq_sig1 = freq_sig1[good_freq1_indices]
            asef_maximum = asef_maximum[good_freq1_indices]
            spsd_maximum = spsd_maximum[good_freq1_indices]
            spsd_weights = spsd_maximum/np.sum(spsd_maximum)
            asef_weights = asef_maximum/np.sum(asef_maximum)
            sampling_counts = sampling_counts[good_freq1_indices]
            sampling_weights = sampling_counts/np.sum(sampling_counts)
            angular_degree = angular_degree[good_freq1_indices]
            order_number = order_number[good_freq1_indices]
            range_maximum = range_maximum[:,good_freq1_indices]
            divisions_maximum = divisions_maximum[:,good_freq1_indices]
            n_freq = len(freq1)   

        # Finally remove any possible frequency identified above the actual l=0 mode.
        tmp_above_radial = np.where(freq1 > freq_radial_chunk)[0]
        tmp_mask = np.zeros(len(freq1),dtype=bool)
        tmp_mask[tmp_above_radial] = True 
        good_freq1_indices = np.arange(len(freq1),dtype=int)[~tmp_mask]
        if len(tmp_above_radial) > 0:
            freq1 = freq1[good_freq1_indices]
            freq_sig1 = freq_sig1[good_freq1_indices]
            asef_maximum = asef_maximum[good_freq1_indices]
            spsd_maximum = spsd_maximum[good_freq1_indices]
            spsd_weights = spsd_maximum/np.sum(spsd_maximum)
            asef_weights = asef_maximum/np.sum(asef_maximum)
            sampling_counts = sampling_counts[good_freq1_indices]
            sampling_weights = sampling_counts/np.sum(sampling_counts)
            angular_degree = angular_degree[good_freq1_indices]
            order_number = order_number[good_freq1_indices]
            range_maximum = range_maximum[:,good_freq1_indices]
            divisions_maximum = divisions_maximum[:,good_freq1_indices]
            n_freq = len(freq1)

        # Update radial and quadrupole indices.
        # Perform a Lorentzian profile fit to the PSD to assess the FWHM of the radial peak.
        # Consider the largest range around the peak. If no radial peak is present, then
        # take as FWHM of the radial mode that of the adjacent dipole mode.
        if n_radial_chunk != 0:
            radial_index = closest(freq_radial_chunk,freq1,index=True)
            quadrupole_index = closest(freq_quadrupole_chunk,freq1,index=True)
            run_subdir = run + '_radial_fwhm'
            run_names = [run_subdir+str(x) for x in range(self.cp.n_fwhm_fit)]

            if not os.path.isfile(self.star_dir/self.cp.pb_subdir/run_subdir/'0'/'peakbagging_computationParameters.txt') or force:
                right_bound = max([range_maximum[1,radial_index],divisions_maximum[1,radial_index]])
                left_bound = max([range_maximum[0,radial_index],divisions_maximum[0,radial_index]])
        
                tmp0 = np.where((freq >= left_bound) & (freq <= right_bound))
                freq0 = freq[tmp0]
                spsd0 = spsd[tmp0]
                bg_peak = np.mean(bg_level_local[tmp0])
                response_peak = np.mean((np.sin(np.pi/2. * freq0/nyq) / (np.pi/2. * freq0/nyq))**2)
                amplitude_radial = np.sqrt((abs(spsd_maximum[radial_index] - bg_peak)/response_peak)*np.pi*fwhm_radial*self.cp.fwhm_magnification_factor_radial)
                freq_prior_radial = [range_maximum[0,radial_index],range_maximum[1,radial_index]]
                amplitude_prior_radial = [0.0,amplitude_radial]
                data_freq_boundaries = [left_bound,right_bound]
                
                data_range_filenames = np.zeros(self.cp.n_fwhm_fit,dtype='U200')
                prior_filenames = np.zeros(self.cp.n_fwhm_fit,dtype='U200')
  
                # Make sure to enlarge FWHM prior if fit fails.
                flag_computation_completed = 0
                iterations = 1
        
                left_fwhm = self.cp.fwhm_lower_bound

                if self.cp.print_on_screen:
                    print('\n Fitting Linewidth of radial mode of the chunk.')

                while flag_computation_completed != 1:
                    # Allow for a very narrow FWHM
                    fwhm_prior = [left_fwhm,fwhm_radial*self.cp.fwhm_magnification_factor_radial*iterations]
                    if best_dnu < self.cp.dnu_threshold:
                        boundaries = [freq_prior_radial,amplitude_prior_radial,fwhm_prior]
                    else:
                        amplitude_quadrupole = np.sqrt((abs(spsd_maximum[quadrupole_index] - bg_peak)/response_peak)*np.pi*fwhm_radial*self.cp.fwhm_magnification_factor_radial)
                        freq_prior_quadrupole = [range_maximum[0,quadrupole_index],range_maximum[1,quadrupole_index]]
                        amplitude_prior_quadrupole = [0.0,amplitude_quadrupole]
                        boundaries = [freq_prior_quadrupole,amplitude_prior_quadrupole,fwhm_prior,freq_prior_radial,amplitude_prior_radial,fwhm_prior]

                    for k in range(0,self.cp.n_fwhm_fit):
                        data_range_filenames[k] = self.star_dir/self.cp.pb_subdir/('frequencyRange_' + run_names[k] + '.txt')
                        prior_filenames[k] = self.star_dir/self.cp.pb_subdir/(self.cp.prior_filename + '_' + run_names[k] + '.txt')

                        diamonds.write_data_range(data_range_filenames[k],data_freq_boundaries)
                        diamonds.write_uniform_prior(prior_filenames[k],boundaries)
            
                    peakbagging_parameters = {'subdir':          self.cp.pb_subdir,
                                              'run':             run_names,
                                              'background':      background_name,    
                                              'fwhm':            -1.0,
                                              'filename_run':    run_subdir}          

                    flag_computation_completed_array = diamonds.run_peakbagging(self.catalog_id,self.star_id,peakbagging_parameters, 2, 0, 0, self.cp.dp_pb, self.cp.diamonds_path, self.cp.n_threads, self.cp.prior_filename)
                    flag_computation_completed = max(flag_computation_completed_array)

                    iterations += 1

                if self.cp.save_test_files != 1:
                    for k in range(0, self.cp.n_fwhm_fit):
                        os.remove(data_range_filenames[k])
                        os.remove(prior_filenames[k])
                        

            fwhm_radial_fit_array = np.zeros(self.cp.n_fwhm_fit)
            for k in range(0, self.cp.n_fwhm_fit):
                filenames = glob.glob(str(self.star_dir/self.cp.pb_subdir/run_names[k]/'peakbagging_parameter0*txt'))
        
                if len(filenames)==3:
                    fwhm_parameter = '002'
                else:
                    fwhm_parameter = '005'
        
                # Read the sampled FWHM of the radial mode.
                par_fwhm0 = np.loadtxt(self.star_dir/self.cp.pb_subdir/run_names[k]/('peakbagging_parameter' + fwhm_parameter + '.txt'))
                # Load the posterior samples to compute Bayesian mean estimate for FWHM_0
                post = np.loadtxt(self.star_dir/self.cp.pb_subdir/run_names[k]/'peakbagging_posteriorDistribution.txt')
                # Compute parameter estimate, using a weighted average
                post /= max(post)
                fwhm_radial_fit_array[k] = np.sum(par_fwhm0*post)/np.sum(post)
    

            fwhm_radial_fit = np.median(fwhm_radial_fit_array)
    
            if self.cp.print_on_screen:
                print(' FWHM (l=0) = {:.3f} muHz \n'.format(fwhm_radial_fit))


        # Define a detection probability array for each frequency peak. Set it to -99, meaning all
        # frequency peaks are not tested (i.e. considered detected at the:ning, with p_detection = 100 %).
        # In the peak blending test, save whether the peak has a blending or not (blending_profile_flag either 1 or 0).
        # In the sinc^2 test (Sinc^2 profile vs. Lorentzian profile), save which profile is better to model (sinc_profile_flag either 1 or 0).
        # Define a rotation probability array, specifying the probability that the peak is split
        # by rotation. Assume a default no rotation test performed (set to -99), meaning that rotation is assumed not detected 
        # at the:ning (with p_rotation = 0 %). 
        # Define a duplet probability array, specifying whether the peak is a duplet or not, and assume no duplicity at 
        # the:ning (set to -99), meaning that the peak is a single peak (p_duplet = 0 %). 

        detection_probability = np.zeros(len(freq1)) - 99.0
        rotation_probability = np.zeros(len(freq1)) - 99.0
        duplicity_probability = np.zeros(len(freq1)) - 99.0
        sinc_profile_flag = np.zeros(len(freq1),dtype=int)
        blending_profile_flag = np.zeros(len(freq1),dtype=int)

        if force:
            delete_filenames = glob.glob(str(self.star_dir/self.cp.pb_subdir/('detectionProbability_'+ run + '*.txt')))
            if len(delete_filenames) > 0:
                tmp = [os.remove(delete_filename) for delete_filename in delete_filenames]
    
            delete_filenames = glob.glob(str(self.star_dir/self.cp.pb_subdir/('rotationProbability_'+ run + '*.txt')))
            if len(delete_filenames) > 0:
                tmp = [os.remove(delete_filename) for delete_filename in delete_filenames]
 
        # Skip this part if no radial mode is found
        if n_radial_chunk != 0:
            # Peak significance and blending test for l=2,0

            # Perform the peak significance test on each frequency peak found to have a critical SNR
            # Evaluate SNR using an estimate for the amplitude of the peak
            # Clear previous peak detection probabilities if fit is forced
            right_bound = range_maximum[1,radial_index]
            left_bound = range_maximum[0,quadrupole_index]
            tmp_freq_peak = np.where((freq <= right_bound) & (freq >= left_bound))[0]
            bg_peak = np.mean(bg_level_local[tmp_freq_peak])
            response_peak = np.mean((np.sin(np.pi/2. * freq[tmp_freq_peak]/nyq) / (np.pi/2. * freq[tmp_freq_peak]/nyq))**2)
            height_ratio_quadrupole = spsd_maximum[quadrupole_index]/(bg_peak*self.cp.height_ratio_threshold)
            height_ratio_radial = spsd_maximum[radial_index]/(bg_peak*self.cp.height_ratio_threshold)

            # Distinguish between RG + late SG stars, and MS + early SG stars (especially hot stars, where linewidths are larger).
            if (best_dnu < self.cp.dnu_threshold):
                # This is the case of either RG or late SG.
                # Here perform a standard significance test only for the l=2,0 pair, but separately
                # for each mode. Assume Lorentzian profiles only.
                # For the quadrupole mode adopt a larger linewidth upper limit for the prior. This is because
                # l=2 modes are generally affected by the presence of l=2 mixed modes in evolved stars.
                fwhm_quadrupole = fwhm_radial_fit*self.cp.fwhm_magnification_factor_quadrupole
                fwhm_radial = fwhm_radial_fit*self.cp.fwhm_magnification_factor_radial
                fwhm = [fwhm_quadrupole,fwhm_radial]
                freq_left = [freq1[quadrupole_index] - freq_sig1[quadrupole_index], freq1[radial_index] - freq_sig1[radial_index]]
                freq_right = [freq1[quadrupole_index] + freq_sig1[quadrupole_index], freq1[radial_index] + freq_sig1[radial_index]]
                amplitude_quadrupole = np.sqrt((abs(spsd_maximum[quadrupole_index] - bg_peak)/response_peak)*np.pi*fwhm_quadrupole)
                amplitude_radial = np.sqrt((abs(spsd_maximum[radial_index] - bg_peak)/response_peak)*np.pi*fwhm_radial)
                amplitude = [amplitude_quadrupole,amplitude_radial]
                height_ratio = [height_ratio_quadrupole,height_ratio_radial]
      
                # For each peak consider the maximum data range between the range and the division    
                left_bound_quadrupole = min([range_maximum[0,quadrupole_index],divisions_maximum[0,quadrupole_index]]) 
                right_bound_quadrupole = max([range_maximum[1,quadrupole_index],divisions_maximum[1,quadrupole_index]]) 
                left_bound_radial = min([range_maximum[0,radial_index],divisions_maximum[0,radial_index]]) 
                right_bound_radial = max([range_maximum[1,radial_index],divisions_maximum[1,radial_index]]) 

                range_quadrupole_radial = [[left_bound_quadrupole,right_bound_quadrupole],[left_bound_radial,right_bound_radial]]
                angular_degree_quadrupole_radial = [angular_degree[quadrupole_index],angular_degree[radial_index]]
                quadrupole_radial_index = [quadrupole_index,radial_index]
        
                test_dirs = np.zeros((2,2),dtype='U150')
                prior_filenames = np.zeros((2,2),dtype='U150')
                data_range_filenames = np.zeros((2,2),dtype='U150')
                probability_filenames = np.zeros(2,dtype='U150')
                run_test = []
                flag_run_test = 0

                for i in range(0, len(height_ratio)):
                    if height_ratio[i] < 1.0:
                        peak_number = str(int(quadrupole_radial_index[i]))
                        test_name = run+'_'+peak_number 
                
                        test_nameA = test_name + 'A'           # Only background
                        test_nameB = test_name + 'B'           # Lorentzian profile and background
                        test_dirA = self.star_dir/self.cp.pb_subdir/test_nameA
                        test_dirB = self.star_dir/self.cp.pb_subdir/test_nameB
                        filenameA = self.star_dir/self.cp.pb_subdir/(self.cp.prior_filename + '_' + test_nameA + '.txt')
                        filenameB = self.star_dir/self.cp.pb_subdir/(self.cp.prior_filename + '_' + test_nameB + '.txt')
                        filename_rangeA = self.star_dir/self.cp.pb_subdir/('frequencyRange_' + test_nameA + '.txt')
                        filename_rangeB = self.star_dir/self.cp.pb_subdir/('frequencyRange_' + test_nameB + '.txt')
                        probability_filename = self.star_dir/self.cp.pb_subdir/('detectionProbability_' + test_name + '.txt')

                        test_dirs[:,i] = [test_dirA,test_dirB]
                        probability_filenames[i] = probability_filename
                        prior_filenames[:,i] = [filenameA,filenameB]
                        data_range_filenames[:,i] = [filename_rangeA,filename_rangeB]

                        if not os.path.isfile(probability_filename):
                            # Set up priors for the peak test profile
                            data_freq_boundaries = [range_quadrupole_radial[i][0],range_quadrupole_radial[i][1]]
                            freq_prior = [freq_left[i],freq_right[i]]
                            amplitude_prior = [0.0,amplitude[i]]
                            left_fwhm = self.cp.fwhm_lower_bound
                            fwhm_prior = [left_fwhm,fwhm[i]]
                            bkg_prior = [self.cp.bkg_prior_lower,self.cp.bkg_prior_upper]
                            boundariesA = [bkg_prior]
                            boundariesB = [freq_prior,amplitude_prior,fwhm_prior,bkg_prior]
                            
                            diamonds.write_data_range(data_range_filenames[0,i],data_freq_boundaries)
                            diamonds.write_data_range(data_range_filenames[1,i],data_freq_boundaries)
                            diamonds.write_uniform_prior(prior_filenames[0,i],boundariesA)
                            diamonds.write_uniform_prior(prior_filenames[1,i],boundariesB)

                            run_test.extend([test_nameA,test_nameB])
                            flag_run_test += 1

                tested_peak_indices = np.where(np.array(height_ratio) < 1.0)[0]
                if len(tested_peak_indices) > 0:
                    if self.cp.print_on_screen:
                        print('\n ---------------------------------------------')
                        print(' Peak Detection Test for l=2,0 pair.' )

                    if flag_run_test > 0:
                        p_BA = np.zeros((self.cp.n_peak_test,len(tested_peak_indices)))
                        
                        for k in range(0, self.cp.n_peak_test):
                            peakbagging_parameters = {'subdir':       self.cp.pb_subdir,
                                                      'run':          run_test,
                                                      'background':   background_name,
                                                      'fwhm':         -1.0,
                                                      'filename_run': test_name}
                
                            flag_computation_completed = diamonds.run_peakbagging(self.catalog_id,self.star_id,peakbagging_parameters,1,0,0, self.cp.dp_pb, self.cp.diamonds_path, self.cp.n_threads, self.cp.prior_filename)
                    
                            if self.cp.print_on_screen:
                                print(' Iteration #{} completed.'.format(k))
                    
                            for i in range(0,len(tested_peak_indices)):
                                ln_evidA = np.loadtxt(test_dirs[0,tested_peak_indices[i]]+'/peakbagging_evidenceInformation.txt',usecols=(0,))
                                ln_evidB = np.loadtxt(test_dirs[1,tested_peak_indices[i]]+'/peakbagging_evidenceInformation.txt',usecols=(0,))                       
                                if (ln_evidA >= ln_evidB):
                                    ln_evidAB = ln_evidA + np.log(1.+np.exp(ln_evidB-ln_evidA)) 
                                else:
                                    ln_evidAB = ln_evidB + np.log(1.+np.exp(ln_evidA-ln_evidB))
                    
                                p_BA[k,i] = np.exp(ln_evidB - ln_evidAB)
                    
                        for i in range(0, len(tested_peak_indices)):
                            with open(probability_filenames[tested_peak_indices[i]],'w') as f:
                                f.write('# Detection probabilities for l={} at {:.2f} muHz\n'.format(angular_degree_quadrupole_radial[tested_peak_indices[i]],freq1[quadrupole_radial_index[tested_peak_indices[i]]])) 
                                f.write('# Each line corresponds to a different run of the same test.\n')
                                f.write('# Col 1: p(BA).\n')
                                for val in p_BA[:,i]:
                                    f.write('{:.4f}\n'.format(val))

                            if self.cp.save_test_files==0:
                                # Remove test folders and prior files if required
                                shutil.rmtree(test_dirs[0,tested_peak_indices[i]])
                                shutil.rmtree(test_dirs[1,tested_peak_indices[i]])
                                os.remove(data_range_filenames[0,tested_peak_indices[i]])
                                os.remove(data_range_filenames[1,tested_peak_indices[i]])
                                os.remove(prior_filenames[0,tested_peak_indices[i]])
                                os.remove(prior_filenames[1,tested_peak_indices[i]])
                 
                    for i in range(0,len(tested_peak_indices)):
                        p_BA = np.loadtxt(probability_filenames[tested_peak_indices[i]])
                        # Consider the maximum probability among the different runs. Approximate up to third decimal digit.
                        max_p_BA = round(max(p_BA),3)

                        if self.cp.print_on_screen:
                            print('\n Peak detection test for l={} at {:.2f} muHz'.format(angular_degree_quadrupole_radial[tested_peak_indices[i]],freq1[quadrupole_radial_index[tested_peak_indices[i]]]))
                            print(' P (One Lorentzian vs Only Background): {} \n'.format(max_p_BA))
                           
                        detection_probability[quadrupole_radial_index[tested_peak_indices[i]]] = max_p_BA
            

                    if self.cp.print_on_screen:
                        print(' Peak Detection Test for l=2,0 pair Completed.' )
                        print(' ---------------------------------------------\n')
            else:
                # This is the case of MS and early SG stars.
                # Here perform a blending test on top of the peak significance test, only for the l=2,0 pair.
                # This is because it is likely that the l=2,0 modes are blended due to the large linewidths in hotter stars.
                # Testing l=2 and l=0 only separately from one another may not represent a realistic condition.
                fwhm_quadrupole = fwhm_radial_fit*self.cp.fwhm_magnification_factor_radial
                fwhm_radial = fwhm_radial_fit*self.cp.fwhm_magnification_factor_radial
                amplitude_quadrupole = np.sqrt((abs(spsd_maximum[quadrupole_index] - bg_peak)/response_peak)*np.pi*fwhm_quadrupole)
                amplitude_radial = np.sqrt((abs(spsd_maximum[radial_index] - bg_peak)/response_peak)*np.pi*fwhm_radial)
                amplitude_max = max([amplitude_quadrupole,amplitude_radial])

                if (height_ratio_quadrupole < 1.0) or (height_ratio_radial < 1.0):
                    peak_number = str(int(radial_index))
                    test_name = run + '_' + peak_number 
                    test_nameA = test_name + 'A'           # Only background
                    test_nameB = test_name + 'B'           # Lorentzian profile and background
                    test_nameD = test_name + 'D'           # Two Lorentzian profiles and background
                    test_dirA = self.star_dir/self.cp.pb_subdir/test_nameA
                    test_dirB = self.star_dir/self.cp.pb_subdir/test_nameB
                    test_dirD = self.star_dir/self.cp.pb_subdir/test_nameD        
                    filenameA = self.star_dir/self.cp.pb_subdir/(self.cp.prior_filename + '_' + test_nameA + '.txt')
                    filenameB = self.star_dir/self.cp.pb_subdir/(self.cp.prior_filename + '_' + test_nameB + '.txt')
                    filenameD = self.star_dir/self.cp.pb_subdir/(self.cp.prior_filename + '_' + test_nameD + '.txt')
                    filename_rangeA = self.star_dir/self.cp.pb_subdir/('frequencyRange_' + test_nameA + '.txt')
                    filename_rangeB = self.star_dir/self.cp.pb_subdir/('frequencyRange_' + test_nameB + '.txt')
                    filename_rangeD = self.star_dir/self.cp.pb_subdir/('frequencyRange_' + test_nameD + '.txt')

                    probability_filename = self.star_dir/self.cp.pb_subdir/('detectionProbability_' + test_name + '.txt')
                    
                    if self.cp.print_on_screen:
                        print('\n ----------------------------------------------------------')
                        print(' Peak Detection and Blending Test for l=2,0 pair.' )

                    if not os.path.isfile(probability_filename):
                        # Set up priors for the peak test profile
                        # For each peak consider the maximum data range between the range and the division 
                        left_bound_quadrupole = min([range_maximum[0,quadrupole_index],divisions_maximum[0,quadrupole_index]])
                        right_bound_radial = max([range_maximum[1,radial_index],divisions_maximum[1,radial_index]]) 
                        data_freq_boundaries = [left_bound_quadrupole,right_bound_radial]
                
                        freq_prior_quadrupole = [freq1[quadrupole_index]-freq_sig1[quadrupole_index],freq1[quadrupole_index]+freq_sig1[quadrupole_index]]
                        freq_prior_radial = [freq1[radial_index]-freq_sig1[radial_index],freq1[radial_index]+freq_sig1[radial_index]]
                        freq_prior = [range_maximum[0,quadrupole_index],range_maximum[1,radial_index]]
                
                        amplitude_prior_quadrupole = [0.0,amplitude_quadrupole]
                        amplitude_prior_radial = [0.0,amplitude_radial]
                        amplitude_prior = [0.0,amplitude_max]
               
                        right_fwhm = fwhm_radial_fit*self.cp.fwhm_magnification_factor_radial
                        left_fwhm = self.cp.fwhm_lower_bound
                        fwhm_prior = [left_fwhm,right_fwhm]
                        bkg_prior = [self.cp.bkg_prior_lower,self.cp.bkg_prior_upper]
                
                        boundariesA = [bkg_prior]
                        boundariesB = [freq_prior,amplitude_prior,fwhm_prior,bkg_prior]
                        boundariesD = [freq_prior_quadrupole,amplitude_prior_quadrupole,fwhm_prior,freq_prior_radial,amplitude_prior_radial,fwhm_prior,bkg_prior]
                
                        diamonds.write_data_range(filename_rangeA,data_freq_boundaries)
                        diamonds.write_data_range(filename_rangeB,data_freq_boundaries)
                        diamonds.write_data_range(filename_rangeD,data_freq_boundaries)
                        diamonds.write_uniform_prior(filenameA,boundariesA)
                        diamonds.write_uniform_prior(filenameB,boundariesB)
                        diamonds.write_uniform_prior(filenameD,boundariesD)
    
                        p_BA = np.zeros(self.cp.n_peak_test)
                        p_DA = np.zeros(self.cp.n_peak_test)
                        p_DB = np.zeros(self.cp.n_peak_test)
                        nu0 = np.zeros(self.cp.n_peak_test)
                
                        for k in range(0, self.cp.n_peak_test):
                            peakbagging_parameters = {'subdir':       self.cp.pb_subdir,
                                                      'run':          [test_nameA,test_nameB,test_nameD],
                                                      'background':   background_name,
                                                      'fwhm':         -1.0,
                                                      'filename_run': test_name}

                            flag_computation_completed = diamonds.run_peakbagging(self.catalog_id,self.star_id,peakbagging_parameters,1,0,0,  self.cp.dp_pb, self.cp.diamonds_path, self.cp.n_threads, self.cp.prior_filename)
                    
                            if self.cp.print_on_screen:
                                print(' Iteration #{} completed.'.format(k))
                    
                            ln_evidA = np.loadtxt(test_dirA/'peakbagging_evidenceInformation.txt', usecols=(0,))
                            ln_evidB = np.loadtxt(test_dirB/'peakbagging_evidenceInformation.txt',usecols=(0,))
                            ln_evidD = np.loadtxt(test_dirD/'peakbagging_evidenceInformation.txt',usecols=(0,))
                        
                            if (ln_evidA >= ln_evidB):
                                ln_evidAB = ln_evidA + np.log(1.+np.exp(ln_evidB-ln_evidA)) 
                            else:
                                ln_evidAB = ln_evidB + np.log(1.+np.exp(ln_evidA-ln_evidB))
                    
                            if (ln_evidA >= ln_evidD):
                                ln_evidAD = ln_evidA + np.log(1.+np.exp(ln_evidD-ln_evidA)) 
                            else:
                                ln_evidAD = ln_evidD + np.log(1.+np.exp(ln_evidA-ln_evidD))
                    
                            if (ln_evidB >= ln_evidD):
                                ln_evidBD = ln_evidB + np.log(1.+np.exp(ln_evidD-ln_evidB)) 
                            else:
                                ln_evidBD = ln_evidD + np.log(1.+np.exp(ln_evidB-ln_evidD))

                            p_BA[k] = np.exp(ln_evidB - ln_evidAB)
                            p_DA[k] = np.exp(ln_evidD - ln_evidAD)
                            p_DB[k] = np.exp(ln_evidD - ln_evidBD)
                    
                            # Read the estimated central frequency from model B.
                            par_nu0 = np.loadtxt(self.star_dir/self.cp.pb_subdir/test_nameB/'peakbagging_parameter000.txt')                    
                            # Load the posterior samples to compute Bayesian mean estimate for nu0                 
                            post = np.loadtxt(self.star_dir/self.cp.pb_subdir/test_nameB/'peakbagging_posteriorDistribution.txt')      
                            # Compute parameter estimate, using a weighted average           
                            post /= max(post)
                            nu0[k] = np.sum(par_nu0*post)/np.sum(post)
                        
                        with open(probability_filename,'w') as f:
                            f.write('# Detection probabilities for l=2,0 duplet at {:.2f} + {:.2f} muHz\n'.format(freq1[quadrupole_index],freq1[radial_index]))
                            f.write('# Each line corresponds to a different run of the same test.\n')
                            f.write('# Col 1: p(BA), Col 2: p(DA), Col 3: p(DB), Col 4: nu0 (microHz).\n')
                            for l in range(0,self.cp.n_peak_test):
                                f.write('{:.4f}  {:10.4f}  {:10.4f}  {:15.4f}\n'.format(p_BA[l],p_DA[l],p_DB[l],nu0[l]))
                            
                        if self.cp.save_test_files==0:
                            # Remove test folders and prior files if required
                            shutil.rmtree(test_dirA)
                            shutil.rmtree(test_dirB)
                            shutil.rmtree(test_dirD)
                            os.remove(filename_rangeA)
                            os.remove(filename_rangeB)
                            os.remove(filename_rangeD)
                            os.remove(filenameA)
                            os.remove(filenameB)
                            os.remove(filenameD)
                
                    else:
                        p_BA,p_DA,p_DB,nu0 = np.loadtxt(probability_filename,unpack=True)
                        # Consider the maximum probabilities among the different runs. Approximate up to third decimal digit.
                        
                    max_p_BA = round(max(p_BA),3)
                    max_p_DA = round(max(p_DA),3)
                    max_p_DB = round(max(p_DB),3)

                    if self.cp.print_on_screen:
                        print('\n Peak detection test for l=2,0 duplet at {:.2f} + {:.2f} muHz\n'.format(freq1[quadrupole_index],freq1[radial_index]))
                        print(' P (One Lorentzian vs Only Background): ',max_p_BA)
                        print(' P (Two Lorentzians vs Only Background): ',max_p_DA)
                        print(' P (Two Lorentzians vs One Lorentzian): \n',max_p_DB)
                
                    if max_p_DB > 0.5:
                        # Here we have a peak blending. The two peaks are treated together and are considered
                        # detected only if their p_DA >= 0.993.
                        detection_probability[quadrupole_index] = max_p_DA
                        detection_probability[radial_index] = max_p_DA
                        blending_profile_flag[quadrupole_index] = 1
                        blending_profile_flag[radial_index] = 1
                    else:
                        detection_probability[quadrupole_index] = 0.0
                        detection_probability[radial_index] = 0.0
                
                        if (max_p_BA >= self.cp.detection_probability_threshold):
                            # Here we have no peak blending, and only one peak can be considered. 
                            # Find which one between l=2 and l=0 has to be deeemd significant.
                            avg_nu0 = np.mean(nu0)
                            detected_index = closest(avg_nu0,[freq1[quadrupole_index],freq1[radial_index]],index=True)
                            detection_probability[detected_index+quadrupole_index] = max_p_BA
                
                    if self.cp.print_on_screen:
                        print(' Peak Detection and Blending Test for l=2,0 pair Completed.') 
                        print(' ----------------------------------------------------------\n')
                        

        # --------------------------------------------------------------
        # Dipole mode(s) and Octupole mode
        # --------------------------------------------------------------

        # At this stage, check for the presence of dipole mode(s), and subsequently l=3.
        # Distinguish red giant stars from subgiants and main sequence (White et al. 2011). 
        # If the star is hot (F-type, Teff > 6300 K), its Dnu will be smaller as compared to G-type star in same relative evolutionary stage.
        # Then it is less likely to have an F-type star that is also a subgiant, hence showing mixed modes.
        # For this reason, also incorporate cooler stars with higher Dnu as possible subgiants (see e.g. Appourchaux et al. 2012).
        # If no dipole modes were identified from the global fit, still give it a try in the chunk modality
        # and check whether there could be a potential dipole mode. This is assuming that l=0 was found!
        # Start by defining the l=3 search region
        if best_dnu < self.cp.dnu_rg:
            octupole_freq_asymp = freq_radial_chunk - best_dnu/2. - d03
            octupole_freq_lower = freq_radial_chunk - best_dnu/2. - d03*self.cp.d03_upper_scaling_factor
            octupole_freq_upper = freq_radial_chunk - best_dnu/2. - d03*self.cp.d03_lower_scaling_factor
        else:
            # Adopt approximated asymptotic relation by Bedding & Kjeldsen 2003 to locate l=3 region in less evolved stars.
            if (best_dnu < self.cp.dnu_sg and teff < self.cp.teff_sg):
                d02_upper_scaling_factor = self.cp.d02_upper_scaling_factor_sg
                d02_lower_scaling_factor = self.cp.d02_lower_scaling_factor_sg
            else:
                d02_upper_scaling_factor = self.cp.d02_upper_scaling_factor_ms
                d02_lower_scaling_factor = self.cp.d02_lower_scaling_factor_ms
    
            octupole_freq_lower = freq_radial_chunk - best_dnu/2. - d02_upper_scaling_factor*actual_d02
            octupole_freq_upper = freq_radial_chunk - best_dnu/2. - d02_lower_scaling_factor*actual_d02
            octupole_freq_asymp = (octupole_freq_lower + octupole_freq_upper)/2.

        if (best_dnu < self.cp.dnu_sg) & (teff < self.cp.teff_sg):
            dipole_indices = np.where(angular_degree==1)[0]
        else:
            # For MS stars, only check the modes below nu0 - Dnu/2 + freq_dipole_sig, and above the lower limit of the l=3 region. 
            # All remaining ones are flagged as undetected. This will speed up the mode identification process.
            if n_radial_chunk != 0:
                upper_dipole_limit = freq_radial_chunk - best_dnu/2.0 + freq_sig_radial
            else:
                upper_dipole_limit = freq_radial_chunk - best_dnu/2.0 + freq_sig_dipole
    
            if upper_dipole_limit > freq_radial_chunk - best_dnu/4.0:
                upper_dipole_limit = freq_radial_chunk - best_dnu/4.0
    
            dipole_indices = np.where((freq1 < upper_dipole_limit) & (freq1 >= octupole_freq_lower) & (angular_degree==1))[0]
            bad_dipole_indices = np.where(((freq1 >= upper_dipole_limit) | (freq1 < octupole_freq_lower)) & (angular_degree==1))[0]
    
            if len(bad_dipole_indices) > 0:
                detection_probability[bad_dipole_indices] = 0.0

        if len(dipole_indices) > 0:
            # Peak significance test for candidate l=1 mode(s)
            
            # Perform a peak significance test in order to identify all the significant peaks marked as l=1.
            # If the star is a RG, add up the test for sinc^2 profile to account for possible unresolved mixed modes.
            height_ratio_array = np.zeros(len(dipole_indices))
            test_dirs = np.zeros((3,len(dipole_indices)),dtype='U200')
            prior_filenames = np.zeros((3,len(dipole_indices)),dtype='U200')
            data_range_filenames = np.zeros((3,len(dipole_indices)),dtype='U200')
            probability_filenames = np.zeros(len(dipole_indices),dtype='U200')
            run_test = []
            flag_run_test = 0

            for i in range(0, len(dipole_indices)):
                # Take smallest range possible that is available to build the frequency prior of the given peak. This will avoid that the peak may result not
                # significant if the parameter space is too large.
                # For each peak consider the maximum data range between the range and the division  
                freq_index = dipole_indices[i]
                left_bound = min([range_maximum[0,freq_index],divisions_maximum[0,freq_index]])
                right_bound = max([range_maximum[1,freq_index],divisions_maximum[1,freq_index]])

                # For RG and SG stars, adopt a narrower FWHM prior for those l=1 modes outside the l=3 region.
                # This will speed up the computation of the peak test by reducing the parameter space.
        
                right_fwhm = fwhm_radial_fit*self.cp.fwhm_magnification_factor_radial*2
                left_fwhm = self.cp.fwhm_lower_bound

                if (best_dnu < self.cp.dnu_sg) & (teff < self.cp.teff_sg):
                    right_fwhm = fwhm_radial_fit*self.cp.fwhm_magnification_factor_radial
            
                    # Allow for a very narrow FWHM accounting for unresolved mixed modes
                    if (freq1[freq_index] >= octupole_freq_upper) or (freq1[freq_index] <= octupole_freq_lower):
                        right_fwhm = fwhm_radial_fit*self.cp.fwhm_magnification_factor_dipole
            
                    # Do not exceed the width given by the frequency range of the peak if the star is a RG, because l=1 modes are narrow. 
                    if best_dnu < self.cp.dnu_rg:
                        range_peak = abs(range_maximum[1,freq_index] - range_maximum[0,freq_index])
                        if right_fwhm > range_peak:
                            right_fwhm = range_peak

                tmp_freq_peak = np.where((freq <= right_bound) & (freq >= left_bound))[0]
                bg_peak = np.mean(bg_level_local[tmp_freq_peak])
                response_peak = np.mean((np.sin(np.pi/2. * freq[tmp_freq_peak]/nyq) / (np.pi/2. * freq[tmp_freq_peak]/nyq))**2)
                amplitude = np.sqrt((abs(max(spsd[tmp_freq_peak]) - bg_peak)/response_peak)*np.pi*right_fwhm)
                height_ratio = spsd_maximum[freq_index]/(bg_peak*self.cp.height_ratio_threshold)
                height_ratio_array[i] = height_ratio
                        
                peak_number = str(int(freq_index))
                test_name = run + '_' + peak_number
                test_nameA = test_name + 'A'           # Only background
                test_nameB = test_name + 'B'           # Background and Lorentzian profile
                test_nameC = test_name + 'C'           # Background and Sinc^2 profile
                test_dirA = self.star_dir/self.cp.pb_subdir/test_nameA
                test_dirB = self.star_dir/self.cp.pb_subdir/test_nameB
                test_dirC = self.star_dir/self.cp.pb_subdir/test_nameC
                probability_filename = self.star_dir/self.cp.pb_subdir/('detectionProbability_' + test_name + '.txt')    
                filenameA = self.star_dir/self.cp.pb_subdir/(self.cp.prior_filename + '_' + test_nameA + '.txt')
                filenameB = self.star_dir/self.cp.pb_subdir/(self.cp.prior_filename + '_' + test_nameB + '.txt')
                filenameC = self.star_dir/self.cp.pb_subdir/(self.cp.prior_filename + '_' + test_nameC + '.txt')

                filename_rangeA = self.star_dir/self.cp.pb_subdir/('frequencyRange_' + test_nameA + '.txt')
                filename_rangeB = self.star_dir/self.cp.pb_subdir/('frequencyRange_' + test_nameB + '.txt')
                filename_rangeC = self.star_dir/self.cp.pb_subdir/('frequencyRange_' + test_nameC + '.txt')
                test_dirs[:,i] = [test_dirA,test_dirB,test_dirC]
                probability_filenames[i] = probability_filename
                prior_filenames[:,i] = [filenameA,filenameB,filenameC]
                data_range_filenames[:,i] = [filename_rangeA,filename_rangeB,filename_rangeC]
        
                if height_ratio < 1.0:
                    if not os.path.isfile(probability_filename):
                        # Set up priors for the peak test profile
                        data_freq_boundaries = [left_bound,right_bound] 
                        freq_prior = [freq1[freq_index] - freq_sig1[freq_index],freq1[freq_index] + freq_sig1[freq_index]]
                        amplitude_prior = [0.0,amplitude]
                        fwhm_prior = [left_fwhm,right_fwhm]
                        bkg_prior = [self.cp.bkg_prior_lower,self.cp.bkg_prior_upper]
                        height_prior = [0.0,spsd_maximum[freq_index]*1.5]
                        sinc_prior = [-1,-1]
               
                        boundariesA = [bkg_prior]
                        boundariesB = [freq_prior,amplitude_prior,fwhm_prior,bkg_prior]
                
                        diamonds.write_data_range(data_range_filenames[0,i],data_freq_boundaries)
                        diamonds.write_data_range(data_range_filenames[1,i],data_freq_boundaries)
                        diamonds.write_uniform_prior(prior_filenames[0,i],boundariesA)
                        diamonds.write_uniform_prior(prior_filenames[1,i],boundariesB)
                
                        if best_dnu < self.cp.dnu_rg:
                            boundariesC = [freq_prior,height_prior,bkg_prior,sinc_prior]
                            diamonds.write_data_range(data_range_filenames[2,i],data_freq_boundaries)
                            diamonds.write_uniform_prior(prior_filenames[2,i],boundariesC)
                            run_test.extend([test_nameA,test_nameB,test_nameC])
                        else:
                            run_test.extend([test_nameA,test_nameB])
                        flag_run_test +=1

            tested_peak_indices = np.where(height_ratio_array < 1.0)[0]
            if len(tested_peak_indices) > 0:
                if self.cp.print_on_screen:
                    print('\n --------------------------------------------------------')
                    print(' Peak Detection Test for candidate l=1 mode(s).')
 
                if flag_run_test > 0:
                    if self.cp.print_on_screen:
                        print(' A total of {} tests must be performed.'.format(len(run_test))) 
            
                    fwhm_detection_fit = np.zeros((self.cp.n_peak_test,len(tested_peak_indices)))
                    p_BA = np.zeros((self.cp.n_peak_test,len(tested_peak_indices)))
            
                    if best_dnu < self.cp.dnu_rg:
                        p_CA = np.zeros((self.cp.n_peak_test,len(tested_peak_indices)))
                        p_CB = np.zeros((self.cp.n_peak_test,len(tested_peak_indices)))
            
                    for k in range(0, self.cp.n_peak_test):
                        peakbagging_parameters = {'subdir':       self.cp.pb_subdir,
                                                  'run':          run_test,
                                                  'background':   background_name,
                                                  'fwhm':         -1.0,
                                                  'filename_run': test_name}
                
                        flag_computation_completed = diamonds.run_peakbagging(self.catalog_id,self.star_id,peakbagging_parameters,1,0,0, self.cp.dp_pb, self.cp.diamonds_path, self.cp.n_threads, self.cp.prior_filename)
                
                        if self.cp.print_on_screen:
                            print(' Iteration #{} completed.'.format(k))
 
                        for i in range (0, len(dipole_indices[tested_peak_indices])):
                            # Compute the estimate for FWHM from model B.
                            par_fwhm = np.loadtxt(test_dirs[1,tested_peak_indices[i]]+'/peakbagging_parameter002.txt')
                            post = np.loadtxt(test_dirs[1,tested_peak_indices[i]]+'/peakbagging_posteriorDistribution.txt')
                            post /= max(post)
                            fwhm_detection_fit[k,i] = np.sum(par_fwhm*post)/np.sum(post)
                        
                            ln_evidA = np.loadtxt(test_dirs[0,tested_peak_indices[i]]+'/peakbagging_evidenceInformation.txt', usecols=(0,))
                            ln_evidB = np.loadtxt(test_dirs[1,tested_peak_indices[i]]+'/peakbagging_evidenceInformation.txt', usecols=(0,))
                    
                            if (ln_evidA >= ln_evidB):
                                ln_evidAB = ln_evidA + np.log(1+np.exp(ln_evidB-ln_evidA)) 
                            else:
                                ln_evidAB = ln_evidB + np.log(1+np.exp(ln_evidA-ln_evidB))
                    
                            p_BA[k,i] = np.exp(ln_evidB - ln_evidAB)
                    
                            if best_dnu < self.cp.dnu_rg:
                                ln_evidC = np.loadtxt(test_dirs[2,tested_peak_indices[i]]+'/peakbagging_evidenceInformation.txt', usecols=(0,))
                        
                                if (ln_evidA >= ln_evidC):
                                    ln_evidAC = ln_evidA + np.log(1+np.exp(ln_evidC-ln_evidA)) 
                                else:
                                    ln_evidAC = ln_evidC + np.log(1+np.exp(ln_evidA-ln_evidC))
                        
                                if (ln_evidB >= ln_evidC):
                                    ln_evidBC = ln_evidB + np.log(1+np.exp(ln_evidC-ln_evidB)) 
                                else:
                                    ln_evidBC = ln_evidC + np.log(1+np.exp(ln_evidB-ln_evidC))
                    
                                p_CA[k,i] = np.exp(ln_evidC - ln_evidAC)
                                p_CB[k,i] = np.exp(ln_evidC - ln_evidBC)
                            
                    for i in range(0, len(dipole_indices[tested_peak_indices])):
                        
                        with open(probability_filenames[tested_peak_indices[i]],'w') as f:
                            f.write('# Detection probabilities for frequency peak {:.3f} muHz\n'.format(freq1[dipole_indices[tested_peak_indices[i]]]))
                            f.write('# Each line corresponds to a different run of the same test.\n')
                            if best_dnu < self.cp.dnu_rg:
                                f.write('# Col 1: p(BA), Col 2: p(CA), Col 3: p(CB), Col 4: FWHM single (microHz).\n')
                                for k in range(0,self.cp.n_peak_test):
                                    f.write('{:.4f}  {:10.4f}  {:10.4f}  {:10.4f}\n'.format(p_BA[k,i],p_CA[k,i],p_CB[k,i],fwhm_detection_fit[k,i]))
                            else:
                                f.write('# Col 1: p(BA), Col 2: FWHM single (microHz).\n')
                                for k in range(0,self.cp.n_peak_test):
                                    f.write('{:.4f}  {:10.4f}\n'.format(p_BA[k,i],fwhm_detection_fit[k,i]))
                                    
                        if self.cp.save_test_files==0:
                            # Remove test folders and prior files if required
                            shutil.rmtree(test_dirs[0,tested_peak_indices[i]])
                            shutil.rmtree(test_dirs[1,tested_peak_indices[i]])
                            os.remove(prior_filenames[0,tested_peak_indices[i]])
                            os.remove(prior_filenames[1,tested_peak_indices[i]])
                            os.remove(data_range_filenames[0,tested_peak_indices[i]])
                            os.remove(data_range_filenames[1,tested_peak_indices[i]])
                            if best_dnu < self.cp.dnu_rg:
                                shutil.rmtree(test_dirs[2,tested_peak_indices[i]])
                                os.remove(prior_filenames[2,tested_peak_indices[i]])
                                os.remove(data_range_filenames[2,tested_peak_indices[i]])     

                for i in range(0, len(dipole_indices[tested_peak_indices])):
                    if best_dnu < self.cp.dnu_rg:
                        p_BA,p_CA,p_CB = np.loadtxt(probability_filenames[tested_peak_indices[i]],usecols=(0,1,2),unpack=True)
                    else:
                        p_BA = np.loadtxt(probability_filenames[tested_peak_indices[i]],usecols=(0,))
            
                    # Consider the maximum probabilities among the different runs. Approximate up to third decimal digit.
                    max_p_BA = round(max(p_BA),3)
                    detection_probability[dipole_indices[tested_peak_indices[i]]] = max_p_BA
            
                    if best_dnu < self.cp.dnu_rg:
                        max_p_CA = round(max(p_CA),3)
                        max_p_CB = round(max(p_CB),3)
                
                        # For a RG star take the maximum probability between the two (Lorentzian and Sinc^2) to deem the peak as 
                        # detected. In this way one makes sure to incorporate the result from the best model possible,
                        # between the two tested.
                        if max_p_CB > 0.5:
                            detection_probability[dipole_indices[tested_peak_indices[i]]] = max_p_CA
                            sinc_profile_flag[dipole_indices[tested_peak_indices[i]]] = 1
                        else:
                            detection_probability[dipole_indices[tested_peak_indices[i]]] = max_p_BA

                    if self.cp.print_on_screen:
                        print('\n Peak detection test for frequency peak #{} at {:.2f} muHz'.format(tested_peak_indices[i],freq1[dipole_indices[tested_peak_indices[i]]]))
                        print(' P (Lorentzian vs Only Background): {}'.format(max_p_BA))
                        if best_dnu < self.cp.dnu_rg:
                            print(' P (Sinc^2 vs Only Background): '.format(max_p_CA))
                            print(' P (Sinc^2 vs Lorentzian): '.format(max_p_CB))
                        print('\n')
            
                if self.cp.print_on_screen:
                    print(' Peak Detection Test for candidate l=1 mode(s) completed.') 
                    print(' --------------------------------------------------------\n')
                    
            # -------------------------------------------------------------
            # Peak rotation and duplet test for all significant candidate l=1 modes
            # -------------------------------------------------------------
            
            # Update the number of dipole modes in the chunk by considering only the significant peaks. 
            # For the significant peaks (if any) perform a rotation test to check whether the peaks are individuals or multiplets. If the rotation
            # test is positive, then also save the rotational splitting associated with the fit, so that the individual
            # frequency components can be reconstructed. For RG stars, add the duplet test, to check whether the peak is actually constituted by
            # two adjacent mixed modes of different g-mode order n_g (i.e. not a rotational splitting effect).
            # Perform the test only if explicitly requested through the input configuring parameter file.
           
            detected_dipole_indices = np.where((angular_degree==1) & ((detection_probability >= self.cp.detection_probability_threshold) | (detection_probability==-99.0)) & (sinc_profile_flag != 1))[0]
            if self.cp.rotation_test_activated:    
                if len(detected_dipole_indices) > 0:
                    n_dipole_chunk = len(detected_dipole_indices)
                    test_dirs = np.zeros((3,len(detected_dipole_indices)),dtype='U200')
                    prior_filenames = np.zeros((3,len(detected_dipole_indices)),dtype='U200')
                    data_range_filenames = np.zeros((3,len(detected_dipole_indices)),dtype='U200')
                    probability_filenames = np.zeros(len(detected_dipole_indices),dtype='U200')
                    run_test = []
                    flag_run_test = 0

                    for j in range(0, len(detected_dipole_indices)):
                        dipole_index = detected_dipole_indices[j]
                        # For each peak consider the maximum data range between the range and the division
                        left_bound = min([range_maximum[0,dipole_index],divisions_maximum[0,dipole_index]])
                        right_bound = max([range_maximum[1,dipole_index],divisions_maximum[1,dipole_index]])
                        peak_number = str(int(dipole_index))
                        test_name = run + '_' + peak_number
                        test_nameE = test_name + 'E'           # Lorentzian profile, fixed background
                        test_nameF = test_name + 'F'           # Lorentzian profile split by rotation (considered as l=1), fixed background
                        test_nameG = test_name + 'G'           # Two Lorentzian profiles, fixed background
                        test_dirE = self.star_dir/self.cp.pb_subdir/test_nameE
                        test_dirF = self.star_dir/self.cp.pb_subdir/test_nameF
                        test_dirG = self.star_dir/self.cp.pb_subdir/test_nameG
                        rotation_probability_filename = self.star_dir/self.cp.pb_subdir/('rotationProbability_' + test_name + '.txt')
                        filenameE = self.star_dir/self.cp.pb_subdir/(self.cp.prior_filename + '_' + test_nameE + '.txt')
                        filenameF = self.star_dir/self.cp.pb_subdir/(self.cp.prior_filename + '_' + test_nameF + '.txt')
                        filenameG = self.star_dir/self.cp.pb_subdir/(self.cp.prior_filename + '_' + test_nameG + '.txt')
                        filename_rangeE = self.star_dir/self.cp.pb_subdir/('frequencyRange_' + test_nameE + '.txt')
                        filename_rangeF = self.star_dir/self.cp.pb_subdir/('frequencyRange_' + test_nameF + '.txt')
                        filename_rangeG = self.star_dir/self.cp.pb_subdir/('frequencyRange_' + test_nameG + '.txt')

                        test_dirs[:,j] = [test_dirE,test_dirF,test_dirG]
                        probability_filenames[j] = rotation_probability_filename
                        prior_filenames[:,j] = [filenameE,filenameF,filenameG]
                        data_range_filenames[:,j] = [filename_rangeE,filename_rangeF,filename_rangeG]
               
                        if not os.path.isfile(rotation_probability_filename):
                            # For RG and SG stars, adopt a narrower FWHM prior for those l=1 modes outside the l=3 region.
                            # This will speed up the computation of the peak test and improve the actual fits of the 
                            # rotational multiplets.
                            # Also, if the star is a MS, make the prior on the frequency centroid as narrow as possible.                    
                            right_fwhm = fwhm_radial_fit*self.cp.fwhm_magnification_factor_radial
                            left_fwhm = self.cp.fwhm_lower_bound

                            freq_prior = [freq1[dipole_index] - freq_sig1[dipole_index],freq1[dipole_index] + freq_sig1[dipole_index]]

                            if (best_dnu < self.cp.dnu_sg) & (teff < self.cp.teff_sg):
                                freq_prior = [range_maximum[0,dipole_index],range_maximum[1,dipole_index]]
                        
                                if (freq1[dipole_index] >= octupole_freq_upper) or (freq1[dipole_index] <= octupole_freq_lower):
                                    if best_dnu < self.cp.dnu_threshold:
                                        right_fwhm = fwhm_radial_fit
                                        # Do not exceed the width given by the frequency range of the peak if the star is a RG or late SG, because l=1 modes are narrow. 
                                        range_peak = abs(range_maximum[1,freq_index] - range_maximum[0,freq_index])
                                        if right_fwhm > range_peak:
                                            right_fwhm = range_peak
                                    else:
                                        right_fwhm = fwhm_radial_fit*self.cp.fwhm_magnification_factor_dipole
                            
                            tmp_peak = np.where((freq >= left_bound) & (freq <= right_bound))[0]
                            bg_peak = np.mean(bg_level_local[tmp_peak])
                            response_peak = np.mean((np.sin(np.pi/2. * freq[tmp_peak]/nyq) / (np.pi/2. * freq[tmp_peak]/nyq))**2)

                            # For the amplitude estimate of the rotational multiplet incorporate the np.sqrt(2) factor to take into account a case with i=90 degrees
                            amplitude = np.sqrt((abs(max(spsd[tmp_peak]) - bg_peak)/response_peak)*np.pi*right_fwhm)    
                            # Set up priors for the peak test profile
                            data_freq_boundaries = [left_bound,right_bound]
                            freq_duplet_prior = [range_maximum[0,dipole_index],range_maximum[1,dipole_index]]
                            duplet_split_prior = [freqbin*2,abs(freq_duplet_prior[1]-freq_duplet_prior[0])]
                            amplitude_prior = [0.0,amplitude]
                            fwhm_prior = [left_fwhm,right_fwhm]

                            # Divide the frequency prior range into nearly three parts
                            rot_split = abs(freq_prior[1] - freq_prior[0])/self.cp.rot_split_scaling
                            rot_split_prior = [freqbin*2,rot_split]
 
                            # If the frequency resolution is comparable to the expected separation among the fine-structure peaks, then skip the test. 
                            if duplet_split_prior[0] >= duplet_split_prior[1]*0.5:
                                continue  
                            if rot_split_prior[0] >= rot_split_prior[1]*0.5:
                                continue
                    
                            cosi_prior = [self.cp.cosi_prior_lower,self.cp.cosi_prior_upper]
                            boundariesE = [freq_prior,amplitude_prior,fwhm_prior]
                            boundariesF = [freq_prior,amplitude_prior,fwhm_prior,rot_split_prior,cosi_prior] 
                            diamonds.write_data_range(data_range_filenames[0,j],data_freq_boundaries)
                            diamonds.write_data_range(data_range_filenames[1,j],data_freq_boundaries)
                            diamonds.write_uniform_prior(prior_filenames[0,j],boundariesE)
                            diamonds.write_uniform_prior(prior_filenames[1,j],boundariesF)
                   
                            if best_dnu < self.cp.dnu_rg:
                                boundariesG = [freq_duplet_prior,amplitude_prior,fwhm_prior,duplet_split_prior,amplitude_prior,fwhm_prior]
                        
                                diamonds.write_data_range(data_range_filenames[2,j],data_freq_boundaries)
                                diamonds.write_uniform_prior(prior_filenames[2,j],boundariesG)
                        
                                run_test.extend([test_nameE,test_nameF,test_nameG])
                            else:
                                run_test.extend([test_nameE,test_nameF])
                            flag_run_test += 1

                    if self.cp.print_on_screen:
                        print('\n ---------------------------------------------')
                        print(' Peak Rotation Test for l=1 mode(s).' )
            
                    if flag_run_test > 0:
                        if self.cp.print_on_screen:
                            print(' A total of {} tests must be performed.'.format(len(run_test))) 
                                
                        p_FE = np.zeros((self.cp.n_peak_test,len(detected_dipole_indices)))
                        fwhm_rotation_fit = np.zeros((self.cp.n_peak_test,len(detected_dipole_indices)))
                        rot_split_fit = np.zeros((self.cp.n_peak_test,len(detected_dipole_indices)))
                        cosi_fit = np.zeros((self.cp.n_peak_test,len(detected_dipole_indices)))
                        central_freq_fit = np.zeros((self.cp.n_peak_test,len(detected_dipole_indices)))

                        if best_dnu < self.cp.dnu_rg:
                            p_GE = np.zeros((self.cp.n_peak_test,len(detected_dipole_indices)))
                            p_GF = np.zeros((self.cp.n_peak_test,len(detected_dipole_indices)))
                            left_freq_fit = np.zeros((self.cp.n_peak_test,len(detected_dipole_indices)))
                            right_freq_fit = np.zeros((self.cp.n_peak_test,len(detected_dipole_indices)))
                            left_fwhm_fit = np.zeros((self.cp.n_peak_test,len(detected_dipole_indices)))
                            right_fwhm_fit = np.zeros((self.cp.n_peak_test,len(detected_dipole_indices)))
            
                        for k in range(0, self.cp.n_peak_test):
                            peakbagging_parameters = {'subdir':       self.cp.pb_subdir,
                                                      'run':          run_test,
                                                      'background':   background_name,
                                                      'fwhm':         -1.0,
                                                      'filename_run': test_name}

                            flag_computation_completed = diamonds.run_peakbagging(self.catalog_id,self.star_id,peakbagging_parameters,1,0,0,self.cp.dp_pb, self.cp.diamonds_path, self.cp.n_threads, self.cp.prior_filename)
                    
                            if self.cp.print_on_screen:
                                print,' Iteration #{} completed.'.format(k) 
                                        
                            for j in range(0, len(detected_dipole_indices)):
                                # Make sure that the given peak has been tested.
                                if os.path.isfile(test_dirs[0,j]+'/peakbagging_evidenceInformation.txt'):
                                    ln_evidE = np.loadtxt(test_dirs[0,j]+'/peakbagging_evidenceInformation.txt', usecols=(0,))
                                    ln_evidF = np.loadtxt(test_dirs[1,j]+'/peakbagging_evidenceInformation.txt', usecols=(0,))
                                    if (ln_evidE >= ln_evidF):
                                        ln_evidEF = ln_evidE + np.log(1.+np.exp(ln_evidF-ln_evidE)) 
                                    else:
                                        ln_evidEF = ln_evidF + np.log(1.+np.exp(ln_evidE-ln_evidF))
                            
                                    p_FE[k,j] = np.exp(ln_evidF - ln_evidEF)

                                    # Compute the estimate for FWHM from model E.
                                    par_fwhm = np.loadtxt(test_dirs[0,j]+'/peakbagging_parameter002.txt')
                                    post = np.loadtxt(test_dirs[0,j]+'/peakbagging_posteriorDistribution.txt')
                                    post /= max(post)
                                    fwhm_rotation_fit[k,j] = np.sum(par_fwhm*post)/np.sum(post)
                            
                                    # Compute the estimates for nu0, rot_split and cosi from model F.
                            
                                    par_central_freq = np.loadtxt(test_dirs[1,j]+'/peakbagging_parameter000.txt')
                                    par_rot_split = np.loadtxt(test_dirs[1,j]+'/peakbagging_parameter003.txt')
                                    par_cosi = np.loadtxt(test_dirs[1,j]+'/peakbagging_parameter004.txt')
                                    post = np.loadtxt(test_dirs[1,j]+'/peakbagging_posteriorDistribution.txt')

                                    # In this case apply an additional trick. Consider only the last 100 nested iterations to get a first
                                    # estimate of the frequency centroid, and compute a standard deviation of the entire sampling. 
                                    # Then recompute the parameter mean by using the sampling located only within 1-sigma of the first estimate.
                                    # This will improve the estimation of the frequency centroid and of the rotational splitting in case
                                    # multiple solutions have been found during the nested sampling process.
                                    post /= max(post)
                                    rot_split_fit[k,j] = np.sum(par_rot_split*post)/np.sum(post)
                                    cosi_fit[k,j] = np.sum(par_cosi*post)/np.sum(post)
                                    central_freq_fit[k,j] = np.sum(par_central_freq*post)/np.sum(post)
                            
                                    if best_dnu < self.cp.dnu_rg:
                                        ln_evidG = np.loadtxt(test_dirs[2,j]+'/peakbagging_evidenceInformation.txt', usecols=(0,))
                                        if (ln_evidG >= ln_evidF):
                                            ln_evidGF = ln_evidG + np.log(1.+np.exp(ln_evidF-ln_evidG)) 
                                        else:
                                            ln_evidGF = ln_evidF + np.log(1.+np.exp(ln_evidG-ln_evidF))

                                        if (ln_evidG >= ln_evidE):
                                            ln_evidGE = ln_evidG + np.log(1.+np.exp(ln_evidE-ln_evidG)) 
                                        else:
                                            ln_evidGE = ln_evidE + np.log(1.+np.exp(ln_evidG-ln_evidE))
                                
                                        p_GE[k,j] = np.exp(ln_evidG - ln_evidGE)
                                        p_GF[k,j] = np.exp(ln_evidG - ln_evidGF)
                            
                                        # Compute the estimates for the duplet frequencies and FWHM from model G
                                        par_left_freq = np.loadtxt(test_dirs[2,j]+'/peakbagging_parameter000.txt')
                                        par_duplet_split = np.loadtxt(test_dirs[2,j]+'/peakbagging_parameter003.txt')
                                        par_left_fwhm = np.loadtxt(test_dirs[2,j]+'/peakbagging_parameter002.txt')
                                        par_right_fwhm = np.loadtxt(test_dirs[2,j]+'/peakbagging_parameter005.txt')
                                        post = np.loadtxt(test_dirs[2,j]+'/peakbagging_posteriorDistribution.txt')
                                        post /= max(post)
                                        left_freq_fit[k,j] = np.sum(par_left_freq*post)/np.sum(post)
                                        right_freq_fit[k,j] = np.sum(par_duplet_split*post)/np.sum(post) + left_freq_fit[k,j]
                                        left_fwhm_fit[k,j] = np.sum(par_left_fwhm*post)/np.sum(post)
                                        right_fwhm_fit[k,j] = np.sum(par_right_fwhm*post)/np.sum(post)
                            
                        for j in range(0, len(detected_dipole_indices)):
                            with open(probability_filenames[j],'w') as f:
                                f.write('# Rotation probabilities for frequency peak {:.2f} muHz\n'.format(freq1[detected_dipole_indices[j]]))
                                f.write('# Each line corresponds to a different run of the same test.\n')
                                if best_dnu < self.cp.dnu_rg:
                                    f.write('# Col 1: p(FE), Col 2: p(GE), Col 3: p(GF), Col 4: FWHM single (microHz), Col 5: Left freq duplet (microHz) \n')
                                    f.write('# Col 6: Right freq duplet (microHz), Col 7: Left FWHM (microHz), Col 8: Right FWHM (microHz) \n')
                                    f.write('# Col 9: Central freq (microHz), Col 10: rotational splitting (microHz), Col 11: cos i \n')
                                    for k in range(0,self.cp.n_peak_test):
                                        f.write('{:.4f}  {:10.4f}  {:10.4f}  {:10.4f}  {:10.4f}  {:10.4f}  {:10.4f}  {:10.4f}  {:10.4f}  {:10.4f}  {:10.4f}\n'.format(p_FE[k,j], p_GE[k,j], p_GF[k,j], fwhm_rotation_fit[k,j], left_freq_fit[k,j], right_freq_fit[k,j], left_fwhm_fit[k,j], right_fwhm_fit[k,j], central_freq_fit[k,j], rot_split_fit[k,j], cosi_fit[k,j]))
                                else:
                                    f.write('# Col 1: p(FE), Col 2: FWHM single (microHz), Col 3: Central freq (microHz), Col 4: rotational splitting (microHz), Col 5: cos i \n')
                                    for k in range(0,self.cp.n_peak_test):
                                        f.write('{:.4f}  {:10.4f}  {:10.4f}  {:10.4f}  {:10.4f}\n'.format(p_FE[k,j], fwhm_rotation_fit[k,j], central_freq_fit[k,j], rot_split_fit[k,j], cosi_fit[k,j]))
                                    
                    
                            if self.cp.save_test_files==0:
                                # Remove test folders and prior files if required
                                if os.path.isdir(test_dirs[0,j]):
                                    shutil.rmtree(test_dirs[0,j])
                                    shutil.rmtree(test_dirs[1,j])
                                    os.remove(data_range_filenames[0,j])
                                    os.remove(data_range_filenames[1,j])
                                    os.remove(prior_filenames[0,j])
                                    os.remove(prior_filenames[1,j])
                            
                                    if best_dnu < self.cp.dnu_rg:
                                        shutil.rmtree(test_dirs[2,j])
                                        os.remove(data_range_filenames[2,j])
                                        os.remove(prior_filenames[2,j])
                            
                    for j in range(0, len(detected_dipole_indices)):
                        if os.path.isfile(probability_filenames[j]):
                            if best_dnu < self.cp.dnu_rg:
                                p_FE,p_GE,p_GF = np.loadtxt(probability_filenames[j], usecols=(0,1,2), unpack=True)
                                max_p_GE = round(max(p_GE),3)
                                max_p_GF = round(max(p_GF),3)
                                duplicity_probability[detected_dipole_indices[j]] = max_p_GE
                            else:
                                p_FE = np.loadtxt(probability_filenames[j], usecols=(0,))

                            # Consider the maximum probabilities among the different runs. Approximate up to second decimal digit.
                            
                            max_p_FE = round(max(p_FE),3)
                            rotation_probability[detected_dipole_indices[j]] = max_p_FE
                    
                            if self.cp.print_on_screen:
                                print('\n Peak rotation and duplet test for frequency peak #{} at {:.2f} muHz'.format(detected_dipole_indices[j],freq1[detected_dipole_indices[j]]))
                                print(' P (Rotation vs No rotation): {}'.format(max_p_FE))
                    
                                if best_dnu < self.cp.dnu_rg:
                                    print(' P (Duplet vs No Rotation): {}'.format(max_p_GE))
                                    print(' P (Duplet vs Rotation): {}'.format(max_p_GF))
                                print()
            
                    if self.cp.print_on_screen:
                        print('\n Peak Rotation Test for l=1 mode(s) completed.')
                        print(' ---------------------------------------------\n')

            
                else:
                    n_dipole_chunk = 0
            else:
                if self.cp.print_on_screen:
                    print('\n --------------------------------------------------------')
                    print(' Peak Rotation Test Disabled.')
                    print(' --------------------------------------------------------\n')

         
            # If at least one candidate l=1 mode is detected, make sure that the n_dipole_chunk is activated in case it is not (i.e. if there is no dipole
            # mode found from the global modality)
    
            if (len(detected_dipole_indices) > 0) & (n_dipole_chunk==0): 
                n_dipole_chunk = 1
                freq_dipole = freq_radial_chunk - best_dnu/2. - d01
        else:
            # If here, then no dipole peaks were found in the chunk from the multi-modal fit. 
            n_dipole_chunk = 0

        largest_octupole_fwhm = 0.0
        
        if n_dipole_chunk != 0:
            # -------------------------------------------------------------------
            # CASE 1: One or more significant dipole modes are found in the chunk 
            # -------------------------------------------------------------------
    
            # Evolutionary stage: RG, SG, MS
            # In general, search for an l=3 only if more than one l=1 mode is detected. This is to give priority to l=1 mode identification over l=3.
            # If only one l=1 mode is detected, perform the l=3 search only in the case of a RG star or a depressed dipole star (which can be a SG too). 
    
            if (len(detected_dipole_indices) > 1) or ((len(detected_dipole_indices)==1) & (best_dnu <= self.cp.dnu_rg)) or (flag_depressed_dipole==1):
                # Check whether an l=3 is present. To do so attempt to find at least one local maximum in the l=3 region as 
                # computed from the asymptotic relation. The search region depends on the evolutionary stage of the star.
                # Apply the additional condition that if the mode is flagged as a sinc^2 profile, then it is excluded from the candidate octupole mode list.

                detected_octupole_index = np.where((freq1 <= octupole_freq_upper) & (freq1 >= octupole_freq_lower) & ((detection_probability >= self.cp.detection_probability_threshold) | (detection_probability==-99.0)) & (sinc_profile_flag != 1))[0]
        
                if len(detected_octupole_index) > 0:
                    # In this case there is at least one significant peak in this l=3 region.
                    # Therefore check whether one of the l=3 candidates is a true l=3 mode, or it is favoring the scenario of a l=1 mixed mode.
                    # If the sinc^2 is not favored, then verify that the detection favors a single Lorentzian peak (i.e. both the rotational splitting and the duplicity are not present). 
                    # Finally use the FWHM from the single Lorentzian peak test of the candidate l=3 peak to assess whether it is an actual l=3 or a mixed mode peak.
                    # Do so for each peak found in the l=3 range from the asymptotic relation.
                    
                    # Obtain the largest FWHM of the set of modes. This will be used as discriminant in case more than one candidate l=3 is present.
                    fwhm_octupole_fit = np.zeros(len(detected_octupole_index))
                    run_subdir = run + '_octupole_fwhm'

                    if self.cp.print_on_screen:
                        print('\n Fitting Linewidth of candidate octupole modes of the chunk.')
                        
                    for j in range(0, len(detected_octupole_index)):
                        octupole_index = detected_octupole_index[j]
                        peak_number = str(int(octupole_index))
                        run_names = [run_subdir+str(x)+'_'+peak_number for x in range(self.cp.n_fwhm_fit)]
           
                        if not os.path.isfile(self.star_dir/self.cp.pb_subdir/run_names[j]/'peakbagging_computationParameters.txt') or force:
                            right_bound = max([range_maximum[1,octupole_index],divisions_maximum[1,octupole_index]])
                            left_bound = min([range_maximum[0,octupole_index],divisions_maximum[0,octupole_index]])
                            tmp0 = np.where((freq >= left_bound) & (freq <= right_bound))[0]
                            freq0 = freq[tmp0]
                            spsd0 = spsd[tmp0]
                            bg_peak = np.mean(bg_level_local[tmp0])
                            response_peak = np.mean((np.sin(np.pi/2. * freq0/nyq) / (np.pi/2. * freq0/nyq))**2)
                            amplitude_octupole = np.sqrt((abs(spsd_maximum[octupole_index] - bg_peak)/response_peak)*np.pi*fwhm_radial_fit*self.cp.fwhm_magnification_factor_octupole)
                            freq_prior_octupole = [range_maximum[0,octupole_index],range_maximum[1,octupole_index]]
                            amplitude_prior_octupole = [0.0,amplitude_octupole]
                            data_freq_boundaries = [left_bound,right_bound]
                            data_range_filenames = np.zeros(self.cp.n_fwhm_fit,dtype='U200')
                            prior_filenames = np.zeros(self.cp.n_fwhm_fit,dtype='U200')

                            # Allow for a very narrow FWHM
                            left_fwhm = self.cp.fwhm_lower_bound
                            fwhm_prior = [left_fwhm,fwhm_radial_fit*self.cp.fwhm_magnification_factor_octupole]
                            boundaries = [freq_prior_octupole,amplitude_prior_octupole,fwhm_prior]
                            
                            for k in range(0, self.cp.n_fwhm_fit):
                                data_range_filenames[k] = self.star_dir/self.cp.pb_subdir/('frequencyRange_' + run_names[k] + '.txt')
                                prior_filenames[k] = self.star_dir/self.cp.pb_subdir/(self.cp.prior_filename + '_' + run_names[k] + '.txt')
                                diamonds.write_data_range(data_range_filenames[k],data_freq_boundaries)
                                diamonds.write_uniform_prior(prior_filenames[k],boundaries)
                    

                            peakbagging_parameters = {'subdir':          self.cp.pb_subdir,
                                                      'run':             run_names,
                                                      'background':      background_name,
                                                      'fwhm':            -1.0,
                                                      'filename_run':    run_subdir}

                            flag_computation_completed_array = diamonds.run_peakbagging(self.catalog_id,self.star_id,peakbagging_parameters,2,0,0, self.cp.dp_pb, self.cp.diamonds_path, self.cp.n_threads, self.cp.prior_filename)

                            if not self.cp.save_test_files:
                                for k in range(0, self.cp.n_fwhm_fit):
                                    os.remove(data_range_filenames[k])
                                    os.remove(prior_filenames[k])
                        
                    
                        fwhm_octupole_fit_array = np.zeros(self.cp.n_fwhm_fit)

                        for k in range(0, self.cp.n_fwhm_fit):
                            fwhm_parameter = '002'
                            # Read the sampled FWHM of the octupole mode.
                            par_fwhm0 = np.loadtxt(self.star_dir/self.cp.pb_subdir/run_names[k]/('peakbagging_parameter' + fwhm_parameter + '.txt'))
                            # Load the posterior samples to compute Bayesian mean estimate for FWHM_0
                            post = np.loadtxt(self.star_dir/self.cp.pb_subdir/run_names[k]/'peakbagging_posteriorDistribution.txt')
                            # Compute parameter estimate, using a weighted average
                            post /= max(post)
                            fwhm_octupole_fit_array[k] = np.sum(par_fwhm0*post)/np.sum(post)
                

                        # Save the largest FWHM of the set as well as the corresponding absolute frequency index
                        fwhm_octupole_fit[j] = np.median(fwhm_octupole_fit_array)
                        if fwhm_octupole_fit[j] >= largest_octupole_fwhm:
                            largest_octupole_fwhm = fwhm_octupole_fit[j]
                            best_octupole_index = octupole_index

                    if self.cp.print_on_screen:
                        print(' Maximum FWHM (candidate l=3) = {:.3f} muHz \n'.format(largest_octupole_fwhm))

                    # Now proceed by performing the octupole mode test on the peak that has the largest FWHM if it has been detected.
                    peak_number = str(int(best_octupole_index))
                    test_name = run + '_' + peak_number
        
                    # If not, now verify whether the peak is split by rotation. If the rotation test is activated, given that all peaks considered here are
                    # significant, by default the rotation test is performed. In the case of a RG, check also for the duplicity.
                    flag_octupole_fwhm_test = 0

                    if self.cp.rotation_test_activated: 
                        rotation_probability_filename = self.star_dir/self.cp.pb_subdir/('rotationProbability_' + test_name + '.txt')
                
                        if os.path.isfile(rotation_probability_filename)==1:
                            if best_dnu < self.cp.dnu_rg:
                                # This is the case of RG stars
                                 p_FE,p_GE,p_GF,freq_left,freq_right,fwhm_left,fwhm_right = np.loadtxt(rotation_probability_filename, usecols=(0,1,2,4,5,6,7),unpack=True)
                    
                                 max_p_FE = round(max(p_FE),3)
                                 max_p_GE = round(max(p_GE),3)
                                 max_p_GF = round(max(p_GF),3)
                    
                                 if (max_p_FE >= self.cp.rotation_probability_threshold) or (max_p_GE >= self.cp.duplicity_probability_threshold):
                                     # Here either rotation or duplicity (or both) were detected.
                                     # Then check whether the duplicity is detected over rotation.
                            
                                     if (max_p_GF > 0.5):
                                         # Here the peak is considered as a duplet. In this case verify whether any of the two peaks of the duplet has a
                                         # FWHM comparable to that of the adjacent l=0 mode. If this is true, then consider such peak as a
                                         # l=3, and split it up from the other peak, which is therefore flagged as a l=1 mixed mode only at the end of the CHUNK modality.
                                
                                         fwhm_duplet = [np.mean(fwhm_left),np.mean(fwhm_right)]

                                         if (fwhm_duplet[0] >= fwhm_radial_fit*self.cp.fwhm_octupole_radial_fraction) or (fwhm_duplet[1] >= fwhm_radial_fit*self.cp.fwhm_octupole_radial_fraction):
                                             # At least one of the peaks in the duplet has a FWHM comparable to that of the radial mode.
                                             # Then consider the mode a potential octupole.
                                             # Check whether the frequency resolution is enough to be able to distinguish between a l=1 mixed mode and a l=3.
                                             # If the peak will pass the test, it will be split up later.
                                             flag_octupole_fwhm_test = 1
                                 else:
                                     # No rotation and no duplicity were detected. Then the peak is considered as a single Lorentzian profile.
                                     # Finally check its FWHM against that of the adjacent l=0 mode.
                                     # If the fitted linewidth is comparable to that of the adjacent radial mode, then consider the mode a potential octupole.
                                     # Check whether the frequency resolution is enough to be able to distinguish between a l=1 mixed mode and a l=3.
                                     flag_octupole_fwhm_test = 1
                        
                            else:
                                # This is the case of SG and MS stars
                                p_FE = np.loadtxt(rotation_probability_filename, usecols=(0,))
                                max_p_FE = round(max(p_FE),3)
                        
                                if max_p_FE < self.cp.rotation_probability_threshold:
                                    # Rotation is not detected. Then the peak is considered as a single Lorentzian profile.
                                    # Finally check its FWHM against that of the adjacent l=0 mode.
                                    # If the fitted linewidth is comparable to that of the adjacent radial mode, then consider the mode as a potential octupole.
                                    # Check whether the frequency resolution is enough to be able to distinguish between a l=1 mixed mode and a l=3.

                                    flag_octupole_fwhm_test = 1
                        
                        else:
                            # Here the rotation test was not performed because the frequency resolution is not sufficiently high. 
                            # Then look for an l=3 based solely on the fitted linewidth from the detection tests.
                            flag_octupole_fwhm_test = 1
                
                    else:
                        # Here the rotation test was not performed because deactivated by the user. 
                        # Then look for an l=3 based solely on the fitted linewidth from the detection tests.
                        flag_octupole_fwhm_test = 1
            
                    # If required, assess the FWHM and ASEF of the l=3 candidate against those of the l=0 mode.
                    if flag_octupole_fwhm_test: 
                        # Impose that l=3, if present, has a low ASEF value (< 3/4 ASEF maximum)
                        asef_threshold = (self.cp.dp_isla['max_nested_it'] + self.cp.dp_isla['n_live']) * self.cp.asef_threshold_fraction
                
                        if asef_maximum[best_octupole_index] < asef_threshold:
                            if largest_octupole_fwhm >= fwhm_radial_fit*self.cp.fwhm_octupole_radial_fraction: 
                                angular_degree[best_octupole_index] = 3
                                order_number[best_octupole_index] = enn_radial - 2
                    
            # If the star is a SG, check that the detected l=1 modes are not closer to one another than Dnu/X.
            # First update the detected dipole indices, in case some manipulation has occurred from the inspection of l=3 modes.
            detected_dipole_indices = np.where((angular_degree==1) & ((detection_probability >= self.cp.detection_probability_threshold) | (detection_probability==-99.0)))[0]

            if best_dnu > self.cp.dnu_rg and (len(detected_dipole_indices) > 1):
                for k in range(0, len(detected_dipole_indices)-1):
                    # Check the separation between two consecutive l=1 candidate mixed modes
                    freq_left_dipole = freq1[detected_dipole_indices[k]]
                    freq_right_dipole = freq1[detected_dipole_indices[k+1]]

                    if abs(freq_left_dipole - freq_right_dipole) < best_dnu**2/self.cp.dnu_mixed_modes_separation_scaling/self.cp.dnu_rg:
                        # Since the two peaks are too close, these peaks are likely to originate from a single multiplet. 
                        # Hence take as detected only the one with the largest sampling counts
                        worst_peak_index = np.argmin([sampling_counts[detected_dipole_indices[k]],sampling_counts[detected_dipole_indices[k+1]]])
                        worst_peak_index = detected_dipole_indices[k + worst_peak_index]
                        detection_probability[worst_peak_index] = 0.0
            
            # If the star is a MS or a high-luminosity RGB star, make sure that there is only one dipole mode for this chunk.
            if (((best_dnu >= self.cp.dnu_rg) & (best_dnu < self.cp.dnu_sg) & (teff >= self.cp.teff_sg)) or (best_dnu >= self.cp.dnu_sg) or (best_dnu <= self.cp.dnu_tip)):
                detected_dipole_indices = np.where((angular_degree==1) & ((detection_probability >= self.cp.detection_probability_threshold) | (detection_probability==-99.0)))[0]
                # If more than one dipole mode is found, then pick up the one with the best combination of SPSD maximum, ASEF maximum, sampling counts and 
                # frequency position with respect to the global frequency of the dipole mode.
                if len(detected_dipole_indices) > 1:
                    freq_diff = abs(freq1[detected_dipole_indices] - freq_dipole)
                    freq_weights = 1/freq_diff
                    freq_weights /= np.sum(freq_weights)
                    freq_ww = freq_weights/max(freq_weights)
                    asef_ww = asef_weights[detected_dipole_indices]/max(asef_weights[detected_dipole_indices])
                    spsd_ww = spsd_weights[detected_dipole_indices]/max(spsd_weights[detected_dipole_indices])
                    sampling_ww = np.log(sampling_counts[detected_dipole_indices])/max(np.log(sampling_counts[detected_dipole_indices]))
                    total_ww = self.cp.weight_freq_fraction*freq_ww + self.cp.weight_asef_fraction*asef_ww + self.cp.weight_spsd_fraction*spsd_ww + self.cp.weight_sampling_fraction*sampling_ww
                    total_ww /= np.sum(total_ww)
                    max_ww = max(total_ww)
                    index = np.argmax(total_ww)
                    dipole_index = detected_dipole_indices[index]
                    freq_dipole_chunk = freq1[dipole_index]
                    bad_dipole_indices = np.where(detected_dipole_indices != dipole_index)[0]

                    # Flag the bad modes as undetected to remove them from the list of good frequencies.
                    detection_probability[detected_dipole_indices[bad_dipole_indices]] = 0.0

                    
        # Produce the final list of detected frequency peaks, and attached mode identification, to be stored as output and
        # overplotted on the PSD of the star.

        # First, select only modes that are significant with respect to noise.
        detected_indices = np.where((detection_probability >= self.cp.detection_probability_threshold) | (detection_probability==-99.0))[0]
        local_dp = 0.

        if len(detected_indices) >0 :
            freq1_final = freq1[detected_indices]
            freq_sig1_final = freq_sig1[detected_indices]
            asef_maximum_final = asef_maximum[detected_indices]
            spsd_maximum_final = spsd_maximum[detected_indices]
            sampling_counts_final = sampling_counts[detected_indices]
            angular_degree_final = angular_degree[detected_indices]
            order_number_final = order_number[detected_indices]
            range_maximum_final = range_maximum[:,detected_indices]
            divisions_maximum_final = divisions_maximum[:,detected_indices]
            detection_probability_final = detection_probability[detected_indices]
            rotation_probability_final = rotation_probability[detected_indices]
            duplicity_probability_final = duplicity_probability[detected_indices]
            blending_profile_flag_final = blending_profile_flag[detected_indices]
            sinc_profile_flag_final = sinc_profile_flag[detected_indices]

            # If possible, evaluate a local value for the observed period spacing of dipolar mixed modes
            tmp_mixed_dipole = np.where((angular_degree_final==1) & (freq1_final < freq_quadrupole_chunk))[0]
            n_mixed_dipole = len(tmp_mixed_dipole)
            
            if n_mixed_dipole > 1:
                # At least two dipolar modes are found. Therefore DP can be estimated.
                freq_mixed_dipole = freq1_final[tmp_mixed_dipole]
                dp_array = np.zeros(n_mixed_dipole-1)

                for i in range(0,n_mixed_dipole-1):
                    dp_array[i] = abs(1/freq_mixed_dipole[i] - 1/freq_mixed_dipole[i+1]) * 1e6   # units in seconds
        
                local_dp = np.mean(dp_array)
                
            # Add information about the rotational components, if these will become available.
            azimuthal_number_final = np.zeros(len(detected_indices),dtype=int) - 99
            cosi_final = np.zeros(len(detected_indices)) - 99.0

            # Second, check whether there are modes split by rotation or that are duplets. Apply separate probability thresholds for each case.
            rotation_duplicity_indices = np.where((rotation_probability_final >= self.cp.rotation_probability_threshold) | (duplicity_probability_final >= self.cp.duplicity_probability_threshold))[0]

            if len(rotation_duplicity_indices) > 0:
                # First store the values of the detected frequencies with rotation and/or duplicity in order to split them up and apply 
                # a mode identification to each component.
                freq_local_array = freq1_final[rotation_duplicity_indices]
                for j in range(0, len(rotation_duplicity_indices)):
                    freq_local = freq_local_array[j]
                    local_index = closest(freq_local,freq1_final,index=True)
                    freq_index = closest(freq_local,freq1,index=True)
                    peak_number = str(int(freq_index))
                    test_name = run + '_' + peak_number
                    rotation_probability_filename = self.star_dir/self.cp.pb_subdir/('rotationProbability_' + test_name + '.txt')
                    flag_duplicity_split = 0

                    if os.path.isfile(rotation_probability_filename):
                        if (duplicity_probability_final[local_index] != -99.0):
                            p_FE,p_GE,p_GF,left_freq,right_freq,left_fwhm,right_fwhm,central_freq,rot_split,cosi = np.loadtxt(rotation_probability_filename, usecols=(0,1,2,4,5,6,7,8,9,10),unpack=True)
                    
                            max_p_FE = round(max(p_FE),3)
                            max_p_GE = round(max(p_GE),3)
                            max_p_GF = round(max(p_GF),3)
                        
                            if (max_p_GE >= self.cp.duplicity_probability_threshold) & (max_p_GF > 0.5):
                                flag_duplicity_split = 1
                                freq_duplet = [np.mean(left_freq),np.mean(right_freq)]
                        else:
                            p_FE,central_freq,rot_split,cosi = np.loadtxt(rotation_probability_filename, usecols=(0,2,3,4),unpack=True)
                            max_p_FE = round(max(p_FE),3)
                               
                        external_indices = np.where(freq1_final != freq_local)[0]
                        if flag_duplicity_split==1:
                            # In this case the duplicity is the favored scenario.
                            range_midpoint = (freq_duplet[1] + freq_duplet[0]) / 2.0
                            range_maximum_new = [[range_maximum_final[0,local_index],range_midpoint],[range_midpoint,range_maximum_final[1,local_index]]]
                            divisions_maximum_new = [[divisions_maximum_final[0,local_index],range_midpoint],[range_midpoint,divisions_maximum_final[1,local_index]]]

                            
                            left_left_bound = range_maximum_new[0][0]
                            left_right_bound = range_maximum_new[1][0]
                            right_left_bound = range_maximum_new[0][1]
                            right_right_bound = range_maximum_new[1][1]
                   
                            # Left peak
                            tmp_range = np.where((par0 < left_right_bound) & (par0 >= left_left_bound))[0]
                            tmp_freq_range = np.where((freq < left_right_bound) & (freq >= left_left_bound))[0]
                            tmp_hist_range = np.where((par_hist < left_right_bound) & (par_hist >= left_left_bound))[0]
                            
                            if len(tmp_freq_range) == 0:
                                continue
                        
                            spsd_range = spsd[tmp_freq_range]
                            spsd_maximum_left = max(spsd_range)

                            if len(tmp_range) > 0:
                                par0_range = par0[tmp_range]
                                sampling_counts_left = np.sum(nest_iter[tmp_range])
                                freq_sig_left = np.sqrt(np.sum((par0_range-freq_duplet[0])**2*tmp_range**2)/np.sum(tmp_range**2))
                            else:                    
                                right_bound_new = right_right_bound
                                left_bound_new = left_left_bound
                                tmp_range = np.where((par0 < right_bound_new) & (par0 >= left_bound_new))[0]
                                par0_range = par0[tmp_range]
                                sampling_counts_left = np.sum(nest_iter[tmp_range])/2.
                                freq_sig_left = np.sqrt(np.sum((par0_range-freq_duplet[0])**2*tmp_range**2)/np.sum(tmp_range**2))/np.sqrt(2.)

                            if len(tmp_hist_range) > 0:
                                asef_maximum_left = max(asef_hist[tmp_hist_range])
                            else:
                                right_bound_new = right_right_bound
                                left_bound_new = left_left_bound
                                tmp_hist_range = np.where((par_hist < right_bound_new) & (par_hist >= left_bound_new))[0]
                                asef_maximum_left = max(asef_hist[tmp_hist_range])
                    
                            # Right peak
                            tmp_range = np.where((par0 < right_right_bound) & (par0 >= right_left_bound))[0]
                            tmp_freq_range = np.where((freq < right_right_bound) & (freq >= right_left_bound))[0]
                            tmp_hist_range = np.where((par_hist < right_right_bound) & (par_hist >= right_left_bound))[0]
                    
                            if len(tmp_freq_range)==0:
                                continue
                            spsd_range = spsd[tmp_freq_range]
                            spsd_maximum_right = max(spsd_range)
                    
                            if len(tmp_range) > 0:
                                par0_range = par0[tmp_range]
                                sampling_counts_right = np.sum(nest_iter[tmp_range])
                                freq_sig_right = np.sqrt(np.sum((par0_range-freq_duplet[1])**2*tmp_range**2)/np.sum(tmp_range**2))
                            else:
                                right_bound_new = right_right_bound
                                left_bound_new = left_left_bound
                                tmp_range = np.where((par0 < right_bound_new) & (par0 >= left_bound_new))[0]
                                par0_range = par0[tmp_range]
                                sampling_counts_right = np.sum(nest_iter[tmp_range])/2.
                                freq_sig_right = np.sqrt(np.sum((par0_range-freq_duplet[1])**2*tmp_range**2)/np.sum(tmp_range**2))/np.sqrt(2.)
                                
                            if len(tmp_hist_range) > 0:
                                asef_maximum_right = max(asef_hist[tmp_hist_range])
                            else:
                                right_bound_new = right_right_bound
                                left_bound_new = left_left_bound
                                tmp_hist_range = np.where((par_hist < right_bound_new) & (par_hist >= left_bound_new))[0]
                                asef_maximum_right = max(asef_hist[tmp_hist_range])
                                       
                            # If the peak corresponding to the duplet is flagged as a l = 3 mode, check which of the two peaks is the l=3, and 
                            # flag the other one as a l=1. The l=3 will be the one with the largest FWHM of the two.

                            if angular_degree_final[local_index]==3:
                                max_duplet_fwhm = max([np.mean(left_fwhm),np.mean(right_fwhm)])
                                octupole_index_local = [np.mean(left_fwhm),np.mean(right_fwhm)].index(max_duplet_fwhm)
                                angular_degree_duplet = [1,1]
                                order_number_duplet = [enn_radial-1,enn_radial-1]
                                angular_degree_duplet[octupole_index_local] = 3
                                order_number_duplet[octupole_index_local] = enn_radial-2
                            else:
                                angular_degree_duplet = [1,1]
                                order_number_duplet = [enn_radial-1,enn_radial-1]
                    
                            if len(external_indices) > 0:
                                range_maximum_final = np.hstack((range_maximum_final[:,external_indices],range_maximum_new))
                                divisions_maximum_final = np.hstack((divisions_maximum_final[:,external_indices],divisions_maximum_new))
                                freq1_final = np.append(freq1_final[external_indices],freq_duplet)
                                freq_sig1_final = np.append(freq_sig1_final[external_indices],[freq_sig_left,freq_sig_right])
                                asef_maximum_final = np.append(asef_maximum_final[external_indices],[asef_maximum_left,asef_maximum_right])
                                spsd_maximum_final = np.append(spsd_maximum_final[external_indices],[spsd_maximum_left,spsd_maximum_right])
                                sampling_counts_final = np.append(sampling_counts_final[external_indices],[sampling_counts_left,sampling_counts_right])

                                angular_degree_final = np.append(angular_degree_final[external_indices],[angular_degree_duplet])
                                order_number_final = np.append(order_number_final[external_indices],[order_number_duplet])
                                azimuthal_number_final = np.append(azimuthal_number_final[external_indices],[-99,-99])
                                cosi_final = np.append(cosi_final,[-99.0,-99.0])
                        
                                # Assign the duplet the same detection probability, but remove the duplicity and rotation test from it, in order to avoid
                                # that the peak is split up again later on.
                        
                                detection_probability_final = np.append(detection_probability_final[external_indices],np.zeros(2) + detection_probability_final[local_index])
                                rotation_probability_final = np.append(rotation_probability_final[external_indices],np.zeros(2)-99.0)
                                duplicity_probability_final = np.append(duplicity_probability_final[external_indices],np.zeros(2)-99.0)
                                blending_profile_flag_final = np.append(blending_profile_flag_final[external_indices],np.zeros(2,dtype=int))
                                sinc_profile_flag_final = np.append(sinc_profile_flag_final[external_indices],np.zeros(2,dtype=int))
                       
                                # Sort all arrays with increasing frequency order.
                                sorted_indices = np.argsort(freq1_final)
                                range_maximum_final = range_maximum_final[:,sorted_indices]
                                divisions_maximum_final = divisions_maximum_final[:,sorted_indices]
                                freq1_final = freq1_final[sorted_indices]
                                freq_sig1_final = freq_sig1_final[sorted_indices]
                                asef_maximum_final = asef_maximum_final[sorted_indices]
                                spsd_maximum_final = spsd_maximum_final[sorted_indices]
                                sampling_counts_final = sampling_counts_final[sorted_indices]
                                angular_degree_final = angular_degree_final[sorted_indices]
                                order_number_final = order_number_final[sorted_indices]
                                azimuthal_number_final = azimuthal_number_final[sorted_indices]
                                cosi_final = cosi_final[sorted_indices]
                                
                                detection_probability_final = detection_probability_final[sorted_indices]
                                rotation_probability_final = rotation_probability_final[sorted_indices]
                                duplicity_probability_final = duplicity_probability_final[sorted_indices]
                                blending_profile_flag_final = blending_profile_flag_final[sorted_indices]
                                sinc_profile_flag_final = sinc_profile_flag_final[sorted_indices]
                            else:
                                range_maximum_final = range_maximum_new
                                divisions_maximum_final = divisions_maximum_new
                                freq1_final = freq_duplet
                                freq_sig1_final = [freq_sig_left,freq_sig_right]
                                asef_maximum_final = [asef_maximum_left,asef_maximum_right]
                                spsd_maximum_final = [spsd_maximum_left,spsd_maximum_right]
                                sampling_counts_final = [sampling_counts_left,sampling_counts_right]
                                angular_degree_final = angular_degree_duplet
                                order_number_final = order_number_duplet
                                azimuthal_number_final = [-99.0,-99.0]
                                cosi_final = [-99.0,-99.0]
                        
                                # Assign the duplet the same detection probability, but remove the duplicity and rotation test from it, in order to avoid
                                # that the peak is split up again later on.
                        
                                detection_probability_final = np.zeros(2) + detection_probability_final[local_index]
                                rotation_probability_final = np.zeros(2) - 99.0
                                duplicity_probability_final = np.zeros(2) - 99.0
                                blending_profile_flag_final = np.zeros(2,dtype=int)
                                sinc_profile_flag_final = np.zeros(2,dtype=int)
                    
                        else:
                            # In this case rotational splitting is the favored scenario.
                            cosi = np.mean(cosi)
                            freq_triplet = [np.mean(central_freq)-np.mean(rot_split),np.mean(central_freq),np.mean(central_freq)+np.mean(rot_split)]
   
                            range_midpoint_left = (freq_triplet[0] + freq_triplet[1]) / 2.0
                            range_boundary_left = freq_triplet[0] - abs(freq_triplet[1]-freq_triplet[0])/2.
                            range_midpoint_right = (freq_triplet[1] + freq_triplet[2]) / 2.0
                            range_boundary_right = freq_triplet[2] + abs(freq_triplet[2]-freq_triplet[1])/2.
                    
                            range_maximum_new = np.array([[range_boundary_left,range_midpoint_left],[range_midpoint_left,range_midpoint_right],[range_midpoint_right,range_boundary_right]]).T
                            divisions_maximum_new = np.array([[divisions_maximum_final[0,local_index],range_midpoint_left],[range_midpoint_left,range_midpoint_right],[range_midpoint_right,divisions_maximum_final[1,local_index]]]).T
                            left_left_bound = range_maximum_new[0,0]
                            left_right_bound = range_maximum_new[1,0]
                            central_left_bound = range_maximum_new[0,1] 
                            central_right_bound = range_maximum_new[1,1] 
                            right_left_bound = range_maximum_new[0,2]
                            right_right_bound = range_maximum_new[1,2]
                           
                            
                            # For low frequency modes the multi-modal sampling from DIAMONDS may turn out to be insufficient
                            # to split up the mode in three peaks. If this happens, then take as frequency uncertainty the 
                            # uncertainty from the un-split peak divided by sqrt(3).
                            
                            # Left peak
                            tmp_range = np.where((par0 < left_right_bound) & (par0 >= left_left_bound))[0]
                            tmp_freq_range = np.where((freq < left_right_bound) & (freq >= left_left_bound))[0]
                            tmp_hist_range = np.where((par_hist < left_right_bound) & (par_hist >= left_left_bound))[0]
                    
                            if len(tmp_freq_range) == 0:
                                continue
                            spsd_range = spsd[tmp_freq_range]
                            spsd_maximum_left = max(spsd_range)

                            if len(tmp_range) > 0:
                                par0_range = par0[tmp_range]
                                sampling_counts_left = np.sum(nest_iter[tmp_range])
                                freq_sig_left = np.sqrt(np.sum((par0_range-freq_triplet[0])**2*tmp_range**2)/np.sum(tmp_range**2))
                            else:
                                right_bound_new = right_right_bound
                                left_bound_new = left_left_bound
                                tmp_range = np.where((par0 < right_bound_new) & (par0 >= left_bound_new))[0]
                                par0_range = par0[tmp_range]
                                sampling_counts_left = np.sum(nest_iter[tmp_range])/3.
                                freq_sig_left = np.sqrt(np.sum((par0_range-freq_triplet[0])**2*tmp_range**2)/np.sum(tmp_range**2))/np.sqrt(3.)
                    
 
                            if len(tmp_hist_range) > 0:
                                asef_maximum_left = max(asef_hist[tmp_hist_range])
                            else:
                                right_bound_new = right_right_bound
                                left_bound_new = left_left_bound
                                tmp_hist_range = np.where((par_hist < right_bound_new) & (par_hist >= left_bound_new))[0]
                                asef_maximum_left = max(asef_hist[tmp_hist_range])
                    
 
                            # Central peak
                            tmp_range = np.where((par0 < central_right_bound) & (par0 >= central_left_bound))[0]
                            tmp_freq_range = np.where((freq < central_right_bound) & (freq >= central_left_bound))[0]
                            tmp_hist_range = np.where((par_hist < central_right_bound) & (par_hist >= central_left_bound))[0]
                    
                            if len(tmp_freq_range)==0:
                                continue
                            spsd_range = spsd[tmp_freq_range]
                            spsd_maximum_central = max(spsd_range)
                    
                            if len(tmp_range) > 0:
                                par0_range = par0[tmp_range]
                                sampling_counts_central = np.sum(nest_iter[tmp_range])
                                freq_sig_central = np.sqrt(np.sum((par0_range-freq_triplet[1])**2*tmp_range**2)/np.sum(tmp_range**2))
                            else:
                                right_bound_new = right_right_bound
                                left_bound_new = left_left_bound
                                tmp_range = np.where((par0 < right_bound_new) & (par0 >= left_bound_new))[0]
                                par0_range = par0[tmp_range]
                                sampling_counts_central = np.sum(nest_iter[tmp_range])/3.
                                freq_sig_central = np.sqrt(np.sum((par0_range-freq_triplet[1])**2*tmp_range**2)/np.sum(tmp_range**2))/np.sqrt(3.)
                    
                    
                            if len(tmp_hist_range) > 0:
                                asef_maximum_central = max(asef_hist[tmp_hist_range])
                            else:
                                right_bound_new = right_right_bound
                                left_bound_new = left_left_bound
                                tmp_hist_range = np.where((par_hist < right_bound_new) & (par_hist >= left_bound_new))[0]
                                asef_maximum_central = max(asef_hist[tmp_hist_range])
                    
 
                            # Right peak
                            tmp_range = np.where((par0 < right_right_bound) & (par0 >= right_left_bound))[0]
                            tmp_freq_range = np.where((freq < right_right_bound) & (freq >= right_left_bound))[0]
                            tmp_hist_range = np.where((par_hist < right_right_bound) & (par_hist >= right_left_bound))[0]
                            par0_range = par0[tmp_range]
                    
                            if len(tmp_freq_range) == 0:
                                continue
                            spsd_range = spsd[tmp_freq_range]
                            spsd_maximum_right = max(spsd_range)

                            if len(tmp_range) > 0:
                                par0_range = par0[tmp_range]
                                sampling_counts_right = np.sum(nest_iter[tmp_range])
                                freq_sig_right = np.sqrt(np.sum((par0_range-freq_triplet[2])**2*tmp_range**2)/np.sum(tmp_range**2))
                            else:
                                right_bound_new = right_right_bound
                                left_bound_new = left_left_bound
                                tmp_range = np.where((par0 < right_bound_new) & (par0 >= left_bound_new))[0]
                                par0_range = par0[tmp_range]
                                sampling_counts_right = np.sum(nest_iter[tmp_range])/3.
                                freq_sig_right = np.sqrt(np.sum((par0_range-freq_triplet[2])**2*tmp_range**2)/np.sum(tmp_range**2))/np.sqrt(3.)                   
                    
                            if len(tmp_hist_range) > 0:
                                asef_maximum_right = max(asef_hist[tmp_hist_range])
                            else:
                                right_bound_new = right_right_bound
                                left_bound_new = left_left_bound
                                tmp_hist_range = np.where((par_hist < right_bound_new) & (par_hist >= left_bound_new))[0]
                                asef_maximum_right = max(asef_hist[tmp_hist_range])
                    

                            # Check whether there are frequencies other than this rotational multiplet 

                            if len(external_indices) > 0:
                                range_maximum_final = np.hstack((range_maximum_final[:,external_indices],range_maximum_new))
                                divisions_maximum_final = np.hstack((divisions_maximum_final[:,external_indices],divisions_maximum_new))
                                freq1_final = np.append(freq1_final[external_indices],freq_triplet)
                                freq_sig1_final = np.append(freq_sig1_final[external_indices],[freq_sig_left,freq_sig_central,freq_sig_right])
                                asef_maximum_final = np.append(asef_maximum_final[external_indices],[asef_maximum_left,asef_maximum_central,asef_maximum_right])
                                spsd_maximum_final = np.append(spsd_maximum_final[external_indices],[spsd_maximum_left,spsd_maximum_central,spsd_maximum_right])
                                sampling_counts_final = np.append(sampling_counts_final[external_indices],[sampling_counts_left,sampling_counts_central,sampling_counts_right])
                                angular_degree_final = np.append(angular_degree_final[external_indices],[1,1,1])
                                order_number_final = np.append(order_number_final[external_indices],[enn_radial-1,enn_radial-1,enn_radial-1])
                                azimuthal_number_final = np.append(azimuthal_number_final[external_indices],[-1,0,1])
                                cosi_final = np.append(cosi_final[external_indices],[cosi,cosi,cosi])
                                
                                # Assign the duplet the same detection probability, as well as the duplicity and rotation probabilities.
                        
                                detection_probability_final = np.append(detection_probability_final[external_indices],np.zeros(3) + detection_probability_final[local_index])
                                rotation_probability_final = np.append(rotation_probability_final[external_indices],np.zeros(3) + rotation_probability_final[local_index])
                                duplicity_probability_final = np.append(duplicity_probability_final[external_indices],np.zeros(3) + duplicity_probability_final[local_index])
                                blending_profile_flag_final = np.append(blending_profile_flag_final[external_indices],np.zeros(3,dtype=int))
                                sinc_profile_flag_final = np.append(sinc_profile_flag_final[external_indices],np.zeros(3,dtype=int))
                                
                                # Sort all arrays with increasing frequency order.
                                
                                sorted_indices = np.argsort(freq1_final)
                                range_maximum_final = range_maximum_final[:,sorted_indices]
                                divisions_maximum_final = divisions_maximum_final[:,sorted_indices]
                                freq1_final = freq1_final[sorted_indices]
                                freq_sig1_final = freq_sig1_final[sorted_indices]
                                asef_maximum_final = asef_maximum_final[sorted_indices]
                                spsd_maximum_final = spsd_maximum_final[sorted_indices]
                                sampling_counts_final = sampling_counts_final[sorted_indices]
                                angular_degree_final = angular_degree_final[sorted_indices]
                                order_number_final = order_number_final[sorted_indices]
                                azimuthal_number_final = azimuthal_number_final[sorted_indices]
                                cosi_final = cosi_final[sorted_indices]
                                detection_probability_final = detection_probability_final[sorted_indices]
                                rotation_probability_final = rotation_probability_final[sorted_indices]
                                duplicity_probability_final = duplicity_probability_final[sorted_indices]
                                blending_profile_flag_final = blending_profile_flag_final[sorted_indices]
                                sinc_profile_flag_final = sinc_profile_flag_final[sorted_indices]
                            else:
                                range_maximum_final = range_maximum_new
                                divisions_maximum_final = divisions_maximum_new
                                freq1_final = freq_triplet
                                freq_sig1_final = [freq_sig_left,freq_sig_central,freq_sig_right]
                                asef_maximum_final = [asef_maximum_left,asef_maximum_central,asef_maximum_right]
                                spsd_maximum_final = [spsd_maximum_left,spsd_maximum_central,spsd_maximum_right]
                                sampling_counts_final = [sampling_counts_left,sampling_counts_central,sampling_counts_right]
                                angular_degree_final = [1,1,1]
                                order_number_final = [enn_radial-1,enn_radial-1,enn_radial-1]
                                azimuthal_number_final = [-1,0,1]
                                cosi_final = [cosi,cosi,cosi]
                                
                                # Assign the duplet the same detection probability, but remove the duplicity and rotation test from it, in order to avoid
                                # that the peak is split up again later on.
                                
                                detection_probability_final = np.zeros(3) + detection_probability_final[local_index]
                                rotation_probability_final = np.zeros(3) + rotation_probability_final[local_index]
                                duplicity_probability_final = np.zeros(3) + duplicity_probability_final[local_index]
                                blending_profile_flag_final = np.zeros(3,dtype=int)
                                sinc_profile_flag_final = np.zeros(3,dtype=int)
                    

                        # Check if any of the new peaks of the duplet or triplet is falling within the l=2 mixed mode region. If this is the case,
                        # remove the corresponding peak(s) from the final list.
            
                        if flag_quadrupole_found==1: 
                            candidate_mixed = np.where((freq1_final <= freq_quadrupole_chunk + d02/self.cp.d02_scaling_merge_mixed) & (freq1_final >= freq_quadrupole_chunk - d02/self.cp.d02_scaling_merge_mixed) & (freq1_final != freq_quadrupole_chunk))[0]
                            tmp_mask = np.zeros(len(freq1_final),dtype=bool)
                            tmp_mask[candidate_mixed] = True 
                            good_freq = np.arange(len(freq1_final))
               
                            if len(candidate_mixed) > 0 :
                                range_maximum_final = range_maximum_final[:,good_freq]
                                divisions_maximum_final = divisions_maximum_final[:,good_freq]
                                freq1_final = freq1_final[good_freq]
                                freq_sig1_final = freq_sig1_final[good_freq]
                                asef_maximum_final = asef_maximum_final[good_freq]
                                spsd_maximum_final = spsd_maximum_final[good_freq]
                                sampling_counts_final = sampling_counts_final[good_freq]
                                angular_degree_final = angular_degree_final[good_freq]
                                order_number_final = order_number_final[good_freq]
                                azimuthal_number_final = azimuthal_number_final[good_freq]
                                cosi_final = cosi_final[good_freq]
                                
                                detection_probability_final = detection_probability_final[good_freq]
                                rotation_probability_final = rotation_probability_final[good_freq]
                                duplicity_probability_final = duplicity_probability_final[good_freq]
                                blending_profile_flag_final = blending_profile_flag_final[good_freq]
                                sinc_profile_flag_final = sinc_profile_flag_final[good_freq]
                                
            n_freq_final = len(freq1_final)
        else:
            # If here, then it means that no significant modes were detected in the whole chunk.
            n_freq_final = 0
            freq1_final = []
            freq_sig1_final = []
            angular_degree_final = []
            order_number_final = []
            azimuthal_number_final = []
            

        separations = np.zeros(n_freq_final + 1)
        if self.cp.print_on_screen:
            # Update smoothing window using fitted FWHM of the chunk radial mode
            smth_bins = int(fwhm_radial_fit/1.5/freqbin)
            spsd_total = smooth(psd_total, window_len=smth_bins, window='flat')
            tmp_good = np.where((freq_total <= max(par_hist)) & (freq_total >= min(par_hist)))[0]
            spsd = spsd_total[tmp_good]
            freq = freq_total[tmp_good]
            psd = psd_total[tmp_good]

        
        # Print summary infomation of the run
        if self.cp.print_on_screen:
            if n_freq_final != 0:
                print('\nFrequency     sig\t n\t l\t m')
                for i in range(0, n_freq_final):
                    print('{:.4f}  {:10.4f}\t {:d}\t {:d}\t {:d} '.format(freq1_final[i],freq_sig1_final[i],order_number_final[i],angular_degree_final[i],azimuthal_number_final[i]))
                print(' ------------------------------------------\n')
                    

        # Save final outputs
        if self.cp.save_complete_lists:
            # Save total list of frequencies, without selection of significant peaks.
            with open(str(peakbagging_filename_chunk) + run + '_' + self.modality + '.all.txt','w') as f:
                f.write('# Col 1: n, Col 2: l, Col 3: Frequency (microHz), Col 4: 1-sigma Frequency (microHz), \n')
                f.write('# Col 5: Left frequency range (microHz), Col 6: Right frequency range (microHz), \n')
                f.write('# Col 7: Left frequency division (microHz), Col 8: Right frequency division (microHz), \n')
                f.write('# Col 9: ASEF maximum (nested iterations), Col 10: Sampling counts (counts), Col 11: SPSD maximum (ppm^2/microHz), \n')
                f.write('# Col 12: P (detection), Col 13: P (rotation), Col 14: P (duplicity), Col 15: Peak blending flag, Col 16: Sinc^2 profile flag.\n')

                for i in range(0,n_freq):
                    f.write('{:d} {:d}  {:15.5f}  {:12.5f}  {:15.5f}  {:15.5f}  {:15.5f}  {:15.5f}  {:15.2f}  {:15d}  {:15.5e}  {:12.3f}  {:12.3f}  {:10d}  {:10d}\n'.format(order_number[i],angular_degree[i],freq1[i],freq_sig1[i],range_maximum[0,i],range_maximum[1,i],divisions_maximum[0,i],divisions_maximum[1,i],asef_maximum[i],sampling_counts[i],spsd_maximum[i],detection_probability[i],rotation_probability[i],duplicity_probability[i],blending_profile_flag[i],sinc_profile_flag[i]))


        if n_freq_final != 0:
            # Save local asymptotic parameters and the final list of frequencies with peak significance and rotation tests applied.
            with open((str(peakbagging_filename_chunk)) + run + '_' + self.modality + '.txt','w') as f:
                f.write('# Local epsilon, Local d02 (microHz), local DeltaP, FWHM radial mode (microHz), FWHM octupole mode (microHz), FWHM multi-modal fit (microHz)\n')
                f.write('{:.3f}  {:10.3f}  {:10.3f}  {:10.4f}  {:10.4f}  {:10.4f} \n'.format(local_epsi,local_d02,local_dp,fwhm_radial_fit,largest_octupole_fwhm,fit_linewidth))
                
                # Save the individual oscillation frequencies from the global fit, their uncertainties, mode identification    
                f.write('# Col 1: n, Col 2: l, Col 3: m, Col 4: Frequency (microHz), Col 5: 1-sigma Frequency (microHz), \n')
                f.write('# Col 6: Left frequency range (microHz), Col 7: Right frequency range (microHz), \n')
                f.write('# Col 8: Left frequency division (microHz), Col 9: Right frequency division (microHz), \n')
                f.write('# Col 10: ASEF maximum (nested iterations), Col 11: Sampling counts (counts), Col 12: SPSD maximum (ppm^2/microHz), Col 13: cosi \n')
                f.write('# Col 14: P (detection), Col 15: P (rotation), Col 16: P (duplicity), Col 17: Peak blending profile flag, Col 18: Sinc^2 profile flag. \n')

                for i in range(0, n_freq_final):
                    f.write('{:2d} {:5d}  {:8.1f}  {:15.5f}  {:12.5f}  {:15.5f}  {:15.5f}  {:15.5f}  {:15.5f}  {:15.2f}  {:15d}  {:15.5e}  {:12.3f}  {:12.3f}  {:12.3f}  {:12.3f}  {:10d}  {:10d}\n'.format(order_number_final[i],angular_degree_final[i],azimuthal_number_final[i],freq1_final[i],freq_sig1_final[i],range_maximum_final[0,i],range_maximum_final[1,i],divisions_maximum_final[0,i],divisions_maximum_final[1,i],asef_maximum_final[i],int(sampling_counts_final[i]),spsd_maximum_final[i],cosi_final[i],detection_probability_final[i],rotation_probability_final[i],duplicity_probability_final[i],blending_profile_flag_final[i],sinc_profile_flag_final[i]))

        # Values to track
        self.local_epsi[run_number] = local_epsi
        self.local_d02[run_number] = local_d02
        self.local_dp[run_number] = local_dp
        self.low_cut_frequency[run_number] = low_cut_frequency
        
        self.n_radial_chunk[run_number] = n_radial_chunk
        self.fwhm_radial_fit[run_number] = fwhm_radial_fit 
        self.avg_fwhm[run_number] = avg_fwhm
        self.upper_height[run_number] = upper_height
        self.fit_linewidth[run_number] = fit_linewidth
        
        self.n_freqs[run_number] = n_freq_final
        self.orders[run_number] = order_number_final
        self.degrees[run_number] = angular_degree_final
        self.freqs[run_number] = freq1_final
        self.freqs_sig[run_number] = freq_sig1_final
        self.freqs_radial_chunk[run_number] = freq_radial_chunk
        self.octupole_freq_asymp[run_number] = octupole_freq_asymp
        self.octupole_freq_lower[run_number] = octupole_freq_lower
        self.octupole_freq_upper[run_number] = octupole_freq_upper
        
        self.par_hist[run_number] = par_hist
        self.asef_hist[run_number] = asef_hist
        self.asef_bins[run_number] = n_bins
        
                    
        # Save stuff into pickle for later steps...
        if self.cp.save_progress_pickle:
            pickle.dump(self,open(self.star_dir/(self.catalog_id+self.star_id+'_chunk.pickle'),'wb'))

        # Return True if we made it all the way through
        return True
                    

    def make_chunk_plots(self,chunk_number=None):
        """
        Produce all plots related to the CHUNK modality and save as desired.
        """
        plt.style.use(self.cp.famed_path/self.cp.mplstyle)
        
        if chunk_number is not None:
            print('Making plot for chunk:',chunk_number)
            famed_plots.chunk_plot(self,chunk_number)
        else:
            print('Need to provide a chunk number to plot.')
            return
            
        if self.cp.save_png:
            plt.savefig(self.star_dir/self.cp.figs_subdir/(self.catalog_id+self.star_id+'_'+self.cp.isla_subdir+'_'+str(chunk_number)+'_CHUNK.png'))
        if self.cp.save_eps:
            plt.savefig(self.star_dir/self.cp.figs_subdir/(self.catalog_id+self.star_id+'_'+self.cp.isla_subdir+'_'+str(chunk_number)+'_CHUNK.eps'))
        
