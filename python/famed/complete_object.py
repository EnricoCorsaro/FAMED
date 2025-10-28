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

__all__ = ['Complete']

class Complete(FamedStar):
    """
    FamedStar sub-class specific to the running of the CHUNK modality.    

    Parameters
    ----------
    catalog_id : str
        Catalogue ID of the star (e.g. 'KIC' for Kepler).
    star_id : str
        ID of the star as a string (e.g. '0012008916' or '7037405').
    background_run_number : str or int
        Number of the background subfolder that contains the results of
        the background fit.
    load_islands : bool, default: False
        Flag to read the pickled object if `make_islands` or `find_islands` has 
        already been ran.
    """
    def __init__(self, catalog_id, star_id, background_run_number=None):
        FamedStar.__init__(self, catalog_id, star_id, background_run_number=background_run_number)
       
        if self.cp.print_on_screen:
            print(' Loading saved information from the pickled object.\n')

        # Load any variables into attributes that should haven been
        # loaded or created during CHUNK modality. Will ignore new config

        temp = pickle.load(open(self.star_dir/(self.catalog_id+self.star_id+'_chunk.pickle'),'rb'))
        self.__dict__ = temp.__dict__.copy()


    def make_complete(self,run,fit=True):
        """
        Compute a high-dimensional fit with DIAMONDS for each chunk specified.

        This function generates priors and executes a chunk high-dimensional fit in order to retrieve the individual
        oscillation peak parameters (frequencies, amplitudes, linewidiths) and their uncertainties with the
        best precision achievable from the data. This function cannot be executed if both the GLOBAL and CHUNK 
        modalities are not previously performed.

        Parameters
        ----------
        run : int or str
            An integer specifying the number of the chunk for which the sampling must be evaluated. 
            If a number < 0 is provided, all the chunks are taken into account at the same time.
        fit : bool
            A flag to specify whether only prior configuring files must be generated, or also the actual high-dimensional 
            fits with DIAMONDS should be executed. Default is set to execute also the fits.
        """
        run = int(run)
        peakbagging_filename_global = self.star_dir/self.cp.summary_subdir/(self.catalog_id + self.star_id + self.cp.peakbagging_filename_label + self.cp.isla_subdir + '_' + self.cp.global_subdir + '_GLOBAL.txt')
        peakbagging_filename_chunk = self.star_dir/self.cp.summary_subdir/(self.catalog_id + self.star_id + self.cp.peakbagging_filename_label + self.cp.isla_subdir + '_')

        if (run < 0):
            filename_summary = np.sort(glob.glob(str(peakbagging_filename_chunk) + '*' + self.modality + '.txt'))
        else:
            filename_summary = np.sort(glob.glob(str(peakbagging_filename_chunk) + str(run) + '_' + self.modality + '.txt'))

        n_chunks = len(filename_summary)
        bg_names = []
        run_labels = []
        fwhm_values = []

        tail_str = str('_' + self.modality + '.txt')
        len_tail = len(tail_str)
        len_head = len(str(peakbagging_filename_chunk))

        if self.cp.print_on_screen:
            print(' Number of chunks available for high-dimensional fitting: ', n_chunks, '\n')

        self.chunk_number_complete = [None]*n_chunks
        self.enn = [None]*n_chunks
        self.ell = [None]*n_chunks
        self.emm = [None]*n_chunks
        self.p_det = [None]*n_chunks
        self.mode_freq = [None]*n_chunks
        self.mode_amp = [None]*n_chunks
        self.mode_fwhm = [None]*n_chunks
        self.mode_freq_sig_low = [None]*n_chunks
        self.mode_freq_sig_up = [None]*n_chunks
        self.mode_amp_sig_low = [None]*n_chunks
        self.mode_amp_sig_up = [None]*n_chunks
        self.mode_fwhm_sig_low = [None]*n_chunks
        self.mode_fwhm_sig_up = [None]*n_chunks
        self.freq_left_margin = [None]*n_chunks
        self.freq_right_margin = [None]*n_chunks
        self.n_chunks = n_chunks

        best_dnu,n_total_chunks = np.loadtxt(peakbagging_filename_global,max_rows=1,skiprows=1,usecols=(2,7))
        n_total_chunks = int(n_total_chunks)
        freq_left,freq_right = np.loadtxt(peakbagging_filename_global,max_rows=n_total_chunks,skiprows=3,usecols=(1,2),unpack=True)
       
        for i in range(0,n_chunks):
            fwhm_0 = np.loadtxt(filename_summary[i], usecols=(3,), skiprows=1, max_rows=1)
            enn,ell,emm,nu,nu_left,nu_right,nu_left_division,nu_right_division,spsd_max,p_det,profile_flag = np.loadtxt(filename_summary[i],usecols=(0,1,2,3,5,6,7,8,11,13,17), skiprows=7, unpack=True,ndmin=2)
            tmp_detected = np.where((p_det >= self.cp.detection_probability_threshold) | (p_det == -99.0))[0]
        
            if len(tmp_detected) > 0:
                run_subdir = filename_summary[i][len_head:-len_tail]
                run_labels = np.append(run_labels, run_subdir)
                bg_names = np.append(bg_names, self.bgp['name'])
                fwhm_values = np.append(fwhm_values, -1.0)

                enn = enn[tmp_detected]
                ell = ell[tmp_detected]
                emm = emm[tmp_detected]
                nu = nu[tmp_detected]
                nu_left = nu_left[tmp_detected]
                nu_right = nu_right[tmp_detected]
                nu_left_division = nu_left_division[tmp_detected]
                nu_right_division = nu_right_division[tmp_detected]
                spsd_max = spsd_max[tmp_detected]
                profile_flag = profile_flag[tmp_detected]

                amp = np.sqrt(spsd_max*fwhm_0*np.pi)*np.sqrt(2.0)
       
                if (best_dnu < self.cp.dnu_rg) or ((best_dnu < self.cp.dnu_sg) & (teff < self.cp.teff_sg)):
                    freq_left_margin = freq_left[int(run_subdir)] - best_dnu*self.cp.dnu_overlap_fraction_rg
                else:
                    freq_left_margin = freq_left[int(run_subdir)] - best_dnu*self.cp.dnu_overlap_fraction_ms
                
                freq_right_margin = freq_right[int(run_subdir)]

                if ((self.cp.print_on_screen) and (run >=0)):
                    print(" Generating prior hyper-parameters for the high-dimensional fit with DIAMONDS for chunk #:", run_subdir)

                with open(self.star_dir/self.cp.complete_subdir/(self.cp.prior_filename + '_' + run_subdir + '.txt'),'w') as f:
            
                    for jj in range(0,len(nu)):
                        # Write the prior hyper-parameters file containing all the peaks to be fit within this chunk
                        f.write('{:.5f}  \t{:.5f}\n'.format(nu_left[jj],nu_right[jj]))
                        f.write('{:.5f}  \t{:.5f}\n'.format(amp[jj]*0.1,amp[jj]*1.4))
                        
                        if (profile_flag[jj] == 1):
                            f.write('{:.5f}  \t{:.5f}\n'.format(self.cp.fwhm_lower_bound,fwhm_0*0.4))
                        else:
                            f.write('{:.5f}  \t{:.5f}\n'.format(self.cp.fwhm_lower_bound,fwhm_0*1.2))

                        f.write('\n')

                with open(self.star_dir/self.cp.complete_subdir/('frequencyRange_' + run_subdir + '.txt'),'w') as f:
                    f.write('{:.5f}\n'.format(freq_left_margin))
                    f.write('{:.5f}\n'.format(freq_right_margin))

                self.freq_left_margin[i] = freq_left_margin
                self.freq_right_margin[i] = freq_right_margin
                self.chunk_number_complete[i] = int(run_subdir)
                self.enn[i] = enn
                self.ell[i] = ell
                self.emm[i] = emm
                self.p_det[i] = p_det

        if self.cp.save_progress_pickle:
            self.modality = 'COMPLETE'
            pickle.dump(self,open(self.star_dir/(self.catalog_id+self.star_id+'_complete.pickle'),'wb'))

        if fit:
            if n_chunks == 1:
                run_labels = run_labels[0]
                bg_names = bg_names[0]
                fwhm_values = fwhm_values[0]

            if self.cp.print_on_screen:
                print(" Performing the high-dimensional fit with DIAMONDS for #: ", n_chunks, " chunks.\n")
                print(" This operation may take up to several minutes...")

            peakbagging_parameters = { 'subdir':     self.cp.complete_subdir,
                                        'run':        run_labels,
                                        'background': bg_names,
                                        'fwhm':       fwhm_values
                                        }

            flag_computation_completed = diamonds.run_peakbagging(self.catalog_id,self.star_id,peakbagging_parameters,0,0,0,self.cp.dp_pb,self.cp.diamonds_path, self.cp.n_threads, self.cp.prior_filename,merge=True)

        return self.n_chunks, self.chunk_number_complete

    def prepare_results(self,run):
        """
        Prepare the output solution from the a high-dimensional fit with DIAMONDS for each chunk specified.

        This function generates summary ASCII files contaning frequency, amplitude, and linewidth for each
        oscillation mode (and related 1-sigma uncertainties), along with its corresponding radial order, 
        angular degree, and azimuthal number, and with the associated detection probability. 
        This function delivers the output even if the DIAMONDS fit did not converge by producing a resulting
        parameter summary file.

        Parameters
        ----------
        run : int or str
            An integer specifying the number of the chunk for which the sampling must be evaluated. 
            If a number < 0 is provided, all the chunks are taken into account at the same time.
        """
        run = int(run)

        if (run >= 0):
            index = np.where(np.array(self.chunk_number_complete) == run)[0]
            
            if (len(index) > 0):
                peakbagging_filename_complete = self.star_dir/self.cp.complete_subdir/str(run)
                peakbagging_filename_chunk_complete = self.star_dir/self.cp.summary_subdir/(self.catalog_id + self.star_id + self.cp.peakbagging_filename_label + self.cp.complete_subdir + '_')
                     
                # Check if the given run exists. If not, exit from the routine
                if not os.path.isfile(peakbagging_filename_complete/'peakbagging_computationParameters.txt'):
                    if self.cp.print_on_screen:
                        print(" The high-dimensional fit for this chunk was not executed. Quitting program.")
                    return False
                else:
                    if not os.path.isfile(peakbagging_filename_complete/'peakbagging_parameterSummary.txt'):
                        if self.cp.print_on_screen:
                            print(' CHUNK #:', str(run), 'PeakBagging fit did not produce Summary file.')
                            print(' Using sampling evolution and posterior values to compute results.')
                            
                            parameter_filenames = np.sort(glob.glob(str(peakbagging_filename_complete) + '/peakbagging_parameter0*.txt'))
                            n_params = parameter_filenames.size
                            params = np.zeros(n_params)
                            n_modes = int(n_params/3)
                            lower_error = np.zeros(n_params)
                            posterior = np.loadtxt(str(peakbagging_filename_complete) + '/peakbagging_posteriorDistribution.txt')
                            posterior /= posterior.max()

                            for par in range(0,n_params):
                                if par < 10:
                                    parstr = '0' + str(par)
                                else:
                                    parstr = str(par)

                                par_sampling = np.loadtxt(str(peakbagging_filename_complete) + '/peakbagging_parameter0' + parstr + '.txt')
                                params[par] = ((posterior*par_sampling).sum())/posterior.sum()
                                lower_error[par] = np.sqrt(((np.square(par_sampling - params[par])*posterior).sum())/posterior.sum())
                                upper_error = lower_error
                    else:
                        if self.cp.print_on_screen:
                            print(' CHUNK #:', str(run), 'Using parameterSummary file from PeakBagging fit.')

                        params,lowerpar,upperpar = np.loadtxt(str(peakbagging_filename_complete) + '/peakbagging_parameterSummary.txt',unpack=True,usecols=(1,4,5))   # Median value of the free parameter
                        lower_error = params - lowerpar
                        upper_error = upperpar - params
                        n_params = len(params)
                        n_modes = int(n_params/3)

                    index = int(index)
                    enn = self.enn[index]
                    ell = self.ell[index]
                    emm = self.emm[index]
                    p_det = self.p_det[index]

                    freq = params[0:n_params:3]
                    amp = params[1:n_params:3]
                    fwhm = params[2:n_params:3]
                    freq_sig_low = lower_error[0:n_params:3]
                    freq_sig_up = upper_error[0:n_params:3] 
                    amp_sig_low = lower_error[1:n_params:3]
                    amp_sig_up = upper_error[1:n_params:3] 
                    fwhm_sig_low = lower_error[2:n_params:3]
                    fwhm_sig_up = upper_error[2:n_params:3]

                    print(" Number of modes that will be printed:", n_modes, "\n")

                    self.mode_freq[index] = freq
                    self.mode_amp[index] = amp
                    self.mode_fwhm[index] = fwhm
                    self.mode_freq_sig_low[index] = freq_sig_low
                    self.mode_freq_sig_up[index] = freq_sig_up
                    self.mode_amp_sig_low[index] = amp_sig_low
                    self.mode_amp_sig_up[index] = amp_sig_up
                    self.mode_fwhm_sig_low[index] = fwhm_sig_low
                    self.mode_fwhm_sig_up[index] = fwhm_sig_up

                    # Save peakbagging mode parameters with peak significance tests applied.
                    with open(str(peakbagging_filename_chunk_complete) + str(run) + '_' + self.modality + '.txt','w') as f:
                        # Save the individual oscillation frequencies from the global fit, their uncertainties, mode identification    
                        f.write('# Col 1: n, Col 2: l, Col 3: m, Col 4: Frequency (uHz), Col 5: Lower 1-sigma Frequency (uHz), \n')
                        f.write('# Col 6: Upper 1-sigma Frequency (uHz), Col 7: Amplitude (ppm), \n')
                        f.write('# Col 8: Lower 1-sigma Amplitude (ppm), Col 9: Upper 1-sigma Amplitude (ppm), \n')
                        f.write('# Col 10: Linewidth (uHz), Col 11: Lower 1-sigma Linewidth (uHz), Col 12: Upper 1-sigma Linewidth (uHz), Col 13: P (detection) \n')

                        for i in range(0, n_modes):
                            f.write('{:2d} {:5d}  {:8.1f}  {:15.5f}  {:12.5f}  {:12.5f}  {:15.5f}  {:12.5f}  {:12.5f}  {:15.5f}  {:12.5f}  {:12.5f}  {:12.3f}\n'.format(int(enn[i]),int(ell[i]),emm[i],freq[i],freq_sig_low[i],freq_sig_up[i],amp[i],amp_sig_low[i],amp_sig_up[i],fwhm[i],fwhm_sig_low[i],fwhm_sig_up[i],p_det[i]))
            else:
                if self.cp.print_on_screen:
                    print(" The requested chunk is not in the list of fitted chunks. Quitting program.")
                    return False

        if self.cp.save_progress_pickle:
            # Update pickle with recovered fit information, if applicable
            pickle.dump(self,open(self.star_dir/(self.catalog_id+self.star_id+'_complete.pickle'),'wb'))

        return True

    def make_complete_plots(self,run):
        """
        Produce all the plots for each chunk that show the high-dimensional fit performed on the stellar PSD
        and save them as output according to the requested format.

        Parameters
        ----------
        run : int or str
            An integer specifying the number of the chunk for which the sampling must be evaluated. 
            If a number < 0 is provided, all the chunks are taken into account at the same time.
        """
      
        plt.style.use(self.cp.famed_path/self.cp.mplstyle)
        print(' Making plot for chunk:', run)
        famed_plots.complete_plot(self,run)
        
        if self.cp.save_png:
            plt.savefig(self.star_dir/self.cp.figs_subdir/(self.catalog_id+self.star_id+'_'+self.cp.complete_subdir+'_'+str(run)+'_COMPLETE.png'))
        if self.cp.save_eps:
            plt.savefig(self.star_dir/self.cp.figs_subdir/(self.catalog_id+self.star_id+'_'+self.cp.complete_subdir+'_'+str(run)+'_COMPLETE.eps')) 



