pro setup_computation
; -------------------------------------------------------------------------------------------------------
; Author:  Enrico Corsaro
; e-mail:  enrico.corsaro@inaf.it
; Date:    February 2019
; Place:   Catania, Italy
; Purpose: This procedure sets up all the configuring parameters needed to run the FAMED pipeline.
;          It is not depending on the star that is selected, but requires that the input file
;          famed_configuring_parameters.txt is available. It also sets up all parameters related to
;          the plotting of the results into output devices X-term and EPS.
; -------------------------------------------------------------------------------------------------------
COMMON CONFIG,cp
COMMON STAR,info
COMMON GRAPHIC,pp,lp,sp,ppe,lpe
COMMON DIAMONDS,dp

readcol,'famed_configuring_parameters.txt',par_name,par_value,format='A,A',comment='#',/silent
root_path = (par_value(where(par_name eq 'root_path')))[0]
global_string = (par_value(where(par_name eq 'global_subdir')))[0]
save_png = (par_value(where(par_name eq 'save_png')))[0]
save_eps = (par_value(where(par_name eq 'save_eps')))[0]

if save_eps eq 1 then begin
    if save_png eq 1 then save_png = 0
endif


; -------------------------------------------------------------------------------------------------------
; Setup input and output directories and filenames for computation 
; -------------------------------------------------------------------------------------------------------
peakbagging_data_dir = root_path + 'PeakBagging/data/'
peakbagging_results_dir =  root_path + 'PeakBagging/results/'
background_data_dir = root_path + 'Background/data/'
background_results_dir = root_path + 'Background/results/'
peakbagging_filename_label = '_peakbagging_'
local_configuring_parameters_filename = 'famed_configuring_parameters_'

info = {   root_path:                             root_path,                                                                       $
           peakbagging_path:                      root_path + 'PeakBagging/build',                                                 $
           asymptotic_path:                       root_path + 'Asymptotic/build',                                                  $
           configuring_parameters_filename:       'famed_configuring_parameters.txt',                                              $
           local_configuring_parameters_filename: local_configuring_parameters_filename,                                           $
           famed_path:                            (par_value(where(par_name eq 'famed_path')))[0],                                 $
           isla_subdir:                           (par_value(where(par_name eq 'isla_subdir')))[0],                                $
           pb_subdir:                             (par_value(where(par_name eq 'pb_subdir')))[0],                                  $
           as_subdir:                             (par_value(where(par_name eq 'as_subdir')))[0],                                  $
           figs_subdir:                           (par_value(where(par_name eq 'figs_subdir')))[0],                                $
           summary_subdir:                        (par_value(where(par_name eq 'summary_subdir')))[0],                             $
           global_subdir:                         global_string,                                                                   $
           save_complete_lists:                   (par_value(where(par_name eq 'save_complete_lists')))[0],                        $
           save_test_files:                       (par_value(where(par_name eq 'save_test_files')))[0],                            $
           save_png:                              save_png,                                                                        $
           save_eps:                              save_eps,                                                                        $
           print_on_screen:                       (par_value(where(par_name eq 'print_on_screen')))[0],                            $
           background_run_number:                 (par_value(where(par_name eq 'background_run_number')))[0],                      $
           n_threads:                             (par_value(where(par_name eq 'n_threads')))[0],                                  $
           logo_filename:                         (par_value(where(par_name eq 'logo_filename')))[0],                              $
           prior_filename:                        (par_value(where(par_name eq 'prior_filename')))[0],                             $
           peakbagging_data_dir:                  peakbagging_data_dir,                                                            $
           background_data_dir:                   background_data_dir,                                                             $
           peakbagging_results_dir:               peakbagging_results_dir,                                                         $
           background_results_dir:                background_results_dir,                                                          $
           external_background_results_dir:       (par_value(where(par_name eq 'external_background_results_dir')))[0],            $
           external_background_filename_suffix:   (par_value(where(par_name eq 'external_background_filename_suffix')))[0],        $
           peakbagging_filename_label:            peakbagging_filename_label                                                       $
       }

; -------------------------------------------------------------------------------------------------------
; Setup configuring parameters of the pipeline as given from an input file.
; -------------------------------------------------------------------------------------------------------
cp = { teff_sun:                   (float(par_value(where(par_name eq 'teff_sun'))))[0],           $
       dnu_sun:                    (float(par_value(where(par_name eq 'dnu_sun'))))[0],            $
       numax_sun:                  (float(par_value(where(par_name eq 'numax_sun'))))[0],          $
       ; Multi-modal sampling production
       n_sigma_envelope:           (float(par_value(where(par_name eq 'n_sigma_envelope'))))[0],                   $
       n_sigma_envelope_cl:        (float(par_value(where(par_name eq 'n_sigma_envelope_cl'))))[0],                $
       n_sigma_envelope_tip:       (float(par_value(where(par_name eq 'n_sigma_envelope_tip'))))[0],               $
       n_dnu_envelope:             (float(par_value(where(par_name eq 'n_dnu_envelope'))))[0],                     $
       dnu_tip:                    (float(par_value(where(par_name eq 'dnu_tip'))))[0],                            $
       dnu_cl:                     (float(par_value(where(par_name eq 'dnu_cl'))))[0],                             $
       dnu_cl2:                    (float(par_value(where(par_name eq 'dnu_cl2'))))[0],                            $
       dnu_rg:                     (float(par_value(where(par_name eq 'dnu_rg'))))[0],                             $
       dnu_sg:                     (float(par_value(where(par_name eq 'dnu_sg'))))[0],                             $
       teff_sg:                    (float(par_value(where(par_name eq 'teff_sg'))))[0],                            $ 
       fwhm_global_scaling:            (float(par_value(where(par_name eq 'fwhm_global_scaling'))))[0],            $
       fwhm_global_scaling_tip:        (float(par_value(where(par_name eq 'fwhm_global_scaling_tip'))))[0],        $
       fwhm_chunk_scaling_ms:          (float(par_value(where(par_name eq 'fwhm_chunk_scaling_ms'))))[0],          $
       fwhm_chunk_scaling_sg:          (float(par_value(where(par_name eq 'fwhm_chunk_scaling_sg'))))[0],          $
       fwhm_chunk_scaling_rg:          (float(par_value(where(par_name eq 'fwhm_chunk_scaling_rg'))))[0],          $
       fwhm_chunk_scaling_cl:          (float(par_value(where(par_name eq 'fwhm_chunk_scaling_cl'))))[0],          $
       fwhm_chunk_scaling_tip:         (float(par_value(where(par_name eq 'fwhm_chunk_scaling_tip'))))[0],         $
       dnu_overlap_fraction_ms:        (float(par_value(where(par_name eq 'dnu_overlap_fraction_ms'))))[0],        $
       dnu_overlap_fraction_rg:        (float(par_value(where(par_name eq 'dnu_overlap_fraction_rg'))))[0],        $
       threshold_asef_global:          (float(par_value(where(par_name eq 'threshold_asef_global'))))[0],          $
       threshold_asef_chunk_ms:        (float(par_value(where(par_name eq 'threshold_asef_chunk_ms'))))[0],        $
       threshold_asef_chunk_sg:        (float(par_value(where(par_name eq 'threshold_asef_chunk_sg'))))[0],        $
       threshold_asef_chunk_rg:        (float(par_value(where(par_name eq 'threshold_asef_chunk_rg'))))[0],        $
       skim_frequency_tolerance_ms:    (float(par_value(where(par_name eq 'skim_frequency_tolerance_ms'))))[0],    $
       skim_frequency_tolerance_sg:    (float(par_value(where(par_name eq 'skim_frequency_tolerance_sg'))))[0],    $
       skim_frequency_tolerance_rg:    (float(par_value(where(par_name eq 'skim_frequency_tolerance_rg'))))[0],    $
       skim_frequency_tolerance_tip:   (float(par_value(where(par_name eq 'skim_frequency_tolerance_tip'))))[0],   $
       ; Multi-modal sampling analysis
       min_sep_scaling_global:     (float(par_value(where(par_name eq 'min_sep_scaling_global'))))[0],                     $
       min_sep_scaling_chunk_ms:   (float(par_value(where(par_name eq 'min_sep_scaling_chunk_ms'))))[0],                   $
       min_sep_scaling_chunk_rg:   (float(par_value(where(par_name eq 'min_sep_scaling_chunk_rg'))))[0],                   $
       min_n_bins:                 (float(par_value(where(par_name eq 'min_n_bins'))))[0],                                 $
       n_bins_acf_scaling:         (float(par_value(where(par_name eq 'n_bins_acf_scaling'))))[0],                         $
       dnu_acf_range_side:         (float(par_value(where(par_name eq 'dnu_acf_range_side'))))[0],                         $
       min_bin_separation:         (fix(par_value(where(par_name eq 'min_bin_separation'))))[0],                           $
       n_bins_max_fraction:        (float(par_value(where(par_name eq 'n_bins_max_fraction'))))[0],                        $
       drop_tolerance_global:      (float(par_value(where(par_name eq 'drop_tolerance_global'))))[0],                      $
       drop_tolerance_chunk_ms:       (float(par_value(where(par_name eq 'drop_tolerance_chunk_ms'))))[0],                 $
       drop_tolerance_chunk_rg:       (float(par_value(where(par_name eq 'drop_tolerance_chunk_rg'))))[0],                 $
       max_iterations_frequency:   (fix(par_value(where(par_name eq 'max_iterations_frequency'))))[0],                     $
       max_sigma_range:            (float(par_value(where(par_name eq 'max_sigma_range'))))[0],                            $
       min_sigma_range:            (float(par_value(where(par_name eq 'min_sigma_range'))))[0],                            $
       max_skim_iterations_global:     (fix(par_value(where(par_name eq 'max_skim_iterations_global'))))[0],               $
       alpha_radial_universal:         (float(par_value(where(par_name eq 'alpha_radial_universal'))))[0],                 $
       correct_radial_frequencies:     (float(par_value(where(par_name eq 'correct_radial_frequencies'))))[0],             $
       save_asymptotic_radial:         (float(par_value(where(par_name eq 'save_asymptotic_radial'))))[0],                 $
       separations_dnu_tolerance_rg:       (float(par_value(where(par_name eq 'separations_dnu_tolerance_rg'))))[0],       $
       separations_dnu_tolerance_ms:       (float(par_value(where(par_name eq 'separations_dnu_tolerance_ms'))))[0],       $
       weight_freq_fraction:               (float(par_value(where(par_name eq 'weight_freq_fraction'))))[0],               $
       weight_freq_fraction_enhanced:      (float(par_value(where(par_name eq 'weight_freq_fraction_enhanced'))))[0],      $
       weight_asef_fraction:               (float(par_value(where(par_name eq 'weight_asef_fraction'))))[0],               $
       weight_spsd_fraction:               (float(par_value(where(par_name eq 'weight_spsd_fraction'))))[0],               $
       weight_sampling_fraction:           (float(par_value(where(par_name eq 'weight_sampling_fraction'))))[0],           $
       upper_limit_freq_radial:            (float(par_value(where(par_name eq 'upper_limit_freq_radial'))))[0],            $
       plot_weights_radial:                (fix(par_value(where(par_name eq 'plot_weights_radial'))))[0],                  $
       max_ratio_search_radial:            (float(par_value(where(par_name eq 'max_ratio_search_radial'))))[0],            $
       previous_radial_range_fraction:     (float(par_value(where(par_name eq 'previous_radial_range_fraction'))))[0],     $
       sampling_counts_fraction:           (float(par_value(where(par_name eq 'sampling_counts_fraction'))))[0],           $
       asef_saturation_fraction:           (float(par_value(where(par_name eq 'asef_saturation_fraction'))))[0],           $
       asef_threshold_scaling_radial:      (float(par_value(where(par_name eq 'asef_threshold_scaling_radial'))))[0],      $
       dnu_lower_cut_fraction:             (float(par_value(where(par_name eq 'dnu_lower_cut_fraction'))))[0],             $
       low_cut_frequency:                  (float(par_value(where(par_name eq 'low_cut_frequency'))))[0],                  $
       d02_scaling_merge_mixed:            (float(par_value(where(par_name eq 'd02_scaling_merge_mixed'))))[0],            $
       d02_factor_search_range:            (float(par_value(where(par_name eq 'd02_factor_search_range'))))[0],            $
       d02_prior_lower_duplet_fit:         (float(par_value(where(par_name eq 'd02_prior_lower_duplet_fit'))))[0],         $
       d02_prior_upper_duplet_fit:         (float(par_value(where(par_name eq 'd02_prior_upper_duplet_fit'))))[0],         $
       d03_upper_scaling_factor:           (float(par_value(where(par_name eq 'd03_upper_scaling_factor'))))[0],           $
       d03_lower_scaling_factor:           (float(par_value(where(par_name eq 'd03_lower_scaling_factor'))))[0],           $
       d02_upper_scaling_factor_ms:        (float(par_value(where(par_name eq 'd02_upper_scaling_factor_ms'))))[0],        $
       d02_lower_scaling_factor_ms:        (float(par_value(where(par_name eq 'd02_lower_scaling_factor_ms'))))[0],        $
       d02_upper_scaling_factor_sg:        (float(par_value(where(par_name eq 'd02_upper_scaling_factor_sg'))))[0],        $
       d02_lower_scaling_factor_sg:        (float(par_value(where(par_name eq 'd02_lower_scaling_factor_sg'))))[0],        $
       smoothing_fwhm_factor_rg:           (float(par_value(where(par_name eq 'smoothing_fwhm_factor_rg'))))[0],           $
       ; Sliding-pattern fit
       asef_threshold_fraction:               (float(par_value(where(par_name eq 'asef_threshold_fraction'))))[0],               $
       n_central_orders_side:                 (float(par_value(where(par_name eq 'n_central_orders_side'))))[0],                 $
       depressed_dipole_fraction:             (float(par_value(where(par_name eq 'depressed_dipole_fraction'))))[0],             $
       dnu_ridge_threshold:                   (float(par_value(where(par_name eq 'dnu_ridge_threshold'))))[0],                   $
       dnu_echelle_threshold:                 (float(par_value(where(par_name eq 'dnu_echelle_threshold'))))[0],                 $
       epsilon_threshold:                     (float(par_value(where(par_name eq 'epsilon_threshold'))))[0],                     $
       dnu_lower_threshold_epsilon:           (float(par_value(where(par_name eq 'dnu_lower_threshold_epsilon'))))[0],           $
       n_sliding_test:                        (float(par_value(where(par_name eq 'n_sliding_test'))))[0],                        $
       input_radial_freq_reference:           (float(par_value(where(par_name eq 'input_radial_freq_reference'))))[0],           $
       force_epsilon_dnu_value:               (float(par_value(where(par_name eq 'force_epsilon_dnu_value'))))[0],               $
       n_orders_side_prior_ms:                (float(par_value(where(par_name eq 'n_orders_side_prior_ms'))))[0],                $
       n_orders_side_prior_sg:                (float(par_value(where(par_name eq 'n_orders_side_prior_sg'))))[0],                $
       n_orders_side_prior_rg:                (float(par_value(where(par_name eq 'n_orders_side_prior_rg'))))[0],                $
       dnu_prior_lower_fraction:              (float(par_value(where(par_name eq 'dnu_prior_lower_fraction'))))[0],              $
       dnu_prior_upper_fraction:              (float(par_value(where(par_name eq 'dnu_prior_upper_fraction'))))[0],              $
       d02_prior_lower_ms:                    (float(par_value(where(par_name eq 'd02_prior_lower_ms'))))[0],                    $
       d02_prior_upper_ms:                    (float(par_value(where(par_name eq 'd02_prior_upper_ms'))))[0],                    $
       d02_prior_lower_sg:                    (float(par_value(where(par_name eq 'd02_prior_lower_sg'))))[0],                    $
       d02_prior_upper_sg:                    (float(par_value(where(par_name eq 'd02_prior_upper_sg'))))[0],                    $
       d01_prior_lower_ms:                    (float(par_value(where(par_name eq 'd01_prior_lower_ms'))))[0],                    $
       d01_prior_upper_ms:                    (float(par_value(where(par_name eq 'd01_prior_upper_ms'))))[0],                    $
       d13_prior_lower_ms:                    (float(par_value(where(par_name eq 'd13_prior_lower_ms'))))[0],                    $
       rot_split_prior_lower_ms:              (float(par_value(where(par_name eq 'rot_split_prior_lower_ms'))))[0],              $
       rot_split_prior_upper_ms:              (float(par_value(where(par_name eq 'rot_split_prior_upper_ms'))))[0],              $
       cosi_prior_lower:                      (float(par_value(where(par_name eq 'cosi_prior_lower'))))[0],                      $
       cosi_prior_upper:                      (float(par_value(where(par_name eq 'cosi_prior_upper'))))[0],                      $
       quadrupole_radial_height_ratio_ms:     (float(par_value(where(par_name eq 'quadrupole_radial_height_ratio_ms'))))[0],     $
       quadrupole_radial_height_ratio_sg:     (float(par_value(where(par_name eq 'quadrupole_radial_height_ratio_sg'))))[0],     $
       quadrupole_radial_height_ratio_rg:     (float(par_value(where(par_name eq 'quadrupole_radial_height_ratio_rg'))))[0],     $
       dipole_radial_height_ratio_ms:         (float(par_value(where(par_name eq 'dipole_radial_height_ratio_ms'))))[0],         $
       dipole_radial_height_ratio_sg:         (float(par_value(where(par_name eq 'dipole_radial_height_ratio_sg'))))[0],         $
       dipole_radial_height_ratio_rg:         (float(par_value(where(par_name eq 'dipole_radial_height_ratio_rg'))))[0],         $
       octupole_radial_height_ratio:          (float(par_value(where(par_name eq 'octupole_radial_height_ratio'))))[0],          $
       dipole_radial_fwhm_ratio_ms:           (float(par_value(where(par_name eq 'dipole_radial_fwhm_ratio_ms'))))[0],           $
       dipole_radial_fwhm_ratio_sg:           (float(par_value(where(par_name eq 'dipole_radial_fwhm_ratio_sg'))))[0],           $
       dipole_radial_fwhm_ratio_rg:           (float(par_value(where(par_name eq 'dipole_radial_fwhm_ratio_rg'))))[0],           $
       upper_epsilon_rg_slope:                (float(par_value(where(par_name eq 'upper_epsilon_rg_slope'))))[0],                $
       upper_epsilon_rg_offset:               (float(par_value(where(par_name eq 'upper_epsilon_rg_offset'))))[0],               $
       lower_epsilon_rg_slope:                (float(par_value(where(par_name eq 'lower_epsilon_rg_slope'))))[0],                $
       lower_epsilon_rg_offset:               (float(par_value(where(par_name eq 'lower_epsilon_rg_offset'))))[0],               $
       ; Asymptotic fit
       dnu_prior_lower_fraction_as:           (float(par_value(where(par_name eq 'dnu_prior_lower_fraction_as'))))[0],           $
       dnu_prior_upper_fraction_as:           (float(par_value(where(par_name eq 'dnu_prior_upper_fraction_as'))))[0],           $
       alpha_prior_lower_as:                  (float(par_value(where(par_name eq 'alpha_prior_lower_as'))))[0],                  $
       alpha_prior_upper_as:                  (float(par_value(where(par_name eq 'alpha_prior_upper_as'))))[0],                  $
       epsi_prior_lower_fraction_as:          (float(par_value(where(par_name eq 'epsi_prior_lower_fraction_as'))))[0],          $
       epsi_prior_upper_fraction_as:          (float(par_value(where(par_name eq 'epsi_prior_upper_fraction_as'))))[0],          $
       ; Peak testing
       n_fwhm_fit:                             (fix(par_value(where(par_name eq 'n_fwhm_fit'))))[0],                                 $
       n_peak_test:                            (fix(par_value(where(par_name eq 'n_peak_test'))))[0],                                $
       detection_probability_threshold:        (float(par_value(where(par_name eq 'detection_probability_threshold'))))[0],          $
       rotation_probability_threshold:         (float(par_value(where(par_name eq 'rotation_probability_threshold'))))[0],           $
       duplicity_probability_threshold:        (float(par_value(where(par_name eq 'duplicity_probability_threshold'))))[0],          $
       rotation_test_activated:                (fix(par_value(where(par_name eq 'rotation_test_activated'))))[0],                    $
       height_ratio_threshold:                 (float(par_value(where(par_name eq 'height_ratio_threshold'))))[0],                   $
       bkg_prior_lower:                        (float(par_value(where(par_name eq 'bkg_prior_lower'))))[0],                          $
       bkg_prior_upper:                        (float(par_value(where(par_name eq 'bkg_prior_upper'))))[0],                          $
       rot_split_scaling:                      (float(par_value(where(par_name eq 'rot_split_scaling'))))[0],                        $
       fwhm_lower_bound:                       (float(par_value(where(par_name eq 'fwhm_lower_bound'))))[0],                         $
       fwhm_magnification_factor_radial:       (float(par_value(where(par_name eq 'fwhm_magnification_factor_radial'))))[0],         $
       fwhm_magnification_factor_dipole:       (float(par_value(where(par_name eq 'fwhm_magnification_factor_dipole'))))[0],         $
       fwhm_magnification_factor_quadrupole:   (float(par_value(where(par_name eq 'fwhm_magnification_factor_quadrupole'))))[0],     $
       fwhm_magnification_factor_octupole:     (float(par_value(where(par_name eq 'fwhm_magnification_factor_octupole'))))[0],       $
       fwhm_octupole_radial_fraction:          (float(par_value(where(par_name eq 'fwhm_octupole_radial_fraction'))))[0],            $
       dnu_mixed_modes_separation_scaling:     (float(par_value(where(par_name eq 'dnu_mixed_modes_separation_scaling'))))[0],       $
       ; Asteroseismic relations
       d03_offset:                 (float(par_value(where(par_name eq 'd03_offset'))))[0],             $
       d03_slope:                  (float(par_value(where(par_name eq 'd03_slope'))))[0],              $
       d01_mass_offset:            (float(par_value(where(par_name eq 'd01_mass_offset'))))[0],        $
       d01_mass_slope:             (float(par_value(where(par_name eq 'd01_mass_slope'))))[0],         $
       d01_offset:                 (float(par_value(where(par_name eq 'd01_offset'))))[0],             $
       d02_mass_offset:            (float(par_value(where(par_name eq 'd02_mass_offset'))))[0],        $
       d02_mass_slope:             (float(par_value(where(par_name eq 'd02_mass_slope'))))[0],         $
       d02_offset:                 (float(par_value(where(par_name eq 'd02_offset'))))[0],             $
       d02_unique_offset:          (float(par_value(where(par_name eq 'd02_unique_offset'))))[0],      $
       d02_unique_slope:           (float(par_value(where(par_name eq 'd02_unique_slope'))))[0],       $
       dnu_threshold:              (float(par_value(where(par_name eq 'dnu_threshold'))))[0],          $
       numax_threshold:            (float(par_value(where(par_name eq 'numax_threshold'))))[0],        $
       numax_coeff_high:           (float(par_value(where(par_name eq 'numax_coeff_high'))))[0],       $
       numax_exponent_high:        (float(par_value(where(par_name eq 'numax_exponent_high'))))[0],    $
       numax_coeff_low:            (float(par_value(where(par_name eq 'numax_coeff_low'))))[0],        $
       numax_exponent_low:         (float(par_value(where(par_name eq 'numax_exponent_low'))))[0],     $
       epsilon_offset:             (float(par_value(where(par_name eq 'epsilon_offset'))))[0],         $
       epsilon_slope:              (float(par_value(where(par_name eq 'epsilon_slope'))))[0]           $
     }


; -------------------------------------------------------------------------------------------------------
; Set up panel box positions for plotting output window, distinguishing between chunk and global modality.
; Plotting parameters based on whether the current graphic device is a X-window or EPS.
; -------------------------------------------------------------------------------------------------------
if info.save_eps eq 0 then begin
    ; Chunk modality
    chunk_pmulti = [0,3,1]
    chunk_position_sampling = [0.055,0.09,0.33,0.99]      ; Sampling
    chunk_position_asef = [0.4,0.09,0.98,0.45]             ; ASEF
    chunk_position_psd = [0.4,0.46,0.98,0.86]             ; PSD

    ; Global modality
    global_pmulti = [0,6,1]
    global_position_inset = [0.81,0.62,0.96,0.73]         ; PSD inset at nuMax
    global_position_epsi = [0.04,0.72,0.19,0.99]          ; epsilon-Teff
    global_position_sampling = [0.055,0.09,0.33,0.65]     ; Sampling
    global_position_asef = [0.4,0.09,0.98,0.45]            ; ASEF
    global_position_psd = [0.4,0.46,0.98,0.86]            ; PSD
    global_position_acf_dnu = [0.23,0.72,0.33,0.99]       ; Dnu ACF
    
    ; Echelle modality
    echelle_pmulti = [0,1,1]
    position_echelle = [0.08,0.09,0.88,0.86]              ; PSD échelle
    position_echelle_bar = [0.90,0.09,0.93,0.86]          ; PSD colorbar for échelle
endif else begin
    ; Chunk modality
    chunk_pmulti = [0,3,1]
    chunk_position_sampling = [0.055,0.06,0.33,0.99]      ; Sampling
    chunk_position_mixed = [0.04,0.68,0.33,0.99]          ; Mixed-modes pattern solution
    chunk_position_asef = [0.4,0.06,0.98,0.45]             ; ASEF
    chunk_position_psd = [0.4,0.46,0.98,0.86]             ; PSD
    
    ; Global modality
    global_pmulti = [0,6,1]
    global_position_inset = [0.81,0.62,0.96,0.73]         ; PSD inset at nuMax
    global_position_epsi = [0.04,0.68,0.19,0.99]          ; epsilon-Teff
    global_position_sampling = [0.055,0.06,0.33,0.62]     ; Sampling
    global_position_asef = [0.4,0.06,0.98,0.45]            ; ASEF
    global_position_psd = [0.4,0.46,0.98,0.86]            ; PSD
    global_position_acf_dnu = [0.23,0.68,0.33,0.99]       ; Dnu ACF
    
    ; Echelle modality
    echelle_pmulti = [0,1,1]
    position_echelle = [0.11,0.09,0.88,0.86]              ; PSD échelle
    position_echelle_bar = [0.90,0.09,0.93,0.86]          ; PSD colorbar for échelle
endelse

chunk_pp = {   position_sampling:      chunk_position_sampling,        $
               position_asef:          chunk_position_asef,            $
               position_psd:           chunk_position_psd,             $
               pmulti:                 chunk_pmulti                    $
           }

global_pp = {  position_inset:       global_position_inset,       $
               position_acf_dnu:     global_position_acf_dnu,     $
               position_epsi:        global_position_epsi,        $
               position_sampling:    global_position_sampling,    $
               position_asef:        global_position_asef,        $
               position_psd:         global_position_psd,         $
               pmulti:               global_pmulti                $
            }

dimensions = get_screen_size(resolution=res)

if info.save_eps eq 0 then begin
    color_line =         0   
    color_envelope =     255 
    color_hist =         160 
    color_interval =     255 
    thick =              1  
    xthick =             1   
    ythick =             1   
    charthick =          1.5 
    xcharsize =          1.0 
    ycharsize =          1.0 
    charsize =           3.0 
    label_charsize =     1.5 
    label_charthick =    2.0 
    label_color =        200 
    psd_smth_thick =     2.0 
    psd_smth_color =     100
    bg_level_color =     160 
    psd_smth_inset_color =   110 
    background_psd =     50
    background_inset =   0  
    psd_color_global =   140 
    psd_color_inset =    120 
    psd_color_chunk =    140 
    psd_thick =          2
    tickscolor_global =  255
    tickscolor_chunk =   255
    symsize =            1.3 
    symsize2 =           1.0 
    symsize3 =           1.05
    symsize_mixed =      1.3
    symthickness =       1  
    freq_symcolor =      240
    freq_symborder =     150
    freq_symthick =      2
    cross_color =        120 
    sample_color =       230 
    color_box =          0  
    fit_color =          205
    brewer_colorbar_asef =      25
    thick_border_asef =         2
    chunk_band_color1 =         255
    chunk_band_color2 =         240
    numax_arrow_hsize =         10
    numax_arrow_color =         200
    numax_arrow_thick =         3
    numax_chunk_arrow_color =   90
    freq_asymp_arrow_hsize =    10
    freq_asymp_arrow_color =    207
    freq_asymp_arrow_thick =    3
   
    ; Graphic window size and position relative to screen
    ;xsize =     1700
    ;ysize =     720
    xsize =     dimensions(0)
    ysize =     dimensions(1)*0.7
    xpos =      0
    ypos =      60
endif else begin
    color_line =         0  
    color_envelope =     0 
    color_hist =         200
    color_interval =     0  
    thick =              3  
    xthick =             3  
    ythick =             3  
    charthick =          3  
    xcharsize =          1.3
    ycharsize =          1.3
    charsize =           1.1
    label_charsize =     0.6
    label_charthick =    2.5
    label_color =        60 
    psd_smth_thick =     3 
    psd_smth_color =     100
    bg_level_color =     160 
    psd_smth_inset_color =   208
    background_psd =     50
    background_inset =   140
    psd_color_global =   160
    psd_color_inset =    100
    psd_color_chunk =    30 
    psd_thick =          2
    tickscolor_global =  255
    tickscolor_chunk =   0
    symsize =            0.9
    symsize2 =           0.4
    symsize3 =           0.45
    symsize_mixed =      0.6
    symthickness =       3 
    freq_symcolor =      240
    freq_symborder =     150
    freq_symthick =      2
    cross_color =        80 
    sample_color =       220
    color_box =          255
    fit_color =          40  
    brewer_colorbar_asef =      25
    thick_border_asef =         2
    chunk_band_color1 =         255
    chunk_band_color2 =         240
    numax_arrow_hsize =         180
    numax_arrow_color =         200
    numax_arrow_thick =         5
    numax_chunk_arrow_color =   60
    freq_asymp_arrow_hsize =    180
    freq_asymp_arrow_color =    207
    freq_asymp_arrow_thick =    5
   
    ; Graphic window size and position relative to screen
    xsize =     26
    ysize =     14
    xpos =      0
    ypos =      60
endelse
 
acf_pp = {  thick:         thick,           $
            xthick:        xthick,          $
            ythick:        ythick,          $
            charthick:     charthick,       $
            charsize:      charsize*0.8,    $
            xcharsize:     xcharsize,       $
            ycharsize:     ycharsize,       $
            symsize:       symsize,         $
            symthickness:  symthickness,    $
            cross_color:   cross_color,     $
            sample_color:  sample_color,    $
            fit_color:     fit_color        $
         }

epsi_pp = {  thick:          thick,               $
             xthick:         xthick,              $
             ythick:         ythick,              $
             charthick:      charthick,           $
             charsize:       charsize*0.8,        $
             xcharsize:      xcharsize,           $
             ycharsize:      ycharsize,           $
             symsize:        symsize,             $
             symthickness:   symthickness,        $
             cross_color:    cross_color,         $
             sample_color:   sample_color,        $
             color_box:      color_box,           $
             label_charsize: label_charsize,      $
             fit_color:      fit_color            $
          }

mixed_pp = {  thick:          thick,               $
               xthick:         xthick,              $
               ythick:         ythick,              $
               charthick:      charthick,           $
               charsize:       charsize*0.8,        $
               xcharsize:      xcharsize,           $
               ycharsize:      ycharsize,           $
               symsize:        symsize_mixed,       $
               symthickness:   symthickness,        $
               color_box:      color_box,           $
               cross_color:    cross_color,         $
               label_charsize: label_charsize,      $
               fit_color:      fit_color            $
            }

pp = { color_line:             color_line,             $
       color_envelope:         color_envelope,         $
       color_hist:             color_hist,             $
       color_interval:         color_interval,         $
       thick:                  thick,                  $
       xthick:                 xthick,                 $
       ythick:                 ythick,                 $
       charthick:              charthick,              $
       xcharsize:              xcharsize,              $
       ycharsize:              ycharsize,              $
       charsize:               charsize,               $
       label_charsize:         label_charsize,         $
       label_charthick:        label_charthick,        $
       label_color:            label_color,            $
       psd_smth_thick:         psd_smth_thick,         $
       psd_smth_inset_color:   psd_smth_inset_color,   $
       psd_smth_color:         psd_smth_color,         $
       bg_level_color:         bg_level_color,         $
       background_psd:         background_psd,         $
       background_inset:       background_inset,       $
       psd_color_global:       psd_color_global,       $
       psd_color_inset:        psd_color_inset,        $
       psd_color_chunk:        psd_color_chunk,        $
       psd_thick:              psd_thick,              $
       tickscolor_global:      tickscolor_global,      $
       tickscolor_chunk:       tickscolor_chunk,       $
       symsize:                symsize,                $
       symsize2:               symsize2,               $
       symsize3:               symsize3,               $
       symsize_mixed:          symsize_mixed,          $
       symthickness:           symthickness,           $
       freq_symcolor:          freq_symcolor,          $
       freq_symborder:         freq_symborder,         $
       freq_symthick:          freq_symthick,          $
       cross_color:            cross_color,            $
       sample_color:           sample_color,           $
       color_box:              color_box,              $
       fit_color:              fit_color,              $
       brewer_colorbar_asef:   brewer_colorbar_asef,   $
       thick_border_asef:      thick_border_asef,      $
       chunk_band_color1:      chunk_band_color1,      $
       chunk_band_color2:      chunk_band_color2,      $
       numax_arrow_hsize:      numax_arrow_hsize,      $
       numax_arrow_color:      numax_arrow_color,      $
       numax_arrow_thick:      numax_arrow_thick,      $
       numax_chunk_arrow_color:    numax_chunk_arrow_color,    $
       freq_asymp_arrow_hsize:     freq_asymp_arrow_hsize,     $
       freq_asymp_arrow_color:     freq_asymp_arrow_color,     $
       freq_asymp_arrow_thick:     freq_asymp_arrow_thick,     $
       xsize:                  xsize,                  $
       ysize:                  ysize,                  $
       xpos:                   xpos,                   $
       ypos:                   ypos,                   $
       chunk:                  chunk_pp,               $
       global:                 global_pp,              $
       acf:                    acf_pp,                 $
       epsi:                   epsi_pp,                $
       mixed:                  mixed_pp                $
     }

if info.save_eps eq 0 then begin
    xcharsize =                 0.8
    ycharsize =                 0.8
    charsize =                  2.0
    label_charsize =            1.4
    label_charthick =           2.0
    label_color =               200
    symsize =                   1.8
    symsize2 =                  1.0
    symsize3 =                  1.05
    symthickness =              1
    color_box =                 0
    label_charsize =            charsize*0.5
    fit_color =                 205
    xticklen =                  0.04
    yticklen =                  0.03
    brewer_colorbar_echelle =   30    ;5,17,13
   
    ; Graphic window size and position relative to screen
    ;xsize = 1100
    ;ysize = 720
    xsize =                     dimensions(0)*0.55
    ysize =                     dimensions(1)*0.7
    xpos =                      0
    ypos =                      60
endif else begin
    xcharsize =                 1.0
    ycharsize =                 1.0
    charsize =                  1.0
    label_charsize =            1.0
    label_charthick =           2.0
    label_color =               60
    symsize =                   0.8
    symsize2 =                  0.4
    symsize3 =                  0.45
    symthickness =              2
    color_box =                 255
    label_charsize =            charsize*0.6
    fit_color =                 40
    xticklen =                  0.04
    yticklen =                  0.03
    brewer_colorbar_echelle =   30
   
    ; Graphic window size and position relative to screen
    xsize =                     18
    ysize =                     14
    xpos =                      0
    ypos =                      60
endelse

ppe = { position_echelle:          position_echelle,              $
        position_echelle_bar:      position_echelle_bar,          $
        echelle_pmulti:            echelle_pmulti,                $
        thick:                     thick,                         $
        xthick:                    xthick,                        $
        ythick:                    ythick,                        $
        charthick:                 charthick,                     $
        xcharsize:                 xcharsize,                     $
        ycharsize:                 ycharsize,                     $
        charsize:                  charsize,                      $
        label_charsize:            label_charsize,                $
        label_charthick:           label_charthick,               $
        label_color:               label_color,                   $
        symsize:                   symsize,                       $
        symsize2:                  symsize2,                      $
        symsize3:                  symsize3,                      $
        symthickness:              symthickness,                  $
        color_box:                 color_box,                     $
        xticklen:                  xticklen,                      $
        yticklen:                  yticklen,                      $
        brewer_colorbar_echelle:   brewer_colorbar_echelle,       $
        xsize:                     xsize,                         $
        ysize:                     ysize,                         $
        xpos:                      xpos,                          $
        ypos:                      ypos                           $
      }

; -------------------------------------------------------------------------------------------------------
; Setup plotting parameters for output label summaries based on whether the current graphic device is 
; a X-window or EPS. This is used for GLOBAL and CHUNK modalities.
; -------------------------------------------------------------------------------------------------------
if info.save_eps eq 0 then begin
    block1_1 = [0.34,0.96]
    block1_2 = [0.34,0.935]
    block1_3 = [0.34,0.91]
    block1_4 = [0.34,0.885]
    block2_1 = [0.425,0.96]
    block2_2 = [0.425,0.935]
    block2_3 = [0.425,0.91]
    block2_4 = [0.425,0.885]
    block2_5 = [0.475,0.91]
    block2_6 = [0.495,0.885]
    block2_7 = [0.53,0.935]
    block3_1 = [0.59,0.96]
    block3_2 = [0.59,0.935]
    block3_3 = [0.59,0.91]
    block3_4 = [0.59,0.885]
    block4_1 = [0.74,0.96]
    block4_2 = [0.74,0.935]
    block4_3 = [0.74,0.91]
    block4_4 = [0.74,0.885]
    block4_5 = [0.80,0.885]
    block1_color = 200 
    block2_color = 160
    block3_color = 100
    summary_charthick_title = 2.0
    summary_charthick = 1.5
    summary_charsize_title = 1.5
    summary_charsize = 1.5
endif else begin
    block1_1 = [0.34,0.96]
    block1_2 = [0.34,0.935]
    block1_3 = [0.34,0.91]
    block1_4 = [0.34,0.885]
    block2_1 = [0.455,0.96]
    block2_2 = [0.455,0.935]
    block2_3 = [0.455,0.91]
    block2_4 = [0.455,0.885]
    block2_5 = [0.515,0.91]
    block2_6 = [0.545,0.885]
    block2_7 = [0.58,0.935]
    block3_1 = [0.66,0.96]
    block3_2 = [0.66,0.935]
    block3_3 = [0.66,0.91]
    block3_4 = [0.66,0.885]
    block4_1 = [0.85,0.96]
    block4_2 = [0.85,0.935]
    block4_3 = [0.85,0.91]
    block4_4 = [0.85,0.885]
    block4_5 = [0.92,0.885]
    block1_color = 250 
    block2_color = 80
    block3_color = 40
    summary_charthick_title = 4.0
    summary_charthick = 3.0
    summary_charsize_title = 0.8
    summary_charsize = 0.8
endelse

lp = {  block1_1:                  block1_1,                   $
        block1_2:                  block1_2,                   $
        block1_3:                  block1_3,                   $
        block1_4:                  block1_4,                   $
        block2_1:                  block2_1,                   $
        block2_2:                  block2_2,                   $
        block2_3:                  block2_3,                   $
        block2_4:                  block2_4,                   $
        block2_5:                  block2_5,                   $
        block2_6:                  block2_6,                   $
        block2_7:                  block2_7,                   $
        block3_1:                  block3_1,                   $
        block3_2:                  block3_2,                   $
        block3_3:                  block3_3,                   $
        block3_4:                  block3_4,                   $
        block4_1:                  block4_1,                   $
        block4_2:                  block4_2,                   $
        block4_3:                  block4_3,                   $
        block4_4:                  block4_4,                   $
        block4_5:                  block4_5,                   $
        block1_color:              block1_color,               $
        block2_color:              block2_color,               $
        block3_color:              block3_color,               $
        summary_charthick_title:   summary_charthick_title,    $
        summary_charthick:         summary_charthick,          $
        summary_charsize_title:    summary_charsize_title,     $
        summary_charsize:          summary_charsize            $
     }
 

; -------------------------------------------------------------------------------------------------------
; Setup plotting parameters for output label summaries based on whether the current graphic device is 
; a X-window or EPS. This is used for the ECHELLE modality.
; -------------------------------------------------------------------------------------------------------
if info.save_eps eq 0 then begin
    block1_1 = [0.02,0.96]
    block1_2 = [0.02,0.935]
    block1_3 = [0.02,0.91]
    block1_4 = [0.02,0.885]
    block2_1 = [0.19,0.96]
    block2_2 = [0.19,0.935]
    block2_3 = [0.19,0.91]
    block2_4 = [0.19,0.885]
    block2_5 = [0.28,0.885]
    block3_1 = [0.36,0.96]
    block3_2 = [0.36,0.935]
    block3_3 = [0.36,0.91]
    block3_4 = [0.36,0.885]
    block4_1 = [0.49,0.96]
    block4_2 = [0.49,0.935]
    block4_3 = [0.49,0.91]
    block4_4 = [0.49,0.885]
    block5_1 = [0.62,0.96]
    block5_2 = [0.62,0.935]
    block5_3 = [0.62,0.91]
    block5_4 = [0.62,0.885]
    block1_color = 200 
    block2_color = 160
    block3_color = 100
    block4_color = 210
    block5_color = 255
    summary_charthick_title = 2.0
    summary_charthick = 1.5
    summary_charsize_title = 1.5
    summary_charsize = 1.5
endif else begin
    block1_1 = [0.02,0.96]
    block1_2 = [0.02,0.935]
    block1_3 = [0.02,0.91]
    block1_4 = [0.02,0.885]
    block2_1 = [0.23,0.96]
    block2_2 = [0.23,0.935]
    block2_3 = [0.23,0.91]
    block2_4 = [0.23,0.885]
    block2_5 = [0.34,0.885]
    block3_1 = [0.44,0.96]
    block3_2 = [0.44,0.935]
    block3_3 = [0.44,0.91]
    block3_4 = [0.44,0.885]
    block4_1 = [0.61,0.96]
    block4_2 = [0.61,0.935]
    block4_3 = [0.61,0.91]
    block4_4 = [0.61,0.885]
    block5_1 = [0.78,0.96]
    block5_2 = [0.78,0.935]
    block5_3 = [0.78,0.91]
    block5_4 = [0.78,0.885]
    block1_color = 250 
    block2_color = 80
    block3_color = 40
    block4_color = 210
    block5_color = 0
    summary_charthick_title = 4.0
    summary_charthick = 3.0
    summary_charsize_title = 0.8
    summary_charsize = 0.8
endelse

lpe = { block1_1:                  block1_1,                   $
        block1_2:                  block1_2,                   $
        block1_3:                  block1_3,                   $
        block1_4:                  block1_4,                   $
        block2_1:                  block2_1,                   $
        block2_2:                  block2_2,                   $
        block2_3:                  block2_3,                   $
        block2_4:                  block2_4,                   $
        block2_5:                  block2_5,                   $
        block3_1:                  block3_1,                   $
        block3_2:                  block3_2,                   $
        block3_3:                  block3_3,                   $
        block3_4:                  block3_4,                   $
        block4_1:                  block4_1,                   $
        block4_2:                  block4_2,                   $
        block4_3:                  block4_3,                   $
        block4_4:                  block4_4,                   $
        block5_1:                  block5_1,                   $
        block5_2:                  block5_2,                   $
        block5_3:                  block5_3,                   $
        block5_4:                  block5_4,                   $
        block1_color:              block1_color,               $
        block2_color:              block2_color,               $
        block3_color:              block3_color,               $
        block4_color:              block4_color,               $
        block5_color:              block5_color,               $
        summary_charthick_title:   summary_charthick_title,    $
        summary_charthick:         summary_charthick,          $
        summary_charsize_title:    summary_charsize_title,     $
        summary_charsize:          summary_charsize            $
     }


; -------------------------------------------------------------------------------------------------------
; Setup configuring parameters for DIAMONDS fits
; -------------------------------------------------------------------------------------------------------
dp_isla = { n_live:                           2000,   $
            n_live_end:                       2000,   $
            max_draw_attempts:                50000,  $
            n_initial_it:                     1000,   $
            n_it_same_clust:                  50,     $
            enlarg_fraction:                  2.00,   $
            shrinking_rate:                   0.0,    $
            termination_factor:               1.0,    $
            max_nested_it:                    8000,   $
            min_ncluster:                     1,      $
            max_ncluster:                     10,     $
            max_ncluster_duplet:              3       $
          }

dp_slid = { n_live:                           2000,   $
          n_live_end:                         2000,   $
          max_draw_attempts:                  50000,  $
          n_initial_it:                       1000,   $
          n_it_same_clust:                    50,     $
          enlarg_fraction:                    2.00,   $
          shrinking_rate:                     0.0,    $
          termination_factor:                 1.0,    $
          max_nested_it:                      4000,   $
          min_ncluster:                       1,      $
          max_ncluster:                       3       $
        }

dp_pb = { n_live:                             500,    $
          n_live_end:                         500,    $
          n_live_test:                        1000,   $
          n_live_end_test:                    1000,   $
          max_draw_attempts:                  50000,  $
          n_initial_it:                       1000,   $
          n_it_same_clust:                    50,     $
          enlarg_fraction:                    2.00,   $
          shrinking_rate:                     0.0,    $
          termination_factor:                 1.0,    $
          termination_factor_test:            0.1,    $
          max_nested_it:                      0,      $
          min_ncluster:                       1,      $
          max_ncluster:                       4       $
        }

dp = { isla:   dp_isla,    $
       slid:   dp_slid,    $
       pb:     dp_pb       $
     }


; -------------------------------------------------------------------------------------------------------
; Setup strings to be used for output plotting axis labeling
; -------------------------------------------------------------------------------------------------------
alpha_letter = "141B       ;"
beta_letter = "142B        ;"
gamma_letter = "103B       ;"
delta_letter = "104B       ;"
deltasmall_letter = "144B      ;" 
epsi_letter = "145B    ;"
mu_letter = "154B     ;" 
nu_letter = "155B      ;"
pi_letter = "120B     ;"
tau_letter = "163B     ;"
copyright_letter = "251B       ;"
alpha_str = '!4'+string(alpha_letter)+'!X'
beta_str = '!4'+string(beta_letter)+'!X'
gamma_str = '!4'+string(gamma_letter)+'!X'
epsi_str = '!4'+string(epsi_letter)+'!X'
mu_str = '!4'+string(mu_letter)+'!X'
nu_str = '!4'+string(nu_letter)+'!X'
tau_str = '!4'+string(tau_letter)+'!X'
copyright_str = '!3' + String(copyright_letter) + '!X'
dnu_str = '!4'+string(delta_letter)+string(nu_letter)+'!X'
dp1_str = '!4'+string(delta_letter)+string(pi_letter)+'!X!D1!n'
dp_str = '!4'+string(delta_letter)+'!XP!D1!n'
d01_str = '!4'+string(deltasmall_letter)+string(nu_letter)+'!X!D01!n'
d02_str = '!4'+string(deltasmall_letter)+string(nu_letter)+'!X!D02!n'
d03_str = '!4'+string(deltasmall_letter)+string(nu_letter)+'!X!D03!n'
freq_unit_str = mu_str+'!3Hz'
psd_unit_str = 'ppm!U2!n/'+mu_str+'!3Hz'

sp = {  alpha_str:             alpha_str,          $
        beta_str:              beta_str,           $
        gamma_str:             gamma_str,          $
        epsi_str:              epsi_str,           $
        mu_str:                mu_str,             $
        nu_str:                nu_str,             $
        tau_str:               tau_str,            $
        copyright_str:         copyright_str,      $
        dnu_str:               dnu_str,            $
        d01_str:               d01_str,            $
        d02_str:               d02_str,            $
        d03_str:               d03_str,            $
        dp1_str:               dp1_str,            $
        dp_str:                dp_str,             $
        freq_unit_str:         freq_unit_str,      $
        psd_unit_str:          psd_unit_str        $
     }
end
