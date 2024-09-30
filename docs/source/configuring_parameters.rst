.. _configuring_parameters:

Configuring Parameters
======================
The FAMED pipeline requires an input file containing a list of 170 configuring parameters. This file is labeled ``famed_configuring_parameters.txt`` and by default it is located under the ``FAMED/idl/`` folder. The parameters initialized inside this file can all be changed, depending on the needs, allowing the code to be highly configurable and flexible to be adapted to different conditions. 

In this section we provide a summary description of each of the configuring parameters currently adopted by the pipeline. The default values are already provided with the package content. These values have been obtained from calibration sessions over hundreds of stars, both from mock and real samples, meaning that they are expected to be relatively accurate. We therefore discourage the user from modifying the configuring parameters that are the result of a calibration phase, unless particular circumstances may occur, such as the analysis of a challenging target. Any change, if done, should always go into the direction of tuning values with respect to the default ones in order to avoid abrupt alteration of the pipeline functionalities. 

The configuring parameters are divided into subgroups, each one related to a class of tasks of the pipeline.

Global paths and filenames
^^^^^^^^^^^^^^^^^^^^^^^^^^
* ``root_path`` 
YOUR_LOCAL_ROOT_PATH_HERE/

This root path will be used to access and execute the DIAMONDS-based codes, as well as to set the paths for storing the output files. If you are not using the shell script installation, you will need to replace the name YOUR_LOCAL_ROOT_PATH_HERE with the actual path containing the DIAMONDS, Background, PeakBagging, and Asymptotic folders (see the Installation section for more details). Make sure to put the final slash `/`. 

* ``famed_path`` 
YOUR_LOCAL_ROOT_PATH_HERE/FAMED/idl

The path where you have stored the FAMED pipeline. If you are not using the shell script installation, you will need to replace the name YOUR_LOCAL_ROOT_PATH_HERE with the actual path containing the FAMED package. This path has to refer to the specific folder that contains either the IDL or Python code (in the example shown, we are using the idl folder). This path is used by FAMED to actually run the pipeline. This path is by default within the ``root_path`` parameter but the user may decide to use a different location.

* ``isla_subdir``
The name of the folder that contains all the output products created by the PeakBagging code when running in the multi-modal modalities (see Corsaro et al. 2020). These outputs include configuring parameters to run DIAMONDS, prior hyper-parameter files and fit products of parameter sampling.

* ``pb_subdir``
The name of the folder that contains all the output products created by the PeakBagging code when running in the uni-modal modalities, which therefore comprise all the peak testing suite (see Corsaro et al. 2020). These outputs include configuring parameters to run DIAMONDS, prior hyper-parameter files and fit products of parameter sampling, parameter estimates, and Bayesian evidences for model comparison.

* ``as_subdir``
The name of the folder that contains all the output products created by the Asymptotic code (see Corsaro et al. 2020). These outputs include configuring parameters to run DIAMONDS, prior hyper-parameter files and fit products of parameter sampling.

* ``figs_subdir``
The name of the folder that countains all the output figures of the computation, in either PNG or EPS format. There is an overall figure for the GLOBAL modality, and as many figures as the number of chunks identified for the CHUNK modality. Each overall figure comprises a plot of the nested sampling, the ASEF, the PSD with mode identification overlaid, and in the case of the GLOBAL modality also the *ACF*:math:`^2` for :math:`\Delta\nu`, and the :math:`\epsilon` - *T*:subscript:`eff` or :math:`\epsilon` - :math:`\Delta\nu` diagram (see Corsaro et al. 2020).

* ``summary_subdir``
The name of the folder that contains all the output ASCII files of the computation. These files store all the useful information that is obtained by FAMED, in particular the oscillation frequencies and their mode identification :math:`(n,\ell,m)`, the peak significance probabilities (for detection, rotation, and duplicity), flags for peak blending and peak *sinc*:math:`^2` profiles, other information in relation to the analysis of the ASEF, such as sampling counts, ASEF maxima, smoothed PSD maxima, as well as asymptotic parameters of :math:`\Delta\nu`, :math:`\epsilon`, :math:`\delta\nu_{02}`, and where available also :math:`\Delta P`.

* ``global_subdir``
This is the name of the sub-folder that contains the results for the GLOBAL modality multi-modal sampling, as well as a filename suffix for priors and results related to the GLOBAL modality.

* ``save_complete_lists``
Allows the user storing as output the complete lists of frequencies identified from the ASEF. These additional lists are generated for both GLOBAL and CHUNK modalities. These lists can be useful in case of post-production inspection, or as a reference to quantify the skimming process applied by the pipeline on the initial frequency set. The output files are in ASCII format and are stored in the summary_subdir folder. The default value is 0, meaning that this option is deactivated. Set it to 1 to activate.

* ``save_png``
Allows the user storing as output the plots generated by FAMED with a PNG image format. The default value is 1, meaning that this option is activated. To deactivate it, set it to 0.

* ``save_eps``
Allows the user storing as output the plots generated by FAMED with an EPS format. The default value is 0, meaning that this option is deactivated. To activate it, set it to 1. Note that if both ``save_png`` and ``save_eps`` are set active, the EPS format will be the one adopted by the pipeline. 

* ``plot_total_solution``
Generates a plot of the PSD of the star with overlaid the frequencies and mode identification extracted during the analysis of the CHUNK module. The plot of the PSD comprises all the frequency range explored during the peak bagging and can be considered as a joint plot of the individual chunk plots generated by the CHUNK module. It can be useful for diagnostic purposes, e.g. to inspect how the detailed mode identification was carried out and to check whether the extracted frequencies are reliable. The default value is set to 1, meaning that the plot is generated. 

* ``print_on_screen``
Set this keyword to print information from the analysis throughout the computation directly on the terminal window. Activating this option will also allow to produce the output plots. The default value is 1, meaning that this option is activated. To deactivate it, set it to 0.

* ``prior_filename``
This is the prefix name of the filename that contains the prior hyper-parameters used by DIAMONDS-related codes.

* ``background_run_number``
The sub-folder label contaning the output results from the Background fit. This option is adopted only if a background fit using the Background code extension of DIAMONDS has been used. This run number refers to the file system format adopted by the Background code. See also the tutorial provided on the `Background code GitHub repository <https://github.com/EnricoCorsaro/Background>`_ for more details.

* ``n_threads``
An integer specifying the total number of CPU threads to parallelize the computation of the multi-modal chunk peak bagging, of the sliding pattern fit, and of the peak test and peak FWHM fits. This number should be set to the maximum allowed by the user depending on the available resources of the system. It has to be a number > 1. Each thread will automatically be used to spawn an individual computation with DIAMONDS where more independent ones are to be performed at the same time.

* ``logo_filename`` 
YOUR_LOCAL_ROOT_PATH_HERE/FAMED/docs/figures/FAMED_LOGO_BLACK_small.jpeg

The filename, including the full path, of the FAMED logo. By default when downloading the package it is placed inside FAMED/docs/figures, but the user may decide to place the file in any other place. Change the name YOUR_LOCAL_ROOT_PATH_HERE according to your system if you are not using the shell script installation.

* ``background_data_dir``
The full path specifying the folder containing the stellar Power Spectral Density used for the background fit. The default value is set to -99, meaning that the data are stored inside the ``Background/data/`` folder of the DIAMONDS+Background code. As an example you can use this parameter in case you have your data stored in a separate hard drive.

* ``background_results_dir``
The full path specifying the folder containing the results of the background fit obtained with the DIAMONDS+Background code. The default value is set to -99, meaning that the results are stored inside the ``Background/results/`` folder of the DIAMONDS+Background code. As an example you can use this parameter in case you have background fit results stored in a separate hard drive.

* ``external_background_results_dir``
The full path of the folder containing the file of the parameters that configure the background model. This folder is used only if there is no result available from a background fit using the Background code extension of DIAMONDS. The default value is -99, but it has to be changed if the user wants to supply a result from a fitting code different than DIAMONDS-Background. When specifying a directory path, make sure that this path ends with the slash symbol `/`.

* ``external_background_filename_suffix``
The suffix that is attached to Catalog_ID + Star_ID of the star being analyzed, for the file that contains the background parameters. This file is located inside ``external_background_results_dir``. The file has a 1-column ASCII format, where the first line contains the background model name, among those available in the PeakBagging code extension of DIAMONDS, and the remaining lines contain each one a parameter estimate for the background model that is adopted, with the same order as those defined in the Background model (see `this page <https://famed.readthedocs.io/en/latest/background_models.html>`_ for a list of the available background model names and of their free parameters). This keyword is adopted only if the previous one is different than its default value.

Solar reference parameters for asteroseismic scaling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* ``teff_sun`` 
The solar :math:`T_\mathrm{eff}`, set to 5777 K

* ``dnu_sun``
The solar :math:`\Delta\nu`, set to 134.9 :math:`\mu\mbox{Hz}`

* ``numax_sun``
The solar :math:`\nu_{max}`, set to 3150 :math:`\mu\mbox{Hz}`

Multi-modal sampling production
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* ``n_sigma_envelope``
The number (float) of standard deviations from the background fit of the Gaussian envelope, to set the total frequency range of power spectrum of the star on each side of :math:`\nu_{max}`. The default value is 4.5 to accomodate MS stars with a very broad oscillation spectrum and 15-20 different radial orders observed.

* ``n_sigma_envelope_cl``
Similar keyword as ``n_sigma_envelope`` but being used for stars having :math:`\Delta\nu_\mathrm{AGB} < \Delta\nu \leq \Delta\nu_\mathrm{CL2}`. The default value is 2.5, significantly smaller than that of less evolved stars.

* ``n_sigma_envelope_tip``
Similar keyword as ``n_sigma_envelope`` but being used for stars having :math:`\Delta\nu \leq \Delta\nu_\mathrm{AGB}`. The default value is 1.2 because of the small extent of the Gaussian envelope in very evolved stars.

* ``n_dnu_envelope``
The maximum number (float) of times that :math:`\Delta\nu` is repeated to set the frequency range of the power spectrum on each side of :math:`\nu_{max}`. This parameter is superimposing on the previous ones if the frequency range obtained from the sigma envelope is larger than that set by this parameter. It can be useful for stars with low SNR data, where sigma envelope could result in being significantly larger than the frequency range that is actually covered by observable oscillations in the power spectrum. The default value is 5.0.

* ``dnu_tip``
A threshold value of :math:`\Delta\nu`, set to 1.55 :math:`\mu\mbox{Hz}` according to the findings by Kallinger et al. (2012), which identifies a regime of stars evolved toward the RGB tip (:math:`\Delta\nu \leq \Delta\nu_\mathrm{tip}`).

* ``dnu_agb``
A threshold value of :math:`\Delta\nu`, set to 3.15 :math:`\mu\mbox{Hz}` according to the findings by Kallinger et al. (2012), which identifies stars that are potentially AGB or evolved toward the RGB tip (:math:`\Delta\nu_\mathrm{tip} < \Delta\nu \leq \Delta\nu_\mathrm{AGB}`), allowing to separate them from stars that are in the red clump (:math:`\Delta\nu > \Delta\nu_\mathrm{AGB}`).

* ``dnu_cl``
A threshold value of :math:`\Delta\nu`, set to 5.0 :math:`\mu\mbox{Hz}` according to the findings by Kallinger et al. (2012), which separates the primary and secondary red clump phases.

* ``dnu_cl2``
A threshold value of :math:`\Delta\nu`, set to 9.0 :math:`\mu\mbox{Hz}` according to the findings by Kallinger et al. (2012), which separates stars that are either in the secondary clump or RGB from stars belonging to the early RGB phase.

* ``dnu_rg``
A threshold value of :math:`\Delta\nu`, set to 15.0 :math:`\mu\mbox{Hz}` (see Corsaro et al. 2020), which identifies stars evolving along the low-luminosity RGB (:math:`\Delta\nu_\mathrm{CL2} < \Delta\nu \leq \Delta\nu_\mathrm{RG}`).

* ``dnu_sg``
A threshold value of :math:`\Delta\nu`, set to 90.0 :math:`\mu\mbox{Hz}` according to Appourchaux et al. (2012), which identifies stars evolving along the subgiant branch (:math:`\Delta\nu_\mathrm{RG} < \Delta\nu < \Delta\nu_\mathrm{SG}` and :math:`T_\mathrm{eff} < T_\mathrm{eff,SG}`, see the parameter ``teff_sg``). Stars in the MS phase of stellar evolution could also appear in this regime, in which case a cross-check about the presence of avoided crossing is necessary to be able to classify the star as a subgiant.

* ``teff_sg``
A threshold value of :math:`T_\mathrm{eff}`, set to 6350 K (see Corsaro et al. 2020), to distinguish hot F-type stars from G-type stars based on the impact of the stellar effective temperature on the width of the radial oscillation modes. This parameter is used in conjunction with ``dnu_sg`` to distinguish subgiant stars from hot MS stars when :math:`\Delta\nu > \Delta\nu_\mathrm{RG}`.

* ``fwhm_global_scaling``
The scaling value applied on the FWHM prediction to set the linewidth of the Lorentzian profile of the islands peak bagging model adopted in the GLOBAL modality. The default value is set to 1, meaning that no scaling is taking place.

* ``fwhm_global_scaling_tip``
Similar keyword as ``fwhm_global_scaling`` but for the specific case of stars evolved toward the tip of the RGB, or into early AGB (:math:`\Delta\nu \leq \Delta\nu_\mathrm{AGB}`). The default value is set to 10 to accomodate the very narrow oscillation mode linewidths found in these cool stars.

* ``fwhm_chunk_scaling_ms``
Similar keyword as ``fwhm_global_scaling`` but for the CHUNK modality and the case of MS stars, having either :math:`\Delta\nu_\mathrm{RG} < \Delta\nu \leq \Delta\nu_\mathrm{SG}` and :math:`T_\mathrm{eff} \geq T_\mathrm{eff,SG}` or :math:`\Delta\nu \geq \Delta\nu_\mathrm{SG}`. The default value is set to 10 for improving the resolution on the multi-modal sampling due to the very large oscillation mode linewidths found in these stars.

* ``fwhm_chunk_scaling_early_sg``
Similar keyword as ``fwhm_chunk_scaling_ms`` but for the CHUNK modality and the case of early subgiant stars, having :math:`\Delta\nu_\mathrm{RG} < \Delta\nu < \Delta\nu_\mathrm{SG}` and :math:`T_\mathrm{eff}` towards the upper limit :math:`T_\mathrm{eff,SG}`. This value is set to 1.0, implying an almost negligible scaling to better accomodate the large mode linewidths founds in hotter subgiants, but it increases linearly with :math:`T_\mathrm{eff} as the temperature decreases, up to the liming value imposed by ``fwhm_chunk_scaling_late_sg`` (see below).

* ``fwhm_chunk_scaling_late_sg``
Similar keyword as ``fwhm_chunk_scaling_early_sg`` but for the case of late (hence cool) subgiant stars, close to the transition to the low RGB, having :math:`\Delta\nu_\mathrm{RG} < \Delta\nu < \Delta\nu_\mathrm{SG}` and :math:`T_\mathrm{eff}` towards the lower limit of 4500 K. This value is set to 3.0, allowing to obtain a finer resolution on the spectrum due to the significantly narrower mode linewidths found for these stars as compared to those of the early subgiants. Again this scaling varies linearly with :math:`T_\mathrm{eff}` and it decreases as the temperature increases, up to the limiting value imposed by ``fwhm_chunk_scaling_early_sg``.

* ``fwhm_chunk_scaling_rg``
Similar keyword as ``fwhm_chunk_scaling_ms`` but for the case of low-luminosity RGB stars, having :math:`\Delta\nu_\mathrm{CL2} < \Delta\nu \leq \Delta\nu_\mathrm{RG}`. The default value is set to 5 to accomodate the presence of very narrow oscillation peaks arising from dipolar mixed modes.

* ``fwhm_chunk_scaling_cl``
Similar keyword as ``fwhm_chunk_scaling_ms`` but for the case of primary and secondary clump stars, having :math:`\Delta\nu_\mathrm{AGB} < \Delta\nu \leq \Delta\nu_\mathrm{CL2}`. The default value is set to 8 to accomodate the presence of very narrow oscillation peaks arising from dipolar mixed modes and the large mode linewidths caused by the hotter temperatures of clump stars relative to RGB stars.

* ``fwhm_chunk_scaling_tip``
Similar keyword as ``fwhm_chunk_scaling_ms`` but for the case of stars evolving toward the tip of the RGB, or into early AGB, having :math:`\Delta\nu \leq \Delta\nu_\mathrm{AGB}`. The default value is set to 12 to accomodate the very narrow oscillation mode linewidths found in these cool stars.

* ``dnu_overlap_fraction_ms``
The fraction of :math:`\Delta\nu` by which each chunk identified by GLOBAL is extended on the left side of the range for MS stars, having either :math:`\Delta\nu_\mathrm{RG} < \Delta\nu \leq \Delta\nu_\mathrm{SG}` and :math:`T_\mathrm{eff} \geq T_\mathrm{eff,SG}` or :math:`\Delta\nu \geq \Delta\nu_\mathrm{SG}`. The default value is set to 0.15.

* ``dnu_overlap_fraction_rg``
The fraction of :math:`\Delta\nu \geq \Delta\nu_\mathrm{SG}` by which each chunk identified by GLOBAL is extended on the left side of the range for all stars not classified as MS stars, having either :math:`\Delta\nu < \Delta\nu_\mathrm{SG}` and :math:`T_\mathrm{eff} < T_\mathrm{eff,SG}`. This is set by default to 0.25, thus larger than ``dnu_overlap_fraction_ms``, to allow for more pronounced curvature effects in evolved stars and for accomodating their more complex oscillation mode patterns.

* ``threshold_asef_global``
The fraction of ASEF amplitude with respect to the global ASEF maximum, that is necessary to consider to locate a local maximum using the hill-climbing algorithm. This amplitude refers to the total amount of amplitude found in a rising phase of a local maximum, i.e. starting from its closest local minimum. This parameter is used only within the GLOBAL module of the pipeline. The default value is set to 0.01.

* ``threshold_asef_chunk_ms``
Similar description as for ``threshold_asef_global`` but for MS stars only, having either :math:`\Delta\nu_\mathrm{RG} < \Delta\nu \leq \Delta\nu_\mathrm{SG}` and :math:`T_\mathrm{eff} \geq T_\mathrm{eff,SG}` or :math:`\Delta\nu \geq \Delta\nu_\mathrm{SG}`. This parameter is used within the CHUNK module of the pipeline. The default value is set to 0.03.

* ``threshold_asef_chunk_sg``
Similar description as for ``threshold_asef_global`` but for subgiant stars, having :math:`\Delta\nu_\mathrm{RG} < \Delta\nu < \Delta\nu_\mathrm{SG}` and :math:`T_\mathrm{eff} < T_\mathrm{eff,SG}`. This parameter is used within the CHUNK module of the pipeline. The default value is set to 0.05, which is larger than that used in MS stars because the stellar spectra appear to be more complex.

* ``threshold_asef_chunk_rg``
Similar description as for ``threshold_asef_global`` but for red giant stars, having :math:`\Delta\nu_\mathrm{AGB} \leq \Delta\nu < \Delta\nu_\mathrm{RG}`. This parameter is used within the CHUNK module of the pipeline. The default value is set to 0.05, which similarly as for subgiant stars, is larger than that used in MS stars because the stellar spectra appear to be more complex.

* ``threshold_asef_chunk_tip``
Similar description as for ``threshold_asef_global`` but for red giant stars evolved toward the RGB tip and for early AGB stars, having :math:`\Delta\nu < \Delta\nu_\mathrm{AGB}`. This parameter is used within the CHUNK module of the pipeline. The default value is set to 0.03, which similarly as for main sequence stars is smaller to increase the chances of detecting peaks with a low number of databins.

* ``skim_frequency_tolerance_ms``
The fraction in units of :math:`\Delta\nu` by which the frequencies extracted from the multi-modal sampling in the GLOBAL module should match within the predictions from the asymptotic pattern of p modes. When the extracted frequencies exceed this fraction with respect to their asymptotic prediction, they are excluded from the final sample. This parameter regulates the skimming process of radial and dipole modes in the case of MS stars, having either :math:`\Delta\nu_\mathrm{RG} < \Delta\nu \leq \Delta\nu_\mathrm{SG}` and :math:`T_\mathrm{eff} \geq T_\mathrm{eff,SG}` or :math:`\Delta\nu \geq \Delta\nu_\mathrm{SG}`. See Corsaro et al. (2020) for more details. The default value is set to 0.18.

* ``skim_frequency_tolerance_sg``
Similar description as for ``skim_frequency_tolerance_ms`` but for subgiant stars, having :math:`\Delta\nu_\mathrm{RG} < \Delta\nu < \Delta\nu_\mathrm{SG}` and :math:`T_\mathrm{eff} < T_\mathrm{eff,SG}`. The default value is set to 0.20 to accomodate for a more complex pattern of the oscillation modes found in these stars.

* ``skim_frequency_tolerance_rg``
Similar description as for ``skim_frequency_tolerance_ms`` but for red giant stars, having :math:`\Delta\nu_\mathrm{AGB} < \Delta\nu \leq \Delta\nu_\mathrm{RG}`, thus excluding high-luminosity RGB stars. The default value is set to 0.25 for overcoming potential issues arising from pronounced curvature effects and from the complex stellar spectra found in these stars.

* ``skim_frequency_tolerance_tip``
Similar description as for ``skim_frequency_tolerance_ms`` but for high-luminosity RGB stars (including also early AGB), having :math:`\Delta\nu \leq \Delta\nu_\mathrm{AGB}`. The default value is set to 0.20 because of the simpler structure of the oscillation spectra found in these stars with respect to stars with larger value of :math:`\Delta\nu`.

Multi-modal sampling analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* ``min_sep_scaling_global``
This is the number of times :math:`\Delta\nu` is divded by in order to define the minimum separation in frequency that is used to evaluate the total number of bins in the ASEF of the GLOBAL module (see also the parameter ``min_n_bins`` below). This is based on an estimate of the minimum frequency width required to obtain a good sampling of the actual frequency peaks obtained in the multi-modal fits. The default value is set to 2.0.

* ``min_sep_scaling_chunk_ms`` 
Similar as for ``min_sep_scaling_global`` but in the case of the ASEF for the CHUNK module, and for stars less evolved than red giants, having :math:`\Delta\nu \geq \Delta\nu_\mathrm{RG}`. The scaling is here referred to the small frequency separation :math:`\delta\nu_\mathrm{02}`, where this latter one is obtained from the :math:`\Delta\nu`-:math:`\delta\nu_\mathrm{02}` relation by Huber et al. (2010). The default value is set to 2.0.

* ``min_sep_scaling_chunk_rg``
Similar as for ``min_sep_scaling_global`` but in the case of the ASEF for the CHUNK module, and for stars evolved into red giants, having :math:`\Delta\nu < \Delta\nu_\mathrm{RG}`. The default value is set to 20.0 to provide a finer ASEF that is able to detect the presence of very narrow mixed modes in the spectrum. This value is likely to be updated in future versions of the pipeline by adopting a more accurate estimate exploiting the relation of the mode linewidth with the frequency in the stellar PSD.

* ``min_n_bins``
This is the number of ASEF bins one wants to represent the minimum separation in frequency evaluated using the previous parameter(s). This parameter is used to define the actual binwidth of the ASEF. The default value is set to 6. The total number of ASEF bins is thus computed automatically as the total frequency range of the stellar PSD being analyzed divided by the obtained binwidth.

* ``n_bins_acf_scaling``
The value of scaling to define the total number of bins when the ASEF is computed in high-resolution. The default value is 10, meaning that for a given number of total nested iterations :math:`N_\mathrm{nest}`, a number of ASEF bins equal to :math:`N_\mathrm{nest}/`n_bins_acf_scaling is used. This corresponds to 800 bins in the default setup of the pipeline.

* ``dnu_acf_range_side``
This parameter defines the fraction of :math:`\Delta\nu` being varied with respect to the value of :math:`\Delta\nu` given by the :math:`\Delta\nu-\nu_\mathrm{max}` relation by Huber et al. (2011) and it is used to set up the search frequency range for :math:`\Delta\nu` in the :math:`\mbox{ACF}^2` in the GLOBAL module. The default value is 0.3, meaning a 30% variation of :math:`\Delta\nu` on each side.

* ``min_bin_separation``
This is the minimum number of ASEF bins to consider for having two local maxima separated from one another. If the hill-climbing algorithm identifies two local maxima in two different bins that are within less than this number from one another, then only the highest amplitude one will be considered. The default value is set to 2.

* ``n_bins_max_fraction``
This parameter controls the maximum extent of ASEF bins that is allowed for a single frequency range :math:`r_\mathrm{a,b}`, in the CHUNK module. It is a number by which the total number of ASEF bins of the chunk is divided. The default value is set to 8. This is useful especially for MS stars where a big local maximum can be produced, resulting in ranges that cover a large part of the frequency range inspected, and therefore producing unreliably large frequency uncertainties.

* ``drop_tolerance_global``
This is the fraction of the ASEF amplitude of the local maximum being inspected within the GLOBAL module. This value is considered as a tolerance to obtain the frequency ranges :math:`r_\mathrm{a,b}` around the local maximum considered. This parameter is useful for cases where the local maxima produced in the ASEF are spread over multiple ASEF bins, with the ASEF in each bin showing small fluctuations, or for cases where the ASEF outside the local maximum is quite flat, which makes it difficult to properly define a frequency range. The default value is set to 0.01.

* ``drop_tolerance_chunk_ms``
Similar as for ``drop_tolerance_global`` but in the case of the CHUNK module and for stars in the MS, having :math:`\Delta\nu > \Delta\nu_\mathrm{SG}`. The default value is 0.06 to allow for bigger fluctuations as the result of broader peak structures in the stellar PSD.

* ``drop_tolerance_chunk_rg``
Similar as for ``drop_tolerance_global`` but in the case of the CHUNK module and for evolved stars, having :math:`\Delta\nu \leq \Delta\nu_\mathrm{SG}`. The default value is 0.03, smaller than that of MS stars because stronger fluctuations are expected in more evolved stars as the result of the presence of narrower oscillation modes.

* ``max_iterations_frequency``
The maximum number of iterations to improve the evaluation of the estimated frequency :math:`\nu_{f,i}` and uncertanties :math:`\sigma_{f,i}`. The default value is set to 2, meaning that once :math:`r_\mathrm{a,b}` are evaluated for the first time, the frequency estimates and uncertanties are recomputed once more.

* ``max_sigma_range``
Within the loop done with the previous parameter, the frequency ranges :math:`r_\mathrm{a,b}` can be adjusted by exploiting the newly computed frequency uncertainty :math:`\sigma_{f,i}`. This parameter sets the maximum number of :math:`\sigma_{f,i}` to be considered to redefine the frequency ranges on each side of the estimated frequency :math:`\nu_{f,i}`. This parameter is used only if the current frequency range exceeds this upper bound imposed by the new frequency uncertainty. The default value is set to 2.

* ``min_sigma_range``
Similar as for ``max_sigma_range`` but now defining the minimum number of :math:`\sigma_{f,i}` to be considered to redefine the frequency ranges on each side of the estimated frequency :math:`\nu_{f,i}`. This parameter is used only if the current frequency range is smaller than this lower bound imposed by the new frequency uncertainty. The default value is set to 1.

* ``max_skim_iterations_global``
The maximum number of iterations to skim the set of estimated frequencies :math:`\nu_{f,i}` obtained within the GLOBAL module. See Corsaro et al. (2020) for more details on the skimming process. The default value is set to 2.

* ``alpha_radial_universal``
The curvature term of the asymptotic relation of radial modes as measured by Mosser et al. (2011) in their universal pattern. The default value is 0.008.

* ``correct_radial_frequencies``
A flag specifying whether or not the radial mode frequencies obtained from the GLOBAL module should be corrected by the shift with respect to the central radial mode frequency obtained from the sliding pattern fit (see Corsaro et al. 2020 for more details). This correction generally improves the position of the global radial mode frequencies, thus allowing a better decomposition of the PSD into chunks. The default value is 1, meaning that the correction is applied.

* ``save_asymptotic_radial``
A flag specifying whether or not the global radial mode frequencies obtained within the GLOBAL module are in the end replaced by their asymptotic predictions. This may turn out to provide more stable estimates of the radial mode frequencies to be used for decomposing the PSD into chunks, especially toward the tails of the Gaussian envelope of the oscillations. In these regions the frequencies extracted from the ASEF and the subsequent mode identification can be more affected by the presence of noise peaks, by mixed modes, and by the curvature effects of the asymptotic pattern of p modes. This correction is mostly useful for stars with low signal-to-noise ratio. The default value is set to 0, meaning that the asymptotic predictions of the global radial modes are not used as a reference to identify the chunks.

* ``separations_dnu_tolerance_rg``
When defining the separation s_n to the right side of each global radial mode frequency, an additional fraction of :math:`\Delta\nu`, set by this parameter, is added to the value of the asymptotic prediction for that radial mode. This is used in the GLOBAL module for stars in the SG and RG regime, having :math:`\Delta\nu < \Delta\nu_\mathrm{SG}` and :math:`T_\mathrm{eff} < T_\mathrm{eff,SG}`. The default value is set to 0.25.

* ``separations_dnu_tolerance_ms``
Similar as for the parameter ``separations_dnu_tolerance_rg`` but for stars classified as MS, having :math:`\Delta\nu_\mathrm{RG} < \Delta\nu \leq \Delta\nu_\mathrm{SG}` and :math:`T_\mathrm{eff} \geq T_\mathrm{eff,SG}` or :math:`\Delta\nu \geq \Delta\nu_\mathrm{SG}`. The default value is set 0.20, lower than that of SG and RGs because of the less pronounced curvature effects in the asymptotic pattern.

* ``weight_freq_fraction``
This parameters controls the weight assigned to the frequency difference between the global radial mode frequency and the candidate chunk radial mode frequency. It allows obtaining an overall weight when combined with other information (see below), to identify the frequency peak that best corresponds to an expected radial mode frequency inside a chunk. It is being used in the CHUNK module, in case a global radial mode frequency was found from the GLOBAL module. The default value is 1.0. 

* ``weight_freq_fraction_enhanced``
This parameter replaces ``weight_freq_fraction`` in the case a chunk radial mode frequency is found from a neighboring chunk. The default value is 1.5. It assigns a significantly larger weight than in the case of ``weight_freq_fraction`` because the neighboring chunk radial mode frequency is assumed to be accurate.

* ``weight_asef_fraction``
Similar description as for ``weight_freq_fraction`` but referring to the maximum ASEF value of the frequency peak being considered. The larger the ASEF maximum, the more it will impact on the overal weight. The default value is 0.1. 

* ``weight_asef_integral_fraction``
Similar description as for ``weight_freq_fraction`` but referring to the integral value of the ASEF for the frequency peak being considered and evaluated within its maximum range. The larger the ASEF integral, the more it will impact on the overal weight. The default value is 0.8. 

* ``weight_spsd_fraction``
Similar description as for ``weight_freq_fraction`` but referring to the maximum value of smoothed PSD of the frequency peak being considered. The larger the smoothed PSD maximum, the more it will impact on the overal weight. The default value is 0.1. 

* ``weight_sampling_fraction``
Similar description as for ``weight_freq_fraction`` but referring to the number of sampling points (in logarithmic units) that are contained within the frequency peak being considered. The larger the number of sampling points, the more the peak will be prominent, and therefore it will impact more on the overal weight. The default value is 0.2. 

* ``upper_limit_freq_radial``
The limiting upper frequency of the chunk in the CHUNK module. The user has the possibility to input this frequency to force the chunk frequency range up to a specific value. The default value is 0, meaning that the limiting frequency is evaluated automatically during the analysis.

* ``plot_weights_radial``
Allows plotting the weights for candidate radial mode peaks as a function of the frequency extracted inside the chunk in the CHUNK module. Can be useful for debugging purposes. The default value is 0, meaning that no plotting is produced.

* ``plot_weights_dipole``
Allows plotting the weights for candidate dipole mode peaks as a function of the frequency extracted inside the chunk in the CHUNK module. Can be useful for debugging purposes. The default value is 0, meaning that no plotting is produced.

* ``threshold_search_radial_asef_integral``
In the search for the proper chunk radial mode frequency in the CHUNK module, this parameter sets a lower limit threshold on the weight of the integral of the ASEF for the peak under inspection. This parameter is used to check whether the selected peak is not the adjacent quadrupole mode. If the weight of the new peak under inspection exceeds the fraction of weight imposed by this threshold (with respect to the weight of the peak initially selected as a radial mode), then the peak can be considered for further analysis as a candidate radial mode peak. The default value is 0.90.

* ``threshold_search_radial_asef_maximum``
n the search for the proper chunk radial mode frequency in the CHUNK module, this parameter sets a lower limit threshold on the weight of the maximum of the ASEF for the peak under inspection. This parameter is used to check whether the selected peak is not the adjacent quadrupole mode. If the weight of the new peak under inspection exceeds the fraction of weight imposed by this threshold (with respect to the weight of the peak initially selected as a radial mode), then the peak can be considered for further analysis as a candidate radial peak. The default value is 0.90.

* ``previous_radial_range_fraction``
The fraction of frequency range of the chunk that is used to locate the first (potential) frequency candidates of the same chunk. These frequencies are used to check whether a potential radial mode from the previous chunk is still present, so that it can be excluded from further analysis. This is used only for chunks having a mean frequency above :math:`\nu_\mathrm{max}`. The default value is 8.

* ``sampling_counts_fraction``
This is the fraction of sampling counts (i.e. number of sampling points) with respect to those of the actual chunk radial mode that has to be verified in order to adjust the limiting lower frequency of the chunk. This is because in the Gaussian tail corresponding to frequencies above :math:`\nu_\mathrm{max}`, one usually expects the radial modes of chunks with smaller mean frequency to be more prominent then their higher frequency neighbors. If this condition is verified, the current chunk radial mode may also be improved. The default value is 1.3.

* ``asef_threshold_scaling_radial``
If a radial mode from a previous chunk is identified within the range defined by the previous parameter, then this parameter is used to define a threshold on the ASEF maximum of the frequency peak of the new candidate radial mode of the chunk. If the ASEF maximum of the new peak is above this threshold, then the chunk radial mode is updated. The default value is set to 3.

* ``dnu_lower_cut_fraction``
A safety condition on the limiting lower frequency of the chunk in the CHUNK module. This limiting lower frequency can never exceed the frequency obtained by subtracting a fraction of :math:`\Delta\nu` defined by this parameter to the frequency of the chunk radial mode. The default value is 0.75, meaning that the chunk frequency range has to extend to at least 75% of :math:`\Delta\nu` below the chunk radial mode frequency.

* ``low_cut_frequency``
The input value for the limiting lower frequency of the chunk in the CHUNK module. If specified by the user, it will force the range to the input value. The default value is 0, meaning that the limiting frequency is evaluated automatically in the analysis.

* ``d02_scaling_merge_mixed``
A parameter specifying the scaling to :math:`\delta\nu_{02}` to define a range around the chunk quadrupole frequency that is used to merge potential quadrupole mixed modes into a single quadrupole mode. This is used in the CHUNK module for RG stars, having :math:`\Delta\nu < \Delta\nu_\mathrm{RG}`. The default value is 2, meaning that the range is defined as :math:`\left[ \nu_2 - \delta\nu_{02}/2, \nu_2 + \delta\nu_{02}/2 \right]`, with :math:`\nu_2` the chunk quadrupole frequency.

* ``d02_factor_search_range``
This parameter is used to compute the upper prior bound for the free parameter :math:`\delta\nu_{02}` in the duplet fit for MS and SG stars. This is used in the CHUNK module in case other solutions for :math:`\delta\nu_{02}` are already available from other chunks. The default value is 1.3.

* ``d02_fraction_prior_lower_duplet_fit``
The fraction to set the lower prior bound for the free parameter :math:`\delta\nu_{02}` in the duplet fit for MS and SG stars. This is used in the CHUNK module, where the default value is 0.5 starting from the median :math:`\delta\nu_{02}` evaluated from the available chunks.

* ``d02_fraction_prior_upper_duplet_fit``
The fraction to set the upper prior bound for the free parameter :math:`\delta\nu_{02}` in the duplet fit for MS and SG stars. This is used in the CHUNK module, where the default value is 1.3 starting from the median :math:`\delta\nu_{02}`, and it is applied only if such median is computed from other existing chunks.

* ``d02_prior_upper_duplet_fit``
The upper prior bound for the free parameter :math:`\delta\nu_{02}` in the duplet fit for MS and SG stars. This is used in the CHUNK module as an input value to force the prior boundary to a specific number. It is set by default to 0, meaning that the parameter ``d02_factor_search_range`` is defining the upper prior bound instead of this one.

* ``d03_upper_scaling_factor``
In case of stars classified as RGs, having :math:`\Delta\nu < \Delta\nu_\mathrm{RG}`, this parameter is the scaling value to be multiplied to the asymptotic prediction :math:`\delta\nu_{03}`. It is used to define a lower bound for the search range of the octupole modes in the CHUNK module (see also the next parameter). The default value is 1.12. This scaling can be directly adopted on a prediction for :math:`\delta\nu_{03}` because the relation :math:`\Delta\nu`-:math:`\delta\nu_{03}` is rather tight for RG stars.

* ``d03_lower_scaling_factor``
Similar as the previous parameter but defining the upper bound for the search range of the octupole modes. The default value is 0.84. The octupole mode search range is thus defined as :math:`\left[ \nu_0 - \Delta\nu/2 - 1.12 \delta\nu_{03}, \nu_0 - \Delta\nu/2 - 0.84 \delta\nu_{03} \right]`.

* ``d02_upper_scaling_factor_ms``
In case of stars classified as MS, having either :math:`\Delta\nu_\mathrm{RG} < \Delta\nu \leq \Delta\nu_\mathrm{SG}` and :math:`T_\mathrm{eff} \geq T_\mathrm{eff,SG}` or :math:`\Delta\nu \geq \Delta\nu_\mathrm{SG}`, this parameter is the scaling value to the asymptotic prediction for :math:`\delta\nu_{03}` as a function of :math:`\delta\nu_{02}`, i.e. :math:`\delta\nu_{03} \simeq 2\delta\nu_{02}` (see Corsaro et al. 2020, Eq. 26). It is used to define a lower bound for the search range of the octupole modes in the CHUNK module (see also the next parameter). The default value is 4.0.

* ``d02_lower_scaling_factor_ms``
Similar as the previous parameter but defining the upper bound for the search range of the octupole modes in MS stars. The default value is 1.5. The octupole mode search range is thus defined as :math:`\left[ \nu_0 - \Delta\nu/2 - 4.0 \delta\nu_{02}, \nu_0 - \Delta\nu/2 - 1.5 \delta\nu_{02} \right]`.

* ``d02_upper_scaling_factor_sg``
Analogous to the parameter ``d02_upper_scaling_factor_ms`` but for stars classified as SGs, having :math:`\Delta\nu_\mathrm{RG} < \Delta\nu \leq \Delta\nu_\mathrm{SG}` and :math:`T_\mathrm{eff} \geq T_\mathrm{eff,SG}`. The default value is 3.0. This value is smaller than that used for MS stars because the relation :math:`\Delta\nu`-:math:`\delta\nu_{02}` is more constrained for SG stars.

* ``d02_lower_scaling_factor_sg``
Analogous to the parameter ``d02_lower_scaling_factor_ms`` but for stars classified as SGs, having :math:`\Delta\nu_\mathrm{RG} < \Delta\nu \leq \Delta\nu_\mathrm{SG}` and :math:`T_\mathrm{eff} \geq T_\mathrm{eff,SG}`. The default value is 1.8.

* ``smoothing_fwhm_factor_rg``
A smoothing factor to apply to the FWHM adopted to perform the chunk multi-modal sampling. This parameter is used only for RG stars, having :math:`\Delta\nu < \Delta\nu_\mathrm{RG}`. This is used in the CHUNK module to obtain the smoothed PSD. The default value is 3.

Sliding-pattern model fitting
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* ``asef_threshold_fraction``
A threshold value with respect to the global maximum ASEF found in the GLOBAL module analysis. If a sufficient number of peaks extracted from the multi-modal sampling has an ASEF maximum lower than the value obtained from this threshold, then the star is considered as a potential depressed dipole star during the sliding-pattern model fitting (see also the two following parameters). The default value is 0.50.

* ``n_central_orders_side``
The number of radial orders (expressed in units of :math:`\Delta\nu` on each side of :math:`\nu_\mathrm{max}` that are used to define the control frequency range to setup the sliding-pattern model in the GLOBAL module. The default value is 2.0, meaning that four radial orders are considered.

* ``depressed_dipole_fraction``
The fraction of frequency peaks found in the frequency range defined by the parameter ``n_central_orders_side``, having an ASEF maximum lower than the threshold imposed by the parameter ``asef_threshold_fraction``. The default value is 0.5, meaning that in order to classify the star as a potential depressed dipole star, at least half of the ASEF peaks considered must exhibit a low ASEF maximum. 

* ``dnu_ridge_threshold``
A threshold in units of :math:`\Delta\nu` that sets the minimum separation that can be considered for pairs of frequencies (even and odd) computed out of the frequencies found using the parameter ``n_central_orders_side``. This parameter is used to understand whether the frequencies group division into even and odd frequencies is adequate or a double step procedure may be needed (see also Corsaro et al. 2020 for more details). The default value is 0.75.

* ``dnu_echelle_threshold``
A threshold in percentage of :math:`\Delta\nu` that sets the maximum deviation that can be found in an echelle ridge computed out of the frequencies found using the parameter ``n_central_orders_side``. This helps in classifying the star based on how much each frequency is deviating from an expected regular asymptotic pattern for p modes. The default value is 6, meaning that if the maximum deviation found exceeds a 6% variation in :math:`\Delta\nu`, then the star is flagged as an evolved star containing mixed modes. This is useful for setting up the sliding-pattern model as a function of different evolutionary stages.

* ``n_sliding_test``
The number of times the sliding-pattern fit has to be repeated. This improves the reliability of the final outcome for more challenging targets. The default value is 5.

* ``input_radial_freq_reference``
The input value in :math:`\mu\mbox{Hz}` specifying the frequency of the reference radial mode of the star. This frequency is then used to evaluate :math:`\epsilon`. The default value is 0, meaning that no input frequency is used, which is instead obtained from the sliding-pattern fit.

* ``force_epsilon_dnu_value``
A flag specifying that :math:`\epsilon` has to be obtained from the :math:`\epsilon`-:math:`\Delta\nu` relation calibrated by Corsaro et al. (2012)b instead of being evaluated through the sliding-pattern fit. This provides a reasonable estimate only for stars that are evolving along the RGB, having :math:`\Delta\nu \leq \Delta\nu_\mathrm{thresh}`. The default value is set to 0, meaning that this condition is not applied. Forcing this parameter to 1 can be useful if the sliding pattern fit keeps failing in obtaining a correct global mode identification despite other attempts have been made (e.g. see the parameter ``remove_dipole_peak``). However, caution should be used in general when forcing this condition. This is because for red giants that are not RGB (e.g. RC or 2nd RC stars) the :math:`\epsilon` value from the :math:`\epsilon`-:math:`\Delta\nu` relation is likely higher than the actual value of the star, thus yielding a potential incorrect global mode identification. When activating this keyword we recommend double-checking the numerical result with the visual solution of the GLOBAL summary plot in order to validate any mode identification that is obtained by the pipeline. This keyword is automatically not activated if the keyword ``input_radial_freq_reference`` is in use.

* ``n_orders_side_prior_ms``
The number of radial orders on each side of :math:`\nu_\mathrm{max}` that are used to set up the prior boundaries on the central frequency :math:`\nu_0` of the sliding-pattern model within the GLOBAL module. This is used for stars classified as MS, having :math:`\Delta\nu > \Delta\nu_\mathrm{thresh}`, without depressed dipole modes, and without the presence of mixed modes. The default value is 1.5, meaning that the prior covers a total of 3 radial orders.

* ``n_orders_side_prior_sg``
Similar as for the parameter ``n_orders_side_prior_ms`` but used for early SG stars, having :math:`\Delta\nu > \Delta\nu_\mathrm{thresh}` and either depressed dipole modes or with the presence of mixed modes. The default value is 1.5.

* ``n_orders_side_prior_rg``
Similar as for the parameter ``n_orders_side_prior_ms`` but used for late SG and RG stars, having :math:`\Delta\nu \leq \Delta\nu_\mathrm{thresh}`. The default value is 1.5.

* ``remove_dipole_peak``
Activate this keyword to remove the :math:`\ell = 1` peak from the sliding-pattern model for stars having :math:`\Delta\nu < \Delta\nu_\mathrm{thresh}`. This can be useful in some circumstances, for example when the :math:`\ell = 1` mode region is particularly confusing (i.e. crowded), especially if the PSD has a low SNR, or if the star has a very prominent, single, :math:`\ell = 1` peak, which can be confused with a :math:`\ell = 0` peak. This keyword can thus be activated to improve the fit and avoid cases where the pipeline could end up in obtaining a swapped global mode identification (:math:`\ell = 1` identified as :math:`\ell = 0` and viceversa). The keyword is set to 0 by default, meaning that the :math:`\ell = 1` peak is included in the model but kept fixed to the position of a pure :math:`\ell = 1` p mode.

* ``dnu_prior_lower_fraction``
The fraction of :math:`\Delta\nu` with respect to the asymptotic fit value obtained in the GLOBAL modality, used to set up the uniform prior lower bound on :math:`\Delta\nu` for the sliding-pattern model. The default value is 0.96.

* ``dnu_prior_upper_fraction``
Similar as for the parameter ``dnu_prior_lower_fraction`` but now definining the uniform prior upper bound on the parameter :math:`\Delta\nu` of the sliding-pattern model. The default value is 1.04.

* ``d02_prior_lower_ms``
The value of :math:`\delta\nu_{02}`, expressed in :math:`\mu\mbox{Hz}`, used to define the uniform prior lower bound for the sliding-pattern model of MS stars. The default value is 
1.5 :math:`\mu\mbox{Hz}`.

* ``d02_prior_upper_ms``
Similar as for the parameter ``d02_prior_lower_ms`` but now specifying the uniform prior upper bound. The default value is 13 :math:`\mu\mbox{Hz}`. This value is obtained from the findings by White et al. (2011) and Lund et al. (2017).

* ``d02_prior_lower_sg``
Similar as for the parameter ``d02_prior_lower_ms`` but used for SG stars. The default value is 1.5 :math:`\mu\mbox{Hz}`.

* ``d02_prior_upper_sg``
Similar as for the parameter ``d02_prior_lower_sg`` but now specifying the uniform prior upper bound. The default value is 7 :math:`\mu\mbox{Hz}`. This value stems from the findings by White et al. (2011).

* ``d01_prior_lower_ms``
The value of :math:`\delta\nu_{01}`, expressed in :math:`\mu\mbox{Hz}`, used to define the uniform prior lower bound for the sliding-pattern model of MS stars. The default value is 
0.1 :math:`\mu\mbox{Hz}`.

* ``d01_prior_upper_ms``
Similar as for the parameter ``d01_prior_lower_ms`` but now specifying the uniform prior upper bound. The default value is 7 :math:`\mu\mbox{Hz}`. 

* ``d13_prior_lower_ms``
The value of :math:`\delta\nu_{13}`, expressed in :math:`\mu\mbox{Hz}`, used to define the uniform prior lower bound for the sliding-pattern model of MS stars. The default value is 
0.1 :math:`\mu\mbox{Hz}`.

* ``rot_split_prior_lower_ms``
The value of :math:`\delta\nu_\mathrm{rot}`, expressed in :math:`\mu\mbox{Hz}`, used to define the uniform prior lower bound for the sliding-pattern model of MS stars. The default value is 0.1 :math:`\mu\mbox{Hz}`.

* ``rot_split_prior_upper_ms``
Similar as for the parameter ``rot_split_prior_lower_ms`` but now specifying the uniform prior upper bound. The default value is 7 :math:`\mu\mbox{Hz}`.

* ``cosi_prior_lower``
The value of :math:`\cos i`, with i the inclination angle of the rotation axis of the star, used to define the uniform prior lower bound for the sliding-pattern model of stars in any evolutionary stage. The default value is 0, meaning that the minimum value allowed corresponds to an inclination angle of 90 degrees, i.e. the star is seen edge on.

* ``cosi_prior_upper``
Similar as for the parameter ``cosi_prior_lower`` but now specifying the uniform prior upper bound. The default value is 1.0, corresponding to an inclination angle of 0 degrees, i.e. the star is seen pole on.

* ``quadrupole_radial_height_ratio_ms``
The visibility (or heights ratio) of the quadrupole modes with respect to the radial modes, used for the sliding-pattern model of MS stars. The default value is 0.62, obtained following the result from Lund et al. (2017). 

* ``quadrupole_radial_height_ratio_sg``
Similar as for the parameter ``quadrupole_radial_height_ratio_ms`` but used for SG stars. The default value is 0.62, the same as for MS stars.

* ``quadrupole_radial_height_ratio_rg``
Similar as for the parameter ``quadrupole_radial_height_ratio_ms`` but used for RG stars. The default value is 0.8 for incorporating potential broadening from the presence of quadrupole mixed modes (see also Corsaro et al. 2020).

* ``dipole_radial_height_ratio_ms``
The visibility (or heights ratio) of the dipole modes with respect to the radial modes, used for the sliding-pattern model of MS stars. The default value is 1.5, obtained following the result from Lund et al. (2017). 

* ``dipole_radial_height_ratio_sg``
Similar as for the parameter ``dipole_radial_height_ratio_ms`` but used for SG stars. The default value is 1.5, the same as for MS stars.

* ``dipole_radial_height_ratio_rg``
Similar as for the parameter ``dipole_radial_height_ratio_ms`` but used for RG stars. The default value is 0.7, significantly smaller than that adopted for MS and SG stars, because the dipole mode power is spread over a large frequency range and because the sliding-pattern model utilizes only one dipole mode peak per radial order.

* ``octupole_radial_height_ratio``
The visibility (or heights ratio) of the octupole modes with respect to the radial modes, used for the sliding-pattern model for stars in any evolutionary stage. The default value is 0.07, obtained following the result from Lund et al. (2017). 

* ``dipole_radial_fwhm_ratio_ms``
The FWHM magnification factor of dipole modes with respect to the FWHM of radial modes for MS stars used in the sliding-pattern model. The default value is 1, meaning that in this case the linewidth of dipole modes is the same as that of the radial modes.

* ``dipole_radial_fwhm_ratio_sg``
Similar as for the parameter ``dipole_radial_fwhm_ratio_ms`` but used for SG stars. The default value is 1, the same as for MS stars.

* ``dipole_radial_fwhm_ratio_rg``
Similar as for the parameter ``dipole_radial_fwhm_ratio_ms`` but used for RG stars. The default value is 4, significantly larger than that used for less evolved stars, to allow compensating the varying frequency position of the dipole modes in such stars.

* ``upper_epsilon_rgb_slope``
The slope :math:`a` for the linear relation :math:`\epsilon_\mathrm{upper} = a \log \Delta\nu + b`, which provides the upper limit for :math:`\epsilon` in stars with :math:`\Delta\nu_\mathrm{AGB} \leq \Delta\nu \leq \Delta\nu_\mathrm{thresh}`, where :math:`\Delta\nu_\mathrm{AGB}` is set by the configuring parameter ``dnu_agb``. This is used as a control check for the sliding pattern fit to avoid wrong inferences of the central radial mode position. It is particularly useful for low-resolution datasets of intermediate-mass red giant stars, where the unresolved dipole mode region can be confused with a radial mode. The default value is set to 0.253694 as obtained from a linear fit to the sample presented by Kallinger et al. (2012).

* ``upper_epsilon_rgb_offset``
The offset :math:`b` for the linear relation presented for the configuring parameter ``upper_epsilon_rgb_slope``. The default value is set to 0.76. 

* ``lower_epsilon_cl_slope``
The slope :math:`a` for the linear relation :math:`\epsilon_\mathrm{lower} = a \log \Delta\nu + b`, which provides the lower limit for :math:`\epsilon` in stars with :math:`\Delta\nu_\mathrm{AGB} \leq \Delta\nu \leq \Delta\nu_\mathrm{CL2}`, where :math:`\Delta\nu_\mathrm{CL2}` is set by the configuring parameter ``dnu_cl2``. This is used as a control check for the sliding pattern fit to avoid wrong inferences of the central radial mode position. It is particularly useful for low-resolution datasets of intermediate-mass red giant stars, where the unresolved dipole mode region can be confused with a radial mode. The default value is set to 0.480525527 as obtained from a linear fit to the sample presented by Kallinger et al. (2012).

* ``lower_epsilon_cl_offset``
The offset :math:`b` for the linear relation presented for the configuring parameter ``lower_epsilon_cl_slope``. The default value is set to -0.14604031. 

* ``upper_epsilon_evolved_rgb_slope``
The slope :math:`a` for the linear relation :math:`\epsilon_\mathrm{upper} = a \log \Delta\nu + b`, which provides the upper limit for :math:`\epsilon` in stars with :math:`\Delta\nu_\mathrm{tip} \leq \Delta\nu < \Delta\nu_\mathrm{AGB}`, where :math:`\Delta\nu_\mathrm{AGB}` is set by the configuring parameter ``dnu_agb``. This is used as a control check for the sliding pattern fit to avoid wrong inferences of the central radial mode position. It is particularly useful for low-resolution datasets of evolved RGB stars, where the unresolved dipole mode region can be confused with a radial mode. The default value is set to 0.25354086 as obtained from a linear fit to the sample presented by Kallinger et al. (2012).

* ``upper_epsilon_evolved_rgb_offset``
The offset :math:`b` for the linear relation presented for the configuring parameter ``upper_epsilon_evolved_rgb_slope``. The default value is set to 0.78433261. 

* ``lower_epsilon_evolved_rgb_slope``
The slope :math:`a` for the linear relation :math:`\epsilon_\mathrm{lower} = a \log \Delta\nu + b`, which provides the lower limit for :math:`\epsilon` in stars with :math:`\Delta\nu < \Delta\nu_\mathrm{AGB}`, where :math:`\Delta\nu_\mathrm{AGB}` is set by the configuring parameter ``dnu_agb``. This is used as a control check for the sliding pattern fit to avoid wrong inferences of the central radial mode position. It is particularly useful for low-resolution datasets of evolved RGB stars, where the unresolved dipole mode region can be confused with a radial mode. The default value is set to 0.44393910 as obtained from a linear fit to the sample presented by Kallinger et al. (2012). This relation is also used as an indication of the transition between evolved RGB and early AGB stars.

* ``lower_epsilon_evolved_rgb_offset``
The offset :math:`b` for the linear relation presented for the configuring parameter ``lower_epsilon_evolved_rgb_slope``. The default value is set to 0.35913168. 

* ``epsilon_division_rgb_cl_slope``
The slope :math:`a` for the linear relation :math:`\epsilon_\mathrm{division} = a \log \Delta\nu + b`, which provides the lower limit for :math:`\epsilon` in stars with :math:`\Delta\nu_\mathrm{AGB} \leq \Delta\nu \leq \Delta\nu_\mathrm{thresh}`, where :math:`\Delta\nu_\mathrm{AGB}` is set by the configuring parameter ``dnu_agb``. This is used to provide an indication of the evolutionary stage of the star. If the star has :math:`\epsilon < \epsilon_\mathrm{division}` then it is classified as a clump star (either CL or CL2 depending on the range in :math:`\Delta\nu`), otherwise it is classified as a RGB. The default value is set to 0.29988410 as obtained from a linear fit to the sample presented by Kallinger et al. (2012). 

* ``epsilon_division_rgb_cl_offset``
The offset :math:`b` for the linear relation presented for the configuring parameter ``epsilon_rgb_cl_division_slope``. The default value is set to 0.50399553. 

* ``lower_epsilon_agb``
The limiting lower value for :math:`\epsilon` in the regime of early AGB star. The default value is set to 0.3 as suggested from the sample presented by Kallinger et al. (2012). 

Asymptotic code fitting
^^^^^^^^^^^^^^^^^^^^^^^
* ``dnu_prior_lower_fraction_as``
The fraction of :math:`\Delta\nu` with respect to the asymptotic fit value obtained in the GLOBAL modality, which is used for computing the uniform prior lower bound on the free parameter :math:`\Delta\nu` of the asymptotic pattern for radial modes in the Asymptotic code extension of DIAMONDS. The default value is 0.95.

* ``dnu_prior_upper_fraction_as``
Similar as for the parameter ``dnu_prior_lower_fraction_as`` but now defining the uniform prior upper bound on the free parameter :math:`\Delta\nu`. The default value is 1.05.

* ``alpha_prior_lower_as``
The value of the curvature term of the asymptotic relation for p modes, :math:`\alpha`, used to set the uniform prior lower bound of the corresponding free parameter in the Asymptotic code extension of DIAMONDS during the analysis performed in the GLOBAL module. The default value is -0.02, meaning that sign changing curvatures with respect to :math:`\nu_\mathrm{max}` are also allowed.

* ``alpha_prior_upper_as``
Similar as for the parameter ``alpha_prior_lower_as`` but now defining the uniform prior upper bound on the free parameter :math:`\alpha`.

* ``epsi_prior_lower_fraction_as``
The fraction of the phase term of the asymptotic relation of p modes, :math:`\epsilon`, which is used to set up the uniform prior lower bound on the corresponding free parameter in the Asymptotic code extension of DIAMONDS during the analysis performed in the GLOBAL module. The default value is 0.9.

* ``epsi_prior_upper_fraction_as``
Similar as for the parameter ``epsi_prior_lower_fraction_as`` but now defining the uniform prior upper bound on the free parameter :math:`\epsilon`. The default value is 1.1.

Peak testing
^^^^^^^^^^^^
* ``n_fwhm_fit``
The number of times the fit to estimate the FWHM of a peak within the CHUNK module is repeated. This is used to evaluate both the FWHM of the radial mode and of the potential octupole mode of the chunk. If ``n_fwhm_fit`` is smaller or equal to the number of available CPU threads, which is set in the parameter ``n_threads``, then the fits are all run in parallel. The default value is 3.

* ``n_peak_test``
The number of times a peak fit for the peak significance test is repeated. Each repetition is not run in parallel with the others. However, all the peak fits within a single repetition are run in parallel, with a maximum number of parallel processes set according to the parameter ``n_threads``. The default value is 2. 

* ``dipole_search_activated``
A flag specifying whether dipole and octupole modes should be considered in the CHUNK analysis. The default value is set to 1, meaning that the pipeline will automatically perform a search and detection test (including rotation and duplicity tests if activated too) for candidate :math:`\ell = 1` modes, as well as the search for the :math:`\ell = 3` mode of the chunk. If only :math:`\ell = 0,2` modes are required within the analysis, e.g. for modelling purporses, then dipole and octupole modes can be discarded. Setting this flag to 0 will significantly speed up the overall computation.

* ``save_test_files``
A flag specifying whether the output files to run the PeakBagging code for the peak testing are to be stored. The default value is 0, meaning that these files are not kept in memory. This option is only recommended for debugging purposes, or for particular testings, since the large number of peak fits produced during the peak testing phase typically yields a large amount of output files and therefore of occupied memory space.

* ``detection_probability_threshold``
The value of the threshold, as a probability, used for the peak significance tests. The default value is set to 0.993, corresponding to a strong evidence condition according to the Jeffreys' scale of strength (see Trotta 2008).

* ``rotation_probability_threshold``
The value of the threshold, as a probability, used for the peak rotation tests. The default value is set to 0.75, corresponding to a weak evidence condition according to the Jeffreys' scale of strength (see Trotta 2008).

* ``duplicity_probability_threshold``
The value of the threshold, as a probability, used for the peak duplicity tests. The default value is set to 0.75, corresponding to a weak evidence condition according to the Jeffreys' scale of strength (see Trotta 2008).

* ``rotation_test_activated``
A flag specifying whether the rotation test has to be activated. The default value is set to 0, meaning that no rotation is inspected in the individual peaks. We typically recommend to activate this option when analysis evolved stars, in order to get the most of the peaks identified already during the analysis performed by the CHUNK module.

* ``height_ratio_threshold``
The minimum threshold on the ratio between the height of the peak to be tested and the level of the background that is required to automatically deem the peak as significant. The default value is set to 10 according to the findings by Corsaro et al. (2015)b, meaning that the height of the peak has to be at least 10 times that of the background in order to skip the peak significance test and consider the peak as automatically detected.

* ``bkg_prior_lower``
The uniform prior lower bound on the free parameter :math:`\sigma_\mathrm{noise}` of the peak fitting models :math:`\mathcal{M}_\mathrm{A}`, :math:`\mathcal{M}_\mathrm{B}`, :math:`\mathcal{M}_\mathrm{C}`, and :math:`\mathcal{M}_\mathrm{D}` (see Corsaro et al. 2020 for more details). The default value is 0.95.

* ``bkg_prior_upper``
Similar as for the parameter ``bkg_prior_lower`` but now defining the uniform prior upper bound on the free parameter :math:`\sigma_\mathrm{noise}`. The default value is 1.05.

* ``rot_split_scaling``
The parameter :math:`\theta` defined in Table A.2 of Corsaro et al. (2020) that specifies the rescaling to apply to the uniform prior upper bound on the free parameter :math:`\delta\nu_\mathrm{rot}` of the peak rotation model :math:`\mathcal{M}_\mathrm{F}`. The default value is 2.8.

* ``fwhm_lower_bound``
The minimum allowed FWHM for the peak fitting that is used as a uniform prior lower bound in all the peak testing models that incorporate the estimation of the peak FWHM. The default value is set to :math:`10^{-4} \mu\mbox{Hz}`.

* ``fwhm_magnification_factor_radial``
A magnification factor to set up the uniform prior upper bound for the fitting of the FWHM of radial modes. The default value is 3.

* ``fwhm_magnification_factor_dipole``
A magnification factor to set up the uniform prior upper bound for the fitting of the FWHM of dipole modes. The default value is 1.5

* ``fwhm_magnification_factor_quadrupole``
A magnification factor to set up the uniform prior upper bound for the fitting of the FWHM of quadrupole modes. The default value is 4.

* ``fwhm_magnification_factor_octupole``
A magnification factor to set up the uniform prior upper bound for the fitting of the FWHM of octupole modes. The default value is 1.5

* ``asef_octupole_fraction``
A threshold on the value of the ASEF maximum of the peak to be tested as an octupole mode. The peak that is tested can be considered as candidate octupole mode only if its ASEF maximum does not exceed the value imposed by this threshold. The default value is set to 0.75.

* ``fwhm_octupole_radial_fraction``
The minimum ratio between the FWHM of an octupole mode and that of a radial mode of the same chunk. If this threshold is verified, then the mode is likely considered as a real octupole mode, following the octupole peak test presented by Corsaro et al. (2020). The default value is 0.2.

* ``dnu_mixed_modes_separation_scaling``
A scaling parameter to :math:`\Delta\nu` that is used to avoid that two neighboring peaks found to be both significant happen to be selected together during the analysis performed by the CHUNK module. This improves the reliability of the results for difficult stars where the fluctuations from the noise and the large mode linewidths can give rise to multiple peak-like structures in the dipole mode region. Since for early SGs mixed modes are greatly separated in frequency, we do not expect to have multiple detection within a small fraction of :math:`\Delta\nu` in the same chunk. The default value is 16.

Asteroseismic relations
^^^^^^^^^^^^^^^^^^^^^^^
* ``d03_offset``
The constant offset :math:`\beta_{03}` of the law :math:`\delta\nu_{03} = \alpha_{03} \Delta\nu + \beta_{03}` as calibrated by Huber et al. (2010). The asymptotic parameter :math:`\delta\nu_{03}` refers to the distance between the midpoint of two adjacent radial modes and the position of the octupole mode contained between the two radial modes. The default value is 0.16.

* ``d03_slope``
The multiplication factor (slope) :math:`\alpha_{03}` of the law :math:`\delta\nu_{03} = \alpha_{03} \Delta\nu + \beta_{03}`, as calibrated by Huber et al. (2010). The default value is 0.282.

* ``d01_mass_offset``
The offset :math:`b_{01}` of the relation :math:`\alpha_{01} = a_{01} \frac{M}{M_{\odot}} + b_{01}` connecting stellar mass to the multiplication factor :math:`\alpha_{01}` of the law :math:`\delta\nu_{01} = \alpha_{01} \Delta\nu + \beta_{01}`, as calibrated by Corsaro et al. (2012)b. The default value is set to -0.073.

* ``d01_mass_slope``
The multiplication factor (slope) :math:`a_{01}` of the relation :math:`\alpha_{01} = a_{01} \left( \frac{M}{M_{\odot}} \right) + b_{01}` connecting stellar mass to the multiplication factor :math:`\alpha_{01}` of the law :math:`\delta\nu_{01} = \alpha_{01} \Delta\nu + \beta_{01}`, as calibrated by Corsaro et al. (2012)b. The default value is set to 0.044.

* ``d01_offset``
The constant offset :math:`\beta_{01}` of the law :math:`\delta\nu_{01} = \alpha_{01} \Delta\nu + \beta_{01}`, as calibrated by Corsaro et al. (2012)b. The default value is set to -0.063.

* ``d02_mass_offset``
The offset :math:`b_{02}` of the relation :math:`\alpha_{02} = a_{02} \left( \frac{M}{M_{\odot}} \right) + b_{02}` connecting stellar mass to the multiplication factor :math:`\alpha_{02}` of the law :math:`\delta\nu_{02} = \alpha_{02} \Delta\nu + \beta_{02}`, as calibrated by Corsaro et al. (2012)b. The default value is set to 0.138.

* ``d02_mass_slope``
The multiplication factor (slope) :math:`a_{02}` of the relation :math:`\alpha_{02} = a_{02} \left( \frac{M}{M_{\odot}} \right) + b_{02}` connecting stellar mass to the multiplication factor :math:`\alpha_{02}` of the law :math:`\delta\nu_{02} = \alpha_{02} \Delta\nu + \beta_{02}`, as calibrated by Corsaro et al. (2012)b. The default value is set to -0.014.

* ``d02_offset``
The constant offset :math:`\beta_{02}` of the law :math:`\delta\nu_{02} = \alpha_{02} \Delta\nu + \beta_{02}`, as calibrated by Corsaro et al. (2012)b. The default value is set to 0.035.

* ``d02_unique_offset``
The constant offset :math:`\beta_{02}` of the law :math:`\delta\nu_{02} = \alpha_{02} \Delta\nu + \beta_{02}`, as calibrated by Huber et al. (2010). This law has no dependency on stellar mass. The default value is set to 0.047.

* ``d02_unique_slope``
The constant multiplication factor (slope) :math:`\alpha_{02}` of the law :math:`\delta\nu_{02} = \alpha_{02} \Delta\nu + \beta_{02}`, as calibrated by Huber et al. (2010). The default value is set to 0.121.

* ``dnu_threshold``
A threshold in :math:`\Delta\nu`, known as :math:`\Delta\nu_\mathrm{thresh}`, to distinguish between the regimes of MS and early SGs and that of late SG and RG stars. The default value is set to 20 :math:`\mu\mbox{Hz}`. This is only an indicative separation that is used by the pipeline to optimize the analysis and has not to be considered as a physically meaningful value for discriminating the evolutionary stage of stars. Based on this threshold, :math:`\epsilon` is verified in either the :math:`\epsilon`-:math:`\Delta\nu` (for :math:`\Delta\nu <= \Delta\nu_\mathrm{thresh}`) or the :math:`\epsilon` - *T*:subscript:`eff` (for :math:`\Delta\nu > \Delta\nu_\mathrm{thresh}`) diagram.

* ``numax_threshold``
A threshold in :math:`\nu_\mathrm{max}`, known as :math:`\nu_\mathrm{max,thresh}`, to distinguish different asteroseismic relations for the evaluation of :math:`\Delta\nu` and of the FWHM of the radial modes (see Corsaro et al. 2020 for more details). The default value is set to 300 :math:`\mu\mbox{Hz}`.

* ``numax_coeff_high``
The multiplication factor of the :math:`\nu_\mathrm{max}`-:math:`\Delta\nu` relation calibrated by Huber et al. (2011) for :math:`\nu_\mathrm{max} > 300 \mu\mbox{Hz}`. The default value is 0.22.

* ``numax_exponent_high``
The exponent of the :math:`\nu_\mathrm{max}`-:math:`\Delta\nu` relation calibrated by Huber et al. (2011) for :math:`\nu_\mathrm{max} > 300 \mu\mbox{Hz}`. The default value is 0.797.

* ``numax_coeff_low``
Similar to the parameter ``numax_coeff_high`` but for stars having :math:`\nu_\mathrm{max} \leq 300 \mu\mbox{Hz}`. The default value is 0.267.

* ``numax_exponent_low``
Similar to the parameter ``numax_coeff_low`` but for stars having :math:`\nu_\mathrm{max} \leq 300 \mu\mbox{Hz}`. The default value is 0.760.

* ``epsilon_offset``
The offset of the :math:`\epsilon`-:math:`\Delta\nu` relation calibrated by Corsaro et al. (2012)b. The default value is 0.601.

* ``epsilon_slope``
The slope of the :math:`\epsilon`-:math:`\Delta\nu` relation calibrated by Corsaro et al. (2012)b. The default value is 0.632.
