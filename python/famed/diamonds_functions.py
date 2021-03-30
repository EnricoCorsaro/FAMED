import numpy as np
import os
import subprocess

from . import asteroseismic_functions as astero

__all__ = ['write_uniform_prior',
           'write_data_range',
           'get_background',
           'set_peakbagging',
           'run_peakbagging',
           'run_asymptotic']

def write_uniform_prior(filename, boundaries):
    """
    Write a uniform prior hyper-parameter ASCII file to be used by DIAMONDS.

    Parameters
    ----------
    filename : str
        Full path for output file.
    boundaries : numpy ndarray, list
        A Nx2 set of lower and upper boundaries to write to file.
    """
    header="""
# Hyper parameters used for setting up uniform priors.
# Each line corresponds to a different free parameter (coordinate).
# Column #1: Minima (lower boundaries)
# Column #2: Maxima (upper boundaries) \n
"""
    np.savetxt(filename, boundaries, fmt='%.3f',header=header)

def write_data_range(filename, boundaries):
    """
    Write data boundaries for a fit to input ASCII file to be used by DIAMONDS.

    Parameters
    ----------
    filename : str
        Full path for output file.
    boundaries : array-like
        Lower and upper data range values to write to file.
    """
    np.savetxt(filename, [boundaries],fmt='%.5f')

def get_background(catalog_id, star_id, background_results_dir, background_run_number='00'):
    """
    Retrieve the background model if previously fitted using DIAMONDS.

    This function retrieves the background name adopted during the background 
    fit with DIAMONDS, as well as its resulting fitted parameters. Currently, 
    only background fits with DIAMONDS are supported.
    
    Parameters
    ---------
    catalog_id : str
        Catalogue ID of the star (e.g. 'KIC' for Kepler ot 'TIC' for TESS).
    star_id : str
        ID of the star as a string (e.g. '0012008916' or '7037405').
    background_results_dir : str
        Full path to location of the DIAMONDS Background results for this star.
    background_run_number : str, default: '00'
        Specific subdirectory of the background model to load. 

    Returns
    -------
    bgp : dict
        A dictionary ``{'parameters': bg_par, 'name': bg_name}``, where
        *parameters* is the array of values for the background parameters and 
        *name* is the type the background model (e.g. 'TwoHarvey' or 
        'ThreeHarveyColor')
    """
    # Read background model fitted parameters
    bg_par = np.loadtxt(background_results_dir/(catalog_id + star_id)/str(background_run_number).zfill(2)/'background_parameterSummary.txt',usecols=(0,))

    # Read background model name
    if os.path.isfile(background_results_dir/(catalog_id + star_id)/str(background_run_number).zfill(2)/'background_computationParameters.txt'):
        config = np.loadtxt(background_results_dir/(catalog_id + star_id)/str(background_run_number).zfill(2)/'background_computationParameters.txt',dtype='str')
    else:
        config = np.loadtxt(background_results_dir/(catalog_id + star_id)/str(background_run_number).zfill(2)/'background_configuringParameters.txt',dtype='str')
    bg_name = config[-2]

    # Apply the following in case the background name is not listed among the
    # computation parameters of the Background fit (e.g. an older version of
    # DIAMONDS was used for the Background fit).
    bg_name_list = np.array(['Flat','Original','OneHarvey','OneHarveyColor','TwoHarvey','TwoHarveyColor','ThreeHarvey','ThreeHarveyColor'])
    bg_name_npar = np.array([1,3,3,5,5,7,7,9])
    bg_name_flag = np.array([1,0,0,1,0,1,0,1])

    tmp_match = np.where(bg_name_list == bg_name)[0]
    if len(tmp_match) < 1:
        n_bg_par = len(bg_par)-3
        tmp_match2 = np.where((bg_name_npar == n_bg_par) & (bg_name_flag == 1))[0]
        bg_name = bg_name_list[tmp_match2[0]]

    bgp = {'parameters':    bg_par,
           'name':          bg_name}

    return bgp


def set_peakbagging(catalog_id, star_id, bgp, diamonds_path, dnu_cl=9, dnu_tip=3.2, n_dnu_envelope=4.5, n_sigma_envelope=4.5, n_sigma_envelope_cl=2.5, n_sigma_envelope_tip=1.2, numax_threshold=300, numax_coeff_low=0.267, numax_coeff_high=0.22, numax_exponent_low=0.76, numax_exponent_high=0.797):
    """
    Verify or create the proper folder structure for DIAMONDS PeakBagging.

    Verifies that the peakbagging folder is properly set up for running DIAMONDS
    PeakBagging fits. If the setup is missing, this function will prepare it for
    you, based on the computation of a background model with DIAMONDS as a 
    previous step. Currently, only background fits with DIAMONDS are supported.

    Parameters
    ----------
    catalog_id : str
        Catalogue ID of the star (e.g. 'KIC' for Kepler).
    star_id : str
        ID of the star as a string (e.g. '0012008916' or '7037405').
    bgp : dict
        A dictionary ``{'parameters': bg_par, 'name': bg_name}``, where
        *parameters* is the array of values for the background parameter and
        *name* is the type the background model (e.g. 'TwoHarvey' or 
        'ThreeHarveyColor')
    diamonds_path : str
        Path to the local DIAMONDS installation folder that contains the 
        PeakBagging and Background packages.    
    dnu_cl : float, default: 9
        Threshold value of `dnu` that identifies that stars are in the clump.
    dnu_tip : float, default: 3.2
        Threshold value of `dnu` that separates clump stars from evolved RGB and
        early AGB stars.
    n_dnu_envelope : float, default: 4.5
        The maximum number of times `dnu` to set the envelope width.
    n_sigma_envelope : float, defualt: 4.5
        The number of times the standard deviation of the Gaussian mode envelope
        to set the total frequency range on either side of `numax`.
    n_sigma_envelope_cl : float, default: 2.5
        The number of times the standard deviation of the Gaussian mode envelope
        to set the total frequency range on either side of `numax` for stars 
        identified as clump stars.
    n_sigma_envelope_tip: float, default: 1.2
        The number of times the standard deviation of the Gaussian mode envelope
        to set the total frequency range on either side of `numax` for stars 
        identified as being near the tip of the RGB.
    numax_threshold : float, default: 300
        Threshold in `numax` which roughly separates RGB stars from subgiant
        and main sequence stars. 
    numax_coeff_low : float, default: 0.267 
        Power law coefficient used for stars below *numax_threshold*.
    numax_coeff_high : float, default: 0.22 
        Power law coefficient used for stars above *numax_threshold*.
    numax_exponent_low : float, default: 0.76
        Power law exponent used for stars below *numax_threshold*.
    numax_exponent_high : float, default0.797
        Power law exponent used for stars above *numax_threshold*.
    """
    peakbagging_results_dir = diamonds_path/'PeakBagging'/'results'
    peakbagging_data_dir  = diamonds_path/'PeakBagging'/'data'
    background_results_dir  = diamonds_path/'Background'/'results'
    background_data_dir  = diamonds_path/'Background'/'data'
    bg_par = bgp['parameters']

    if not os.path.isfile(peakbagging_data_dir/(catalog_id + star_id + '.txt')):
        freq,psd = np.loadtxt(background_data_dir/(catalog_id + star_id + '.txt')).T

        # Trim the global PSD of the star, used for the background fit, to a
        # narrower frequency region centered around nuMax
        numax = bg_par[len(bg_par)-2]
        sig_env = bg_par[len(bg_par)-1]
        dnu = astero.compute_scaling_dnu(numax, numax_threshold, numax_coeff_low, numax_coeff_high, numax_exponent_low, numax_exponent_high)

        if dnu <= dnu_cl:
            if dnu <= dnu_tip:
                width_factor = n_sigma_envelope_tip
            else:
                width_factor = n_sigma_envelope_cl
        else:
            width_factor = n_sigma_envelope

        if n_dnu_envelope*dnu < sig_env*width_factor:
            env_width = sig_env*width_factor
        else:
            env_width = n_dnu_envelope*dnu

        lower_bound = numax - env_width
        upper_bound = numax + env_width

        tmp_clipping = np.where((freq >= lower_bound) & (freq <= upper_bound))[0]

        freq_pb = freq[tmp_clipping]
        psd_pb = psd[tmp_clipping]

        filename = peakbagging_data_dir/(catalog_id + star_id + '.txt')
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        with open(filename, 'w') as f:
            for i in range(0,len(freq_pb)):
                f.write('{}    {}\n'.format(freq_pb[i],psd_pb[i]))

    if not os.path.isfile(peakbagging_results_dir/(catalog_id + star_id)/'backgroundParameters.txt'):
        filename = peakbagging_results_dir/(catalog_id + star_id)/'backgroundParameters.txt'
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        with open(filename, 'w') as f:
            f.write("""
# Configuring parameters for background model of {}{}
# These parameters are the median values of the result derived
# from the Background code based on DIAMONDS.\n
            """.format (catalog_id,star_id))
            for par in bg_par[:-3]:
                f.write('{}\n'.format(par))

    if not os.path.isfile(peakbagging_results_dir/(catalog_id + star_id)/'gaussianEnvelopeParameters.txt'):
        with open(peakbagging_results_dir/(catalog_id + star_id)/'gaussianEnvelopeParameters.txt','w') as f:
            f.write("""
# Parameters of the Gaussian envelope fit for background model of {}{}
# Row #1: height
# Row #2: numax
# Row #3: sigma_env\n
            """.format(catalog_id,star_id))
            for par in bg_par[-3:]:
                f.write('{}\n'.format(par))

    if not os.path.isfile(peakbagging_results_dir/(catalog_id + star_id)/'NyquistFrequency.txt'):
        subprocess.call(('cp {} {}'.format(background_results_dir/(catalog_id + star_id)/'NyquistFrequency.txt',peakbagging_results_dir/(catalog_id + star_id))),shell=True)


def run_peakbagging(catalog_id, star_id, parameters, flag_peaktest, flag_asymptotic, flag_bglevel, dp, diamonds_path, n_threads=12, prior_filename='prior_hyperParameters', merge=False):
    """
    Execute the DIAMONDS PeakBagging code. 

    This function executes the peakbagging code based on DIAMONDS for different
    kind of fits, as determined by the input flags. These flags differentiate 
    between uni-modal, multi-modal, peak testing, and sliding pattern fits.
    
    Parameters
    ----------
    catalog_id : str
        Catalogue ID of the star (e.g. 'KIC' for Kepler).
    star_id : str
        ID of the star as a string (e.g. '0012008916' or '7037405').
    parameters : dict
        A dictionary with the parameters of the specific run.
    flag_peaktest: int
        Flag used to indicate running as a peak test.
    flag_asymptotic : int
        Flag used to indicate running a sliding pattern fit.
    flag_bglevel : int
        Flag used to write the background fit level out to an ASCII file.
    dp : dict
        Dictionary containing the DIAMONDS configuring parameters. See DIAMONDS
        documentation for full description of parameters.
    diamonds_path : str
        Path to the local DIAMONDS installation folder that contains the 
        PeakBagging package.
    n_threads : int, default: 12
        Number of available threads for parallel processing of multiple runs.
    prior_filename : str, default: 'prior_hyperParameters'
        Root of the prior hyperparameters filename.
    merge : bool, default: False
        Merge 

    Returns
    -------
    flag_run_completed : bool
        Returns ``True`` if PeakBagging ran successfully

    """
    # Setup DIAMONDS configuring parameters
    peakbagging_results_dir = diamonds_path/'PeakBagging'/'results'
    star_dir = peakbagging_results_dir/(catalog_id + star_id)
    nsmc_filename = star_dir/parameters['subdir']/'NSMC_configuringParameters.txt'
    xmeans_filename = star_dir/parameters['subdir']/'Xmeans_configuringParameters.txt'
    fwhm_minimum = np.min(parameters['fwhm'])

    
    # Uni-modal peak bagging or peak testing for FWHM fit
    if ((fwhm_minimum < 0) and (flag_peaktest == 0)) or (flag_peaktest == 2):
        nsmc_keys = ['n_live','n_live_end','max_draw_attempts','n_initial_it','n_it_same_clust','enlarg_fraction','shrinking_rate','termination_factor','max_nested_it']
        nsmc_parameters = [dp.get(key) for key in nsmc_keys ]
        xmeans_parameters = [dp['min_ncluster'], dp['max_ncluster']]

    # Peak testing
    if flag_peaktest == 1:
        nsmc_keys = ['n_live_test','n_live_end_test','max_draw_attempts','n_initial_it','n_it_same_clust','enlarg_fraction','shrinking_rate','termination_factor_test','max_nested_it']
        nsmc_parameters = [dp.get(key) for key in nsmc_keys ]
        xmeans_parameters = [dp['min_ncluster'], dp['max_ncluster']]
        
    # Multi-modal peak bagging
    if (fwhm_minimum >= 0) & (flag_peaktest == 0) & (flag_asymptotic == 0):
        nsmc_keys = ['n_live','n_live_end','max_draw_attempts','n_initial_it','n_it_same_clust','enlarg_fraction','shrinking_rate','termination_factor','max_nested_it']
        nsmc_parameters = [dp.get(key) for key in nsmc_keys ]

        if parameters['duplet']:
            xmeans_parameters = [dp['min_ncluster'], dp['max_ncluster_duplet']]
        else:
            xmeans_parameters = [dp['min_ncluster'], dp['max_ncluster']]

    # Sliding pattern fit
    if flag_asymptotic == 1:
        nsmc_keys = ['n_live','n_live_end','max_draw_attempts','n_initial_it','n_it_same_clust','enlarg_fraction','shrinking_rate','termination_factor','max_nested_it']
        nsmc_parameters = [dp.get(key) for key in nsmc_keys ]
        xmeans_parameters = [dp['min_ncluster'], dp['max_ncluster']]

    with open(nsmc_filename, 'w') as f:
        f.write("\n".join(map(str,nsmc_parameters)))
    with open(xmeans_filename, 'w') as f:
        f.write("\n".join(map(str,xmeans_parameters)))


    # Create subdirectories for each run, if not already present
    n_runs = len(parameters['run'])
    for i in range(0, n_runs):
        directory = star_dir/parameters['subdir']/str(parameters['run'][i])
        if not os.path.exists(directory):
            os.makedirs(directory)

    cwd = os.getcwd()
    os.chdir(diamonds_path/'PeakBagging'/'build')
    
    if n_runs > 1:
        n_chunks = np.int(np.ceil(n_runs*1.0/n_threads))
        start_indices = np.zeros(n_chunks,dtype='int')
        end_indices = np.zeros(n_chunks,dtype='int')
        start_index = 0

        if n_runs >= n_threads:
            # Estimate how many chunks to divide the total run set into
            end_index = n_threads
            start_indices[0] = start_index
            end_indices[0] = end_index

            for j in range(1, n_chunks):
                if j == n_chunks-1:
                    start_index = end_index
                    end_index = n_runs 
                    start_indices[j] = start_index
                    end_indices[j] = end_index
                else:
                    start_index = end_index
                    end_index = end_index + n_threads
                    start_indices[j] = start_index
                    end_indices[j] = end_index
        else:
            start_indices[0] = start_index
            end_indices[0] = n_runs - 1

        for j in range(0, n_chunks):
            start_index = start_indices[j]
            end_index = end_indices[j]

            if merge:
                filename_run =  star_dir/parameters['subdir']/(catalog_id + star_id + '_run_indices.txt')
                with open(filename_run,'w') as f:
                    for run in parameters['run'][start_index:end_index+1]:
                        f.write('{}\n'.format(run))

                filename_bg =  star_dir/parameters['subdir']/(catalog_id + star_id + '_bg_names.txt')
                with open(filename_bg,'w') as f:
                    for bg in parameters['background'][start_index:end_index+1]:
                        f.write('{}\n'.format(bg))

                filename_prior = star_dir/parameters['subdir']/(catalog_id + star_id + '_prior_filenames.txt')
                with open(filename_prior,'w') as f:
                    for i in range(0,end_index-start_index+1):
                        f.write('{}\n'.format(prior_filename))

                filename_fwhm = star_dir/parameters['subdir']/(catalog_id + star_id + '_linewidths.txt')
                with open(filename_fwhm,'w') as f:
                    for fwhm in parameters['fwhm'][start_index:end_index+1]:
                        f.write('{}\n'.format(fwhm))

                output_err_filename = star_dir/parameters['subdir']/(catalog_id + star_id + '_' + parameters.subdir + '_parallel.out')
                command = 'parallel ./peakbagging ::: {} ::: {} ::: {} :::: {} ::::+ {} ::::+ {} ::::+ {} ::: {} ::: {} ::: {}'.format(catalog_id,star_id,parameters['subdir'],filename_run,filename_bg,filename_prior,filename_fwhm,flag_peaktest,flag_asymptotic,flag_bglevel,output_err_filename)
                subprocess.run(command,shell=True,check=True,stdout=open(output_err_filename,'w'),stderr=subprocess.STDOUT)
                os.remove(filename_run)
                os.remove(filename_bg)
                os.remove(filename_prior)
                os.remove(filename_fwhm)
            else:
                with open(star_dir/parameters['subdir']/(parameters['filename_run']+'.txt'),'w') as f:
                    for run in parameters['run'][start_index:end_index+1]:
                        f.write('{}\n'.format(run))

                output_err_filename =  star_dir/parameters['subdir']/(catalog_id + star_id + '_' + parameters['subdir'] + '_'+ parameters['filename_run'] + '_parallel.out')
                command='parallel ./peakbagging ::: {} ::: {} ::: {} :::: {}.txt ::: {} ::: {} ::: {} ::: {} ::: {} ::: {}'.format(catalog_id,star_id,parameters['subdir'],star_dir/parameters['subdir']/parameters['filename_run'],parameters['background'],prior_filename,parameters['fwhm'],flag_peaktest,flag_asymptotic,flag_bglevel,output_err_filename)
                process = subprocess.run(command,shell=True,check=True,stdout=open(output_err_filename,'w'),stderr=subprocess.STDOUT)

    else:
        output_err_filename = star_dir/parameters['subdir']/(catalog_id + star_id + '_' + parameters['subdir'] + '_' + str(parameters['run'][0]) + '.out')
        command = './peakbagging {} {} {} {} {} {} {} {} {} {}'.format(catalog_id,star_id,parameters['subdir'],parameters['run'][0],parameters['background'],prior_filename,parameters['fwhm'],flag_peaktest,flag_asymptotic,flag_bglevel)
        print(command)
        subprocess.run(command,shell=True,check=True,stdout=open(output_err_filename,'w'),stderr=subprocess.STDOUT)
        
    os.chdir(cwd)
    flag_computation_completed = np.zeros(len(parameters['run']),dtype='int')

    for i in range(0, len(parameters['run'])):
        if os.path.isfile(star_dir/parameters['subdir']/str(parameters['run'][i])/'peakbagging_parameterSummary.txt'):
            flag_computation_completed[i] = 1

    return flag_computation_completed


def run_asymptotic(catalog_id, star_id, parameters, numax, ell, dp, diamonds_path, prior_filename='prior_hyperParameters'):
    """
    Execute the DIAMONDS Asymptotic code. 

    This function executes the asymptotic code based on DIAMONDS for different 
    asymptotic relations as given by the input angular degree. Radial (l=0) and
    dipole (l=1) modes are supported. 
    
    Parameters
    ----------
    catalog_id : str
        Catalogue ID of the star (e.g. 'KIC' for Kepler).
    star_id : str
        ID of the star as a string (e.g. '0012008916' or '7037405').
    parameters : dict
        A dictionary with the parameters of the specific run.
    numax : float
        The frequency of maximum oscillation.
    ell : int
        The degree of the modes, either 0 or 1.
    dp : dict
        Dictionary containing the DIAMONDS configuring parameters. See DIAMONDS
        documentation for full description of parameters.
    diamonds_path : str
        Path to the local DIAMONDS installation folder that contains the
        PeakBagging and Asymptotic packages.
    prior_filename : str, default: 'prior_hyperParameters'
        Root of the prior hyperparameters filename.
 
    Returns
    -------
    flag_computation_completed : bool
        Returns True if the Asymptotic code ran successfully

    """
    peakbagging_results_dir = diamonds_path/'PeakBagging'/'results'
    asymptotic_path = diamonds_path/'Asymptotic'/'build'
    star_dir = peakbagging_results_dir/(catalog_id+star_id)
    nsmc_filename = star_dir/parameters['subdir']/'NSMC_configuringParameters.txt'
    xmeans_filename = star_dir/parameters['subdir']/'Xmeans_configuringParameters.txt'
    nsmc_keys = ['n_live','n_live_end','max_draw_attempts','n_initial_it','n_it_same_clust','enlarg_fraction','shrinking_rate','termination_factor','max_nested_it']
    nsmc_parameters = [dp.get(key) for key in nsmc_keys ]
    xmeans_parameters = [dp['min_ncluster'], dp['max_ncluster']]

    if not os.path.isfile(nsmc_filename):
        with open(nsmc_filename,'w') as f:
            for par in nsmc_parameters:
                f.write(str(par)+'/n')

    if not os.path.isfile(xmeans_filename):
        with open(xmeans_filename,'w') as f:
            for par in xmeans_parameters:
                f.write(str(par)+'/n')

    if not os.path.isdir(star_dir/parameters['subdir']/parameters['run']):
        subprocess.call(['mkdir',star_dir/parameters['subdir']/parameters['run']])
    if not os.path.isdir(star_dir/parameters['subdir']/'data'):
        subprocess.call(['mkdir',star_dir/parameters['subdir']/'data'])

    # Add in logger here to say path is changing and what it is changed to etc.
    cwd = os.getcwd()
    os.chdir(asymptotic_path)
    output_err_filename = star_dir/parameters['subdir']/(catalog_id + star_id + '_' + parameters['subdir'] + '_' + parameters['run'] + '.out')



    subprocess.call(['./asymptotic', catalog_id, star_id, parameters['subdir'], 
                     parameters['run'], prior_filename, str(numax), str(ell)],
                    stdout=open(output_err_filename,'w'), stderr=subprocess.STDOUT)

    #subprocess.call(['rm', output_err_filename])
    os.chdir(cwd)
    
    flag_computation_completed = False
    if os.path.isfile(star_dir/parameters['subdir']/parameters['run']/'peakbagging_parameterSummary.txt'): 
        flag_computation_completed = True

    return flag_computation_completed
