pro find_islands_global,catalog_id,star_id,threshold_asef,tolerance,teff,force=force,external=external
; -------------------------------------------------------------------------------------------------------------------
; Author:      Enrico Corsaro
; e-mail:      enrico.corsaro@inaf.it
; Date:        January-March 2019
; Place:       Catania, Italy
; Purpose:     This routine processes the global multi-modal sampling obtained with DIAMONDS in order to identify 
;              meaningful radial orders and a first proxy for the location of dipole and radial mode frequencies.
;              This is also used to obtain a rather accurate estimate of DeltaNu from ACF of the multi-moidal 
;              sampling obtained with DIAMONDS. For evolved stars, an additional fit to locate the central
;              radial mode is performed, here also performed with DIAMONDS to gain computational speed. This is 
;              because the epsilon term from asymptotic relation may change depending on the evolutionary stage of 
;              the starmode is performed, as epsilon from asymptotic relation may change depending on the 
;              evolutionary stage of the star.
; -------------------------------------------------------------------------------------------------------------------

COMMON CONFIG,cp
COMMON STAR,info
COMMON GRAPHIC,pp,lp,sp,lpe
COMMON DIAMONDS,dp

modality = 'GLOBAL'

if ~keyword_set(external) then begin
    setup_computation
endif

star_dir = info.peakbagging_results_dir + catalog_id + star_id + '/'
peakbagging_filename_global = info.peakbagging_results_dir + catalog_id + star_id + '/' + info.summary_subdir + '/' $
                              + catalog_id + star_id + info.peakbagging_filename_label + $
                                info.isla_subdir + '_' + info.global_subdir + '_' + modality + '.txt'

; Copy FAMED configuring parameters used in this run

famed_configuring_parameters_filename_copy = info.peakbagging_results_dir + catalog_id + star_id + '/' + info.summary_subdir + '/' $
    + info.local_configuring_parameters_filename + catalog_id + star_id + $
                                '_' + info.isla_subdir + '_' + info.global_subdir + '_' + modality + '.txt'

spawn,'cp ' + info.configuring_parameters_filename + ' ' + famed_configuring_parameters_filename_copy

; Read sampled frequency from DIAMONDS multi-modal fit

readcol,star_dir + info.isla_subdir + '/' + info.global_subdir + '/peakbagging_parameter000.txt',par0,format='D',/silent

; Read prior height (upper limit)

readcol,star_dir + info.isla_subdir + '/' + info.global_subdir + '/peakbagging_hyperParametersUniform.txt',prior_down,prior_up,format='D,D',/silent
left_bound = prior_down(0)
right_bound = prior_up(0)
upper_height = prior_up(1)

; Read input linewidth and background name of the run

readcol,star_dir + info.isla_subdir + '/' + info.global_subdir + '/peakbagging_computationParameters.txt',config,format='A',/silent
fit_linewidth = config(n_elements(config)-1)
background_name = config(n_elements(config)-2)

; Read star PSD

readcol,info.peakbagging_data_dir + catalog_id + star_id + '.txt',freq,psd,format='D,D',/silent
freqbin = freq(1)-freq(0)

; Read Nyquist frequency

readcol,star_dir + 'NyquistFrequency.txt',nyq,format='D',/silent
nyq = nyq[0]

if info.print_on_screen eq 1 then begin
    print,'---------------------------------------------------'
    print,' Parameter range (microHz): [',strcompress(string(min(par0)),/remove_all),', ', strcompress(string(max(par0)),/remove_all),']'
    print,'---------------------------------------------------'
    print,''
    print,'-------------------------------------------------'
    print,' PSD frequency range (microHz): [',strcompress(string(min(freq)),/remove_all),', ', strcompress(string(max(freq)),/remove_all),']'
    print,'-------------------------------------------------'
    print,''
endif


; -------------------------------------------------------------------------------------------------------------------
; Set up plotting panels for either X term or EPS
; -------------------------------------------------------------------------------------------------------------------

if (info.save_eps ne 0) or (info.save_png ne 0) then begin
    !p.multi = pp.global.pmulti
    position_sampling = pp.global.position_sampling             ; Sampling
    position_asef = pp.global.position_asef                     ; ASEF
    position_psd = pp.global.position_psd                       ; PSD
    
    if info.save_eps eq 0 then begin
        set_plot,'x'
        window,1,xs=pp.xsize,ys=pp.ysize,xpos=pp.xpos,ypos=pp.ypos,title=catalog_id + star_id + ' RUN: ' + info.global_subdir
        device,decomposed=0,true_color=8,retain=2
    endif else begin
        set_plot,'PS'
        filename_star = star_dir + info.figs_subdir + '/' + catalog_id + star_id + '_' + info.isla_subdir + '_' + info.global_subdir + '_' + modality + '.eps'
        device,filename=filename_star,xs=pp.xsize,ys=pp.ysize,/encapsulated,/color,bits=8
    endelse
endif


; -------------------------------------------------------------------------------------------------------------------
; Load the nuMax information
; -------------------------------------------------------------------------------------------------------------------

readcol,star_dir + 'gaussianEnvelopeParameters.txt',gauss_par,format='D',/silent
numax = gauss_par(1)
sig_env = gauss_par(2)
scaling_dnu = compute_scaling_dnu(numax)


; -------------------------------------------------------------------------------------------------------------------
; Set a binwidth proportional to the expected separation between adjacent peaks (normally taken as d02/2 if MS
; or smaller if RG star).
; -------------------------------------------------------------------------------------------------------------------

par_range = max(par0) - min(par0)
min_separation = get_minimum_freq_separation(scaling_dnu,1)
binwidth = 1.d0*min_separation/cp.min_n_bins
n_bins = ceil(par_range/binwidth)
nest_iter = findgen(n_elements(par0))


; -------------------------------------------------------------------------------------------------------------------
; Plot the nested sampling evolution of the frequency parameter as obtained by Diamonds.
; -------------------------------------------------------------------------------------------------------------------

if (info.save_eps ne 0) or (info.save_png ne 0) then begin
    plot_sampling,par0,position_sampling
    arrow,n_elements(par0)-1,numax,n_elements(par0)*0.9,numax,/data,thick=4,/solid,hsize=pp.numax_arrow_hsize,color=pp.numax_chunk_arrow_color
endif


; -------------------------------------------------------------------------------------------------------------------
; Compute an Average Shifted Envelope Function (ASEF) of the distribution of nested iterations
; -------------------------------------------------------------------------------------------------------------------
; Compute the ASEF with high resolution in order to determine DeltaNu with high accuracy

n_bins_acf = round(n_elements(nest_iter)/cp.n_bins_acf_scaling)
ash = compute_asef(par0,nest_iter,n_bins_acf)
par_hist = ash.x
asef_hist = ash.y

if (info.save_eps ne 0) or (info.save_png ne 0) then begin
    acf_dnu = compute_acf_dnu(scaling_dnu,par_hist,asef_hist,/plot)
endif else begin
    acf_dnu = compute_acf_dnu(scaling_dnu,par_hist,asef_hist)
endelse

; Compute the ASEF with input resolution for extracting the local maxima

ash = compute_asef(par0,nest_iter,n_bins)
par_hist = ash.x
asef_hist = ash.y


; -------------------------------------------------------------------------------------------------------------------
; Find an appropriate epsilon based on the input temperature.
; -------------------------------------------------------------------------------------------------------------------

if (info.save_eps ne 0) or (info.save_png ne 0) then begin
    interp_epsi = interpolate_epsilon(teff,acf_dnu,/plot)
endif else begin
    interp_epsi = interpolate_epsilon(teff,acf_dnu)
endelse


; -------------------------------------------------------------------------------------------------------------------
; Find the local maxima using a hill-climbing algorithm on the envelope function
; -------------------------------------------------------------------------------------------------------------------

; Give as input the minimum number of adjacent bins required to consider two local maxima separated

threshold = threshold_asef*max(asef_hist)
index_maximum = hill_climbing(par_hist,asef_hist,threshold,cp.min_bin_separation)
maximum = par_hist(index_maximum)
asef_maximum = asef_hist(index_maximum)


; -------------------------------------------------------------------------------------------------------------------
; Identify the divisions among the local maxima. These will be used to select chunks of parameter range that locate
; an optimal frequency region to perform a fit around the peak.
; Also compute ranges for each maximum such that an ASEF peak is considered up to the points of its tails, i.e. 
; the last ones in the decreasing phase off the peak. Ranges are used to compute the actual frequency estimates
; from the sampling.
; -------------------------------------------------------------------------------------------------------------------

n_maxima = n_elements(maximum)
color_step = floor(255./n_maxima)
color_start = fix(abs(255.-color_step*n_maxima))
range_maximum = dblarr(2,n_maxima)
divisions_maximum = dblarr(2,n_maxima)
get_range_divisions,par_hist,asef_hist,index_maximum,range_maximum,divisions_maximum


; -------------------------------------------------------------------------------------------------------------------
; Based on the identified frequency ranges, estimate the oscillation frequencies using the sampling information 
; obtained by DIAMONDS.
; -------------------------------------------------------------------------------------------------------------------

freq1 = dblarr(n_maxima)
freq_sig1 = dblarr(n_maxima)
sampling_counts = fltarr(n_maxima)
spsd_maximum = fltarr(n_maxima)
n_freq = n_maxima

; Compute a smoothed PSD by some average FWHM to evaluate its values at each local maxima of the ASEF.
; In case of RG, adopt a finer smoothing window to overcome the problem of very narrow mixed modes.

avg_fwhm = mean(get_linewidth(maximum,teff,numax))
if acf_dnu le cp.dnu_tip then begin
    avg_fwhm = fit_linewidth*cp.smoothing_fwhm_factor_rg
endif

smth_bins = avg_fwhm/freqbin
spsd = smooth(psd,smth_bins,/edge_truncate)
tmp_good = where(freq le max(par_hist) and freq ge min(par_hist))
psd = temporary(psd(tmp_good))
freq = temporary(freq(tmp_good))
spsd = temporary(spsd(tmp_good))

; Load the background level of the star

readcol,star_dir + 'backgroundLevel.txt',bg_level,format='x,D',/silent,comment='#'
bg_level_local = bg_level(tmp_good)

; Evaluate SNR of the dataset in the case of global modality

snr = max(spsd)/mean(bg_level_local)
sampled_estimates = evaluate_sampling_frequencies(par0,par_hist,freq,spsd,maximum,range_maximum)
freq1 = sampled_estimates.freq1
freq_sig1 = sampled_estimates.freq_sig1
sampling_counts = sampled_estimates.sampling_counts
spsd_maximum = sampled_estimates.spsd_maximum

if (info.save_eps ne 0) or (info.save_png ne 0) then begin 
    plot_asef,par_hist,asef_hist,maximum,range_maximum,threshold,position_asef
endif

; Define some weights useful for identification of proper frequency peaks during mode identification

sampling_weights = sampling_counts/total(sampling_counts)
spsd_weights = spsd_maximum/total(spsd_maximum)
asef_weights = asef_maximum/total(asef_maximum)

if info.save_complete_lists eq 1 then begin
    ; Save the total list of frequencies, uncertainties, ASEF maxima and sampling counts for reference.
    
    get_lun,lun1
    openw,lun1,info.peakbagging_results_dir + catalog_id + star_id + '/' + info.summary_subdir + '/' $
               + catalog_id + star_id + info.peakbagging_filename_label + 'global.all.txt'
    printf,lun1,'# Frequency (microHz), 1-sigma frequency (microHz), ASEF maximum (iterations), sampling counts',format='(A0)'

    for i=0,n_freq-1 do begin
        printf,lun1,freq1(i),freq_sig1(i),asef_maximum(i),sampling_counts(i),format='(F0.5,F15.5,I10,I15)'
    endfor
    
    free_lun,lun1
endif

if info.print_on_screen eq 1 then begin
    print,' Total number of local maxima found: ',strcompress(string(n_maxima,format='(I)'),/remove_all)
endif


; -------------------------------------------------------------------------------------------------------------------
; Sliding Pattern Analysis
; Obtain epsilon from an échelle value of the l=0 ridge échelle frequency position.
; -------------------------------------------------------------------------------------------------------------------

fit_dnu = acf_dnu
central_freq = numax
run_subdir = 'sliding'
run_names = run_subdir + strcompress(string(indgen(cp.n_sliding_test)),/remove_all)
flag_evolved_star = 0
flag_depressed_dipole = 0

; Obtain asymptotic parameters for evolved stars

ap = get_asymptotic_parameters(numax,acf_dnu,teff)


; Start by verifying whether the star has depressed dipole modes

asef_threshold = (dp.isla.max_nested_it + dp.isla.n_live) * cp.asef_threshold_fraction

if fit_dnu le cp.dnu_threshold then begin
    central_indices = where(freq1 gt central_freq - fit_dnu*cp.n_central_orders_side/2. and freq1 lt central_freq + fit_dnu*cp.n_central_orders_side/2.)
endif else begin
    central_indices = where(freq1 gt central_freq - fit_dnu*cp.n_central_orders_side and freq1 lt central_freq + fit_dnu*cp.n_central_orders_side)
endelse

depressed_indices = where(asef_maximum(central_indices) lt asef_threshold)
n_central_freq = n_elements(central_indices)
n_depressed_freq = n_elements(depressed_indices)

; Set up priors for the sliding-pattern fit

if n_depressed_freq ge n_central_freq*cp.depressed_dipole_fraction then begin
    flag_depressed_dipole = 1
    
    if info.print_on_screen eq 1 then begin
        print,''
        print,' The star is likely to have depressed dipole modes.'
        print,''
    endif
endif

if fit_dnu le cp.dnu_threshold then begin
    flag_evolved_star = 1
    dipole_radial_height_ratio = cp.dipole_radial_height_ratio_rg
    quadrupole_radial_height_ratio = cp.quadrupole_radial_height_ratio_rg
    dipole_radial_fwhm_ratio = cp.dipole_radial_fwhm_ratio_rg
    n_orders_side_prior = cp.n_orders_side_prior_rg
    n_orders_side_data = n_orders_side_prior

    ; Set the range of the dataset and frequency prior for the sliding pattern fit

    data_freq_boundaries = [central_freq - n_orders_side_data*fit_dnu,central_freq + n_orders_side_data*fit_dnu]
    tmp_central = where(freq le central_freq + n_orders_side_prior*fit_dnu and freq ge central_freq - n_orders_side_prior*fit_dnu)
    spsd_central = spsd(tmp_central)
    freq_prior = [central_freq - n_orders_side_prior*fit_dnu, central_freq + n_orders_side_prior*fit_dnu]

    dnu_prior = [fit_dnu*cp.dnu_prior_lower_fraction,fit_dnu*cp.dnu_prior_upper_fraction]
    d02_prior = [ap.d02,ap.d02]
    
    if fit_dnu lt cp.dnu_tip then begin
        d01_prior = [ap.d01,ap.d01]
    endif else begin
        d01_prior = [99.0,99.0]
    endelse
    
    d13_prior = [99.0,99.0]
    rot_split_prior = [0.0,0.0]
    cosi_prior = [0.0,0.0]

    n_orders_model = round(2*n_orders_side_prior)
    if (n_orders_model mod 2) eq 0 then n_orders_model++
endif else begin
    if flag_depressed_dipole eq 1 then begin
        flag_evolved_star = 1
        dipole_radial_height_ratio = 0.0
        dipole_radial_fwhm_ratio = 1.0

        dnu_prior = [fit_dnu*cp.dnu_prior_lower_fraction,fit_dnu*cp.dnu_prior_upper_fraction]
        d02_prior = [1.5,cp.d02_prior_upper_sg]
        d01_prior = [99.0,99.0]
        d13_prior = [99.0,99.0]
        rot_split_prior = [0.0,0.0]
        cosi_prior = [0.0,0.0]

        quadrupole_radial_height_ratio = cp.quadrupole_radial_height_ratio_sg
        n_orders_side_prior = cp.n_orders_side_prior_sg
        n_orders_side_data = n_orders_side_prior*1.33333

        ; Set the range of the dataset and frequency prior for the sliding pattern fit
        
        right_freq_bound = central_freq + n_orders_side_data*fit_dnu
        left_freq_bound = central_freq - n_orders_side_data*fit_dnu
        data_freq_boundaries = [left_freq_bound,right_freq_bound]
        tmp_central = where(freq le right_freq_bound and freq ge left_freq_bound)
        spsd_central = spsd(tmp_central)
        freq_prior = [left_freq_bound,right_freq_bound]

        n_orders_model = round(4*n_orders_side_prior)
        if (n_orders_model mod 2) eq 0 then n_orders_model++
    endif else begin
        ; Evaluate the maximum spread in % of DeltaNu of each ridge (both l=2,0 and l=1) with respect to their corresponding median
        ; value. If such deviation for any of the two ridges is larger than a given threshold, then consider that mixed modes
        ; are present and classify the star as a subgiant.
        
        central_indices = where(freq1 gt central_freq - fit_dnu*cp.n_central_orders_side and freq1 lt central_freq + fit_dnu*cp.n_central_orders_side and $
                                asef_maximum ge asef_threshold)
        
        freq1_right = freq1(central_indices(1:*:2))     ; Odd frequencies
        freq1_left = freq1(central_indices(0:*:2))      ; Even frequencies

        ; Check that the extracted frequencies are approximately varying in steps of DeltaNu

        flag_double_step = 0

        if n_elements(freq1_right) ge 2 then begin
            freq1_right_diff = freq1_right(1:*) - freq1_right(0:n_elements(freq1_right)-2)
            
            if median(freq1_right_diff) lt fit_dnu*cp.dnu_ridge_threshold then begin
                flag_double_step = 1
            endif
        endif

        if flag_double_step eq 0 then begin
            if n_elements(freq1_left) ge 2 then begin
                freq1_left_diff = freq1_left(1:*) - freq1_left(0:n_elements(freq1_left)-2)

                if median(freq1_left_diff) lt fit_dnu*cp.dnu_ridge_threshold then begin
                    flag_double_step = 1
                endif
            endif
        endif

        ; If a double step is required, then consider three ridges (left, central, and right)

        if flag_double_step eq 1 then begin
            freq1_left = freq1(central_indices(0:*:3))
            freq1_right = freq1(central_indices(2:*:3))

            if info.print_on_screen eq 1 then begin
                print,''
                print,' Using double step to check the frequency ridges.'
                print,''
            endif
        endif

        freq1_left_modulo = freq1_left mod fit_dnu
        freq1_right_modulo = freq1_right mod fit_dnu

        tmp = where(freq1_left_modulo gt fit_dnu/2.)
        if tmp(0) ne -1 then begin
            freq1_left_modulo(tmp) = fit_dnu - freq1_left_modulo(tmp)
        endif

        tmp = where(freq1_right_modulo gt fit_dnu/2.)
        if tmp(0) ne -1 then begin
            freq1_right_modulo(tmp) = fit_dnu - freq1_right_modulo(tmp)
        endif

        freq1_left_modulo_median = median(freq1_left_modulo)
        freq1_right_modulo_median = median(freq1_right_modulo)
        left_modulo_median_index = where(freq1_left_modulo eq freq1_left_modulo_median)
        right_modulo_median_index = where(freq1_right_modulo eq freq1_right_modulo_median)
        freq1_left_median = freq1_left(left_modulo_median_index)
        freq1_right_median = freq1_right(right_modulo_median_index)
        deviation_left = abs(freq1_left_modulo - freq1_left_modulo_median) / fit_dnu*100.
        deviation_right = abs(freq1_right_modulo - freq1_right_modulo_median) / fit_dnu*100.
        max_dev_left = max(deviation_left)
        max_dev_right = max(deviation_right)

        if flag_double_step eq 1 then begin
            freq1_central = freq1(central_indices(1:*:3))
            freq1_central_modulo = freq1_central mod fit_dnu
            
            tmp = where(freq1_central_modulo gt fit_dnu/2.)
            if tmp(0) ne -1 then begin
                freq1_central_modulo(tmp) = fit_dnu - freq1_central_modulo(tmp)
            endif

            freq1_central_modulo_median = median(freq1_central_modulo)
            central_modulo_median_index = where(freq1_central_modulo eq freq1_central_modulo_median)
            freq1_central_median = freq1_central(central_modulo_median_index)

            deviation_central = abs(freq1_central_modulo - freq1_central_modulo_median) / fit_dnu*100.
            max_dev_central = max(deviation_central)
            max_dev = max([max_dev_left,max_dev_central,max_dev_right])
        endif else begin
            max_dev = max([max_dev_left,max_dev_right])
        endelse

        if max_dev ge cp.dnu_echelle_threshold then flag_evolved_star = 1

        if flag_evolved_star eq 1 then begin
            if info.print_on_screen eq 1 then begin
                print,''
                print,' The star likely contains modes that have undergone avoided crossings, so it is classified as a subgiant.'
                print,''
            endif

            ; Now find the ridge with the smallest maximum spread and from its modulo median value find epsilon and the value of radial mode frequency
           
            if flag_double_step eq 1 then begin
                freq1_median = [freq1_left_median,freq1_central_median,freq1_right_median]
                freq1_modulo_median = [freq1_left_modulo_median,freq1_central_modulo_median,freq1_right_modulo_median]
                min_max_dev = min([max_dev_left,max_dev_central,max_dev_right],ridge_index)
            endif else begin
                freq1_median = [freq1_left_median,freq1_right_median]
                freq1_modulo_median = [freq1_left_modulo_median,freq1_right_modulo_median]
                min_max_dev = min([max_dev_left,max_dev_right],ridge_index)
            endelse

            freq_radial_echelle = freq1_median(ridge_index)
            dipole_radial_height_ratio = cp.dipole_radial_height_ratio_sg
            quadrupole_radial_height_ratio = cp.quadrupole_radial_height_ratio_sg
            dipole_radial_fwhm_ratio = cp.dipole_radial_fwhm_ratio_sg
            n_orders_side_prior = cp.n_orders_side_prior_sg
            n_orders_side_data = n_orders_side_prior*1.33333

            ; Set the range of the dataset and frequency prior for the sliding pattern fit
            
            right_freq_bound = freq_radial_echelle + n_orders_side_data*fit_dnu
            left_freq_bound = freq_radial_echelle - n_orders_side_data*fit_dnu
            data_freq_boundaries = [left_freq_bound,right_freq_bound]
            tmp_central = where(freq le right_freq_bound and freq ge left_freq_bound)
            spsd_central = spsd(tmp_central)
            freq_prior = [left_freq_bound,right_freq_bound]

            dnu_prior = [fit_dnu,fit_dnu]
            d02_prior = [cp.d02_prior_lower_sg,cp.d02_prior_upper_sg]
            d01_prior = [99.0,99.0]
            d13_prior = [99.0,99.0]
            rot_split_prior = [0.0,0.0]
            cosi_prior = [0.0,0.0]
        endif else begin
            if info.print_on_screen eq 1 then begin
                print,''
                print,' The star could be a main sequence.'
                print,''
            endif
            
            dipole_radial_height_ratio = cp.dipole_radial_height_ratio_ms
            quadrupole_radial_height_ratio = cp.quadrupole_radial_height_ratio_ms
            dipole_radial_fwhm_ratio = cp.dipole_radial_fwhm_ratio_ms
            n_orders_side_prior = cp.n_orders_side_prior_ms
            n_orders_side_data = n_orders_side_prior

            ; Set the range of the dataset and frequency prior for the sliding pattern fit
            
            data_freq_boundaries = [central_freq - n_orders_side_data*fit_dnu, central_freq + n_orders_side_data*fit_dnu]
            tmp_central = where(freq le central_freq + n_orders_side_prior*fit_dnu and freq ge central_freq - n_orders_side_prior*fit_dnu)
            spsd_central = spsd(tmp_central)
            freq_prior = [central_freq - n_orders_side_prior*fit_dnu, central_freq + n_orders_side_prior*fit_dnu]

            dnu_prior = [fit_dnu*cp.dnu_prior_lower_fraction,fit_dnu*cp.dnu_prior_upper_fraction]
            
            d02_prior_upper_ms = cp.d02_prior_upper_ms
            if cp.d02_prior_upper_ms ge fit_dnu/4. then d02_prior_upper_ms = fit_dnu/4.
            d02_prior = [cp.d02_prior_lower_ms,d02_prior_upper_ms]
            
            d01_prior_upper = cp.d01_prior_upper_ms
            if d01_prior_upper gt fit_dnu/4. then d01_prior_upper = fit_dnu/4.
            d01_prior = [cp.d01_prior_lower_ms,d01_prior_upper]

            d13_prior_upper = d01_prior_upper
            if d01_prior_upper lt (3./8.*fit_dnu) then d13_prior_upper = fit_dnu/4. - d01_prior_upper
            d13_prior = [cp.d13_prior_lower_ms,abs(d13_prior_upper)]

            rot_split_prior = [cp.rot_split_prior_lower_ms,cp.rot_split_prior_upper_ms]
            cosi_prior = [cp.cosi_prior_lower,cp.cosi_prior_upper]
        endelse

        n_orders_model = round(4*n_orders_side_prior)
        if (n_orders_model mod 2) eq 0 then n_orders_model++
    endelse        
endelse
    
flag_repeat_sliding_fit = 1
sliding_iteration = 0

while ((flag_repeat_sliding_fit eq 1) and (sliding_iteration le 1)) do begin
    ; Check if the run already exists

    if (file_test(star_dir + info.as_subdir + '/' + run_subdir + '0/peakbagging_parameter000.txt') eq 0 or keyword_set(force) or sliding_iteration gt 0) then begin
        ; Make sure that the number of orders to compute the sliding pattern model is an odd number
        ; This will make the selected range symmetric with respect to nuMax

        print,' Total number of radial orders in the sliding model: ',strcompress(string(n_orders_model,format='(I)'),/remove_all)

        asymp_param = [n_orders_model,dipole_radial_height_ratio, $
            quadrupole_radial_height_ratio,cp.octupole_radial_height_ratio,dipole_radial_fwhm_ratio]

        fwhm_asymp = get_linewidth(numax,teff,numax)
        fwhm_asymp = temporary(fwhm_asymp[0])
        asymp_filename = 'asymptoticParameters'
        filename = star_dir + asymp_filename + '.txt'
        
        get_lun,lun1
        openw,lun1,filename
        printf,lun1,'# Asymptotic parameters to set up the asymptotic pattern model.',format='(A0)'
        printf,lun1,'# Row 1: Norders (spanning orders range around nuMax for asymptotic pattern)',format='(A0)'
        printf,lun1,'# Row 2: l=1/l=0 height',format='(A0)'
        printf,lun1,'# Row 3: l=2/l=0 height',format='(A0)'
        printf,lun1,'# Row 3: l=3/l=0 height',format='(A0)'
        printf,lun1,'# Row 4: l=1/l=0 fwhm',format='(A0)'
        printf,lun1,asymp_param, format = '(F0.3)'
        free_lun,lun1

        height_prior = [max(spsd_central)*0.1, max(spsd_central)*1.4]
        boundaries = [freq_prior, height_prior, dnu_prior, d02_prior, d01_prior, d13_prior, rot_split_prior, cosi_prior]

        ; Prepare frequency range and prior filenames for each run
        
        prior_filenames = strarr(cp.n_sliding_test)
        data_range_filenames = strarr(cp.n_sliding_test)

        for k=0, cp.n_sliding_test-1 do begin
            data_range_filenames(k) = star_dir + info.as_subdir + '/frequencyRange_' + run_names(k) + '.txt'
            prior_filenames(k) = star_dir + info.as_subdir + '/' + info.prior_filename + '_' + run_names(k) + '.txt'
            write_diamonds_data_range,data_range_filenames(k),data_freq_boundaries
            write_diamonds_prior,prior_filenames(k),boundaries
        endfor

        if info.print_on_screen eq 1 then begin
            print,''
            print,' Performing sliding-pattern fit with DIAMONDS.'
            print,''
        endif
     
        peakbagging_parameters = { subdir:          info.as_subdir,     $
                                   run:             run_names,          $
                                   background:      background_name,    $
                                   fwhm:            fwhm_asymp,         $
                                   filename_run:    run_subdir          $    
                                 }

        flag_computation_completed = run_peakbagging(catalog_id,star_id,peakbagging_parameters,0,1,0)

        if info.save_test_files ne 1 then begin
            for k=0, cp.n_sliding_test-1 do begin
                file_delete,data_range_filenames(k)
                file_delete,prior_filenames(k)
            endfor
        endif
    endif else begin
        if info.print_on_screen eq 1 then begin
            print,' Load information from sliding-pattern fit with DIAMONDS.'
        endif
    endelse

    echelle_epsi_array = fltarr(cp.n_sliding_test)
    radial_freq_reference_array = fltarr(cp.n_sliding_test)
    d01_array = fltarr(cp.n_sliding_test)

    for k=0, cp.n_sliding_test-1 do begin
        ; Read sampled frequency from DIAMONDS multi-modal fit for nu0 central
        
        readcol,star_dir + info.as_subdir + '/' + run_names(k) + '/peakbagging_parameter000.txt',par_nu0,format='D',/silent

        ; Read sampled posterior distribution from DIAMONDS multi-modal fit
        
        readcol,star_dir + info.as_subdir + '/' + run_names(k) + '/peakbagging_posteriorDistribution.txt',post,format='D',/silent
       
        post /= max(post) 
        radial_freq_reference = total(par_nu0*post)/total(post)

        if cp.input_radial_freq_reference gt 0 then begin
            radial_freq_reference = cp.input_radial_freq_reference
        endif

        radial_freq_reference_array(k) = radial_freq_reference
        modulo_reference = radial_freq_reference mod fit_dnu
        echelle_epsi = modulo_reference/fit_dnu
        
        if echelle_epsi lt cp.epsilon_threshold and fit_dnu gt cp.dnu_lower_threshold_epsilon then echelle_epsi += 1
        echelle_epsi_array(k) = echelle_epsi

        ; If the star is flagged as a MS star, then retrieve also the value of the small spacing d01 from the sliding fit

        if flag_evolved_star eq 0 then begin
            readcol,star_dir + info.as_subdir + '/' + run_names(k) + '/peakbagging_parameter004.txt',par_d01,format='D',/silent
            d01_array(k) = total(par_d01*post)/total(post)
        endif
    endfor

    median_echelle_epsi = median(echelle_epsi_array)

    if fit_dnu le cp.dnu_threshold then begin
        epsilon_upper_limit = cp.upper_epsilon_rg_slope * alog(fit_dnu) + cp.upper_epsilon_rg_offset
        epsilon_lower_limit = cp.lower_epsilon_rg_slope * alog(fit_dnu) + cp.lower_epsilon_rg_offset
        
        ; The sliding pattern without including the l=1 mode peak may have failed in providing a reliable epsilon.
        ; If this is the case repeat the fit by including a l=1 mode peak, having a position fixed to the p-mode frequency of the
        ; asymptotic relation.

        if (median_echelle_epsi ge epsilon_upper_limit) or ((median_echelle_epsi le epsilon_lower_limit) and (fit_dnu ge cp.dnu_cl2)) then begin
            if info.print_on_screen eq 1 then begin
                if sliding_iteration gt 0 then begin
                    print,' Repeating the sliding-pattern-fit did not solve the issue. Epsilon is likely to be wrong for this star.'
                endif else begin
                    if median_echelle_epsi ge epsilon_upper_limit then begin
                        print,' Repeating the sliding-pattern fit because epsilon exceeds the upper limit for RGs.'
                    endif else begin
                        print,' Repeating the sliding-pattern fit because epsilon exceeds the lower limit for RGs.'
                    endelse
                endelse
            endif

            d01_prior = [ap.d01,ap.d01]
        endif else begin
            ; Epsilon from the sliding-pattern fit has been validated. Therefore exit the while loop.
            flag_repeat_sliding_fit = 0
        endelse
    endif else begin
        ; There is no need to repeat the sliding-pattern fit because the star is not a red giant.
        flag_repeat_sliding_fit = 0
    endelse

    sliding_iteration++
endwhile

median_index = where(echelle_epsi_array eq median_echelle_epsi)
radial_freq_reference = radial_freq_reference_array(median_index)
radial_freq_reference = temporary(radial_freq_reference(0))

if flag_evolved_star eq 0 then begin
    fit_d01 = median(d01_array)
endif else begin
    if d01_prior(0) ne 99.0 then begin
        fit_d01 = ap.d01
    endif else begin
        fit_d01 = 0.
    endelse
endelse

; Control check on epsilon for late subgiants/early RGB

interp_epsi_flag = 0

if cp.force_epsilon_dnu_value eq 1 then begin
    if fit_dnu le cp.dnu_threshold and fit_dnu ge cp.dnu_cl then begin
        ; If the star is an evolved subgiant/early RGB check that epsilon is in agreement with the epsilon-DeltaNu relation
        
        radial_freq_reference2 = freq1(closest(radial_freq_reference,freq1))
        modeid_sliding = get_modeid(radial_freq_reference,fit_dnu,median_echelle_epsi,0,numax,0)
        modeid_epsi_dnu = get_modeid(radial_freq_reference2,fit_dnu,interp_epsi,0,numax,0)
        diff_radial = abs(radial_freq_reference2 - radial_freq_reference)
        
        if (modeid_sliding.degree ne modeid_epsi_dnu.degree) or (diff_radial ge fit_dnu/4.) then begin
            median_echelle_epsi = interp_epsi
            interp_epsi_flag = 1

            if modeid_epsi_dnu.degree eq 1 then begin
                radial_freq_reference = freq1(closest(radial_freq_reference2+fit_dnu/2.,freq1))

                if info.print_on_screen eq 1 then begin
                    print,' Applying correction to epsilon from epsilon-Dnu relation.'
                    print,' '
                endif 
            endif else begin
                radial_freq_reference = radial_freq_reference2
            endelse
        endif
    endif
endif

; Control check on epsilon for F-type stars

if teff ge cp.teff_sg then begin
    ; Due to the confusion arising from the strong blending of the modes for F-type stars, the sliding pattern fit may be unreliable.
    ; In this case check the obtained mode identification against the one using the epsilon-Teff relation 
    
    radial_freq_reference2 = freq1(closest(radial_freq_reference,freq1))
    modeid_sliding = get_modeid(radial_freq_reference,fit_dnu,median_echelle_epsi,fit_d01,numax,0)
    modeid_teff = get_modeid(radial_freq_reference2,fit_dnu,interp_epsi,fit_d01,numax,0)
    diff_radial = abs(radial_freq_reference2 - radial_freq_reference)

    if (modeid_sliding.degree ne modeid_teff.degree) or (diff_radial ge fit_dnu/4.) then begin
        median_echelle_epsi = interp_epsi
        interp_epsi_flag = 1
        
        if modeid_teff.degree eq 1 then begin
            radial_freq_reference = freq1(closest(radial_freq_reference2+fit_dnu/2.,freq1))

            if info.print_on_screen eq 1 then begin
                print,' Applying correction to epsilon from epsilon-Teff relation.'
                print,' '
            endif 
        endif else begin
            radial_freq_reference = radial_freq_reference2
        endelse
    endif
endif

if info.print_on_screen eq 1 then begin
    print,' Epsilon values from the sliding fit: ',echelle_epsi_array
    print,' Final epsilon: ',median_echelle_epsi
    print,' The reference radial mode is at: ' + strcompress(string(radial_freq_reference,format='(F0.3)'),/remove_all) + ' muHz'
    print,''
endif

; -------------------------------------------------------------------------------------------------------------------
; Find possible bad frequencies and remove them from the list. Then compute a large frequency separation.
; Finally obtain a simple mode identification (either l=0 or l=1) for the full set of frequencies.
; -------------------------------------------------------------------------------------------------------------------

fit_epsi = median_echelle_epsi 
fit_alpha = 0.
flag_dnu_fit = 1
iterations = 0

while (flag_dnu_fit eq 1) and (iterations lt cp.max_skim_iterations_global) do begin
    n_freq = n_maxima
    order_number = intarr(n_freq)
    angular_degree = intarr(n_freq)
    
    if info.print_on_screen eq 1 then begin
        print,' Iteration: '+strcompress(string(iterations),/remove_all)
    endif
    
    ; Perform a first mode identification based on the ACF value of DeltaNu
    
    for i=0,n_freq-1 do begin
        mode_id = get_modeid(freq1(i),fit_dnu,fit_epsi,fit_d01,numax,fit_alpha)
        enn = mode_id.order
        ell = mode_id.degree
        order_number(i) = enn
        angular_degree(i) = ell
    endfor
    
    ; Using the obtained mode identification, divide the frequency set into dipole and radial modes.
    ; Select only l=0 modes
    
    tmp_radial = where(angular_degree eq 0,complement=tmp_dipole)
    
    if tmp_radial(0) ne -1 then begin
        n_radial = n_elements(tmp_radial)
        freq1_radial = freq1(tmp_radial)
        freq1_radial_org = freq1(tmp_radial)
        freq_sig1_radial = freq_sig1(tmp_radial)
        order_radial = order_number(tmp_radial)
        asef_maximum_radial = asef_maximum(tmp_radial)
    endif else n_radial = 0
   
    ; Select only l=1 modes
    
    if tmp_dipole(0) ne -1 then begin
        n_dipole = n_elements(tmp_dipole)
        freq1_dipole = freq1(tmp_dipole)
        freq_sig1_dipole = freq_sig1(tmp_dipole)
        order_dipole = order_number(tmp_dipole)
        asef_maximum_dipole = asef_maximum(tmp_dipole)
    endif else n_dipole = 0
    
    ; If requested, apply a small correction to the radial mode frequencies that takes into account
    ; the difference between the reference radial mode from the sliding pattern and the one estimated from the ASEF.
    ; If the sliding-pattern fit was unsuccessful, this difference will be zero.

    if cp.correct_radial_frequencies eq 1 then begin
        delta_nu0 = radial_freq_reference - freq1_radial(closest(radial_freq_reference,freq1_radial))
        ap = get_asymptotic_parameters(numax,fit_dnu,teff)
        
        if delta_nu0 lt ap.d02 then begin
            freq1_radial = freq1_radial + delta_nu0
        endif
    endif

    ; Verify the position of each frequency by comparing it to the expected asymptotic value.
    ; Discard those frequencies that deviate from their asymptotic value by more than a given tolerance in DeltaNu.
    
    if n_radial ne 0 then begin
        ; Select only good l=0 frequencies
        
        good_freq_index_radial = assess_freq_asymptotic(freq1_radial,order_radial,0,fit_dnu,fit_epsi,fit_alpha,fit_d01,numax,tolerance)
       
        if info.print_on_screen eq 1 then begin
            if n_elements(good_freq_index_radial) lt n_radial then begin
                print,' Returning fewer frequencies than input for l=0. Input: '+strcompress(string(n_radial),/remove_all)+ ', Actual: '+  $
                strcompress(string(n_elements(good_freq_index_radial)),/remove_all)
            endif
        endif
        
        freq1_radial = temporary(freq1_radial(good_freq_index_radial))
        freq1_radial_org = temporary(freq1_radial_org(good_freq_index_radial))
        freq_sig1_radial = temporary(freq_sig1_radial(good_freq_index_radial))
        order_radial = temporary(order_radial(good_freq_index_radial))
        asef_maximum_radial = temporary(asef_maximum_radial(good_freq_index_radial))
    endif
    
    ; Collect all frequencies after removing the bad ones
    
    n_dipole = n_elements(freq1_dipole)
    n_radial = n_elements(freq1_radial)
    n_freq = n_dipole + n_radial
    
    if n_dipole ne 0 and n_radial ne 0 then begin
        freq1_final = [freq1_radial,freq1_dipole]
        freq_sig1_final = [freq_sig1_radial,freq_sig1_dipole]
        angular_degree = [intarr(n_radial),intarr(n_dipole)+1]
        order_number = [order_radial,order_dipole]
        asef_maximum_final = [asef_maximum_radial,asef_maximum_dipole] 
    endif else if n_dipole eq 0 and n_radial ne 0 then begin
        freq1_final = [freq1_radial]
        freq_sig1_final = [freq_sig1_radial]
        angular_degree = [intarr(n_radial)]
        order_number = [order_radial]
        asef_maximum_final = [asef_maximum_radial] 
    endif
        
    sorted_index = sort(freq1_final)
    freq1_final = temporary(freq1_final(sorted_index))
    freq_sig1_final = temporary(freq_sig1_final(sorted_index))
    angular_degree = temporary(angular_degree(sorted_index))
    order_number = temporary(order_number(sorted_index))
    asef_maximum_final = temporary(asef_maximum_final(sorted_index))
    
    ; Compute optimal value for DeltaNu from individual frequencies only if at least three modes are present.
    ; Otherwise, keep as best Dnu and epsilon, those from the ACF and Dnu-epsilon diagram.

    flag_dnu_fit = 0
    if n_radial ge 2 then begin
        ; Check if the run already exists
        
        run_subdir = 'radial_global'
        filename = star_dir + info.as_subdir + '/' + info.prior_filename + '_' + run_subdir + '.txt'
        dnu_prior = [fit_dnu*cp.dnu_prior_lower_fraction_as,fit_dnu*cp.dnu_prior_upper_fraction_as]

        ; Use epsilon fixed from the sliding pattern fit, if this one was successful
        
        epsi_prior = [fit_epsi,fit_epsi]

        if n_radial eq 2 or n_radial eq 3 then begin
            ; Only DeltaNu can be estimated as a free parameter. Then fix alpha to 0 and epsilon to
            ; its former value.
            
            alpha_prior = [0.0,0.0]
        endif else begin
            alpha_prior = [cp.alpha_prior_lower_as,cp.alpha_prior_upper_as] 
           
            ; If the sliding pattern fit was not successful, i.e. epsilon comes from an interpolated relation, then fit it within a small prior range. 
            
            if interp_epsi_flag eq 1 then begin
                epsi_prior = [fit_epsi*cp.epsi_prior_lower_fraction_as,fit_epsi*cp.epsi_prior_upper_fraction_as]
            endif  
        endelse

        boundaries = [dnu_prior,epsi_prior,alpha_prior]
        write_diamonds_prior,filename,boundaries
       
        data = transpose([[order_radial],[freq1_radial],[freq_sig1_radial]])
        data_filename = star_dir + info.as_subdir + '/data/' + run_subdir + '.txt'
        
        get_lun, lun1
        openw, lun1, data_filename
        printf, lun1, data, format = '(F0.1,F12.4,F12.4)'
        free_lun, lun1
        
        run_parameters = { subdir:     info.as_subdir,     $
                           run:        run_subdir          $
                         }

        flag_computation_completed = run_asymptotic(catalog_id,star_id,run_parameters,numax,0)
       
        spawn,'ls -1 '+ star_dir + info.as_subdir + '/' + run_subdir + '/asymptotic_parameter0*.txt',parameter_filenames
        n_par = n_elements(parameter_filenames) 
        
        readcol,parameter_filenames(0),par,format='D',/silent
        readcol,star_dir + info.as_subdir + '/' + run_subdir + '/asymptotic_posteriorDistribution*.txt',post,format='D',/silent
        
        post /= max(post)
        best_dnu = total(par*post)/total(post)
       
        if n_par eq 1 then begin
            best_epsi = epsi_prior(0)
            best_alpha = alpha_prior(0)
        endif
        
        if n_par eq 2 then begin
            readcol,parameter_filenames(1),par,format='D',/silent
            best_epsi = epsi_prior(0)
            best_alpha = total(par*post)/total(post)
        endif
        
        if n_par eq 3 then begin
            readcol,parameter_filenames(1),par,format='D',/silent
            best_epsi = total(par*post)/total(post)
            readcol,parameter_filenames(2),par,format='D',/silent
            best_alpha = total(par*post)/total(post)
        endif
        
        ; Update local values for asymptotic parameters
       
        fit_dnu = best_dnu
        fit_epsi = best_epsi
        fit_alpha = best_alpha
        flag_dnu_fit = 1
        iterations++
    endif else begin
        ; Not enough radial mode frequencies were found to perform an asymptotic fit
        
        best_dnu = acf_dnu
        best_epsi = fit_epsi
        best_alpha = cp.alpha_radial_universal
    endelse
endwhile

; Save asymptotic radial mode frequencies

freq_radial_asymptotic = best_dnu*(best_epsi + order_radial + best_alpha/2.*(order_radial - numax/best_dnu)^2)

; Select only l=1 modes

tmp_dipole = where(angular_degree eq 1)

if tmp_dipole(0) ne -1 then begin
    n_dipole = n_elements(tmp_dipole)
    freq1_dipole = freq1_final(tmp_dipole)
endif else n_dipole = 0

n_freq = n_elements(freq1_final)

if (info.save_eps ne 0) or (info.save_png ne 0) then begin 
    for i=0, n_freq-1 do begin
        loadct,39,/silent
        xyouts,freq1_final(i)*0.997,max(asef_hist)*1.17,strcompress(string(angular_degree(i)),/remove_all), $
            charsize=pp.label_charsize,charthick=pp.charthick+1,color=pp.label_color
        loadct,0,/silent
        oploterror,[freq1_final(i)],[max(asef_hist)*1.10],[freq_sig1_final(i)],[0],errthick=2,errcolor=100
        loadct,39,/silent
        plotsym,0,pp.symsize2,/fill,color=pp.freq_symcolor
        oplot,[freq1_final(i)],[max(asef_hist)*1.10],psym=8
        loadct,0,/silent
        plotsym,0,pp.symsize3,color=pp.freq_symborder,thick=pp.freq_symthick
        oplot,[freq1_final(i)],[max(asef_hist)*1.10],psym=8
    endfor
endif


; -------------------------------------------------------------------------------------------------------------------
; Calculate radial order positions to divide the PSD into chunks.
; Compute radial mode positions from universal pattern for RG stars (Mosser et al. 2011, A&A, 525, L9).
; This will improve the trimming of chunks with more accurate frequency positions.
; -------------------------------------------------------------------------------------------------------------------
; Check if there is an additional dipolar frequency after the last radial mode.
; If so, add one more chunk to include the last dipole separately.

n_chunks = max(order_number) - min(order_number) + 1
unique_order_number = indgen(n_chunks) + min(order_number)

if max(freq1_final) gt freq1_radial(n_radial-1) then begin
    n_chunks = n_chunks + 1
    unique_order_number = [unique_order_number,max(unique_order_number)+1]
endif

if min(freq1_final) lt freq1_radial(0) then begin
    n_chunks = n_chunks - 1
    unique_order_number = temporary(unique_order_number(1:*))
endif

n_separations = n_chunks + 1
separations = fltarr(n_separations)
freq_asymptotic_separation = best_dnu*(best_epsi + unique_order_number + best_alpha/2.*(unique_order_number - numax/best_dnu)^2)

for i=1, n_separations-2 do begin
    if (best_dnu lt cp.dnu_sg and teff lt cp.teff_sg) then begin
        separations(i) = freq_asymptotic_separation(i-1) + best_dnu*cp.separations_dnu_tolerance_rg
    endif else begin
        separations(i) = freq_asymptotic_separation(i-1) + best_dnu*cp.separations_dnu_tolerance_ms
    endelse
endfor

if separations(1) - best_dnu lt min(par_hist) then begin
    separations(0) = min(par_hist)
endif else begin
    separations(0) = separations(1) - best_dnu
endelse

if separations(n_separations-2) + best_dnu gt max(par_hist) then begin
    separations(n_separations-1) = max(par_hist)
endif else begin
    separations(n_separations-1) = separations(n_separations-2) + best_dnu
endelse


; -------------------------------------------------------------------------------------------------------------------
; Plot the original PSD and a smoothed version of it by the mean linewidth found.
; Overplot an inset of the PSD if a global approach is used, to show the detail of the mode identification.
; -------------------------------------------------------------------------------------------------------------------
if (info.save_eps ne 0) or (info.save_png ne 0) then begin
    parameters = { min_freq:             min(par_hist),          $
                   max_freq:             max(par_hist),          $
                   freq_list:            freq1_final,            $
                   freq_sig_list:        freq_sig1_final,        $
                   freq_radial_asymp:    freq_radial_asymptotic, $
                   freq_radial_org:      freq1_radial_org,       $
                   degree:               angular_degree,         $
                   order:                order_number,           $
                   separations:          separations,            $
                   dnu:                  best_dnu,               $
                   numax:                numax,                  $
                   max_psd:              max(psd)                $
                 }

    plot_psd,freq,psd,spsd,bg_level_local,parameters,position_psd,1
endif

if cp.save_asymptotic_radial eq 1 then begin
    ; If required by the user, replace extracted l=0 frequencies with those from asymptotic pattern.
    ; This will improve characterization of chunk frequencies in evolved stars because the asymptotic 
    ; frequencies for radial modes are in general more accurate than the global frequencies extracted
    ; from the ASEF.
   
    freq1_final = [freq1_dipole,freq_radial_asymptotic]
    tmp_sort_final = sort(freq1_final)
    freq1_final = temporary(freq1_final(tmp_sort_final))
endif


; -------------------------------------------------------------------------------------------------------------------
; Plot logo of the pipeline and summary information of the run
; -------------------------------------------------------------------------------------------------------------------
if info.print_on_screen eq 1 then begin
    print,' ------------------------------------------'
    print,'    Frequency       sig       n       l'
    
    for i=0, n_freq-1 do begin
        print,freq1_final(i),freq_sig1_final(i),order_number(i),angular_degree(i)
    endfor
    
    print,' ------------------------------------------'
    print,' '
endif

if (info.save_eps ne 0) or (info.save_png ne 0) then begin
    parameters = { catalog_id:     catalog_id,             $
                   star_id:        star_id,                $
                   run:            info.global_subdir,     $
                   modality:       modality,               $
                   numax:          numax,                  $
                   teff:           teff,                   $
                   best_dnu:       best_dnu,               $
                   best_alpha:     best_alpha,             $
                   best_epsi:      best_epsi,              $
                   acf_dnu:        acf_dnu,                $
                   snr:            snr,                    $
                   fit_linewidth:  fit_linewidth,          $
                   avg_fwhm:       avg_fwhm,               $
                   upper_height:   upper_height,           $
                   threshold_asef: threshold_asef,         $
                   n_bins:         n_bins,                 $
                   tolerance:      tolerance,              $
                   n_freq:         n_freq,                 $
                   n_chunks:       n_chunks                $
                 }

    plot_summary,parameters,1
endif

if (info.save_eps ne 0) or (info.save_png ne 0) then begin
    if info.save_eps eq 0 then begin
        read_jpeg,info.logo_filename,image
        TV, image, 0.855,0.865,TRUE = 1,/normal
    endif else begin
        xyouts,0.987,0.2,sp.copyright_str+' FAMED',orientation=-90,charsize=lp.summary_charsize,charthick=lp.summary_charthick,/normal,color=250
    endelse
endif

if info.save_eps eq 1 then begin
    device,/close
    spawn,'open ' + filename_star
    set_plot,'x'
endif
   
if info.save_png eq 1 then begin
   write_png, star_dir + info.figs_subdir + '/' + catalog_id + star_id + '_' + info.isla_subdir + '_' + info.global_subdir + '_' + modality + '.PNG', TVRD(/TRUE)
endif


; -------------------------------------------------------------------------------------------------------------------
; Save final outputs
; -------------------------------------------------------------------------------------------------------------------
; Evaluate empirical linewidth for each oscillation frequency

linewidth = get_linewidth(freq1_final,teff,numax)

; Save the value of DeltaNu from ACF and the value of epsilon from the epsilon diagram

get_lun,lun1
openw,lun1,peakbagging_filename_global
printf,lun1,'# nuMax (microHz), DeltaNu_ACF (microHz), DeltaNu_fit (microHz), epsilon, alpha, Teff (K), N_chunks, Flag depressed dipole',format='(A0)'
printf,lun1,numax,acf_dnu,best_dnu,best_epsi,best_alpha,teff,n_chunks,flag_depressed_dipole,format='(F0.4,F10.4,F10.4,F10.4,F10.4,F10.1,I6,I6)'

; Save the frequency positions of each chunk identified with the global fit

printf,lun1,'# Chunk index, start and end frequency values for each chunk (one per line), SNR',format='(A0)'

for i=0, n_separations-2 do begin
    tmp_chunk = where(freq le separations(i+1) and freq ge separations(i))
    spsd_chunk = spsd(tmp_chunk)
    bg_level_chunk = bg_level_local(tmp_chunk)
    snr_chunk = max(spsd_chunk)/mean(bg_level_chunk)
    printf,lun1,i,separations(i),separations(i+1),snr_chunk,format='(I0,F10.3,F15.3,F15.2)'
endfor

; Save the individual oscillation frequencies from the global fit, their uncertainties, mode identification and associated FWHM from Ball+18

printf,lun1,'# n, l, frequency (microHz), 1-sigma frequency (microHz), ASEF maximum (iterations), FWHM from predictions (microHz)',format='(A0)'

for i=0,n_freq-1 do begin
    printf,lun1,order_number(i),angular_degree(i),freq1_final(i),freq_sig1_final(i),asef_maximum_final(i),linewidth(i),format='(I0,I5,F15.5,F15.5,I10,F15.5)'
endfor

free_lun,lun1
end
