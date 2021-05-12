pro find_islands_chunk,catalog_id,star_id,run,threshold_asef,teff,force=force
; -------------------------------------------------------------------------------------------------------------------
; Author:      Enrico Corsaro
; e-mail:      enrico.corsaro@inaf.it
; Date:        August 2019
; Place:       Catania, Italy
; Purpose:     This routine processes the chunk multi-modal sampling obtained with DIAMONDS in order to identify 
;              meaningful frequency peaks related to oscillation modes. It also applies a tag to each selected 
;              frequency peak, corresponding to the expected mode identification, up to the angular degree level.
; -------------------------------------------------------------------------------------------------------------------
COMMON CONFIG,cp
COMMON STAR,info
COMMON GRAPHIC,pp,lp,sp,lpe
COMMON DIAMONDS,dp

modality = 'CHUNK'

setup_computation
run = strcompress(string(run),/remove_all)
star_dir = info.peakbagging_results_dir + catalog_id + star_id + '/'
peakbagging_filename_global = info.peakbagging_results_dir + catalog_id + star_id + '/' + info.summary_subdir + '/' $
                              + catalog_id + star_id + info.peakbagging_filename_label + $
                                info.isla_subdir + '_' + info.global_subdir + '_' + 'GLOBAL.txt'

peakbagging_filename_chunk = info.peakbagging_results_dir + catalog_id + star_id + '/' + info.summary_subdir + '/'  $
                             + catalog_id + star_id + info.peakbagging_filename_label + '_' + info.isla_subdir + '_'

; Copy FAMED configuring parameters used in this run

famed_configuring_parameters_filename_copy = info.peakbagging_results_dir + catalog_id + star_id + '/' + info.summary_subdir + '/' $
    + info.local_configuring_parameters_filename + catalog_id + star_id + $
                                '_' + info.isla_subdir + '_' + run + '_' + modality + '.txt'

spawn,'cp ' + info.configuring_parameters_filename + ' ' + famed_configuring_parameters_filename_copy

; Read sampled frequency from DIAMONDS multi-modal fit

readcol,star_dir + info.isla_subdir + '/' + run + '/peakbagging_parameter000.txt',par0,format='D',/silent

; Read prior height (upper limit) and frequency range

readcol,star_dir + info.isla_subdir + '/' + run + '/peakbagging_hyperParametersUniform.txt',prior_down,prior_up,format='D,D',/silent

freq_left_chunk = prior_down(0)
freq_right_chunk = prior_up(0)
upper_height = prior_up(1)

; Read input linewidth and background name of the run

readcol,star_dir + info.isla_subdir + '/' + run + '/peakbagging_computationParameters.txt',config,format='A',/silent
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
    print,' PSD frequency range (microHz): [',strcompress(string(min(freq)),/remove_all),', ',strcompress(string(max(freq)),/remove_all),']'
    print,'-------------------------------------------------'
    print,''
endif


; -------------------------------------------------------------------------------------------------------------------
; Load the information available from the global fit
; -------------------------------------------------------------------------------------------------------------------

; Load asymptotic parameters

readcol,peakbagging_filename_global,acf_dnu,best_dnu,best_epsi,best_alpha,teff,n_chunks,flag_depressed_dipole,format='x,F,F,F,F,F,I,I',numline=2,comment='#',/silent

acf_dnu = acf_dnu(0)
best_dnu = best_dnu(0)
best_epsi = best_epsi(0)
best_alpha = best_alpha(0)
teff = teff(0)
n_chunks = n_chunks(0)
flag_depressed_dipole = flag_depressed_dipole(0)

; Load the SNR for the given chunk

readcol,peakbagging_filename_global,chunk_number,freq_left,freq_right,snr_chunks,format='I,D,D,F',numline=n_chunks,skipline=3,comment='#',/silent

run_number = fix(run)
snr = snr_chunks(run_number)
left_bound = freq_left(run_number)
right_bound = freq_right(run_number)

; Load global frequencies

readcol,peakbagging_filename_global,enn_global,ell_global,freq_global,freq_sig_global,fwhm_global,     $
format='I,I,F,F,x,F',comment='#',skipline=3+n_chunks,/silent

; Select only global frequencies relevant for the chunk

tmp_chunk = where(freq_global le right_bound and freq_global ge left_bound)
if tmp_chunk(0) ne -1 then begin
    freq_global = freq_global(tmp_chunk)
    freq_sig_global = freq_sig_global(tmp_chunk)
    enn_global = enn_global(tmp_chunk)
    ell_global = ell_global(tmp_chunk)
    fwhm_global = fwhm_global(tmp_chunk)
endif else begin
    print,'This chunk is empty according to global fit. Skip chunk.'
    return
endelse

; Select the radial mode. If the global modality was ran correctly, there should be only one radial mode, but there could be none as well.
; If more than one radial mode is present, then pick up the one closest to the right frequency bound of the chunk.

tmp_radial = where(ell_global eq 0)
if tmp_radial(0) ne -1 then begin
    if n_elements(tmp_radial) eq 1 then begin
        enn_radial = enn_global(tmp_radial)
        enn_radial = enn_radial(0)
        freq_radial = freq_global(tmp_radial)
        freq_radial = freq_radial(0)
        freq_sig_radial = freq_sig_global(tmp_radial)
        freq_sig_radial = freq_sig_radial(0)
        fwhm_radial = fwhm_global(tmp_radial)
        fwhm_radial = fwhm_radial(0)
    endif else begin
        freq_radial = freq_global(tmp_radial)
        good_radial_index = closest(freq_radial,right_bound)
        good_radial_index = good_radial_index(0) + min(tmp_radial)
        enn_radial = enn_global(good_radial_index)
        enn_radial = enn_radial(0)
        freq_radial = freq_global(good_radial_index)
        freq_radial = freq_radial(0)
        freq_sig_radial = freq_sig_global(good_radial_index)
        freq_sig_radial = freq_sig_radial(0)
        fwhm_radial = fwhm_global(good_radial_index)
        fwhm_radial = fwhm_radial(0)
    endelse
    n_radial_chunk = n_elements(freq_radial)
endif else begin
    n_radial_chunk = 0
endelse

; Select the dipole mode(s). There might be none as well, depending if it was selected during the global fit.

tmp_dipole = where(ell_global eq 1)
if tmp_dipole(0) ne -1 then begin
    ; If more than one dipole mode frequency is selected, take the one closest to nu0 - best_dnu/2.

    if n_elements(tmp_dipole) gt 1 then begin
        freq_dipole = freq_global(tmp_dipole)
        
        if n_radial_chunk eq 0 then begin
            ; If the radial mode is not present, take the dipole frequency closest to the central frequency of the chunk

            dipole_global_index = closest((left_bound + right_bound)/2.,freq_dipole)
        endif else begin
            dipole_global_index = closest(freq_radial - best_dnu/2.,freq_dipole)
        endelse

        best_dipole_index = dipole_global_index + min(tmp_dipole)
        freq_dipole = freq_global(best_dipole_index)
        freq_dipole = freq_dipole(0)
        enn_dipole = enn_global(best_dipole_index)
        enn_dipole = enn_dipole(0)
        freq_sig_dipole = freq_sig_global(best_dipole_index)
        freq_sig_dipole = freq_sig_dipole(0)
        fwhm_dipole = fwhm_global(best_dipole_index)
        fwhm_dipole = fwhm_dipole(0)
    endif else begin
        enn_dipole = enn_global(tmp_dipole)
        enn_dipole = enn_dipole(0)
        freq_dipole = freq_global(tmp_dipole)
        freq_dipole = freq_dipole(0)
        freq_sig_dipole = freq_sig_global(tmp_dipole)
        freq_sig_dipole = freq_sig_dipole(0)
        fwhm_dipole = fwhm_global(tmp_dipole)
        fwhm_dipole = fwhm_dipole(0)
    endelse
    n_dipole_chunk = n_elements(freq_dipole)
endif else begin
    n_dipole_chunk = 0
endelse

; -------------------------------------------------------------------------------------------------------------------
; Set up plotting panels for either X term or EPS
; -------------------------------------------------------------------------------------------------------------------

if info.print_on_screen eq 1 then begin
    !p.multi = pp.chunk.pmulti
    position_sampling = pp.chunk.position_sampling       ; Sampling
    position_asef = pp.chunk.position_asef               ; ASEF
    position_psd = pp.chunk.position_psd                 ; PSD
    
    if info.save_eps eq 0 then begin
        set_plot,'x'
        window,1,xs=pp.xsize,ys=pp.ysize,xpos=pp.xpos,ypos=pp.ypos,title=catalog_id + star_id + ' RUN: ' + run
        device,decomposed=0,true_color=8,retain=2
    endif else begin
        set_plot,'PS'
        filename_star = star_dir + info.figs_subdir + '/' + catalog_id + star_id + '_' + info.isla_subdir + '_' + run + '_' + modality + '.eps'
        device,filename=filename_star,xs=pp.xsize,ys=pp.ysize,/encapsulated,/color,bits=8
    endelse
endif


; -------------------------------------------------------------------------------------------------------------------
; Load the nuMax information as obtained from a previous background fit
; -------------------------------------------------------------------------------------------------------------------

readcol,star_dir + 'gaussianEnvelopeParameters.txt',gauss_par,format='D',/silent
numax = gauss_par(1)
scaling_dnu = compute_scaling_dnu(numax)


; -------------------------------------------------------------------------------------------------------------------
; Set a binwidth proportional to the expected separation between adjacent peaks (normally taken as d02/2) if MS
; or smaller if RG star.
; -------------------------------------------------------------------------------------------------------------------
par_range = max(par0) - min(par0)
min_separation = get_minimum_freq_separation(best_dnu,0)
binwidth = 1.d0*min_separation/cp.min_n_bins
n_bins = ceil(par_range/binwidth)
nest_iter = findgen(n_elements(par0))


; -------------------------------------------------------------------------------------------------------------------
; Plot the nested sampling evolution of the frequency parameter as obtained by Diamonds.
; -------------------------------------------------------------------------------------------------------------------
if info.print_on_screen eq 1 then begin
    plot_sampling,par0,position_sampling
endif


; -------------------------------------------------------------------------------------------------------------------
; Compute an Average Shifted Histogram (ASH) of the distribution of nested iterations
; -------------------------------------------------------------------------------------------------------------------
; Compute the ASEF with input resolution for extracting the local maxima

ash = compute_asef(par0,nest_iter,n_bins)
par_hist = ash.x
asef_hist = ash.y


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
get_range_divisions,par_hist,asef_hist,index_maximum,range_maximum,divisions_maximum,/chunk


; -------------------------------------------------------------------------------------------------------------------
; Based on the identified frequency ranges, estimate the oscillation frequencies using the sampling information 
; obtained by DIAMONDS.
; -------------------------------------------------------------------------------------------------------------------
n_freq = n_maxima

; Compute a smoothed PSD by some average FWHM to evaluate its values at each local maxima of the ASEF.
; In case of RG, adopt a finer smoothing window to overcome the problem of very narrow mixed modes.

avg_fwhm = mean(get_linewidth(maximum,teff,numax))
if best_dnu le cp.dnu_rg then begin
    avg_fwhm = fit_linewidth*cp.smoothing_fwhm_factor_rg
endif 

smth_bins = avg_fwhm/freqbin
spsd = smooth(psd,smth_bins,/edge_truncate)
tmp_good = where(freq le max(par_hist) and freq ge min(par_hist))
psd_total = psd
freq_total = freq
psd = psd(tmp_good)
freq = freq(tmp_good)
spsd = spsd(tmp_good)

; Load the background level of the star

readcol,star_dir + 'backgroundLevel.txt',bg_level,format='x,D',/silent,comment='#'

bg_level_local = bg_level(tmp_good)
sampled_estimates = evaluate_sampling_frequencies(par0,par_hist,freq,spsd,maximum,range_maximum)
freq1 = sampled_estimates.freq1
freq_sig1 = sampled_estimates.freq_sig1
sampling_counts = sampled_estimates.sampling_counts
spsd_maximum = sampled_estimates.spsd_maximum


if info.print_on_screen eq 1 then begin
    plot_asef,par_hist,asef_hist,maximum,range_maximum,threshold,position_asef

    for i=0, n_elements(freq_global)-1 do begin
        oplot,[freq_global(i),freq_global(i)],[0,max(asef_hist)*1.4],linestyle=2,color=80,thick=3
    endfor
endif

; Define some weights useful for identification of proper frequency peaks during mode identification

sampling_weights = sampling_counts/total(sampling_counts)
spsd_weights = spsd_maximum/total(spsd_maximum)
asef_weights = asef_maximum/total(asef_maximum)

if info.print_on_screen eq 1 then begin
   print,' Total number of local maxima found: ',n_maxima
endif

; Define some average asymptotic frequency spacings useful for the computation

ap = get_asymptotic_parameters(numax,best_dnu,teff)
angular_degree = intarr(n_freq) + 1
d02 = ap.d02
d01 = ap.d01
d03 = ap.d03

spawn,'ls -1 ' + peakbagging_filename_chunk + '*' + modality + '.txt',filename_summary,/stderr
str_length = strlen(filename_summary(0))

; If possible, update the d02 spacing from existing chunk outputs

flag_median_d02_active = 0
if strmid(filename_summary(0),str_length-3,3) eq 'txt' then begin
    d02_array = fltarr(n_elements(filename_summary))
    for jj=0, n_elements(filename_summary)-1 do begin
        readcol,filename_summary(jj),d02_chunk,format='x,D',numline=2,comment='#',/silent
        d02_array(jj) = d02_chunk
    endfor

    ; Remove zeros from possible solutions (e.g. if a chunk had no l=2 mode detected)

    tmp_zero = where(d02_array eq 0.,complement=tmp_nonzero)
    if tmp_zero(0) ne -1 then begin
        d02_array = d02_array(tmp_nonzero)
    endif

    median_d02 = median(d02_array)
    max_d02 = max(d02_array)

    if (max_d02 lt d02 or max_d02 eq 0) then max_d02 = d02
    if median_d02 eq 0 then median_d02 = d02
    flag_median_d02_active = 1
endif else begin
    median_d02 = d02
    max_d02 = d02
endelse


; -------------------------------------------------------------------------------------------------------------------
; -------------------------------------------------------------------------------------------------------------------
; Raial mode and Quadrupole mode
; -------------------------------------------------------------------------------------------------------------------
; -------------------------------------------------------------------------------------------------------------------
; Start by obtaining a more accurate position of the current radial mode by using the (possible) radial mode of the
; previous (or next) chunk. The previous (or next) radial mode will correspond to a chunk with a higher SNR (also assuming that the chunks
; are being calculated by decreasing SNR order). If no value could be found, double the step.

previous_run_number = strcompress(string(run_number - 1,format='(I)'),/remove_all)
next_run_number = strcompress(string(run_number + 1,format='(I)'),/remove_all)
previous_run_number2 = strcompress(string(run_number - 2,format='(I)'),/remove_all)
next_run_number2 = strcompress(string(run_number + 2,format='(I)'),/remove_all)

flag_previous_radial_mode_found = 0
flag_next_radial_mode_found = 0

reference_central_freq = (freq_left_chunk + freq_right_chunk)/2.

if reference_central_freq ge numax then begin
    if previous_run_number ge 0 then begin
        if file_test(peakbagging_filename_chunk + previous_run_number + '_' + modality + '.txt') ne 0 then begin
            readcol,peakbagging_filename_chunk + previous_run_number + '_' + modality + '.txt',enn_previous_chunk,ell_previous_chunk,freq_previous_chunk,freq_sig_previous_chunk,   $
               format='I,I,x,D,D',comment='#',skipline=3,/silent
            tmp_previous_radial = where(ell_previous_chunk eq 0)
           
            ; Check whether the selected chunk contains a l=0 mode, otherwise move to the one before it.
            
            if tmp_previous_radial(0) ne -1 then begin
                if info.print_on_screen eq 1 then begin
                    print,''
                    print,' Resuming l=0 frequency from previous chunk...'
                    print,''
                endif
                
                freq_previous_radial = freq_previous_chunk(tmp_previous_radial)
                freq_sig_previous_radial = freq_sig_previous_chunk(tmp_previous_radial)
                enn_previous_radial = enn_previous_chunk(tmp_previous_radial)
                enn_radial = enn_previous_radial + 1
                enn_radial = enn_radial(0)
                freq_previous_radial = freq_previous_radial(0)
                freq_sig_previous_radial = freq_sig_previous_radial(0)
                flag_previous_radial_mode_found = 1
               
                
                ; Set new global radial mode frequency
                 
                freq_radial = freq_previous_radial + best_dnu*(1.0 + best_alpha*(enn_radial - 0.5 - numax/best_dnu))
                fwhm_radial = get_linewidth(freq_radial,teff,numax)
               
                if n_radial_chunk eq 0 then begin
                    if freq_radial le max(par0) then begin
                        n_radial_chunk = 1
                        freq_sig_radial = freq_sig_dipole
                    endif
                endif
               
                if n_dipole_chunk eq 0 then begin
                    freq_dipole = freq_radial - best_dnu/2. - d01
                    freq_sig_dipole = freq_sig_radial
                    n_dipole_chunk = 1
                endif
            endif else begin
                if previous_run_number2 ge 0 then begin
                    if file_test(peakbagging_filename_chunk + previous_run_number2 + '_' + modality + '.txt') ne 0 then begin
                        readcol,peakbagging_filename_chunk + previous_run_number2 + '_' + modality + '.txt',enn_previous_chunk,ell_previous_chunk,freq_previous_chunk,freq_sig_previous_chunk,   $
                           format='I,I,x,D,D',comment='#',skipline=3,/silent
                        tmp_previous_radial = where(ell_previous_chunk eq 0)
                        
                        if tmp_previous_radial(0) ne -1 then begin
                            if info.print_on_screen eq 1 then begin
                                print,''
                                print,' Resuming l=0 frequency from second previous chunk...'
                                print,''
                            endif
                            
                            freq_previous_radial = freq_previous_chunk(tmp_previous_radial)
                            freq_sig_previous_radial = freq_sig_previous_chunk(tmp_previous_radial)
                            enn_previous_radial = enn_previous_chunk(tmp_previous_radial)
                            enn_radial = enn_previous_radial + 2
                            enn_radial = enn_radial(0)
                            freq_previous_radial = freq_previous_radial(0)
                            freq_sig_previous_radial = freq_sig_previous_radial(0)
                            flag_previous_radial_mode_found = 2
                           
                            
                            ; Set new global radial mode frequency
                             
                            freq_radial = freq_previous_radial + 2*best_dnu*(1.0 + best_alpha*((2*enn_radial-1)/2. - 0.5 - numax/best_dnu))
                            fwhm_radial = get_linewidth(freq_radial,teff,numax)
                           
                            if n_radial_chunk eq 0 then begin
                                if freq_radial le max(par0) then begin
                                    n_radial_chunk = 1
                                    freq_sig_radial = freq_sig_dipole
                                endif
                            endif
                           
                            if n_dipole_chunk eq 0 then begin
                                freq_dipole = freq_radial - best_dnu/2. - d01
                                freq_sig_dipole = freq_sig_radial
                                n_dipole_chunk = 1
                            endif
                        endif
                    endif
                endif
            endelse
        endif
    endif
endif else begin
    if next_run_number le n_chunks-1 then begin
        if file_test(peakbagging_filename_chunk + next_run_number + '_' + modality + '.txt') ne 0 then begin
            readcol,peakbagging_filename_chunk + next_run_number + '_' + modality + '.txt',enn_next_chunk,ell_next_chunk,freq_next_chunk,freq_sig_next_chunk,   $
               format='I,I,x,D,D',comment='#',skipline=3,/silent
            tmp_next_radial = where(ell_next_chunk eq 0)
           
            ; Check whether the selected chunk contains a l=0 mode, otherwise move to the one next to it.

            if tmp_next_radial(0) ne -1 then begin
                if info.print_on_screen eq 1 then begin
                    print,''
                    print,' Resuming l=0 frequency from next chunk...'
                    print,''
                endif
                
                freq_next_radial = freq_next_chunk(tmp_next_radial)
                freq_sig_next_radial = freq_sig_next_chunk(tmp_next_radial)
                enn_next_radial = enn_next_chunk(tmp_next_radial)
                enn_next_radial = enn_next_radial(0)
                enn_radial = enn_next_radial - 1
                enn_radial = enn_radial(0)
                freq_next_radial = freq_next_radial(0)
                freq_sig_next_radial = freq_sig_next_radial(0)
                flag_next_radial_mode_found = 1
               
               
                ; Set new global radial mode frequency
                
                freq_radial = freq_next_radial - best_dnu*(1.0 + best_alpha*(enn_next_radial - 0.5 - numax/best_dnu))
                fwhm_radial = get_linewidth(freq_radial,teff,numax)
                
                if n_radial_chunk eq 0 then begin
                    if freq_radial ge min(par0) then begin
                        n_radial_chunk = 1
                        freq_sig_radial = freq_sig_dipole
                    endif
                endif
               
                if n_dipole_chunk eq 0 then begin
                    freq_dipole = freq_radial - best_dnu/2. - d01
                    freq_sig_dipole = freq_sig_radial
                    n_dipole_chunk = 1
                endif
            endif else begin
                if next_run_number2 le n_chunks-1 then begin
                    if file_test(peakbagging_filename_chunk + next_run_number2 + '_' + modality + '.txt') ne 0 then begin
                        readcol,peakbagging_filename_chunk + next_run_number2 + '_' + modality + '.txt',enn_next_chunk,ell_next_chunk,freq_next_chunk,freq_sig_next_chunk,   $
                           format='I,I,x,D,D',comment='#',skipline=3,/silent
                        tmp_next_radial = where(ell_next_chunk eq 0)
                        
                        if tmp_next_radial(0) ne -1 then begin
                            if info.print_on_screen eq 1 then begin
                                print,''
                                print,' Resuming l=0 frequency from second next chunk...'
                                print,''
                            endif

                            freq_next_radial = freq_next_chunk(tmp_next_radial)
                            freq_sig_next_radial = freq_sig_next_chunk(tmp_next_radial)
                            enn_next_radial = enn_next_chunk(tmp_next_radial)
                            enn_next_radial = enn_next_radial(0)
                            enn_radial = enn_next_radial - 2
                            enn_radial = enn_radial(0)
                            freq_next_radial = freq_next_radial(0)
                            freq_sig_next_radial = freq_sig_next_radial(0)
                            flag_next_radial_mode_found = 2
                           
                           
                            ; Set new global radial mode frequency
                            
                            freq_radial = freq_next_radial - 2*best_dnu*(1.0 + best_alpha*((2*enn_next_radial-1)/2. - 0.5 - numax/best_dnu))
                            fwhm_radial = get_linewidth(freq_radial,teff,numax)
                            
                            if n_radial_chunk eq 0 then begin
                                if freq_radial ge min(par0) then begin
                                    n_radial_chunk = 1
                                    freq_sig_radial = freq_sig_dipole
                                endif
                            endif
                           
                            if n_dipole_chunk eq 0 then begin
                                freq_dipole = freq_radial - best_dnu/2. - d01
                                freq_sig_dipole = freq_sig_radial
                                n_dipole_chunk = 1
                            endif
                        endif
                    endif
                endif
            endelse
        endif
    endif
endelse

flag_quadrupole_found = 0

if n_radial_chunk ne 0 then begin
    ; -------------------------------------------------------------------------------------------------------------------
    ; CASE 1: One radial mode is found in the chunk from the global modality
    ; -------------------------------------------------------------------------------------------------------------------
    
    ; ---------------------------------
    ; Radial mode
    ; ---------------------------------
    ; There must be only one reference radial mode in the chunk for the routine to work properly. 
    ; If the global modality was run successfully, this should not constitute a problem.
    
    order_number = intarr(n_freq) + enn_radial - 1 
   
    ; Select the frequency corresponding to the closest and highest (in counts) peak to the global l=0. 
    ; This is the most likely one to be the correct radial mode frequency because of its large linewidth as
    ; compared to other dipole modes (in the case of RG) or spurious peaks. 
    ; This applies only if more than one frequency is present.
    ; Possible contamination by l=2 mode should be avoided.
  
    freq_diff = abs(freq1 - freq_radial)
    freq_weights = 1.d0/freq_diff
    freq_weights /= total(freq_weights)
    freq_ww = freq_weights/max(freq_weights)
    asef_ww = asef_weights/max(asef_weights)
    spsd_ww = spsd_weights/max(spsd_weights)
    sampling_ww = alog(sampling_counts)/max(alog(sampling_counts))

    ; If we have a radial mode solution from previous (or next) chunks, then make the frequency position of the global radial mode a much stronger constraint
    
    weight_freq_fraction = cp.weight_freq_fraction
    if flag_previous_radial_mode_found ne 0 or flag_next_radial_mode_found ne 0 then begin
        weight_freq_fraction = cp.weight_freq_fraction_enhanced
    endif

    total_ww = weight_freq_fraction*freq_ww + cp.weight_asef_fraction*asef_ww + cp.weight_spsd_fraction*spsd_ww + cp.weight_sampling_fraction*sampling_ww
    total_ww /= total(total_ww)

    upper_limit_freq_radial = freq_radial + max_d02*cp.d02_factor_search_range

    ; Make sure that the upper limit frequency for radial mode search is not larger than chunk upper limit frequency
    
    if upper_limit_freq_radial gt right_bound then upper_limit_freq_radial = right_bound

    ; Force input value if required 

    if cp.upper_limit_freq_radial ne 0 then begin
        upper_limit_freq_radial = cp.upper_limit_freq_radial
    endif

    if info.print_on_screen eq 1 then begin
        print,' Upper limit frequency for l=0: ' + strcompress(string(upper_limit_freq_radial),/remove_all) + ' muHz'
    endif

    additional_freq_sig_radial = 0.

    if flag_previous_radial_mode_found ne 0 then begin
        additional_freq_sig_radial = freq_sig_previous_radial
    endif 

    if flag_next_radial_mode_found ne 0 then begin
        additional_freq_sig_radial = freq_sig_next_radial
    endif 

    if n_dipole_chunk ne 0 then begin
        tmp_radial = where(freq1 ge (freq_radial - freq_sig_radial - additional_freq_sig_radial - d02) and freq1 lt upper_limit_freq_radial and $
            freq1 gt freq_dipole)
    endif else begin
        tmp_radial = where(freq1 ge (freq_radial - freq_sig_radial - additional_freq_sig_radial - d02) and freq1 lt upper_limit_freq_radial)
    endelse
    
    if tmp_radial(0) ne -1 then begin
        total_ww_radial = total_ww(tmp_radial)
        max_ww = max(total_ww_radial,index)
        radial_index = index + min(tmp_radial)
        min_ww = min(total_ww)
    endif else begin
        print,' The l=0 mode could not be located based on the position from the global fit. '
        print,' This could be caused by the low SNR of the chunk.'
        return
    endelse
    
    ; Control check for assessing radial mode frequency peak identification
    if cp.plot_weights_radial eq 1 then begin
        loadct,39
        !p.multi=[0,2,1]
        window,2
        plot,freq1,sampling_ww,psym=-6,yr=[0,1],xr=[min(par_hist),max(par_hist)],xtitle='!3Frequency ('+ sp.freq_unit_str + ')',ytitle='Weights',   $
            xticklen=0.02,yticklen=0.03,xthick=pp.xthick,ythick=pp.ythick,charthick=pp.charthick
        oplot,freq1,asef_ww,color=200,psym=-6
        oplot,freq1,freq_ww,color=250,psym=-6
        oplot,freq1,spsd_ww,color=160,psym=-6
        plot,freq1,total_ww,psym=-6,xr=[min(par_hist),max(par_hist)],yr=[0,max(total_ww)],/nodata,xtitle='!3Frequency ('+ sp.freq_unit_str + ')',   $
            ytitle='Total Weight',xticklen=0.02,yticklen=0.03,xthick=pp.xthick,ythick=pp.ythick,charthick=pp.charthick
        oplot,freq1,total_ww,psym=-6,color=90
        oplot,[freq1(radial_index),freq1(radial_index)],[0,max(total_ww)*1.2],color=205,thick=3,linestyle=2
        return
    endif
    
    ; Check if the selected peak is not the adjacent l=2 by assessing the total weight of the subsequent modes.
    
    if radial_index lt n_freq-1 then begin
        candidate_radial_index = where(freq1 gt freq1(radial_index) and freq1 lt upper_limit_freq_radial)
        if candidate_radial_index(0) ne -1 then begin
            for kk=0, n_elements(candidate_radial_index)-1 do begin
                local_index = candidate_radial_index(kk)
                if total_ww(local_index) gt (abs(max_ww - min_ww)*cp.max_ratio_search_radial + min_ww) then begin
                    radial_index = local_index
                endif
            endfor
        endif
    endif
    
    freq_radial_chunk = freq1(radial_index)
    freq_sig_radial_chunk = freq_sig1(radial_index)
    low_cut_frequency = freq_radial_chunk - best_dnu*(1.0 + best_alpha*(enn_radial - 0.5 - numax/best_dnu)) + freq_sig_radial/2.0 
  
    ; Make sure that the previous l=0 mode is not included, in case the chunk belongs to the right-end tail of the oscillation envelope.
    ; First verify that if a radial mode has been found from previous chunks (to the left side), the low cut frequency is not including it.

    if flag_previous_radial_mode_found eq 1 then begin
        if low_cut_frequency lt freq_previous_radial + freq_sig_radial/2. then begin
            low_cut_frequency = freq_previous_radial + freq_sig_radial/2.
        endif
    endif

    if flag_previous_radial_mode_found eq 2 then begin
        if low_cut_frequency lt freq_previous_radial + best_dnu*(1.0 + best_alpha*(enn_radial-1 - 0.5 - numax/best_dnu)) + freq_sig_radial/2. then begin
            low_cut_frequency = freq_previous_radial + best_dnu*(1.0 + best_alpha*(enn_radial-1 - 0.5 - numax/best_dnu)) + freq_sig_radial/2.
        endif
    endif

    ; However, the l=0 mode from the previous chunk may turnout to be wrong (e.g. if confused with an adjacent l=1 mixed mode). 
    ; To verify this do some additional check to test whether the selected low cut frequency can be reliable. 

    if freq_radial_chunk gt numax then begin
        first_chunk_indices = where(freq1 lt min(par_hist) + (max(par_hist) - min(par_hist))/cp.previous_radial_range_fraction)
        
        if first_chunk_indices(0) ne -1 then begin
            first_chunk_sampling_counts = sampling_counts(first_chunk_indices)
            first_chunk_asef = asef_maximum(first_chunk_indices)
            
            ; Do a first check on the sampling counts
    
            max_first_chunk_sampling_counts = max(first_chunk_sampling_counts,previous_radial_mode_index)
            if sampling_counts(radial_index) lt max_first_chunk_sampling_counts*cp.sampling_counts_fraction then begin
                ; Here update the lower cutting frequency of the chunk, as well as the radial mode index
    
                low_cut_frequency2 = freq1(previous_radial_mode_index) + freq_sig_radial_chunk
                if low_cut_frequency2 gt low_cut_frequency then low_cut_frequency = low_cut_frequency2

                radial_index_new = closest(freq1(previous_radial_mode_index) + best_dnu*(1.0 + best_alpha*(enn_radial - 0.5 - numax/best_dnu)),freq1)
                if sampling_counts(radial_index_new) ge sampling_counts(radial_index) then begin
                    radial_index = radial_index_new
                    freq_radial_chunk = freq1(radial_index)
                endif
            endif

            ; Do a second check on the ASEF maximum
    
            max_first_chunk_asef = max(first_chunk_asef,previous_radial_mode_index)
            if max_first_chunk_asef*cp.asef_saturation_fraction ge asef_maximum(radial_index) then begin
                ; Here update the lower cutting frequency of the chunk, as well as the radial mode index, if adequate.
    
                low_cut_frequency2 = freq1(previous_radial_mode_index) + freq_sig_radial_chunk
                if low_cut_frequency2 gt low_cut_frequency then low_cut_frequency = low_cut_frequency2

                radial_index_new = closest(freq1(previous_radial_mode_index) + best_dnu*(1.0 + best_alpha*(enn_radial - 0.5 - numax/best_dnu)),freq1)
                asef_threshold = (dp.isla.max_nested_it + dp.isla.n_live)/cp.asef_threshold_scaling_radial 
                
                if asef_maximum(radial_index_new) gt asef_threshold then begin
                    radial_index = radial_index_new
                    freq_radial_chunk = freq1(radial_index)
                endif
            endif
        endif
    endif

    ; Apply a final check on the upper bound for the low cut frequency

    if low_cut_frequency gt freq_radial_chunk - best_dnu*cp.dnu_lower_cut_fraction then low_cut_frequency = freq_radial_chunk - best_dnu*cp.dnu_lower_cut_fraction

    angular_degree(radial_index) = 0
    order_number(radial_index) = enn_radial

    ; Compute a local value for epsilon
    
    local_epsi = (freq_radial_chunk - best_dnu * enn_radial)/best_dnu

    ; ---------------------------------
    ; Quadrupole mode
    ; ---------------------------------
    ; Identify the quadrupole mode by assuming a small spacing d02.
    ; Compute more reliable values for d02 from high SNR chunks, if available.
    ; Mostly useful for the case of MS stars.
    
    freq_radial_chunk_org = freq_radial_chunk
    freq_sig_radial_chunk_org = freq_sig_radial_chunk

    flag_quadrupole_found = 1
    flag_duplet_fit = 0

    if best_dnu lt cp.dnu_rg then begin
        ; ---------------------------------
        ; Evolutionary stage: RG
        ; ---------------------------------
        quadrupole_freq_asymp = freq_radial_chunk - d02
        quadrupole_index = closest(quadrupole_freq_asymp,freq1)
        if quadrupole_index ne radial_index then begin
            freq_quadrupole_chunk = freq1(quadrupole_index)
            freq_sig_quadrupole_chunk = freq_sig1(quadrupole_index)
           
            ;freq_quadrupole_chunk_org = freq_quadrupole_chunk
            ;freq_sig_quadrupole_chunk_org = freq_sig_quadrupole_chunk
            range_maximum_quadrupole = range_maximum(*,quadrupole_index)
            range_maximum_radial = range_maximum(*,radial_index)

            angular_degree(quadrupole_index) = 2
            order_number(quadrupole_index) = enn_radial - 1

            ; Find the l=1 mixed modes closest to l=2 on each side of the peak.
            ; If the mixed mode falls within d02/2 from the l=2 frequency position then merge it with the l=2
   
            candidate_mixed = where(freq1 le freq_quadrupole_chunk + d02/cp.d02_scaling_merge_mixed and $
                                    freq1 ge freq_quadrupole_chunk - d02/cp.d02_scaling_merge_mixed and freq1 ne freq_quadrupole_chunk and $
                                    freq1 lt freq_radial_chunk, $
                                    complement=good_freq)
    
            if candidate_mixed(0) ne -1 then begin
                ; Remove the merged l=1 frequencies from the frequency set
    
                sampling_counts(quadrupole_index) = sampling_counts(quadrupole_index) + total(sampling_counts(candidate_mixed))
                merged_freq = freq1(candidate_mixed)
                   
                new_lower_bound = range_maximum(0, min(candidate_mixed))
                if new_lower_bound lt range_maximum(0,quadrupole_index) then begin
                    range_maximum(0,quadrupole_index) = new_lower_bound
                endif
                
                new_upper_bound = range_maximum(1, max(candidate_mixed))
                if new_upper_bound gt range_maximum(1,quadrupole_index) then begin
                    range_maximum(1,quadrupole_index) = new_upper_bound
                endif
               
                new_lower_bound = divisions_maximum(0, min(candidate_mixed))
                if new_lower_bound lt divisions_maximum(0,quadrupole_index) then begin
                    divisions_maximum(0,quadrupole_index) = new_lower_bound
                endif
                
                new_upper_bound = divisions_maximum(1, max(candidate_mixed))
                if new_upper_bound gt divisions_maximum(1,quadrupole_index) then begin
                    divisions_maximum(1,quadrupole_index) = new_upper_bound
                endif
                
                if good_freq(0) ne -1 then begin
                    ; Deprecate one range and division between the lowest and highest merged frequencies
                    
                    range_maximum = range_maximum(*,good_freq)
                    divisions_maximum = divisions_maximum(*,good_freq)
                   
                    freq1 = freq1(good_freq)
                    freq_sig1 = freq_sig1(good_freq)
                    angular_degree = angular_degree(good_freq)
                    order_number = order_number(good_freq)
                    sampling_counts = sampling_counts(good_freq)
                    sampling_weights = sampling_counts/total(sampling_counts)
                    asef_maximum = asef_maximum(good_freq)
                    spsd_maximum = spsd_maximum(good_freq)
                    spsd_weights = spsd_maximum/total(spsd_maximum)
                    asef_weights = asef_maximum/total(asef_maximum)
                    n_freq = n_elements(freq1)
                endif
                
                quadrupole_index = where(angular_degree eq 2)

                upper_bound = range_maximum(1,quadrupole_index)
                lower_bound = range_maximum(0,quadrupole_index)
                tmp_range = where(par0 lt upper_bound(0) and par0 ge lower_bound(0))
                par0_range = par0(tmp_range)
                
                freq_quadrupole_chunk = total(par0_range*tmp_range^2)/total(tmp_range^2)
                freq_sig_quadrupole_chunk = sqrt(total((par0_range-freq_quadrupole_chunk)^2*tmp_range^2)/total(tmp_range^2))

                freq1(quadrupole_index) = freq_quadrupole_chunk
                freq_sig1(quadrupole_index) = freq_sig_quadrupole_chunk
            endif

        endif else begin 
            ; No quadrupole has been found.
    
            flag_quadrupole_found = 0
        endelse
    endif else begin
        ; ---------------------------------
        ; Evolutionary stage: MS, SG
        ; ---------------------------------
        ; Estimate the position of the l=2 and l=0 peaks. 
        ; Consider as an upper limit for d02 given by the d02 from RG scaling in the case of stars in the subgiant regime. 
        ; This is because here d02_RG  constitutes an upper limit for the actual d02 (see White et al. 2011). 
        ; For MS stars, take as an upper limit for d02 the value Dnu/4 in order to incorporate also cases with very large d02.
        
        if flag_median_d02_active eq 1 then begin 
            max_local_d02 = median_d02*cp.d02_factor_search_range
        endif else begin
            if best_dnu le cp.dnu_sg then begin
                ; Here the star is considered a SG

                max_local_d02 = max_d02
            endif else begin
                ; Here the star is considered a MS

                max_local_d02 = best_dnu/4.
            endelse
        endelse

        ; Set up and run the multi-modal fit for the double peak, after checking whether the run already exists

        run_subdir = run + 'A'
        flag_duplet_fit = 1 
       
        if (file_test(star_dir + info.isla_subdir + '/' + run_subdir + '/peakbagging_computationParameters.txt') eq 0 or keyword_set(force)) then begin
            ; Set up prior boundaries
        
            freq_prior = [freq_radial_chunk - d02,upper_limit_freq_radial]
            height_prior = [prior_down(1),prior_up(1)]

            if cp.d02_prior_upper_duplet_fit ne 0 then begin
                max_local_d02 = cp.d02_prior_upper_duplet_fit
            endif

            d02_prior = [cp.d02_prior_lower_duplet_fit,max_local_d02]
            boundaries = [freq_prior,height_prior,d02_prior]

            filename = star_dir + info.isla_subdir + '/' + info.prior_filename + '_' + run_subdir + '.txt'
            write_diamonds_prior,filename,boundaries

            if info.print_on_screen eq 1 then begin
                print,' Performing l=2,0 duplet fit with DIAMONDS...'
            endif
           
            peakbagging_parameters = { subdir:       info.isla_subdir,  $
                                       run:          run_subdir,        $
                                       background:   background_name,   $
                                       fwhm:         get_linewidth(freq_radial_chunk,teff,numax),     $
                                       duplet:       1                  $
                                     }

            flag_computation_completed = run_peakbagging(catalog_id,star_id,peakbagging_parameters,0,0,0)
        endif else begin
            if info.print_on_screen eq 1 then begin
                print,' Load information from l=2,0 duplet fit with DIAMONDS.'
            endif
        endelse

        ; Read sampled frequency from DIAMONDS multi-modal fit for nu_0
        
        readcol,star_dir + info.isla_subdir + '/' + run_subdir + '/peakbagging_parameter000.txt',par_nu0,format='D',/silent
        
        ; Read sampled frequency from DIAMONDS multi-modal fit for d02
        
        readcol,star_dir + info.isla_subdir + '/' + run_subdir + '/peakbagging_parameter002.txt',par_d02,format='D',/silent
      
        ; Read posterior distribution from DIAMONDS multi-modal fit
        
        readcol,star_dir + info.isla_subdir + '/' + run_subdir + '/peakbagging_posteriorDistribution.txt',post,format='D',/silent

        ; Compute parameter estimates, using a weighted mean from the posterior
      
        post /= max(post) 
        nest_iter = findgen(n_elements(par_nu0)) 
        nu0 = total(post*par_nu0)/total(post)
        nu0_sig = sqrt(total(post*(par_nu0 - nu0)^2)/total(post))
        d02 = total(post*par_d02)/total(post)
        d02_sig = sqrt(total(post*(par_d02 - d02)^2)/total(post))

        freq_quadrupole_chunk = nu0 - d02
        freq_sig_quadrupole_chunk = sqrt(nu0_sig^2 + d02_sig^2)
        freq_radial_chunk = nu0
        freq_sig_radial_chunk = nu0_sig

        ; Now make frequency uncertainties consistent in definition with the remainder of the extracted set of frequencies
        
        midpoint_02 = (freq_radial_chunk + freq_quadrupole_chunk)/2.
        upper_bound_quadrupole = midpoint_02
        best_quadrupole_index = max(where(range_maximum(0,*) le freq_quadrupole_chunk))
        lower_bound_quadrupole = range_maximum(0,best_quadrupole_index)
        lower_bound_quadrupole = lower_bound_quadrupole(0)
       
        ; Put safe conditions on the lower bound for l=2
        
        if lower_bound_quadrupole ge upper_bound_quadrupole then begin
            print,' Adjusting lower frequency bound for l=2 due to swap with upper bound.'
            lower_bound_quadrupole = freq_quadrupole_chunk - d02
        endif

        if lower_bound_quadrupole lt freq_quadrupole_chunk - d02 then begin
            print,' Adjusting lower frequency bound for l=2 due to large range.'
            lower_bound_quadrupole = freq_quadrupole_chunk - d02
        endif

        lower_bound_radial = midpoint_02
        best_radial_index = max(where(range_maximum(0,*) le freq_radial_chunk))
        upper_bound_radial = range_maximum(1,best_radial_index)
        upper_bound_radial = upper_bound_radial(0)
        
        ; Put safe conditions on the lower bound for l=0
        
        if lower_bound_radial ge upper_bound_radial then begin
            print,' Adjusting upper frequency bound for l=0 due to swap with lower bound.'
            upper_bound_radial += abs(upper_bound_radial - lower_bound_radial)*2 
        endif

        ; Compute l=2 frequency uncertainty
        
        tmp_range = where(par0 lt upper_bound_quadrupole and par0 ge lower_bound_quadrupole)
        tmp_freq_range = where(freq lt upper_bound_quadrupole and freq ge lower_bound_quadrupole)
        tmp_hist_range = where(par_hist lt upper_bound_quadrupole and par_hist ge lower_bound_quadrupole)
        
        if tmp_range(0) eq -1 then begin
            print,' Not enough sampling points in the specified frequency range for l=2. Quitting program.'
            return 
        endif 
        
        par0_range = par0(tmp_range)
        spsd_range = spsd(tmp_freq_range)
        spsd_maximum_quadrupole = max(spsd_range)
        asef_maximum_quadrupole = max(asef_hist(tmp_hist_range))
        sampling_counts_quadrupole = total(nest_iter(tmp_range))
        freq_sig_quadrupole_chunk = sqrt(total((par0_range-freq_quadrupole_chunk)^2*tmp_range^2)/total(tmp_range^2))
       
        ; Compute l=0 frequency uncertainty
        
        tmp_range = where(par0 lt upper_bound_radial and par0 ge lower_bound_radial)
        tmp_freq_range = where(freq lt upper_bound_radial and freq ge lower_bound_radial)
        tmp_hist_range = where(par_hist lt upper_bound_radial and par_hist ge lower_bound_radial)
        
        if tmp_range(0) eq -1 then begin
            print,' Not enough sampling points in the specified frequency range for l=0. Quitting program.'
            return 
        endif

        iteration = 0
        while tmp_hist_range(0) eq -1 do begin
            best_radial_index++

            if best_radial_index le n_elements(freq1) then begin
                upper_bound_radial = range_maximum(iteration mod 2,best_radial_index)
                upper_bound_radial = upper_bound_radial(0)
                tmp_hist_range = where(par_hist lt upper_bound_radial and par_hist ge lower_bound_radial)
            endif else begin
                print,'Not enough bins in the specified ASEF range for l=0. Quitting program.'
                return 
            endelse

            iteration++
        endwhile

        par0_range = par0(tmp_range)
        spsd_range = spsd(tmp_freq_range)
        spsd_maximum_radial = max(spsd_range)
        asef_maximum_radial = max(asef_hist(tmp_hist_range))
        sampling_counts_radial = total(nest_iter(tmp_range))
        freq_sig_radial_chunk = sqrt(total((par0_range-freq_radial_chunk)^2*tmp_range^2)/total(tmp_range^2))

        ; Check whether there are old frequencies (from local maxima) inside this region of the l=2,0 duplet
        
        tmp_inside = where(freq1 ge lower_bound_quadrupole and freq1 le upper_bound_radial,complement=tmp_outside)
        tmp_radial_outside = where(angular_degree(tmp_outside) eq 0)
        if tmp_radial_outside(0) ne -1 then begin
            angular_degree(tmp_radial_outside) = 1
            order_number(tmp_radial_outside)--
        endif

        if tmp_inside(0) ne -1 and tmp_outside(0) ne -1 then begin
            ; If at least one frequency is found inside the l=2,0 range, first remove it (them) from the sample of frequencies
            
            range_maximum = range_maximum(*,tmp_outside)
            divisions_maximum = divisions_maximum(*,tmp_outside)
               
            freq1 = freq1(tmp_outside)
            freq_sig1 = freq_sig1(tmp_outside)
            angular_degree = angular_degree(tmp_outside)
            order_number = order_number(tmp_outside)
            sampling_counts = sampling_counts(tmp_outside)
            sampling_weights = sampling_counts/total(sampling_counts)
            asef_maximum = asef_maximum(tmp_outside)
            spsd_maximum = spsd_maximum(tmp_outside)
        endif

        ; Then add up the new l=2,0 frequencies to the list

        if tmp_outside(0) ne -1 then begin
            ; In this case there are other frequencies outside the l=2,0 range, so take them into account

            range_maximum = [[range_maximum],[lower_bound_quadrupole,upper_bound_quadrupole],[lower_bound_radial,upper_bound_radial]]
            divisions_maximum = [[divisions_maximum],[lower_bound_quadrupole,upper_bound_quadrupole],[lower_bound_radial,upper_bound_radial]]
           
            freq1 = [freq1,freq_quadrupole_chunk,freq_radial_chunk]
            freq_sig1 = [freq_sig1,freq_sig_quadrupole_chunk,freq_sig_radial_chunk]
            angular_degree = [angular_degree,2,0]
            order_number = [order_number,enn_radial-1,enn_radial]
            sampling_counts = [sampling_counts,sampling_counts_quadrupole,sampling_counts_radial]
            sampling_weights = sampling_counts/total(sampling_counts)
            asef_maximum = [asef_maximum,asef_maximum_quadrupole,asef_maximum_radial]
            spsd_maximum = [spsd_maximum,spsd_maximum_quadrupole,spsd_maximum_radial]

            ; Resort by increasing frequency order (useful if there are frequencies above l=0)
            
            tmp_sort = sort(freq1)
            range_maximum = range_maximum(*,tmp_sort)
            divisions_maximum = divisions_maximum(*,tmp_sort)
            
            freq1 = freq1(tmp_sort)
            freq_sig1 = freq_sig1(tmp_sort)
            angular_degree = angular_degree(tmp_sort)
            order_number = order_number(tmp_sort)
            sampling_counts = sampling_counts(tmp_sort)
            sampling_weights = sampling_counts/total(sampling_counts)
            asef_maximum = asef_maximum(tmp_sort)
            spsd_maximum = spsd_maximum(tmp_sort)
            spsd_weights = spsd_maximum/total(spsd_maximum)
            asef_weights = asef_maximum/total(asef_maximum)
            n_freq = n_elements(freq1)
        endif else begin
            ; Only the candidate l=2,0 mode frequencies are available

            range_maximum = [[lower_bound_quadrupole,upper_bound_quadrupole],[lower_bound_radial,upper_bound_radial]]
            divisions_maximum = [[lower_bound_quadrupole,upper_bound_quadrupole],[lower_bound_radial,upper_bound_radial]]
           
            freq1 = [freq_quadrupole_chunk,freq_radial_chunk]
            freq_sig1 = [freq_sig_quadrupole_chunk,freq_sig_radial_chunk]
            angular_degree = [2,0]
            order_number = [enn_radial-1,enn_radial]
            sampling_counts = [sampling_counts_quadrupole,sampling_counts_radial]
            sampling_weights = sampling_counts/total(sampling_counts)
            asef_maximum = [asef_maximum_quadrupole,asef_maximum_radial]
            spsd_maximum = [spsd_maximum_quadrupole,spsd_maximum_radial]
            spsd_weights = spsd_maximum/total(spsd_maximum)
            asef_weights = asef_maximum/total(asef_maximum)
            n_freq = n_elements(freq1)
        endelse

        ; Obtain the radial and quadrupole mode indices from frequency array
        
        quadrupole_index = where(angular_degree eq 2)
        quadrupole_index = quadrupole_index(0)
        radial_index = where(angular_degree eq 0)
        radial_index = radial_index(0)
    endelse

    ; Compute local small spacing d02 for this chunk.
   
    if flag_quadrupole_found ne 0 then begin 
        local_d02 = freq_radial_chunk - freq_quadrupole_chunk
    endif else begin
        local_d02 = d02
    endelse

    actual_d02 = local_d02
    
    if flag_median_d02_active eq 1 then begin
        if median_d02 gt local_d02 then actual_d02 = median_d02
    endif

    ; Update the lower frequency limit for this chunk

    low_cut_frequency2 = freq_radial_chunk - best_dnu*(1.0 + best_alpha*(enn_radial - 0.5 - numax/best_dnu)) + freq_sig_radial
    low_cut_frequency3 = freq_radial_chunk - best_dnu + freq_sig_radial
    if low_cut_frequency2 gt low_cut_frequency then low_cut_frequency = low_cut_frequency2
    if low_cut_frequency3 gt low_cut_frequency then begin
        tmp_lower_order = where(freq1 le low_cut_frequency3)
        if tmp_lower_order(0) ne -1 then begin
            low_cut_frequency = low_cut_frequency3
        endif
    endif

;save,filename='PSM06550_2_duplet_fit.idl',par_hist,asef_hist,range_maximum,index_maximum,freq1,freq_sig1,freq_radial,freq_radial_chunk,freq_quadrupole_chunk,$
;    freq,psd,spsd,bg_level_local,order_number,angular_degree,freq_radial_chunk_org,freq_sig_radial_chunk_org    
endif else begin
     ; -------------------------------------------------------------------------------------------------------------------
     ; CASE 2: No radial mode is found in the chunk from the global modality
     ; -------------------------------------------------------------------------------------------------------------------
     ; Compute the estimated radial mode frequency from the asymptotic relation for radial modes.
     ; Assume no quadrupole mode in this case.

     enn_radial = enn_dipole + 1
     freq_radial_chunk = asymptotic_relation_radial(enn_radial,numax,[best_dnu,best_epsi,best_alpha])
     fwhm_radial_fit = fwhm_dipole
     order_number = intarr(n_freq) + enn_dipole

     local_epsi = 0
     local_d02 = 0
     actual_d02 = d02

     freq_quadrupole_chunk = freq_radial_chunk - d02
     low_cut_frequency = freq_radial_chunk - best_dnu + freq_sig_dipole
endelse

; Force input value for lower limit frequency if required

if cp.low_cut_frequency ne 0 then begin
    low_cut_frequency = cp.low_cut_frequency
endif

if info.print_on_screen eq 1 then begin
    print,' Lower limit for chunk frequency range: '+strcompress(string(low_cut_frequency,format='(F0.3)'),/remove_all)+ ' muHz'
    print,''
endif


; Now check that the previous order l=0 frequency was not selected as a potential dipole mode by the ASEF.
; If it is selected, then remove it from the sample of identified frequencies.
; This could happen for RG and SG stars, where a chunk overlapping with the previous one is used during the computation of the fits.

tmp_below_radial_double = where(freq1 le low_cut_frequency,complement=good_freq1_indices)
if tmp_below_radial_double(0) ne -1 then begin
    freq1 = freq1(good_freq1_indices)
    freq_sig1 = freq_sig1(good_freq1_indices)
    asef_maximum = asef_maximum(good_freq1_indices)
    spsd_maximum = spsd_maximum(good_freq1_indices)
    spsd_weights = spsd_maximum/total(spsd_maximum)
    asef_weights = asef_maximum/total(asef_maximum)
    sampling_counts = sampling_counts(good_freq1_indices)
    sampling_weights = sampling_counts/total(sampling_counts)
    angular_degree = angular_degree(good_freq1_indices)
    order_number = order_number(good_freq1_indices)
    range_maximum = range_maximum(*,good_freq1_indices)
    divisions_maximum = divisions_maximum(*,good_freq1_indices)
    n_freq = n_elements(freq1)   
endif

; Finally remove any possible frequency identified above the actual l=0 mode.

tmp_above_radial = where(freq1 gt freq_radial_chunk,complement=good_freq1_indices)
if tmp_above_radial(0) ne -1 then begin
    freq1 = freq1(good_freq1_indices)
    freq_sig1 = freq_sig1(good_freq1_indices)
    asef_maximum = asef_maximum(good_freq1_indices)
    spsd_maximum = spsd_maximum(good_freq1_indices)
    spsd_weights = spsd_maximum/total(spsd_maximum)
    asef_weights = asef_maximum/total(asef_maximum)
    sampling_counts = sampling_counts(good_freq1_indices)
    sampling_weights = sampling_counts/total(sampling_counts)
    angular_degree = angular_degree(good_freq1_indices)
    order_number = order_number(good_freq1_indices)
    range_maximum = range_maximum(*,good_freq1_indices)
    divisions_maximum = divisions_maximum(*,good_freq1_indices)
    n_freq = n_elements(freq1)
endif
  
; Update radial and quadrupole indices.
; Perform a Lorentzian profile fit to the PSD to assess the FWHM of the radial peak.
; Consider the largest range around the peak. If no radial peak is present, then
; take as FWHM of the radial mode that of the adjacent dipole mode.

if n_radial_chunk ne 0 then begin
    radial_index = closest(freq_radial_chunk,freq1)
    quadrupole_index = closest(freq_quadrupole_chunk,freq1)
    run_subdir = run + '_radial_fwhm'
    run_names = run_subdir + strcompress(string(indgen(cp.n_fwhm_fit)),/remove_all) 
    
    if (file_test(star_dir + info.pb_subdir + '/' + run_subdir + '0/peakbagging_computationParameters.txt') eq 0 or keyword_set(force)) then begin
        right_bound = max([range_maximum(1,radial_index),divisions_maximum(1,radial_index)])
        left_bound = max([range_maximum(0,radial_index),divisions_maximum(0,radial_index)])
        
        tmp0 = where(freq ge left_bound and freq le right_bound)
        freq0 = freq(tmp0)
        spsd0 = spsd(tmp0)
        bg_peak = mean(bg_level_local(tmp0))
        response_peak = mean((sin(!pi/2. * freq0/nyq) / (!pi/2. * freq0/nyq))^2)
        amplitude_radial = sqrt((abs(spsd_maximum(radial_index) - bg_peak)/response_peak)*!pi*fwhm_radial*cp.fwhm_magnification_factor_radial)
        freq_prior_radial = [range_maximum(0,radial_index),range_maximum(1,radial_index)]
        amplitude_prior_radial = [0.0,amplitude_radial]
        data_freq_boundaries = [left_bound,right_bound]

        data_range_filenames = strarr(cp.n_fwhm_fit)
        prior_filenames = strarr(cp.n_fwhm_fit)
  
        ; Make sure to enlarge FWHM prior if fit fails.
        
        flag_computation_completed = 0
        iterations = 1
        
        left_fwhm = cp.fwhm_lower_bound

        if info.print_on_screen eq 1 then begin
            print,''
            print,' Fitting Linewidth of radial mode of the chunk.'
        endif

        while flag_computation_completed ne 1 do begin
            ; Allow for a very narrow FWHM
        
            fwhm_prior = [left_fwhm,fwhm_radial*cp.fwhm_magnification_factor_radial*iterations]
            
            if best_dnu lt cp.dnu_threshold then begin
                boundaries = [freq_prior_radial,amplitude_prior_radial,fwhm_prior]
            endif else begin
                amplitude_quadrupole = sqrt((abs(spsd_maximum(quadrupole_index) - bg_peak)/response_peak)*!pi*fwhm_radial*cp.fwhm_magnification_factor_radial)
                freq_prior_quadrupole = [range_maximum(0,quadrupole_index),range_maximum(1,quadrupole_index)]
                amplitude_prior_quadrupole = [0.0,amplitude_quadrupole]
                boundaries = [freq_prior_quadrupole,amplitude_prior_quadrupole,fwhm_prior,freq_prior_radial,amplitude_prior_radial,fwhm_prior]
            endelse

            for k=0, cp.n_fwhm_fit-1 do begin
                data_range_filenames(k) = star_dir + info.pb_subdir + '/frequencyRange_' + run_names(k) + '.txt'
                prior_filenames(k) = star_dir + info.pb_subdir + '/' + info.prior_filename + '_' + run_names(k) + '.txt'
                
                write_diamonds_data_range,data_range_filenames(k),data_freq_boundaries
                write_diamonds_prior,prior_filenames(k),boundaries
            endfor

            peakbagging_parameters = { subdir:          info.pb_subdir,     $
                                       run:             run_names,          $
                                       background:      background_name,    $
                                       fwhm:            -1.0,               $
                                       filename_run:    run_subdir          $
                                     }

            flag_computation_completed_array = run_peakbagging(catalog_id,star_id,peakbagging_parameters,2,0,0)
            flag_computation_completed = max(flag_computation_completed_array)

            iterations++
        endwhile

        if info.save_test_files ne 1 then begin
            for k=0, cp.n_fwhm_fit-1 do begin
                file_delete,data_range_filenames(k)
                file_delete,prior_filenames(k)
            endfor
        endif
    endif

    fwhm_radial_fit_array = fltarr(cp.n_fwhm_fit)

    for k=0, cp.n_fwhm_fit-1 do begin
        spawn,'ls -1 ' + star_dir + info.pb_subdir + '/' + run_names(k) + '/peakbagging_parameter0*txt',filenames
        
        if n_elements(filenames) eq 3 then begin
            fwhm_parameter = '002'
        endif else begin
            fwhm_parameter = '005'
        endelse
        
        ; Read the sampled FWHM of the radial mode.
        
        readcol,star_dir + info.pb_subdir + '/' + run_names(k) + '/peakbagging_parameter' + fwhm_parameter + '.txt',par_fwhm0,format='D',/silent
        
        ; Load the posterior samples to compute Bayesian mean estimate for FWHM_0
        
        readcol,star_dir + info.pb_subdir + '/' + run_names(k) + '/peakbagging_posteriorDistribution.txt',post,format='D',/silent
        
        ; Compute parameter estimate, using a weighted average
        
        post /= max(post)
        fwhm_radial_fit_array(k) = total(par_fwhm0*post)/total(post)
    endfor

    fwhm_radial_fit = median(fwhm_radial_fit_array)
    
    if info.print_on_screen eq 1 then begin
        print,' FWHM (l=0) = '+strcompress(string(fwhm_radial_fit,format='(F0.3)'),/remove_all) + ' muHz'
        print,''
    endif
endif


; Define a detection probability array for each frequency peak. Set it to -99, meaning all
; frequency peaks are not tested (i.e. considered detected at the beginning, with p_detection = 100 %).
; In the peak blending test, save whether the peak has a blending or not (blending_profile_flag either 1 or 0).
; In the sinc^2 test (Sinc^2 profile vs. Lorentzian profile), save which profile is better to model (sinc_profile_flag either 1 or 0).
; Define a rotation probability array, specifying the probability that the peak is split
; by rotation. Assume a default no rotation test performed (set to -99), meaning that rotation is assumed not detected 
; at the beginning (with p_rotation = 0 %). 
; Define a duplet probability array, specifying whether the peak is a duplet or not, and assume no duplicity at 
; the beginning (set to -99), meaning that the peak is a single peak (p_duplet = 0 %). 

detection_probability = fltarr(n_elements(freq1)) - 99.0
rotation_probability = fltarr(n_elements(freq1)) - 99.0
duplicity_probability = fltarr(n_elements(freq1)) - 99.0
sinc_profile_flag = intarr(n_elements(freq1))
blending_profile_flag = intarr(n_elements(freq1))

if keyword_set(force) then begin
    if file_test(star_dir + info.pb_subdir + '/detectionProbability_' + run + '*.txt') ne 0 then begin
        spawn,'ls -1 '+ star_dir + info.pb_subdir + '/detectionProbability_'+ run + '*.txt',delete_filenames
        file_delete,delete_filenames,/recursive
    endif

    if file_test(star_dir + info.pb_subdir + '/rotationProbability_' + run + '*.txt') ne 0 then begin
        spawn,'ls -1 '+ star_dir + info.pb_subdir + '/rotationProbability_' + run + '*.txt',delete_filenames
        file_delete,delete_filenames,/recursive
    endif
endif 

; Skip this part if no radial mode is found

if n_radial_chunk ne 0 then begin
    ; -------------------------------------------------------------------------------------------------------------------
    ; Peak significance and blending test for l=2,0
    ; -------------------------------------------------------------------------------------------------------------------
    ; Perform the peak significance test on each frequency peak found to have a critical SNR
    ; Evaluate SNR using an estimate for the amplitude of the peak
    ; Clear previous peak detection probabilities if fit is forced
    
    right_bound = range_maximum(1,radial_index)
    left_bound = range_maximum(0,quadrupole_index)
    tmp_freq_peak = where(freq le right_bound and freq ge left_bound)
    bg_peak = mean(bg_level_local(tmp_freq_peak))
    response_peak = mean((sin(!pi/2. * freq(tmp_freq_peak)/nyq) / (!pi/2. * freq(tmp_freq_peak)/nyq))^2)
    height_ratio_quadrupole = spsd_maximum(quadrupole_index)/(bg_peak*cp.height_ratio_threshold)
    height_ratio_radial = spsd_maximum(radial_index)/(bg_peak*cp.height_ratio_threshold)

    ; Distinguish between RG + late SG stars, and MS + early SG stars (especially hot stars, where linewidths are larger).

    if (best_dnu lt cp.dnu_threshold) then begin
        ; This is the case of either RG or late SG.
        ; Here perform a standard significance test only for the l=2,0 pair, but separately
        ; for each mode. Assume Lorentzian profiles only.
        ; For the quadrupole mode adopt a larger linewidth upper limit for the prior. This is because
        ; l=2 modes are generally affected by the presence of l=2 mixed modes in evolved stars.
        
        fwhm_quadrupole = fwhm_radial_fit*cp.fwhm_magnification_factor_quadrupole
        fwhm_radial = fwhm_radial_fit*cp.fwhm_magnification_factor_radial
        fwhm = [fwhm_quadrupole,fwhm_radial]
        freq_left = [freq1(quadrupole_index) - freq_sig1(quadrupole_index), freq1(radial_index) - freq_sig1(radial_index)]
        freq_right = [freq1(quadrupole_index) + freq_sig1(quadrupole_index), freq1(radial_index) + freq_sig1(radial_index)]
        amplitude_quadrupole = sqrt((abs(spsd_maximum(quadrupole_index) - bg_peak)/response_peak)*!pi*fwhm_quadrupole)
        amplitude_radial = sqrt((abs(spsd_maximum(radial_index) - bg_peak)/response_peak)*!pi*fwhm_radial)
        amplitude = [amplitude_quadrupole,amplitude_radial]
        height_ratio = [height_ratio_quadrupole,height_ratio_radial]
      
        ; For each peak consider the maximum data range between the range and the division 
        
        left_bound_quadrupole = max([range_maximum(0,quadrupole_index),divisions_maximum(0,quadrupole_index)]) 
        right_bound_quadrupole = max([range_maximum(1,quadrupole_index),divisions_maximum(1,quadrupole_index)]) 
        left_bound_radial = max([range_maximum(0,radial_index),divisions_maximum(0,radial_index)]) 
        right_bound_radial = max([range_maximum(1,radial_index),divisions_maximum(1,radial_index)]) 

        range_quadrupole_radial = [[left_bound_quadrupole,right_bound_quadrupole],[left_bound_radial,right_bound_radial]]
        angular_degree_quadrupole_radial = [angular_degree(quadrupole_index),angular_degree(radial_index)]
        quadrupole_radial_index = [quadrupole_index,radial_index]
        
        test_dirs = strarr(2,2)
        prior_filenames = strarr(2,2)
        data_range_filenames = strarr(2,2)
        probability_filenames = strarr(2)
        run_test = strarr(1)
        flag_run_test = 0

        for i=0, n_elements(height_ratio)-1 do begin
            if height_ratio(i) lt 1.0 then begin
                peak_number = strcompress(string(quadrupole_radial_index(i)),/remove_all)
                test_name = run + '_' + peak_number 
                
                test_nameA = test_name + 'A'           ; Only background
                test_nameB = test_name + 'B'           ; Lorentzian profile and background
                
                test_dirA = star_dir + info.pb_subdir + '/' + test_nameA
                test_dirB = star_dir + info.pb_subdir + '/' + test_nameB
                
                filenameA = star_dir + info.pb_subdir + '/' + info.prior_filename + '_' + test_nameA + '.txt'
                filenameB = star_dir + info.pb_subdir + '/' + info.prior_filename + '_' + test_nameB + '.txt'
                
                filename_rangeA = star_dir + info.pb_subdir + '/frequencyRange_' + test_nameA + '.txt'
                filename_rangeB = star_dir + info.pb_subdir + '/frequencyRange_' + test_nameB + '.txt'
                
                probability_filename = star_dir + info.pb_subdir + '/detectionProbability_' + test_name + '.txt'

                test_dirs(*,i) = [test_dirA,test_dirB]
                probability_filenames(i) = probability_filename
                prior_filenames(*,i) = [filenameA,filenameB]
                data_range_filenames(*,i) = [filename_rangeA,filename_rangeB]

                if file_test(probability_filename) eq 0 then begin
                    ; Set up priors for the peak test profile
                    
                    data_freq_boundaries = [range_quadrupole_radial(0,i),range_quadrupole_radial(1,i)]
                    freq_prior = [freq_left(i),freq_right(i)]
                    amplitude_prior = [0.0,amplitude(i)]
                    left_fwhm = cp.fwhm_lower_bound
                    fwhm_prior = [left_fwhm,fwhm(i)]
                    bkg_prior = [cp.bkg_prior_lower,cp.bkg_prior_upper]
                    boundariesA = [bkg_prior]
                    boundariesB = [freq_prior,amplitude_prior,fwhm_prior,bkg_prior]

                    write_diamonds_data_range,data_range_filenames(0,i),data_freq_boundaries
                    write_diamonds_data_range,data_range_filenames(1,i),data_freq_boundaries
                    write_diamonds_prior,prior_filenames(0,i),boundariesA
                    write_diamonds_prior,prior_filenames(1,i),boundariesB

                    run_test = [run_test,test_nameA,test_nameB]
                    flag_run_test++
                endif
            endif
        endfor

        tested_peak_indices = where(height_ratio lt 1.0)
        
        if tested_peak_indices(0) ne -1 then begin
            if info.print_on_screen eq 1 then begin
                print,' '
                print,' ---------------------------------------------'
                print,' Peak Detection Test for l=2,0 pair.' 
            endif

            if flag_run_test gt 0 then begin
                run_test = run_test(1:*)
                
                p_BA = fltarr(cp.n_peak_test,n_elements(tested_peak_indices))
                        
                for k=0, cp.n_peak_test-1 do begin
                    peakbagging_parameters = { subdir:       info.pb_subdir,              $
                                               run:          run_test,                    $
                                               background:   background_name,             $
                                               fwhm:         -1.0,                        $
                                               filename_run: test_name                    $
                                             }
                
                    flag_computation_completed = run_peakbagging(catalog_id,star_id,peakbagging_parameters,1,0,0)
                    
                    if info.print_on_screen eq 1 then begin
                        print,' Iteration #'+strcompress(string(k),/remove_all)+' completed.' 
                    endif

                    for i=0,n_elements(tested_peak_indices)-1 do begin
                        readcol,test_dirs(0,tested_peak_indices(i)) + '/peakbagging_evidenceInformation.txt', ln_evidA,format='D,x,x',/silent
                        readcol,test_dirs(1,tested_peak_indices(i)) + '/peakbagging_evidenceInformation.txt', ln_evidB,format='D,x,x',/silent
                       
                        if (ln_evidA ge ln_evidB) then begin
                            ln_evidAB = ln_evidA + alog(1.+exp(ln_evidB-ln_evidA)) 
                        endif else begin
                            ln_evidAB = ln_evidB + alog(1.+exp(ln_evidA-ln_evidB))
                        endelse
                        
                        p_BA(k,i) = exp(ln_evidB - ln_evidAB)
                    endfor
                endfor
                
                for i=0,n_elements(tested_peak_indices)-1 do begin
                    get_lun, lun1
                    openw, lun1, probability_filenames(tested_peak_indices(i))
                    printf, lun1, '# Detection probabilities for l=' + strcompress(string(angular_degree_quadrupole_radial(tested_peak_indices(i))),/remove_all) + ' at ' +  $
                            strcompress(string(freq1(quadrupole_radial_index(tested_peak_indices(i))),format='(F0.2)'),/remove_all) + ' muHz', format = '(A0)'
                    printf, lun1, '# Each line corresponds to a different run of the same test.', format = '(A0)'
                    printf, lun1, '# Col 1: p(BA).', format = '(A0)'
                    printf, lun1, p_BA(*,i), format = '(F0.4)'
                    free_lun, lun1

                    if info.save_test_files eq 0 then begin
                        ; Remove test folders and prior files if required
                       
                        file_delete,test_dirs(0,tested_peak_indices(i)),/recursive
                        file_delete,test_dirs(1,tested_peak_indices(i)),/recursive
                        file_delete,data_range_filenames(0,tested_peak_indices(i))
                        file_delete,data_range_filenames(1,tested_peak_indices(i))
                        file_delete,prior_filenames(0,tested_peak_indices(i))
                        file_delete,prior_filenames(1,tested_peak_indices(i))
                    endif
                endfor
            endif
                
            for i=0,n_elements(tested_peak_indices)-1 do begin
                readcol,probability_filenames(tested_peak_indices(i)),p_BA,format='F',/silent,comment='#'
                
                ; Consider the maximum probability among the different runs. Approximate up to third decimal digit.
                
                max_p_BA = float(round(max(p_BA)*1.d3)/1.d3)

                if info.print_on_screen eq 1 then begin
                    print,''
                    print,' Peak detection test for l=' + strcompress(string(angular_degree_quadrupole_radial(tested_peak_indices(i))),/remove_all) + ' at ' +  $
                            strcompress(string(freq1(quadrupole_radial_index(tested_peak_indices(i))),format='(F0.2)'),/remove_all) + ' muHz'
                    print,' P (One Lorentzian vs Only Background): ',max_p_BA
                    print,''
                endif

                detection_probability(quadrupole_radial_index(tested_peak_indices(i))) = max_p_BA
            endfor

            if info.print_on_screen eq 1 then begin
                print,' Peak Detection Test for l=2,0 pair Completed.' 
                print,' ---------------------------------------------'
                print,' '
            endif
        endif    
    endif else begin
        ; This is the case of MS and early SG stars.
        ; Here perform a blending test on top of the peak significance test, only for the l=2,0 pair.
        ; This is because it is likely that the l=2,0 modes are blended due to the large linewidths in hotter stars.
        ; Testing l=2 and l=0 only separately from one another may not represent a realistic condition.
       
        fwhm_quadrupole = fwhm_radial_fit*cp.fwhm_magnification_factor_radial
        fwhm_radial = fwhm_radial_fit*cp.fwhm_magnification_factor_radial
        amplitude_quadrupole = sqrt((abs(spsd_maximum(quadrupole_index) - bg_peak)/response_peak)*!pi*fwhm_quadrupole)
        amplitude_radial = sqrt((abs(spsd_maximum(radial_index) - bg_peak)/response_peak)*!pi*fwhm_radial)
        amplitude_max = max([amplitude_quadrupole,amplitude_radial])

        if (height_ratio_quadrupole lt 1.0) or (height_ratio_radial lt 1.0) then begin
            peak_number = strcompress(string(radial_index),/remove_all)
            test_name = run + '_' + peak_number 
            
            test_nameA = test_name + 'A'           ; Only background
            test_nameB = test_name + 'B'           ; Lorentzian profile and background
            test_nameD = test_name + 'D'           ; Two Lorentzian profiles and background
            
            test_dirA = star_dir + info.pb_subdir + '/' + test_nameA
            test_dirB = star_dir + info.pb_subdir + '/' + test_nameB
            test_dirD = star_dir + info.pb_subdir + '/' + test_nameD
            
            filenameA = star_dir + info.pb_subdir + '/' + info.prior_filename + '_' + test_nameA + '.txt'
            filenameB = star_dir + info.pb_subdir + '/' + info.prior_filename + '_' + test_nameB + '.txt'
            filenameD = star_dir + info.pb_subdir + '/' + info.prior_filename + '_' + test_nameD + '.txt'
            
            filename_rangeA = star_dir + info.pb_subdir + '/frequencyRange_' + test_nameA + '.txt'
            filename_rangeB = star_dir + info.pb_subdir + '/frequencyRange_' + test_nameB + '.txt'
            filename_rangeD = star_dir + info.pb_subdir + '/frequencyRange_' + test_nameD + '.txt'

            probability_filename = star_dir + info.pb_subdir + '/detectionProbability_' + test_name + '.txt'
            
            if info.print_on_screen eq 1 then begin
                print,' '
                print,' ----------------------------------------------------------'
                print,' Peak Detection and Blending Test for l=2,0 pair.' 
            endif

            if file_test(probability_filename) eq 0 then begin
                ; Set up priors for the peak test profile
                ; For each peak consider the maximum data range between the range and the division 
                
                left_bound_quadrupole = max([range_maximum(0,quadrupole_index),divisions_maximum(0,quadrupole_index)])
                right_bound_radial = max([range_maximum(1,radial_index),divisions_maximum(1,radial_index)]) 
                data_freq_boundaries = [left_bound_quadrupole,right_bound_radial]
                
                freq_prior_quadrupole = [freq1(quadrupole_index)-freq_sig1(quadrupole_index),freq1(quadrupole_index)+freq_sig1(quadrupole_index)]
                freq_prior_radial = [freq1(radial_index)-freq_sig1(radial_index),freq1(radial_index)+freq_sig1(radial_index)]
                freq_prior = [range_maximum(0,quadrupole_index),range_maximum(1,radial_index)]
                
                amplitude_prior_quadrupole = [0.0,amplitude_quadrupole]
                amplitude_prior_radial = [0.0,amplitude_radial]
                amplitude_prior = [0.0,amplitude_max]
               
                right_fwhm = fwhm_radial_fit*cp.fwhm_magnification_factor_radial
                left_fwhm = cp.fwhm_lower_bound
                fwhm_prior = [left_fwhm,right_fwhm]
                bkg_prior = [cp.bkg_prior_lower,cp.bkg_prior_upper]
                
                boundariesA = [bkg_prior]
                boundariesB = [freq_prior,amplitude_prior,fwhm_prior,bkg_prior]
                boundariesD = [freq_prior_quadrupole,amplitude_prior_quadrupole,fwhm_prior,freq_prior_radial,amplitude_prior_radial,fwhm_prior,bkg_prior]
                
                write_diamonds_data_range,filename_rangeA,data_freq_boundaries
                write_diamonds_data_range,filename_rangeB,data_freq_boundaries
                write_diamonds_data_range,filename_rangeD,data_freq_boundaries
                write_diamonds_prior,filenameA,boundariesA
                write_diamonds_prior,filenameB,boundariesB
                write_diamonds_prior,filenameD,boundariesD
    
                p_BA = fltarr(cp.n_peak_test)
                p_DA = fltarr(cp.n_peak_test)
                p_DB = fltarr(cp.n_peak_test)
                nu0 = fltarr(cp.n_peak_test)
                
                for k=0, cp.n_peak_test-1 do begin
                    peakbagging_parameters = { subdir:       info.pb_subdir,                       $
                                               run:          [test_nameA,test_nameB,test_nameD],   $
                                               background:   background_name,                      $
                                               fwhm:         -1.0,                                 $
                                               filename_run: test_name                             $
                                             }
                    flag_computation_completed = run_peakbagging(catalog_id,star_id,peakbagging_parameters,1,0,0)
                    
                    if info.print_on_screen eq 1 then begin
                        print,' Iteration #'+strcompress(string(k),/remove_all)+' completed.' 
                    endif
                    
                    readcol, test_dirA + '/peakbagging_evidenceInformation.txt', ln_evidA,format='D,x,x',/silent
                    readcol, test_dirB + '/peakbagging_evidenceInformation.txt', ln_evidB,format='D,x,x',/silent
                    readcol, test_dirD + '/peakbagging_evidenceInformation.txt', ln_evidD,format='D,x,x',/silent
                   
                    if (ln_evidA ge ln_evidB) then begin
                        ln_evidAB = ln_evidA + alog(1.+exp(ln_evidB-ln_evidA)) 
                    endif else begin
                        ln_evidAB = ln_evidB + alog(1.+exp(ln_evidA-ln_evidB))
                    endelse
                    
                    if (ln_evidA ge ln_evidD) then begin
                        ln_evidAD = ln_evidA + alog(1.+exp(ln_evidD-ln_evidA)) 
                    endif else begin
                        ln_evidAD = ln_evidD + alog(1.+exp(ln_evidA-ln_evidD))
                    endelse
                    
                    if (ln_evidB ge ln_evidD) then begin
                        ln_evidBD = ln_evidB + alog(1.+exp(ln_evidD-ln_evidB)) 
                    endif else begin
                        ln_evidBD = ln_evidD + alog(1.+exp(ln_evidB-ln_evidD))
                    endelse

                    p_BA(k) = exp(ln_evidB - ln_evidAB)
                    p_DA(k) = exp(ln_evidD - ln_evidAD)
                    p_DB(k) = exp(ln_evidD - ln_evidBD)
                    
                    ; Read the estimated central frequency from model B.
                    
                    readcol,star_dir + info.pb_subdir + '/' + test_nameB + '/peakbagging_parameter000.txt', par_nu0,format='D',/silent
                    
                    ; Load the posterior samples to compute Bayesian mean estimate for nu0
                    
                    readcol,star_dir + info.pb_subdir + '/' + test_nameB + '/peakbagging_posteriorDistribution.txt',post,format='D',/silent
               
                    ; Compute parameter estimate, using a weighted average
                    
                    post /= max(post)
                    nu0(k) = total(par_nu0*post)/total(post)
                endfor
                
                get_lun, lun1
                openw, lun1, probability_filename
                printf, lun1, '# Detection probabilities for l=2,0 duplet at ' +   $
                        strcompress(string(freq1(quadrupole_index),format='(F0.2)'),/remove_all) + ' + ' + $
                        strcompress(string(freq1(radial_index),format='(F0.2)'),/remove_all) + ' muHz', format = '(A0)'
                printf, lun1, '# Each line corresponds to a different run of the same test.', format = '(A0)'
                printf, lun1, '# Col 1: p(BA), Col 2: p(DA), Col 3: p(DB), Col 4: nu0 (microHz).', format = '(A0)'

                for l=0, cp.n_peak_test-1 do begin
                    printf, lun1, p_BA(l),p_DA(l),p_DB(l),nu0(l), format = '(F0.4,F10.4,F10.4,F15.4)'
                endfor

                free_lun, lun1

                if info.save_test_files eq 0 then begin
                    ; Remove test folders and prior files if required
                    
                    file_delete,test_dirA,/recursive
                    file_delete,test_dirB,/recursive
                    file_delete,test_dirD,/recursive
                    file_delete,filename_rangeA
                    file_delete,filename_rangeB
                    file_delete,filename_rangeD
                    file_delete,filenameA
                    file_delete,filenameB
                    file_delete,filenameD
                endif
            endif else begin
                readcol,probability_filename,p_BA,p_DA,p_DB,nu0,format='F,F,F,F',/silent,comment='#'
            endelse
               
            ; Consider the maximum probabilities among the different runs. Approximate up to third decimal digit.
            max_p_BA = float(round(max(p_BA)*1.d3)/1.d3)
            max_p_DA = float(round(max(p_DA)*1.d3)/1.d3)
            max_p_DB = float(round(max(p_DB)*1.d3)/1.d3)

            if info.print_on_screen eq 1 then begin
                print,''
                print,' Peak detection test for l=2,0 duplet at ' +   $
                        strcompress(string(freq1(quadrupole_index),format='(F0.2)'),/remove_all) + ' + ' + $
                        strcompress(string(freq1(radial_index),format='(F0.2)'),/remove_all) + ' muHz'
                print,' P (One Lorentzian vs Only Background): ',max_p_BA
                print,' P (Two Lorentzians vs Only Background): ',max_p_DA
                print,' P (Two Lorentzians vs One Lorentzian): ',max_p_DB
                print,' '
            endif

            if max_p_DB gt 0.5 then begin
                ; Here we have a peak blending. The two peaks are treated together and are considered
                ; detected only if their p_DA >= 0.993.

                detection_probability(quadrupole_index) = max_p_DA
                detection_probability(radial_index) = max_p_DA
                blending_profile_flag(quadrupole_index) = 1
                blending_profile_flag(radial_index) = 1
            endif else begin
                detection_probability(quadrupole_index) = 0.0
                detection_probability(radial_index) = 0.0
                
                if (max_p_BA ge cp.detection_probability_threshold) then begin
                    ; Here we have no peak blending, and only one peak can be considered. 
                    ; Find which one between l=2 and l=0 has to be deeemd significant.

                    avg_nu0 = mean(nu0)
                    detected_index = closest(avg_nu0,[freq1(quadrupole_index),freq1(radial_index)])
                    detection_probability(detected_index+quadrupole_index) = max_p_BA
                endif
            endelse


            if info.print_on_screen eq 1 then begin
                print,' Peak Detection and Blending Test for l=2,0 pair Completed.' 
                print,' ----------------------------------------------------------'
                print,' '
            endif
        endif
    endelse
endif 


; -------------------------------------------------------------------------------------------------------------------
; -------------------------------------------------------------------------------------------------------------------
; Dipole mode(s) and Octupole mode
; -------------------------------------------------------------------------------------------------------------------
; -------------------------------------------------------------------------------------------------------------------
; At this stage, check for the presence of dipole mode(s), and subsequently l=3.
; Distinguish red giant stars from subgiants and main sequence (White et al. 2011). 
; If the star is hot (F-type, Teff > 6300 K), its Dnu will be smaller as compared to G-type star in same relative evolutionary stage.
; Then it is less likely to have an F-type star that is also a subgiant, hence showing mixed modes.
; For this reason, also incorporate cooler stars with higher Dnu as possible subgiants (see e.g. Appourchaux et al. 2012).
; If no dipole modes were identified from the global fit, still give it a try in the chunk modality
; and check whether there could be a potential dipole mode. This is assuming that l=0 was found!
; Start by defining the l=3 search region

if best_dnu lt cp.dnu_rg then begin
    octupole_freq_asymp = freq_radial_chunk - best_dnu/2. - d03
    octupole_freq_lower = freq_radial_chunk - best_dnu/2. - d03*cp.d03_upper_scaling_factor
    octupole_freq_upper = freq_radial_chunk - best_dnu/2. - d03*cp.d03_lower_scaling_factor
endif else begin
    ; Adopt approximated asymptotic relation by Bedding & Kjeldsen 2003 to locate l=3 region in less evolved stars.
    
    if (best_dnu lt cp.dnu_sg and teff lt cp.teff_sg) then begin
        d02_upper_scaling_factor = cp.d02_upper_scaling_factor_sg
        d02_lower_scaling_factor = cp.d02_lower_scaling_factor_sg
    endif else begin
        d02_upper_scaling_factor = cp.d02_upper_scaling_factor_ms
        d02_lower_scaling_factor = cp.d02_lower_scaling_factor_ms
    endelse
    
    octupole_freq_lower = freq_radial_chunk - best_dnu/2. - d02_upper_scaling_factor*actual_d02
    octupole_freq_upper = freq_radial_chunk - best_dnu/2. - d02_lower_scaling_factor*actual_d02
    octupole_freq_asymp = (octupole_freq_lower + octupole_freq_upper)/2.
endelse

if (best_dnu lt cp.dnu_sg and teff lt cp.teff_sg) then begin
    dipole_indices = where(angular_degree eq 1)
endif else begin
    ; For MS stars, only check the modes below nu0 - Dnu/2 + freq_dipole_sig, and above the lower limit of the l=3 region. 
    ; All remaining ones are flagged as undetected. This will speed up the mode identification process.
    
    if n_radial_chunk ne 0 then begin
        upper_dipole_limit = freq_radial_chunk - best_dnu/2.0 + freq_sig_radial
    endif else begin
        upper_dipole_limit = freq_radial_chunk - best_dnu/2.0 + freq_sig_dipole
    endelse

    if upper_dipole_limit gt freq_radial_chunk - best_dnu/4.0 then begin
        upper_dipole_limit = freq_radial_chunk - best_dnu/4.0
    endif

    dipole_indices = where((freq1 lt upper_dipole_limit) and (freq1 ge octupole_freq_lower) and (angular_degree eq 1))
    bad_dipole_indices = where(((freq1 ge upper_dipole_limit) or freq1 lt octupole_freq_lower) and angular_degree eq 1)
    
    if bad_dipole_indices(0) ne -1 then begin
        detection_probability(bad_dipole_indices) = 0.0
    endif
endelse

if dipole_indices(0) ne -1 then begin
    ; -------------------------------------------------------------------------------------------------------------------
    ; Peak significance test for candidate l=1 mode(s)
    ; -------------------------------------------------------------------------------------------------------------------
    ; Perform a peak significance test in order to identify all the significant peaks marked as l=1.
    ; If the star is a RG, add up the test for sinc^2 profile to account for possible unresolved mixed modes.
    
    height_ratio_array = fltarr(n_elements(dipole_indices))
    test_dirs = strarr(3,n_elements(dipole_indices))
    prior_filenames = strarr(3,n_elements(dipole_indices))
    data_range_filenames = strarr(3,n_elements(dipole_indices))
    probability_filenames = strarr(n_elements(dipole_indices))
    run_test = strarr(1)
    flag_run_test = 0

    for i=0, n_elements(dipole_indices)-1 do begin
        ; Take smallest range possible that is available to build the frequency prior of the given peak. This will avoid that the peak may result not
        ; significant if the parameter space is too large.
        ; For each peak consider the maximum data range between the range and the division  
        
        freq_index = dipole_indices(i)
        left_bound = max([range_maximum(0,freq_index),divisions_maximum(0,freq_index)])
        right_bound = max([range_maximum(1,freq_index),divisions_maximum(1,freq_index)])

        ; For RG and SG stars, adopt a narrower FWHM prior for those l=1 modes outside the l=3 region.
        ; This will speed up the computation of the peak test by reducing the parameter space.
        
        right_fwhm = fwhm_radial_fit*cp.fwhm_magnification_factor_radial*2
        left_fwhm = cp.fwhm_lower_bound

        if (best_dnu lt cp.dnu_sg and teff lt cp.teff_sg) then begin
            right_fwhm = fwhm_radial_fit*cp.fwhm_magnification_factor_radial
            
            ; Allow for a very narrow FWHM accounting for unresolved mixed modes
            
            if (freq1(freq_index) ge octupole_freq_upper) or (freq1(freq_index) le octupole_freq_lower) then begin
                right_fwhm = fwhm_radial_fit*cp.fwhm_magnification_factor_dipole
            endif

            ; Do not exceed the width given by the frequency range of the peak if the star is a RG, because l=1 modes are narrow. 
        
            if best_dnu lt cp.dnu_rg then begin
                range_peak = abs(range_maximum(1,freq_index) - range_maximum(0,freq_index))
                if right_fwhm gt range_peak then right_fwhm = range_peak
            endif
        endif

        tmp_freq_peak = where(freq le right_bound and freq ge left_bound)
        bg_peak = mean(bg_level_local(tmp_freq_peak))
        response_peak = mean((sin(!pi/2. * freq(tmp_freq_peak)/nyq) / (!pi/2. * freq(tmp_freq_peak)/nyq))^2)
        amplitude = sqrt((abs(max(spsd(tmp_freq_peak)) - bg_peak)/response_peak)*!pi*right_fwhm);*2
        height_ratio = spsd_maximum(freq_index)/(bg_peak*cp.height_ratio_threshold)
        height_ratio_array(i) = height_ratio
        
        peak_number = strcompress(string(freq_index),/remove_all)
        test_name = run + '_' + peak_number
        
        test_nameA = test_name + 'A'           ; Only background
        test_nameB = test_name + 'B'           ; Background and Lorentzian profile
        test_nameC = test_name + 'C'           ; Background and Sinc^2 profile
        
        test_dirA = star_dir + info.pb_subdir + '/' + test_nameA
        test_dirB = star_dir + info.pb_subdir + '/' + test_nameB
        test_dirC = star_dir + info.pb_subdir + '/' + test_nameC
        
        probability_filename = star_dir + info.pb_subdir + '/detectionProbability_' + test_name + '.txt'
        
        filenameA = star_dir + info.pb_subdir + '/' + info.prior_filename + '_' + test_nameA + '.txt'
        filenameB = star_dir + info.pb_subdir + '/' + info.prior_filename + '_' + test_nameB + '.txt'
        filenameC = star_dir + info.pb_subdir + '/' + info.prior_filename + '_' + test_nameC + '.txt'
        
        filename_rangeA = star_dir + info.pb_subdir + '/frequencyRange_' + test_nameA + '.txt'
        filename_rangeB = star_dir + info.pb_subdir + '/frequencyRange_' + test_nameB + '.txt'
        filename_rangeC = star_dir + info.pb_subdir + '/frequencyRange_' + test_nameC + '.txt'

        test_dirs(*,i) = [test_dirA,test_dirB,test_dirC]
        probability_filenames(i) = probability_filename
        prior_filenames(*,i) = [filenameA,filenameB,filenameC]
        data_range_filenames(*,i) = [filename_rangeA,filename_rangeB,filename_rangeC]
        
        if height_ratio lt 1.0 then begin
            if file_test(probability_filename) eq 0 then begin
                ; Set up priors for the peak test profile
               
                data_freq_boundaries = [left_bound,right_bound] 
                freq_prior = [freq1(freq_index) - freq_sig1(freq_index),freq1(freq_index) + freq_sig1(freq_index)]
                amplitude_prior = [0.0,amplitude]
                fwhm_prior = [left_fwhm,right_fwhm]
                bkg_prior = [cp.bkg_prior_lower,cp.bkg_prior_upper]
                height_prior = [0.0,spsd_maximum(freq_index)*1.5]
                sinc_prior = [-1,-1]
               
                boundariesA = [bkg_prior]
                boundariesB = [freq_prior,amplitude_prior,fwhm_prior,bkg_prior]
                
                write_diamonds_data_range,data_range_filenames(0,i),data_freq_boundaries
                write_diamonds_data_range,data_range_filenames(1,i),data_freq_boundaries
                write_diamonds_prior,prior_filenames(0,i),boundariesA
                write_diamonds_prior,prior_filenames(1,i),boundariesB
                
                if best_dnu lt cp.dnu_rg then begin
                    boundariesC = [freq_prior,height_prior,bkg_prior,sinc_prior]
                    
                    write_diamonds_data_range,data_range_filenames(2,i),data_freq_boundaries
                    write_diamonds_prior,prior_filenames(2,i),boundariesC
                    
                    run_test = [run_test,test_nameA,test_nameB,test_nameC]
                endif else begin
                    run_test = [run_test,test_nameA,test_nameB]
                endelse
                flag_run_test++
            endif
        endif
    endfor
 
    tested_peak_indices = where(height_ratio_array lt 1.0)
    if tested_peak_indices(0) ne -1 then begin
        if info.print_on_screen eq 1 then begin
            print,' '
            print,' --------------------------------------------------------'
            print,' Peak Detection Test for candidate l=1 mode(s).' 
        endif
 
        if flag_run_test gt 0 then begin
            run_test = run_test(1:*)

            if info.print_on_screen eq 1 then begin
                print,' A total of '+strcompress(string(n_elements(run_test)),/remove_all)+' tests must be performed.' 
            endif
            
            fwhm_detection_fit = fltarr(cp.n_peak_test,n_elements(tested_peak_indices))
            p_BA = fltarr(cp.n_peak_test,n_elements(tested_peak_indices))
            
            if best_dnu lt cp.dnu_rg then begin
                p_CA = fltarr(cp.n_peak_test,n_elements(tested_peak_indices))
                p_CB = fltarr(cp.n_peak_test,n_elements(tested_peak_indices))
            endif
            
            for k=0, cp.n_peak_test-1 do begin
                peakbagging_parameters = { subdir:       info.pb_subdir,        $
                                           run:          run_test,              $
                                           background:   background_name,       $
                                           fwhm:         -1.0,                  $
                                           filename_run: test_name              $
                                         }
                flag_computation_completed = run_peakbagging(catalog_id,star_id,peakbagging_parameters,1,0,0)
                
                if info.print_on_screen eq 1 then begin
                    print,' Iteration #'+strcompress(string(k),/remove_all)+' completed.' 
                endif
 
                for i=0, n_elements(dipole_indices(tested_peak_indices))-1 do begin
                    ; Compute the estimate for FWHM from model B.
                    
                    readcol,test_dirs(1,tested_peak_indices(i)) + '/peakbagging_parameter002.txt', par_fwhm,format='D',/silent
                    readcol,test_dirs(1,tested_peak_indices(i)) + '/peakbagging_posteriorDistribution.txt',post,format='D',/silent
 
                    post /= max(post)
                    fwhm_detection_fit(k,i) = total(par_fwhm*post)/total(post)

                    readcol,test_dirs(0,tested_peak_indices(i)) + '/peakbagging_evidenceInformation.txt', ln_evidA,format='D,x,x',/silent
                    readcol,test_dirs(1,tested_peak_indices(i)) + '/peakbagging_evidenceInformation.txt', ln_evidB,format='D,x,x',/silent                     
                    
                    if (ln_evidA ge ln_evidB) then begin
                        ln_evidAB = ln_evidA + alog(1.d0+exp(ln_evidB-ln_evidA)) 
                    endif else begin
                        ln_evidAB = ln_evidB + alog(1.d0+exp(ln_evidA-ln_evidB))
                    endelse
                    
                    p_BA(k,i) = exp(ln_evidB - ln_evidAB)
                    
                    if best_dnu lt cp.dnu_rg then begin
                        readcol, test_dirs(2,tested_peak_indices(i)) + '/peakbagging_evidenceInformation.txt', ln_evidC,format='D,x,x',/silent
                        
                        if (ln_evidA ge ln_evidC) then begin
                            ln_evidAC = ln_evidA + alog(1.d0+exp(ln_evidC-ln_evidA)) 
                        endif else begin
                            ln_evidAC = ln_evidC + alog(1.d0+exp(ln_evidA-ln_evidC))
                        endelse
                        
                        if (ln_evidB ge ln_evidC) then begin
                            ln_evidBC = ln_evidB + alog(1.d0+exp(ln_evidC-ln_evidB)) 
                        endif else begin
                            ln_evidBC = ln_evidC + alog(1.d0+exp(ln_evidB-ln_evidC))
                        endelse
                    
                        p_CA(k,i) = exp(ln_evidC - ln_evidAC)
                        p_CB(k,i) = exp(ln_evidC - ln_evidBC)
                    endif
                endfor
            endfor     
            
            for i=0, n_elements(dipole_indices(tested_peak_indices))-1 do begin
                get_lun, lun1
                openw, lun1, probability_filenames(tested_peak_indices(i))
                printf, lun1, '# Detection probabilities for frequency peak '+ $
                               strcompress(string(freq1(dipole_indices(tested_peak_indices(i))),format='(F0.3)'),/remove_all) + ' muHz', format = '(A0)'
                printf, lun1, '# Each line corresponds to a different run of the same test.', format = '(A0)'
                
                if best_dnu lt cp.dnu_rg then begin
                    printf, lun1, '# Col 1: p(BA), Col 2: p(CA), Col 3: p(CB), Col 4: FWHM single (microHz)', format = '(A0)'
                    for k=0, cp.n_peak_test-1 do begin
                        printf, lun1, p_BA(k,i),p_CA(k,i),p_CB(k,i),fwhm_detection_fit(k,i), format = '(F0.4,F10.4,F10.4,F10.4)'
                    endfor
                endif else begin
                    printf, lun1, '# Col 1: p(BA), Col 2: FWHM single (microHz).', format = '(A0)'
                    for k=0, cp.n_peak_test-1 do begin
                        printf, lun1, p_BA(k,i),fwhm_detection_fit(k,i), format = '(F0.4,F10.4)'
                    endfor
                endelse
                free_lun, lun1
                    
                if info.save_test_files eq 0 then begin
                    ; Remove test folders and prior files if required
                    
                    file_delete,test_dirs(0,tested_peak_indices(i)),/recursive
                    file_delete,test_dirs(1,tested_peak_indices(i)),/recursive
                    file_delete,prior_filenames(0,tested_peak_indices(i))
                    file_delete,prior_filenames(1,tested_peak_indices(i))
                    file_delete,data_range_filenames(0,tested_peak_indices(i))
                    file_delete,data_range_filenames(1,tested_peak_indices(i))

                    if best_dnu lt cp.dnu_rg then begin
                        file_delete,test_dirs(2,tested_peak_indices(i)),/recursive
                        file_delete,prior_filenames(2,tested_peak_indices(i))
                        file_delete,data_range_filenames(2,tested_peak_indices(i))
                    endif
                endif
            endfor
        endif
            
        for i=0, n_elements(dipole_indices(tested_peak_indices))-1 do begin
            if best_dnu lt cp.dnu_rg then begin
                readcol,probability_filenames(tested_peak_indices(i)),p_BA,p_CA,p_CB,format='F,F,F,x',/silent,comment='#'
            endif else begin
                readcol,probability_filenames(tested_peak_indices(i)),p_BA,format='F,x',/silent,comment='#'
            endelse
            
            ; Consider the maximum probabilities among the different runs. Approximate up to third decimal digit.
            
            max_p_BA = float(round(max(p_BA)*1.d3)/1.d3)
            detection_probability(dipole_indices(tested_peak_indices(i))) = max_p_BA
            
            if best_dnu lt cp.dnu_rg then begin
                max_p_CA = float(round(max(p_CA)*1.d3)/1.d3)
                max_p_CB = float(round(max(p_CB)*1.d3)/1.d3)
                
                ; For a RG star take the maximum probability between the two (Lorentzian and Sinc^2) to deem the peak as 
                ; detected. In this way one makes sure to incorporate the result from the best model possible,
                ; between the two tested.
               
                if max_p_CB gt 0.5 then begin
                    detection_probability(dipole_indices(tested_peak_indices(i))) = max_p_CA
                    sinc_profile_flag(dipole_indices(tested_peak_indices(i))) = 1
                endif else begin
                    detection_probability(dipole_indices(tested_peak_indices(i))) = max_p_BA
                endelse
            endif

            if info.print_on_screen eq 1 then begin
                print,''
                print,' Peak detection test for frequency peak #' + strcompress(string(tested_peak_indices(i)),/remove_all) +' at ' + $
                    strcompress(string(freq1(dipole_indices(tested_peak_indices(i))),format='(F0.2)'),/remove_all) + ' muHz'
                print,' P (Lorentzian vs Only Background): ',max_p_BA

                if best_dnu lt cp.dnu_rg then begin
                    print,' P (Sinc^2 vs Only Background): ',max_p_CA
                    print,' P (Sinc^2 vs Lorentzian): ',max_p_CB
                endif

                print,' '
            endif
        endfor
        
        if info.print_on_screen eq 1 then begin
            print,' Peak Detection Test for candidate l=1 mode(s) completed.' 
            print,' --------------------------------------------------------'
            print,' '
        endif
    endif


    ; -------------------------------------------------------------------------------------------------------------------
    ; Peak rotation and duplet test for all significant candidate l=1 modes
    ; -------------------------------------------------------------------------------------------------------------------
    ; Update the number of dipole modes in the chunk by considering only the significant peaks. 
    ; For the significant peaks (if any) perform a rotation test to check whether the peaks are individuals or multiplets. If the rotation
    ; test is positive, then also save the rotational splitting associated with the fit, so that the individual
    ; frequency components can be reconstructed. For RG stars, add the duplet test, to check whether the peak is actually constituted by
    ; two adjacent mixed modes of different g-mode order n_g (i.e. not a rotational splitting effect).
    ; Perform the test only if explicitly requested through the input configuring parameter file.

    detected_dipole_indices = where(angular_degree eq 1 and (detection_probability ge cp.detection_probability_threshold or detection_probability eq -99.0) and $ 
        sinc_profile_flag ne 1)
  
    if cp.rotation_test_activated eq 1 then begin    
        if detected_dipole_indices(0) ne -1 then begin
            n_dipole_chunk = n_elements(detected_dipole_indices)
            test_dirs = strarr(3,n_elements(detected_dipole_indices))
            prior_filenames = strarr(3,n_elements(detected_dipole_indices))
            data_range_filenames = strarr(3,n_elements(detected_dipole_indices))
            probability_filenames = strarr(n_elements(detected_dipole_indices))
            run_test = strarr(1)
            flag_run_test = 0

            for j=0, n_elements(detected_dipole_indices)-1 do begin
                dipole_index = detected_dipole_indices(j)
               
                ; For each peak consider the maximum data range between the range and the division  

                left_bound = max([range_maximum(0,dipole_index),divisions_maximum(0,dipole_index)])
                right_bound = max([range_maximum(1,dipole_index),divisions_maximum(1,dipole_index)])

                peak_number = strcompress(string(dipole_index),/remove_all)
                test_name = run + '_' + peak_number
                
                test_nameE = test_name + 'E'           ; Lorentzian profile, fixed background
                test_nameF = test_name + 'F'           ; Lorentzian profile split by rotation (considered as l=1), fixed background
                test_nameG = test_name + 'G'           ; Two Lorentzian profiles, fixed background
                test_dirE = star_dir + info.pb_subdir + '/' + test_nameE
                test_dirF = star_dir + info.pb_subdir + '/' + test_nameF
                test_dirG = star_dir + info.pb_subdir + '/' + test_nameG
                
                rotation_probability_filename = star_dir + info.pb_subdir + '/rotationProbability_' + test_name + '.txt'
                
                filenameE = star_dir + info.pb_subdir + '/' + info.prior_filename + '_' + test_nameE + '.txt'
                filenameF = star_dir + info.pb_subdir + '/' + info.prior_filename + '_' + test_nameF + '.txt'
                filenameG = star_dir + info.pb_subdir + '/' + info.prior_filename + '_' + test_nameG + '.txt'

                filename_rangeE = star_dir + info.pb_subdir + '/frequencyRange_' + test_nameE + '.txt'
                filename_rangeF = star_dir + info.pb_subdir + '/frequencyRange_' + test_nameF + '.txt'
                filename_rangeG = star_dir + info.pb_subdir + '/frequencyRange_' + test_nameG + '.txt'

                test_dirs(*,j) = [test_dirE,test_dirF,test_dirG]
                probability_filenames(j) = rotation_probability_filename
                prior_filenames(*,j) = [filenameE,filenameF,filenameG]
                data_range_filenames(*,j) = [filename_rangeE,filename_rangeF,filename_rangeG]
               
                if file_test(rotation_probability_filename) eq 0 then begin
                    ; For RG and SG stars, adopt a narrower FWHM prior for those l=1 modes outside the l=3 region.
                    ; This will speed up the computation of the peak test and improve the actual fits of the 
                    ; rotational multiplets.
                    ; Also, if the star is a MS, make the prior on the frequency centroid as narrow as possible.
                    
                    right_fwhm = fwhm_radial_fit*cp.fwhm_magnification_factor_radial
                    left_fwhm = cp.fwhm_lower_bound

                    freq_prior = [freq1(dipole_index) - freq_sig1(dipole_index),freq1(dipole_index) + freq_sig1(dipole_index)]

                    if (best_dnu lt cp.dnu_sg and teff lt cp.teff_sg) then begin
                        freq_prior = [range_maximum(0,dipole_index),range_maximum(1,dipole_index)]
                        
                        if (freq1(dipole_index) ge octupole_freq_upper) or (freq1(dipole_index) le octupole_freq_lower) then begin
                            if best_dnu lt cp.dnu_threshold then begin
                                right_fwhm = fwhm_radial_fit
                    
                                ; Do not exceed the width given by the frequency range of the peak if the star is a RG or late SG, because l=1 modes are narrow. 
                
                                range_peak = abs(range_maximum(1,freq_index) - range_maximum(0,freq_index))
                                if right_fwhm gt range_peak then right_fwhm = range_peak
                            endif else begin
                                right_fwhm = fwhm_radial_fit*cp.fwhm_magnification_factor_dipole
                            endelse
                        endif
                    endif
                    
                    tmp_peak = where(freq ge left_bound and freq le right_bound)
                    bg_peak = mean(bg_level_local(tmp_peak))
                    response_peak = mean((sin(!pi/2. * freq(tmp_peak)/nyq) / (!pi/2. * freq(tmp_peak)/nyq))^2)

                    ; For the amplitude estimate of the rotational multiplet incorporate the sqrt(2) factor to take into account a case with i=90 degrees
                    
                    amplitude = sqrt((abs(max(spsd(tmp_peak)) - bg_peak)/response_peak)*!pi*right_fwhm)    

                    ; Set up priors for the peak test profile
                    
                    data_freq_boundaries = [left_bound,right_bound]
                    freq_duplet_prior = [range_maximum(0,dipole_index),range_maximum(1,dipole_index)]
                    duplet_split_prior = [freqbin*2,abs(freq_duplet_prior(1)-freq_duplet_prior(0))]

                    amplitude_prior = [0.0,amplitude]
                    fwhm_prior = [left_fwhm,right_fwhm]

                    ; Divide the frequency prior range into nearly three parts

                    rot_split = abs(freq_prior(1) - freq_prior(0))/cp.rot_split_scaling
                    rot_split_prior = [freqbin*2,rot_split]
 
                     ; If the frequency resolution is comparable to the expected separation among the fine-structure peaks, then skip the test. 

                    if duplet_split_prior(0) ge duplet_split_prior(1)*0.5 then begin
                        continue
                    endif
                   
                    if rot_split_prior(0) ge rot_split_prior(1)*0.5 then begin
                        continue
                    endif

                    cosi_prior = [cp.cosi_prior_lower,cp.cosi_prior_upper]
                    
                    boundariesE = [freq_prior,amplitude_prior,fwhm_prior]
                    boundariesF = [freq_prior,amplitude_prior,fwhm_prior,rot_split_prior,cosi_prior]
                   
                    write_diamonds_data_range,data_range_filenames(0,j),data_freq_boundaries
                    write_diamonds_data_range,data_range_filenames(1,j),data_freq_boundaries
                    write_diamonds_prior,prior_filenames(0,j),boundariesE
                    write_diamonds_prior,prior_filenames(1,j),boundariesF
                   
                    if best_dnu lt cp.dnu_rg then begin
                        boundariesG = [freq_duplet_prior,amplitude_prior,fwhm_prior,duplet_split_prior,amplitude_prior,fwhm_prior]
                        
                        write_diamonds_data_range,data_range_filenames(2,j),data_freq_boundaries
                        write_diamonds_prior,prior_filenames(2,j),boundariesG
                        
                        run_test = [run_test,test_nameE,test_nameF,test_nameG]
                    endif else begin
                        run_test = [run_test,test_nameE,test_nameF]
                    endelse
                    flag_run_test++
                endif
            endfor

            if info.print_on_screen eq 1 then begin
                print,' '
                print,' ---------------------------------------------'
                print,' Peak Rotation Test for l=1 mode(s).' 
            endif

            if flag_run_test gt 0 then begin
                run_test = run_test(1:*)
            
                if info.print_on_screen eq 1 then begin
                    print,' A total of '+strcompress(string(n_elements(run_test)),/remove_all)+' tests must be performed.' 
                endif
                
                p_FE = fltarr(cp.n_peak_test,n_elements(detected_dipole_indices))
                fwhm_rotation_fit = fltarr(cp.n_peak_test,n_elements(detected_dipole_indices))
                rot_split_fit = fltarr(cp.n_peak_test,n_elements(detected_dipole_indices))
                cosi_fit = fltarr(cp.n_peak_test,n_elements(detected_dipole_indices))
                central_freq_fit = fltarr(cp.n_peak_test,n_elements(detected_dipole_indices))

                if best_dnu lt cp.dnu_rg then begin
                    p_GE = fltarr(cp.n_peak_test,n_elements(detected_dipole_indices))
                    p_GF = fltarr(cp.n_peak_test,n_elements(detected_dipole_indices))
                    left_freq_fit = fltarr(cp.n_peak_test,n_elements(detected_dipole_indices))
                    right_freq_fit = fltarr(cp.n_peak_test,n_elements(detected_dipole_indices))
                    left_fwhm_fit = fltarr(cp.n_peak_test,n_elements(detected_dipole_indices))
                    right_fwhm_fit = fltarr(cp.n_peak_test,n_elements(detected_dipole_indices))
                endif
            
                for k=0, cp.n_peak_test-1 do begin
                    peakbagging_parameters = { subdir:       info.pb_subdir,   $
                                               run:          run_test,         $
                                               background:   background_name,  $
                                               fwhm:         -1.0,             $
                                               filename_run: test_name         $
                                             }
                    flag_computation_completed = run_peakbagging(catalog_id,star_id,peakbagging_parameters,1,0,0)
                    
                    if info.print_on_screen eq 1 then begin
                        print,' Iteration #'+strcompress(string(k),/remove_all)+' completed.' 
                    endif
                    
                    for j=0, n_elements(detected_dipole_indices)-1 do begin
                        ; Make sure that the given peak has been tested.

                        if file_test(test_dirs(0,j) + '/peakbagging_evidenceInformation.txt') ne 0 then begin
                            readcol,test_dirs(0,j) + '/peakbagging_evidenceInformation.txt',ln_evidE,format='D,x,x',/silent
                            readcol,test_dirs(1,j) + '/peakbagging_evidenceInformation.txt',ln_evidF,format='D,x,x',/silent
                           
                            if (ln_evidE ge ln_evidF) then begin
                                ln_evidEF = ln_evidE + alog(1.+exp(ln_evidF-ln_evidE)) 
                            endif else begin
                                ln_evidEF = ln_evidF + alog(1.+exp(ln_evidE-ln_evidF))
                            endelse
                            
                            p_FE(k,j) = exp(ln_evidF - ln_evidEF)

                            ; Compute the estimate for FWHM from model E.
                            
                            readcol,test_dirs(0,j) + '/peakbagging_parameter002.txt',par_fwhm,format='D',/silent
                            readcol,test_dirs(0,j) + '/peakbagging_posteriorDistribution.txt',post,format='D',/silent
                            
                            post /= max(post)
                            fwhm_rotation_fit(k,j) = total(par_fwhm*post)/total(post)
                            
                            ; Compute the estimates for nu0, rot_split and cosi from model F.
                            
                            readcol,test_dirs(1,j) + '/peakbagging_parameter000.txt',par_central_freq,format='D',/silent
                            readcol,test_dirs(1,j) + '/peakbagging_parameter003.txt',par_rot_split,format='D',/silent
                            readcol,test_dirs(1,j) + '/peakbagging_parameter004.txt',par_cosi,format='D',/silent
                            readcol,test_dirs(1,j) + '/peakbagging_posteriorDistribution.txt',post,format='D',/silent

                            ; In this case apply an additional trick. Consider only the last 100 nested iterations to get a first
                            ; estimate of the frequency centroid, and compute a standard deviation of the entire sampling. 
                            ; Then recompute the parameter mean by using the sampling located only within 1-sigma of the first estimate.
                            ; This will improve the estimation of the frequency centroid and of the rotational splitting in case
                            ; multiple solutions have been found during the nested sampling process.

                            post /= max(post)
                            rot_split_fit(k,j) = total(par_rot_split*post)/total(post)
                            cosi_fit(k,j) = total(par_cosi*post)/total(post)
                            central_freq_fit(k,j) = total(par_central_freq*post)/total(post)
                            
                            if best_dnu lt cp.dnu_rg then begin
                                readcol, test_dirs(2,j) + '/peakbagging_evidenceInformation.txt',ln_evidG,format='D,x,x',/silent
                                
                                if (ln_evidG ge ln_evidF) then begin
                                    ln_evidGF = ln_evidG + alog(1.+exp(ln_evidF-ln_evidG)) 
                                endif else begin
                                    ln_evidGF = ln_evidF + alog(1.+exp(ln_evidG-ln_evidF))
                                endelse
                                
                                if (ln_evidG ge ln_evidE) then begin
                                    ln_evidGE = ln_evidG + alog(1.+exp(ln_evidE-ln_evidG)) 
                                endif else begin
                                    ln_evidGE = ln_evidE + alog(1.+exp(ln_evidG-ln_evidE))
                                endelse
                                
                                p_GE(k,j) = exp(ln_evidG - ln_evidGE)
                                p_GF(k,j) = exp(ln_evidG - ln_evidGF)
                            
                                ; Compute the estimates for the duplet frequencies and FWHM from model G
                                
                                readcol,test_dirs(2,j) + '/peakbagging_parameter000.txt',par_left_freq,format='D',/silent
                                readcol,test_dirs(2,j) + '/peakbagging_parameter003.txt',par_duplet_split,format='D',/silent
                                readcol,test_dirs(2,j) + '/peakbagging_parameter002.txt',par_left_fwhm,format='D',/silent
                                readcol,test_dirs(2,j) + '/peakbagging_parameter005.txt',par_right_fwhm,format='D',/silent
                                readcol,test_dirs(2,j) + '/peakbagging_posteriorDistribution.txt',post,format='D',/silent

                                post /= max(post)

                                left_freq_fit(k,j) = total(par_left_freq*post)/total(post)
                                right_freq_fit(k,j) = total(par_duplet_split*post)/total(post) + left_freq_fit(k,j)
                                left_fwhm_fit(k,j) = total(par_left_fwhm*post)/total(post)
                                right_fwhm_fit(k,j) = total(par_right_fwhm*post)/total(post)
                            endif
                        endif
                    endfor
                endfor
                    
                for j=0, n_elements(detected_dipole_indices)-1 do begin
                    get_lun, lun1
                    openw, lun1, probability_filenames(j)
                    printf, lun1, '# Rotation probabilities for frequency peak '+  $
                                  strcompress(string(freq1(detected_dipole_indices(j)),format='(F0.2)'),/remove_all) + ' muHz', format = '(A0)'
                    printf, lun1, '# Each line corresponds to a different run of the same test.', format = '(A0)'
                    
                    if best_dnu lt cp.dnu_rg then begin
                        printf, lun1, '# Col 1: p(FE), Col 2: p(GE), Col 3: p(GF), Col 4: FWHM single (microHz), Col 5: Left freq duplet (microHz) ',  $
                                format='(A0)'
                        printf, lun1, '# Col 6: Right freq duplet (microHz), Col 7: Left FWHM (microHz), Col 8: Right FWHM (microHz) ',format='(A0)'
                        printf, lun1, '# Col 9: Central freq (microHz), Col 10: rotational splitting (microHz), Col 11: cos i', format = '(A0)'
                    
                        for k=0, cp.n_peak_test-1 do begin
                            printf, lun1, p_FE(k,j), p_GE(k,j), p_GF(k,j), fwhm_rotation_fit(k,j),$
                                          left_freq_fit(k,j), right_freq_fit(k,j), left_fwhm_fit(k,j), right_fwhm_fit(k,j), $
                                          central_freq_fit(k,j), rot_split_fit(k,j), cosi_fit(k,j), $
                                          format = '(F0.4,F10.4,F10.4,F10.4,F10.4,F10.4,F10.4,F10.4,F10.4,F10.4,F10.4)'
                        endfor
                    endif else begin
                        printf, lun1, '# Col 1: p(FE), Col 2: FWHM single (microHz), Col 3: Central freq (microHz), Col 4: rotational splitting (microHz), Col 5: cos i', $
                                format = '(A0)'
                    
                        for k=0, cp.n_peak_test-1 do begin
                            printf, lun1, p_FE(k,j), fwhm_rotation_fit(k,j), central_freq_fit(k,j), rot_split_fit(k,j), cosi_fit(k,j), $
                                          format = '(F0.4,F10.4,F10.4,F10.4,F10.4)'
                        endfor
                    endelse
                    
                    free_lun, lun1
                    
                    if info.save_test_files eq 0 then begin
                        ; Remove test folders and prior files if required
                       
                        if file_test(test_dirs(0,j),/directory) ne 0 then begin 
                            file_delete,test_dirs(0,j),/recursive
                            file_delete,test_dirs(1,j),/recursive
                            file_delete,data_range_filenames(0,j)
                            file_delete,data_range_filenames(1,j)
                            file_delete,prior_filenames(0,j)
                            file_delete,prior_filenames(1,j)
                            
                            if best_dnu lt cp.dnu_rg then begin
                                file_delete,test_dirs(2,j),/recursive
                                file_delete,data_range_filenames(2,j)
                                file_delete,prior_filenames(2,j)
                            endif
                        endif
                    endif
                endfor
            endif
            
            for j=0, n_elements(detected_dipole_indices)-1 do begin
                if file_test(probability_filenames(j)) ne 0 then begin
                    if best_dnu lt cp.dnu_rg then begin
                        readcol,probability_filenames(j),p_FE,p_GE,p_GF,format='F,F,F,x,x,x,x,x,x,x,x',/silent,comment='#'
                        max_p_GE = float(round(max(p_GE)*1.d3)/1.d3)
                        max_p_GF = float(round(max(p_GF)*1.d3)/1.d3)
                        duplicity_probability(detected_dipole_indices(j)) = max_p_GE
                    endif else begin
                        readcol,probability_filenames(j),p_FE,format='F,x,x,x,x',/silent,comment='#'
                    endelse
                    
                    ; Consider the maximum probabilities among the different runs. Approximate up to second decimal digit.
                    
                    max_p_FE = float(round(max(p_FE)*1.d3)/1.d3)
                    rotation_probability(detected_dipole_indices(j)) = max_p_FE
                    
                    if info.print_on_screen eq 1 then begin
                        print,''
                        print,' Peak rotation and duplet test for frequency peak #'+ strcompress(string(detected_dipole_indices(j)),/remove_all) + ' at ' + $
                            strcompress(string(freq1(detected_dipole_indices(j)),format='(F0.2)'),/remove_all) + ' muHz'
                        print,' P (Rotation vs No rotation): ',max_p_FE
                    
                        if best_dnu lt cp.dnu_rg then begin
                            print,' P (Duplet vs No Rotation): ',max_p_GE
                            print,' P (Duplet vs Rotation): ',max_p_GF
                        endif
                        print,' '
                    endif
                endif
            endfor
            
            if info.print_on_screen eq 1 then begin
                print,' '
                print,' Peak Rotation Test for l=1 mode(s) completed.' 
                print,' ---------------------------------------------'
                print,' '
            endif
        endif else begin
            n_dipole_chunk = 0
        endelse
    endif else begin
        if info.print_on_screen eq 1 then begin
            print,' '
            print,' --------------------------------------------------------'
            print,' Peak Rotation Test Disabled.'
            print,' --------------------------------------------------------'
            print,' '
        endif 
    endelse
    
    ; If at least one candidate l=1 mode is detected, make sure that the n_dipole_chunk is activated in case it is not (i.e. if there is no dipole
    ; mode found from the global modality)
    
    if detected_dipole_indices(0) ne -1 and n_dipole_chunk eq 0 then begin 
        n_dipole_chunk = 1
        freq_dipole = freq_radial_chunk - best_dnu/2. - d01
    endif
endif else begin
    ; If here, then no dipole peaks were found in the chunk from the multi-modal fit. 
    
    n_dipole_chunk = 0
endelse

largest_octupole_fwhm = 0.0

if n_dipole_chunk ne 0 then begin
    ; -------------------------------------------------------------------------------------------------------------------
    ; CASE 1: One or more significant dipole modes are found in the chunk 
    ; -------------------------------------------------------------------------------------------------------------------
    
    ; -------------------------------------------------------------------------------------------------------------------
    ; Evolutionary stage: RG, SG, MS
    ; -------------------------------------------------------------------------------------------------------------------
    ; In general, search for an l=3 only if more than one l=1 mode is detected. This is to give priority to l=1 mode identification over l=3.
    ; If only one l=1 mode is detected, perform the l=3 search only in the case of a RG star or a depressed dipole star (which can be a SG too). 
    
    if (n_elements(detected_dipole_indices) gt 1) or (n_elements(detected_dipole_indices) eq 1 and best_dnu le cp.dnu_rg) or (flag_depressed_dipole eq 1) then begin
        ; Check whether an l=3 is present. To do so attempt to find at least one local maximum in the l=3 region as 
        ; computed from the asymptotic relation. The search region depends on the evolutionary stage of the star.
        ; Apply the additional condition that if the mode is flagged as a sinc^2 profile, then it is excluded from the candidate octupole mode list.

        detected_octupole_index = where(freq1 le octupole_freq_upper and freq1 ge octupole_freq_lower and $
        (detection_probability ge cp.detection_probability_threshold or detection_probability eq -99.0) and (sinc_profile_flag ne 1))
        
        if detected_octupole_index(0) ne -1 then begin
            ; In this case there is at least one significant peak in this l=3 region.
            ; Therefore check whether one of the l=3 candidates is a true l=3 mode, or it is favoring the scenario of a l=1 mixed mode.
            ; If the sinc^2 is not favored, then verify that the detection favors a single Lorentzian peak (i.e. both the rotational splitting and the duplicity are not present). 
            ; Finally use the FWHM from the single Lorentzian peak test of the candidate l=3 peak to assess whether it is an actual l=3 or a mixed mode peak.
            ; Do so for each peak found in the l=3 range from the asymptotic relation.
           
            ; Obtain the largest FWHM of the set of modes. This will be used as discriminant in case more than one candidate l=3 is present.
            
            fwhm_octupole_fit = fltarr(n_elements(detected_octupole_index))
            run_subdir = run + '_octupole_fwhm'

            if info.print_on_screen eq 1 then begin
                print,''
                print,' Fitting Linewidth of candidate octupole modes of the chunk.'
            endif

            for j=0, n_elements(detected_octupole_index)-1 do begin
                octupole_index = detected_octupole_index(j)
                peak_number = strcompress(string(octupole_index),/remove_all)
                run_names = run_subdir + strcompress(string(indgen(cp.n_fwhm_fit)),/remove_all) + '_' + peak_number
           
                if (file_test(star_dir + info.pb_subdir + '/' + run_names(j) + '/peakbagging_computationParameters.txt') eq 0 or keyword_set(force)) then begin
                    right_bound = max([range_maximum(1,octupole_index),divisions_maximum(1,octupole_index)])
                    left_bound = max([range_maximum(0,octupole_index),divisions_maximum(0,octupole_index)])
                    
                    tmp0 = where(freq ge left_bound and freq le right_bound)
                    freq0 = freq(tmp0)
                    spsd0 = spsd(tmp0)
                    bg_peak = mean(bg_level_local(tmp0))
                    response_peak = mean((sin(!pi/2. * freq0/nyq) / (!pi/2. * freq0/nyq))^2)
                    amplitude_octupole = sqrt((abs(spsd_maximum(octupole_index) - bg_peak)/response_peak)*!pi*fwhm_radial_fit*cp.fwhm_magnification_factor_octupole)
                    freq_prior_octupole = [range_maximum(0,octupole_index),range_maximum(1,octupole_index)]
                    amplitude_prior_octupole = [0.0,amplitude_octupole]
                    data_freq_boundaries = [left_bound,right_bound]

                    data_range_filenames = strarr(cp.n_fwhm_fit)
                    prior_filenames = strarr(cp.n_fwhm_fit)

                    ; Allow for a very narrow FWHM
                    
                    left_fwhm = cp.fwhm_lower_bound
                    fwhm_prior = [left_fwhm,fwhm_radial_fit*cp.fwhm_magnification_factor_octupole]
                    boundaries = [freq_prior_octupole,amplitude_prior_octupole,fwhm_prior]

                    for k=0, cp.n_fwhm_fit-1 do begin
                        data_range_filenames(k) = star_dir + info.pb_subdir + '/frequencyRange_' + run_names(k) + '.txt'
                        prior_filenames(k) = star_dir + info.pb_subdir + '/' + info.prior_filename + '_' + run_names(k) + '.txt'
                        
                        write_diamonds_data_range,data_range_filenames(k),data_freq_boundaries
                        write_diamonds_prior,prior_filenames(k),boundaries
                    endfor

                    peakbagging_parameters = { subdir:          info.pb_subdir,     $
                                               run:             run_names,          $
                                               background:      background_name,    $
                                               fwhm:            -1.0,               $
                                               filename_run:    run_subdir          $
                                             }

                    flag_computation_completed_array = run_peakbagging(catalog_id,star_id,peakbagging_parameters,2,0,0)

                    if info.save_test_files ne 1 then begin
                        for k=0, cp.n_fwhm_fit-1 do begin
                            file_delete,data_range_filenames(k)
                            file_delete,prior_filenames(k)
                        endfor
                    endif
                endif

                fwhm_octupole_fit_array = fltarr(cp.n_fwhm_fit)

                for k=0, cp.n_fwhm_fit-1 do begin
                    spawn,'ls -1 ' + star_dir + info.pb_subdir + '/' + run_names(k) + '/peakbagging_parameter0*txt',filenames
                    
                    fwhm_parameter = '002'
                    
                    ; Read the sampled FWHM of the octupole mode.
                    
                    readcol,star_dir + info.pb_subdir + '/' + run_names(k) + '/peakbagging_parameter' + fwhm_parameter + '.txt',par_fwhm0,format='D',/silent
                    
                    ; Load the posterior samples to compute Bayesian mean estimate for FWHM_0
                    
                    readcol,star_dir + info.pb_subdir + '/' + run_names(k) + '/peakbagging_posteriorDistribution.txt',post,format='D',/silent
                    
                    ; Compute parameter estimate, using a weighted average
                    
                    post /= max(post)
                    fwhm_octupole_fit_array(k) = total(par_fwhm0*post)/total(post)
                endfor

                ; Save the largest FWHM of the set as well as the corresponding absolute frequency index

                fwhm_octupole_fit(j) = median(fwhm_octupole_fit_array)
                if fwhm_octupole_fit(j) ge largest_octupole_fwhm then begin
                    largest_octupole_fwhm = fwhm_octupole_fit(j)
                    best_octupole_index = octupole_index
                endif
            endfor

            if info.print_on_screen eq 1 then begin
                print,' Maximum FWHM (candidate l=3) = '+strcompress(string(largest_octupole_fwhm,format='(F0.3)'),/remove_all) + ' muHz'
                print,''
            endif
               
            ; Now proceed by performing the octupole mode test on the peak that has the largest FWHM if it has been detected.

            peak_number = strcompress(string(best_octupole_index),/remove_all)
            test_name = run + '_' + peak_number
            
            ; If not, now verify whether the peak is split by rotation. If the rotation test is activated, given that all peaks considered here are
            ; significant, by default the rotation test is performed. In the case of a RG, check also for the duplicity.
           
            flag_octupole_fwhm_test = 0

            if cp.rotation_test_activated eq 1 then begin 
                rotation_probability_filename = star_dir + info.pb_subdir + '/rotationProbability_' + test_name + '.txt'
                
                if file_test(rotation_probability_filename) eq 1 then begin
                    if best_dnu lt cp.dnu_rg then begin
                        ; This is the case of RG stars

                        readcol,rotation_probability_filename,p_FE,p_GE,p_GF,freq_left,freq_right,fwhm_left,fwhm_right,  $
                        format='F,F,F,x,F,F,F,F,x,x,x',/silent,comment='#'
                    
                        max_p_FE = float(round(max(p_FE)*1.d3)/1.d3)
                        max_p_GE = float(round(max(p_GE)*1.d3)/1.d3)
                        max_p_GF = float(round(max(p_GF)*1.d3)/1.d3)
                    
                        if (max_p_FE ge cp.rotation_probability_threshold) or (max_p_GE ge cp.duplicity_probability_threshold) then begin
                            ; Here either rotation or duplicity (or both) were detected.
                            ; Then check whether the duplicity is detected over rotation.
                            
                            if (max_p_GF gt 0.5) then begin
                                ; Here the peak is considered as a duplet. In this case verify whether any of the two peaks of the duplet has a
                                ; FWHM comparable to that of the adjacent l=0 mode. If this is true, then consider such peak as a
                                ; l=3, and split it up from the other peak, which is therefore flagged as a l=1 mixed mode only at the end of the CHUNK modality.
                                
                                fwhm_duplet = [mean(fwhm_left),mean(fwhm_right)]

                                if (fwhm_duplet(0) ge fwhm_radial_fit*cp.fwhm_octupole_radial_fraction) or (fwhm_duplet(1) ge fwhm_radial_fit*cp.fwhm_octupole_radial_fraction) then begin
                                    ; At least one of the peaks in the duplet has a FWHM comparable to that of the radial mode.
                                    ; Then consider the mode a potential octupole.
                                    ; Check whether the frequency resolution is enough to be able to distinguish between a l=1 mixed mode and a l=3.
                                    ; If the peak will pass the test, it will be split up later.
                           
                                    flag_octupole_fwhm_test = 1
                                endif
                            endif
                        endif else begin
                            ; No rotation and no duplicity were detected. Then the peak is considered as a single Lorentzian profile.
                            ; Finally check its FWHM against that of the adjacent l=0 mode.
                            ; If the fitted linewidth is comparable to that of the adjacent radial mode, then consider the mode a potential octupole.
                            ; Check whether the frequency resolution is enough to be able to distinguish between a l=1 mixed mode and a l=3.
                            
                            flag_octupole_fwhm_test = 1
                        endelse
                    endif else begin
                        ; This is the case of SG and MS stars

                        readcol,rotation_probability_filename,p_FE,format='F,x,x,x,x',/silent,comment='#'
                        max_p_FE = float(round(max(p_FE)*1.d3)/1.d3)
                        
                        if max_p_FE lt cp.rotation_probability_threshold then begin
                            ; Rotation is not detected. Then the peak is considered as a single Lorentzian profile.
                            ; Finally check its FWHM against that of the adjacent l=0 mode.
                            ; If the fitted linewidth is comparable to that of the adjacent radial mode, then consider the mode as a potential octupole.
                            ; Check whether the frequency resolution is enough to be able to distinguish between a l=1 mixed mode and a l=3.

                            flag_octupole_fwhm_test = 1
                        endif
                    endelse
                endif else begin
                    ; Here the rotation test was not performed because the frequency resolution is not sufficiently high. 
                    ; Then look for an l=3 based solely on the fitted linewidth from the detection tests.
               
                    flag_octupole_fwhm_test = 1
                endelse
            endif else begin
                ; Here the rotation test was not performed because deactivated by the user. 
                ; Then look for an l=3 based solely on the fitted linewidth from the detection tests.
               
                flag_octupole_fwhm_test = 1
            endelse
      
            ; If required, assess the FWHM and ASEF of the l=3 candidate against those of the l=0 mode. 

            if flag_octupole_fwhm_test eq 1 then begin 
                ; Impose that l=3, if present, has a low ASEF value (< 3/4 ASEF maximum)

                asef_threshold = (dp.isla.max_nested_it + dp.isla.n_live) * cp.asef_threshold_fraction
                
                if asef_maximum(best_octupole_index) lt asef_threshold then begin
                    if largest_octupole_fwhm ge fwhm_radial_fit*cp.fwhm_octupole_radial_fraction then begin 
                        angular_degree(best_octupole_index) = 3
                        order_number(best_octupole_index) = enn_radial - 2
                    endif
                endif
            endif
        endif
    endif

    ; If the star is a SG, check that the detected l=1 modes are not closer to one another than Dnu/X.
    ; First update the detected dipole indices, in case some manipulation has occurred from the inspection of l=3 modes.

    detected_dipole_indices = where(angular_degree eq 1 and (detection_probability ge cp.detection_probability_threshold or detection_probability eq -99.0))

    if best_dnu gt cp.dnu_rg and n_elements(detected_dipole_indices) gt 1 then begin
        for k=0, n_elements(detected_dipole_indices)-2 do begin
            ; Check the separation between two consecutive l=1 candidate mixed modes

            freq_left_dipole = freq1(detected_dipole_indices(k))
            freq_right_dipole = freq1(detected_dipole_indices(k+1))

            if abs(freq_left_dipole - freq_right_dipole) lt best_dnu^2/cp.dnu_mixed_modes_separation_scaling/cp.dnu_rg then begin
                ; Since the two peaks are too close, these peaks are likely to originate from a single multiplet. 
                ; Hence take as detected only the one with the largest sampling counts

                min_peak_sampling_counts = min([sampling_counts(detected_dipole_indices(k)),sampling_counts(detected_dipole_indices(k+1))],worst_peak_index)
                worst_peak_index = detected_dipole_indices(k + worst_peak_index)

                detection_probability(worst_peak_index) = 0.0
            endif
        endfor
    endif

    ; If the star is a MS or a high-luminosity RGB star, make sure that there is only one dipole mode for this chunk.
   
    if ((best_dnu ge cp.dnu_rg and best_dnu lt cp.dnu_sg and teff ge cp.teff_sg) or best_dnu ge cp.dnu_sg or best_dnu le cp.dnu_tip) then begin
        detected_dipole_indices = where(angular_degree eq 1 and (detection_probability ge cp.detection_probability_threshold or detection_probability eq -99.0))
        ; If more than one dipole mode is found, then pick up the one with the best combination of SPSD maximum, ASEF maximum, sampling counts and 
        ; frequency position with respect to the global frequency of the dipole mode.
        
        if n_elements(detected_dipole_indices) gt 1 then begin
            freq_diff = abs(freq1(detected_dipole_indices) - freq_dipole)
            freq_weights = 1.d0/freq_diff
            freq_weights /= total(freq_weights)
            freq_ww = freq_weights/max(freq_weights)
            asef_ww = asef_weights(detected_dipole_indices)/max(asef_weights(detected_dipole_indices))
            spsd_ww = spsd_weights(detected_dipole_indices)/max(spsd_weights(detected_dipole_indices))
            sampling_ww = alog(sampling_counts(detected_dipole_indices))/max(alog(sampling_counts(detected_dipole_indices)))
            total_ww = cp.weight_freq_fraction*freq_ww + cp.weight_asef_fraction*asef_ww + cp.weight_spsd_fraction*spsd_ww + cp.weight_sampling_fraction*sampling_ww
            total_ww /= total(total_ww)
            max_ww = max(total_ww,index)
            dipole_index = detected_dipole_indices(index)
            freq_dipole_chunk = freq1(dipole_index)
            bad_dipole_indices = where(detected_dipole_indices ne dipole_index)

            ; Flag the bad modes as undetected to remove them from the list of good frequencies.
            
            detection_probability(detected_dipole_indices(bad_dipole_indices)) = 0.0
        endif
    endif
endif

if info.print_on_screen eq 1 then begin 
    if n_dipole_chunk ne 0 then begin
        loadct,39,/silent
        if octupole_freq_asymp gt min(par_hist) then begin
            arrow,octupole_freq_asymp,max(asef_hist)*0.9,octupole_freq_asymp,max(asef_hist)*0.8,/data,thick=3,/solid,hsize=pp.numax_arrow_hsize,color=pp.color_interval
        endif
        oplot,[octupole_freq_lower,octupole_freq_upper],[max(asef_hist)*0.75,max(asef_hist)*0.75],color=pp.color_interval,thick=pp.thick_border_asef
        oplot,[octupole_freq_lower,octupole_freq_lower],[max(asef_hist)*0.75,max(asef_hist)*0.72],color=pp.color_interval,thick=pp.thick_border_asef
        oplot,[octupole_freq_upper,octupole_freq_upper],[max(asef_hist)*0.75,max(asef_hist)*0.72],color=pp.color_interval,thick=pp.thick_border_asef
    endif
endif

; -------------------------------------------------------------------------------------------------------------------
; Produce the final list of detected frequency peaks, and attached mode identification, to be stored as output and
; overplotted on the PSD of the star.
; -------------------------------------------------------------------------------------------------------------------
; First, select only modes that are significant with respect to noise.

detected_indices = where(detection_probability ge cp.detection_probability_threshold or detection_probability eq -99.0)
local_dp = 0.

if detected_indices(0) ne -1 then begin
    freq1_final = freq1(detected_indices)
    freq_sig1_final = freq_sig1(detected_indices)
    asef_maximum_final = asef_maximum(detected_indices)
    spsd_maximum_final = spsd_maximum(detected_indices)
    sampling_counts_final = sampling_counts(detected_indices)
    angular_degree_final = angular_degree(detected_indices)
    order_number_final = order_number(detected_indices)
    range_maximum_final = range_maximum(*,detected_indices)
    divisions_maximum_final = divisions_maximum(*,detected_indices)
    detection_probability_final = detection_probability(detected_indices)
    rotation_probability_final = rotation_probability(detected_indices)
    duplicity_probability_final = duplicity_probability(detected_indices)
    blending_profile_flag_final = blending_profile_flag(detected_indices)
    sinc_profile_flag_final = sinc_profile_flag(detected_indices)

    ; If possible, evaluate a local value for the observed period spacing of dipolar mixed modes

    tmp_mixed_dipole = where(angular_degree_final eq 1 and freq1_final lt freq_quadrupole_chunk)
    n_mixed_dipole = n_elements(tmp_mixed_dipole)

    if n_mixed_dipole gt 1 then begin
        ; At least two dipolar modes are found. Therefore DP can be estimated.

        freq_mixed_dipole = freq1_final(tmp_mixed_dipole)
        dp_array = fltarr(n_mixed_dipole-1)

        for i=0,n_mixed_dipole-2 do begin
            dp_array(i) = abs(1./freq_mixed_dipole(i) - 1./freq_mixed_dipole(i+1)) * 1.d6   ; units in seconds
        endfor

        local_dp = mean(dp_array)
    endif

    ; Add information about the rotational components, if these will become available.
    
    azimuthal_number_final = fltarr(n_elements(detected_indices)) - 99.0
    cosi_final = fltarr(n_elements(detected_indices)) - 99.0

    ; Second, check whether there are modes split by rotation or that are duplets. Apply separate probability thresholds for each case.
    
    rotation_duplicity_indices = where(rotation_probability_final ge cp.rotation_probability_threshold or $
                                       duplicity_probability_final ge cp.duplicity_probability_threshold)

    if rotation_duplicity_indices(0) ne -1 then begin
        ; First store the values of the detected frequencies with rotation and/or duplicity in order to split them up and apply 
        ; a mode identification to each component.
    
        freq_local_array = dblarr(n_elements(rotation_duplicity_indices))

        for j=0, n_elements(rotation_duplicity_indices)-1 do begin
            local_index = rotation_duplicity_indices(j)
            freq_local_array(j) = freq1_final(local_index)
        endfor

        for j=0, n_elements(rotation_duplicity_indices)-1 do begin
            freq_local = freq_local_array(j)
            local_index = closest(freq_local,freq1_final)
            freq_index = closest(freq_local,freq1)
            peak_number = strcompress(string(freq_index),/remove_all)
            test_name = run + '_' + peak_number
            rotation_probability_filename = star_dir + info.pb_subdir + '/rotationProbability_' + test_name + '.txt'
            flag_duplicity_split = 0

            if file_test(rotation_probability_filename) eq 1 then begin
                if (duplicity_probability_final(local_index) ne -99.0) then begin
                    readcol,rotation_probability_filename,p_FE,p_GE,p_GF,left_freq,right_freq,left_fwhm,right_fwhm,central_freq,rot_split,cosi,  $
                            format='F,F,F,x,F,F,F,F,F,F,F',/silent,comment='#'
                    
                    max_p_FE = float(round(max(p_FE)*1.d3)/1.d3)
                    max_p_GE = float(round(max(p_GE)*1.d3)/1.d3)
                    max_p_GF = float(round(max(p_GF)*1.d3)/1.d3)
               
                    if (max_p_GE ge cp.duplicity_probability_threshold and max_p_GF gt 0.5) then begin
                        flag_duplicity_split = 1
                        freq_duplet = [mean(left_freq),mean(right_freq)]
                    endif
                endif else begin
                    readcol,rotation_probability_filename,p_FE,central_freq,rot_split,cosi,format='F,x,F,F,F',/silent,comment='#'
                    max_p_FE = float(round(max(p_FE)*1.d3)/1.d3)
                endelse
               
                external_indices = where(freq1_final ne freq_local)

                if flag_duplicity_split eq 1 then begin
                    ; In this case the duplicity is the favored scenario.

                    range_midpoint = (freq_duplet(1) + freq_duplet(0)) / 2.0
                    range_maximum_new = [[range_maximum_final(0,local_index),range_midpoint],           $
                                           [range_midpoint,range_maximum_final(1,local_index)]]
                    divisions_maximum_new = [[divisions_maximum_final(0,local_index),range_midpoint],   $
                                           [range_midpoint,divisions_maximum_final(1,local_index)]]

                    left_left_bound = range_maximum_new(0,0)
                    left_right_bound = range_maximum_new(1,0)
                    right_left_bound = range_maximum_new(0,1)
                    right_right_bound = range_maximum_new(1,1)
                   
                    ; Left peak
                   
                    tmp_range = where(par0 lt left_right_bound and par0 ge left_left_bound)
                    tmp_freq_range = where(freq lt left_right_bound and freq ge left_left_bound)
                    tmp_hist_range = where(par_hist lt left_right_bound and par_hist ge left_left_bound)
                    
                    if tmp_freq_range(0) eq -1 then continue
                    spsd_range = spsd(tmp_freq_range)
                    spsd_maximum_left = max(spsd_range)

                    if tmp_range(0) ne -1 then begin 
                        par0_range = par0(tmp_range)
                        sampling_counts_left = total(nest_iter(tmp_range))
                        freq_sig_left = sqrt(total((par0_range-freq_duplet(0))^2*tmp_range^2)/total(tmp_range^2))
                    endif else begin
                        right_bound_new = right_right_bound
                        left_bound_new = left_left_bound
                        tmp_range = where(par0 lt right_bound_new and par0 ge left_bound_new)
                        par0_range = par0(tmp_range)
                        sampling_counts_left = total(nest_iter(tmp_range))/2.
                        freq_sig_left = sqrt(total((par0_range-freq_duplet(0))^2*tmp_range^2)/total(tmp_range^2))/sqrt(2.)
                    endelse

                    if tmp_hist_range(0) ne -1 then begin
                        asef_maximum_left = max(asef_hist(tmp_hist_range))
                    endif else begin
                        right_bound_new = right_right_bound
                        left_bound_new = left_left_bound
                        tmp_hist_range = where(par_hist lt right_bound_new and par_hist ge left_bound_new)
                        asef_maximum_left = max(asef_hist(tmp_hist_range))
                    endelse

                    ; Right peak
                    
                    tmp_range = where(par0 lt right_right_bound and par0 ge right_left_bound)
                    tmp_freq_range = where(freq lt right_right_bound and freq ge right_left_bound)
                    tmp_hist_range = where(par_hist lt right_right_bound and par_hist ge right_left_bound)
                    
                    if tmp_freq_range(0) eq -1 then continue
                    spsd_range = spsd(tmp_freq_range)
                    spsd_maximum_right = max(spsd_range)
                    
                    if tmp_range(0) ne -1 then begin
                        par0_range = par0(tmp_range)
                        sampling_counts_right = total(nest_iter(tmp_range))
                        freq_sig_right = sqrt(total((par0_range-freq_duplet(1))^2*tmp_range^2)/total(tmp_range^2))
                    endif else begin
                        right_bound_new = right_right_bound
                        left_bound_new = left_left_bound
                        tmp_range = where(par0 lt right_bound_new and par0 ge left_bound_new)
                        par0_range = par0(tmp_range)
                        sampling_counts_right = total(nest_iter(tmp_range))/2.
                        freq_sig_right = sqrt(total((par0_range-freq_duplet(1))^2*tmp_range^2)/total(tmp_range^2))/sqrt(2.)
                    endelse

                    if tmp_hist_range(0) ne -1 then begin
                        asef_maximum_right = max(asef_hist(tmp_hist_range))
                    endif else begin
                        right_bound_new = right_right_bound
                        left_bound_new = left_left_bound
                        tmp_hist_range = where(par_hist lt right_bound_new and par_hist ge left_bound_new)
                        asef_maximum_right = max(asef_hist(tmp_hist_range))
                    endelse
                   
                    ; If the peak corresponding to the duplet is flagged as a l = 3 mode, check which of the two peaks is the l=3, and 
                    ; flag the other one as a l=1. The l=3 will be the one with the largest FWHM of the two.

                    if angular_degree_final(local_index) eq 3 then begin
                        max_duplet_fwhm = max([mean(left_fwhm),mean(right_fwhm)],octupole_index_local)
                        angular_degree_duplet = [1,1]
                        order_number_duplet = [enn_radial-1,enn_radial-1]
                        angular_degree_duplet(octupole_index_local) = 3
                        order_number_duplet(octupole_index_local) = enn_radial-2
                    endif else begin
                        angular_degree_duplet = [1,1]
                        order_number_duplet = [enn_radial-1,enn_radial-1]
                    endelse

                    if external_indices(0) ne -1 then begin
                        range_maximum_final = [[range_maximum_final(*,external_indices)],[range_maximum_new]]
                        divisions_maximum_final = [[divisions_maximum_final(*,external_indices)],[divisions_maximum_new]]
                        freq1_final = [freq1_final(external_indices),freq_duplet]
                        freq_sig1_final = [freq_sig1_final(external_indices),freq_sig_left,freq_sig_right]
                        asef_maximum_final = [asef_maximum_final(external_indices),asef_maximum_left,asef_maximum_right]
                        spsd_maximum_final = [spsd_maximum_final(external_indices),spsd_maximum_left,spsd_maximum_right]
                        sampling_counts_final = [sampling_counts_final(external_indices),sampling_counts_left,sampling_counts_right]

                        angular_degree_final = [angular_degree_final(external_indices),angular_degree_duplet]
                        order_number_final = [order_number_final(external_indices),order_number_duplet]
                        azimuthal_number_final = [azimuthal_number_final(external_indices),-99.0,-99.0]
                        cosi_final = [cosi_final,-99.0,-99.0]
                        
                        ; Assign the duplet the same detection probability, but remove the duplicity and rotation test from it, in order to avoid
                        ; that the peak is split up again later on.
                        
                        detection_probability_final = [detection_probability_final(external_indices),fltarr(2) + detection_probability_final(local_index)]
                        rotation_probability_final = [rotation_probability_final(external_indices),fltarr(2) - 99.0]
                        duplicity_probability_final = [duplicity_probability_final(external_indices),fltarr(2) - 99.0]
                        blending_profile_flag_final = [blending_profile_flag_final(external_indices),intarr(2)]
                        sinc_profile_flag_final = [sinc_profile_flag_final(external_indices),intarr(2)]
                       
                        ; Sort all arrays with increasing frequency order.
                        
                        sorted_indices = sort(freq1_final)
                        range_maximum_final = range_maximum_final(*,sorted_indices)
                        divisions_maximum_final = divisions_maximum_final(*,sorted_indices)
                        freq1_final = freq1_final(sorted_indices)
                        freq_sig1_final = freq_sig1_final(sorted_indices)
                        asef_maximum_final = asef_maximum_final(sorted_indices)
                        spsd_maximum_final = spsd_maximum_final(sorted_indices)
                        sampling_counts_final = sampling_counts_final(sorted_indices)
                        angular_degree_final = angular_degree_final(sorted_indices)
                        order_number_final = order_number_final(sorted_indices)
                        azimuthal_number_final = azimuthal_number_final(sorted_indices)
                        cosi_final = cosi_final(sorted_indices)

                        detection_probability_final = detection_probability_final(sorted_indices)
                        rotation_probability_final = rotation_probability_final(sorted_indices)
                        duplicity_probability_final = duplicity_probability_final(sorted_indices)
                        blending_profile_flag_final = blending_profile_flag_final(sorted_indices)
                        sinc_profile_flag_final = sinc_profile_flag_final(sorted_indices)
                    endif else begin
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
                        
                        ; Assign the duplet the same detection probability, but remove the duplicity and rotation test from it, in order to avoid
                        ; that the peak is split up again later on.
                        
                        detection_probability_final = [fltarr(2) + detection_probability_final(local_index)]
                        rotation_probability_final = [fltarr(2) - 99.0]
                        duplicity_probability_final = [fltarr(2) - 99.0]
                        blending_profile_flag_final = [intarr(2)]
                        sinc_profile_flag_final = [intarr(2)]
                    endelse
                endif else begin
                    ; In this case rotational splitting is the favored scenario.
                 
                    cosi = mean(cosi)
                    freq_triplet = [mean(central_freq)-mean(rot_split),mean(central_freq),mean(central_freq)+mean(rot_split)]
   
                    range_midpoint_left = (freq_triplet(0) + freq_triplet(1)) / 2.0
                    range_boundary_left = freq_triplet(0) - abs(freq_triplet(1)-freq_triplet(0))/2.
                    range_midpoint_right = (freq_triplet(1) + freq_triplet(2)) / 2.0
                    range_boundary_right = freq_triplet(2) + abs(freq_triplet(2)-freq_triplet(1))/2.
                    
                    range_maximum_new = [[range_boundary_left,range_midpoint_left],       $
                                             [range_midpoint_left,range_midpoint_right],                 $
                                             [range_midpoint_right,range_boundary_right]]

                    divisions_maximum_new = [[divisions_maximum_final(0,local_index),range_midpoint_left],       $
                                             [range_midpoint_left,range_midpoint_right],                 $
                                             [range_midpoint_right,divisions_maximum_final(1,local_index)]]

                    left_left_bound = range_maximum_new(0,0)
                    left_right_bound = range_maximum_new(1,0)
                    central_left_bound = range_maximum_new(0,1) 
                    central_right_bound = range_maximum_new(1,1) 
                    right_left_bound = range_maximum_new(0,2)
                    right_right_bound = range_maximum_new(1,2)
 
                    ; For low frequency modes the multi-modal sampling from DIAMONDS may turn out to be insufficient
                    ; to split up the mode in three peaks. If this happens, then take as frequency uncertainty the 
                    ; uncertainty from the un-split peak divided by sqrt(3).
                    
                    ; Left peak
                    
                    tmp_range = where(par0 lt left_right_bound and par0 ge left_left_bound)
                    tmp_freq_range = where(freq lt left_right_bound and freq ge left_left_bound)
                    tmp_hist_range = where(par_hist lt left_right_bound and par_hist ge left_left_bound)
                    
                    if tmp_freq_range(0) eq -1 then continue
                    spsd_range = spsd(tmp_freq_range)
                    spsd_maximum_left = max(spsd_range)

                    if tmp_range(0) ne -1 then begin
                        par0_range = par0(tmp_range)
                        sampling_counts_left = total(nest_iter(tmp_range))
                        freq_sig_left = sqrt(total((par0_range-freq_triplet(0))^2*tmp_range^2)/total(tmp_range^2))
                    endif else begin
                        right_bound_new = right_right_bound
                        left_bound_new = left_left_bound
                        tmp_range = where(par0 lt right_bound_new and par0 ge left_bound_new)
                        par0_range = par0(tmp_range)
                        sampling_counts_left = total(nest_iter(tmp_range))/3.
                        freq_sig_left = sqrt(total((par0_range-freq_triplet(0))^2*tmp_range^2)/total(tmp_range^2))/sqrt(3.)
                    endelse
 
                    if tmp_hist_range(0) ne -1 then begin
                        asef_maximum_left = max(asef_hist(tmp_hist_range))
                    endif else begin
                        right_bound_new = right_right_bound
                        left_bound_new = left_left_bound
                        tmp_hist_range = where(par_hist lt right_bound_new and par_hist ge left_bound_new)
                        asef_maximum_left = max(asef_hist(tmp_hist_range))
                    endelse
 
                    ; Central peak
                    
                    tmp_range = where(par0 lt central_right_bound and par0 ge central_left_bound)
                    tmp_freq_range = where(freq lt central_right_bound and freq ge central_left_bound)
                    tmp_hist_range = where(par_hist lt central_right_bound and par_hist ge central_left_bound)
                    
                    if tmp_freq_range(0) eq -1 then continue
                    spsd_range = spsd(tmp_freq_range)
                    spsd_maximum_central = max(spsd_range)
                    
                    if tmp_range(0) ne -1 then begin
                        par0_range = par0(tmp_range)
                        sampling_counts_central = total(nest_iter(tmp_range))
                        freq_sig_central = sqrt(total((par0_range-freq_triplet(1))^2*tmp_range^2)/total(tmp_range^2))
                    endif else begin
                        right_bound_new = right_right_bound
                        left_bound_new = left_left_bound
                        tmp_range = where(par0 lt right_bound_new and par0 ge left_bound_new)
                        par0_range = par0(tmp_range)
                        sampling_counts_central = total(nest_iter(tmp_range))/3.
                        freq_sig_central = sqrt(total((par0_range-freq_triplet(1))^2*tmp_range^2)/total(tmp_range^2))/sqrt(3.)
                    endelse
                    
                    if tmp_hist_range(0) ne -1 then begin
                        asef_maximum_central = max(asef_hist(tmp_hist_range))
                    endif else begin
                        right_bound_new = right_right_bound
                        left_bound_new = left_left_bound
                        tmp_hist_range = where(par_hist lt right_bound_new and par_hist ge left_bound_new)
                        asef_maximum_central = max(asef_hist(tmp_hist_range))
                    endelse
 
                    ; Right peak
                    
                    tmp_range = where(par0 lt right_right_bound and par0 ge right_left_bound)
                    tmp_freq_range = where(freq lt right_right_bound and freq ge right_left_bound)
                    tmp_hist_range = where(par_hist lt right_right_bound and par_hist ge right_left_bound)
                    par0_range = par0(tmp_range)
                    
                    if tmp_freq_range(0) eq -1 then continue
                    spsd_range = spsd(tmp_freq_range)
                    spsd_maximum_right = max(spsd_range)

                    if tmp_range(0) ne -1 then begin
                        par0_range = par0(tmp_range)
                        sampling_counts_right = total(nest_iter(tmp_range))
                        freq_sig_right = sqrt(total((par0_range-freq_triplet(2))^2*tmp_range^2)/total(tmp_range^2))
                    endif else begin
                        right_bound_new = right_right_bound
                        left_bound_new = left_left_bound
                        tmp_range = where(par0 lt right_bound_new and par0 ge left_bound_new)
                        par0_range = par0(tmp_range)
                        sampling_counts_right = total(nest_iter(tmp_range))/3.
                        freq_sig_right = sqrt(total((par0_range-freq_triplet(2))^2*tmp_range^2)/total(tmp_range^2))/sqrt(3.)
                    endelse
                    
                    if tmp_hist_range(0) ne -1 then begin
                        asef_maximum_right = max(asef_hist(tmp_hist_range))
                    endif else begin
                        right_bound_new = right_right_bound
                        left_bound_new = left_left_bound
                        tmp_hist_range = where(par_hist lt right_bound_new and par_hist ge left_bound_new)
                        asef_maximum_right = max(asef_hist(tmp_hist_range))
                    endelse

                    ; Check whether there are frequencies other than this rotational multiplet 

                    if external_indices(0) ne -1 then begin
                        range_maximum_final = [[range_maximum_final(*,external_indices)],[range_maximum_new]]
                        divisions_maximum_final = [[divisions_maximum_final(*,external_indices)],[divisions_maximum_new]]
                        freq1_final = [freq1_final(external_indices),freq_triplet]
                        freq_sig1_final = [freq_sig1_final(external_indices),freq_sig_left,freq_sig_central,freq_sig_right]
                        asef_maximum_final = [asef_maximum_final(external_indices),asef_maximum_left,asef_maximum_central,asef_maximum_right]
                        spsd_maximum_final = [spsd_maximum_final(external_indices),spsd_maximum_left,spsd_maximum_central,spsd_maximum_right]
                        sampling_counts_final = [sampling_counts_final(external_indices),sampling_counts_left,sampling_counts_central,sampling_counts_right]
                        angular_degree_final = [angular_degree_final(external_indices),1,1,1]
                        order_number_final = [order_number_final(external_indices),enn_radial-1,enn_radial-1,enn_radial-1]
                        azimuthal_number_final = [azimuthal_number_final(external_indices),-1,0,1]
                        cosi_final = [cosi_final(external_indices),cosi,cosi,cosi]
 
                        ; Assign the duplet the same detection probability, as well as the duplicity and rotation probabilities.
                        
                        detection_probability_final = [detection_probability_final(external_indices),fltarr(3) + detection_probability_final(local_index)]
                        rotation_probability_final = [rotation_probability_final(external_indices),fltarr(3) + rotation_probability_final(local_index)]
                        duplicity_probability_final = [duplicity_probability_final(external_indices),fltarr(3) + duplicity_probability_final(local_index)]
                        blending_profile_flag_final = [blending_profile_flag_final(external_indices),intarr(3)]
                        sinc_profile_flag_final = [sinc_profile_flag_final(external_indices),intarr(3)]
 
                        ; Sort all arrays with increasing frequency order.
                        
                        sorted_indices = sort(freq1_final)
                        range_maximum_final = range_maximum_final(*,sorted_indices)
                        divisions_maximum_final = divisions_maximum_final(*,sorted_indices)
                        freq1_final = freq1_final(sorted_indices)
                        freq_sig1_final = freq_sig1_final(sorted_indices)
                        asef_maximum_final = asef_maximum_final(sorted_indices)
                        spsd_maximum_final = spsd_maximum_final(sorted_indices)
                        sampling_counts_final = sampling_counts_final(sorted_indices)
                        angular_degree_final = angular_degree_final(sorted_indices)
                        order_number_final = order_number_final(sorted_indices)
                        azimuthal_number_final = azimuthal_number_final(sorted_indices)
                        cosi_final = cosi_final(sorted_indices)
                        detection_probability_final = detection_probability_final(sorted_indices)
                        rotation_probability_final = rotation_probability_final(sorted_indices)
                        duplicity_probability_final = duplicity_probability_final(sorted_indices)
                        blending_profile_flag_final = blending_profile_flag_final(sorted_indices)
                        sinc_profile_flag_final = sinc_profile_flag_final(sorted_indices)
                    endif else begin
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

                        ; Assign the duplet the same detection probability, but remove the duplicity and rotation test from it, in order to avoid
                        ; that the peak is split up again later on.
                        
                        detection_probability_final = [fltarr(3) + detection_probability_final(local_index)]
                        rotation_probability_final = [fltarr(3) + rotation_probability_final(local_index)]
                        duplicity_probability_final = [fltarr(3) + duplicity_probability_final(local_index)]
                        blending_profile_flag_final = [intarr(3)]
                        sinc_profile_flag_final = [intarr(3)]
                    endelse
                endelse

                ; Check if any of the new peaks of the duplet or triplet is falling within the l=2 mixed mode region. If this is the case,
                ; remove the corresponding peak(s) from the final list.
            
                if flag_quadrupole_found eq 1 then begin 
                    candidate_mixed = where(freq1_final le freq_quadrupole_chunk + d02/cp.d02_scaling_merge_mixed and $
                                            freq1_final ge freq_quadrupole_chunk - d02/cp.d02_scaling_merge_mixed and freq1_final ne freq_quadrupole_chunk, $
                                            complement=good_freq)
               
                    if candidate_mixed(0) ne -1 then begin
                        range_maximum_final = range_maximum_final(*,good_freq)
                        divisions_maximum_final = divisions_maximum_final(*,good_freq)
                        freq1_final = freq1_final(good_freq)
                        freq_sig1_final = freq_sig1_final(good_freq)
                        asef_maximum_final = asef_maximum_final(good_freq)
                        spsd_maximum_final = spsd_maximum_final(good_freq)
                        sampling_counts_final = sampling_counts_final(good_freq)
                        angular_degree_final = angular_degree_final(good_freq)
                        order_number_final = order_number_final(good_freq)
                        azimuthal_number_final = azimuthal_number_final(good_freq)
                        cosi_final = cosi_final(good_freq)
                        
                        detection_probability_final = detection_probability_final(good_freq)
                        rotation_probability_final = rotation_probability_final(good_freq)
                        duplicity_probability_final = duplicity_probability_final(good_freq)
                        blending_profile_flag_final = blending_profile_flag_final(good_freq)
                        sinc_profile_flag_final = sinc_profile_flag_final(good_freq)
                    endif
                endif
            endif
        endfor
    endif

    n_freq_final = n_elements(freq1_final)
endif else begin
    ; If here, then it means that no significant modes were detected in the whole chunk.
    
    n_freq_final = 0
    freq1_final = [0]
    freq_sig1_final = [0]
    angular_degree_final = [-99]
    order_number_final = [-99]
    azimuthal_number_final = [-99]
endelse

if info.print_on_screen eq 1 then begin 
    for i=0, n_freq_final-1 do begin
        loadct,39,/silent
        xyouts,freq1_final(i),max(asef_hist)*1.17,strcompress(string(angular_degree_final(i)),/remove_all),   $
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
; Plot the original PSD and a smoothed version of it by the mean linewidth found.
; Overplot an inset of the PSD if a global approach is used, to show the detail of the mode identification.
; -------------------------------------------------------------------------------------------------------------------
separations = fltarr(n_freq_final + 1)

if info.print_on_screen eq 1 then begin
    ; Update smoothing window using fitted FWHM of the chunk radial mode

    smth_bins = fwhm_radial_fit/1.5/freqbin
    spsd_total = smooth(psd_total,smth_bins,/edge_truncate)
    tmp_good = where(freq_total le max(par_hist) and freq_total ge min(par_hist))
    spsd = spsd_total(tmp_good)
    freq = freq_total(tmp_good)
    psd = psd_total(tmp_good)

    parameters = {  min_freq:       min(par_hist),            $
                    max_freq:       max(par_hist),            $
                    freq_list:      freq1_final,              $
                    freq_sig_list:  freq_sig1_final,          $
                    degree:         angular_degree_final,     $
                    order:          order_number_final,       $
                    azimuthal:      azimuthal_number_final,   $
                    lower_freq:     low_cut_frequency,        $
                    separations:    separations,              $
                    freq_radial:    freq_radial_chunk,        $
                    dnu:            best_dnu,                 $
                    numax:          numax,                    $
                    n_radial:       n_radial_chunk,           $
                    max_psd:        max(psd_total)            $
                 }
        
    plot_psd,freq,psd,spsd,bg_level_local,parameters,position_psd,0
endif


; -------------------------------------------------------------------------------------------------------------------
; Plot logo of the pipeline and summary information of the run
; -------------------------------------------------------------------------------------------------------------------
if info.print_on_screen eq 1 then begin
    if n_freq_final ne 0 then begin
        print,' ------------------------------------------'
        print,'       Frequency      sig             n        l      m'
        
        for i=0, n_freq_final-1 do begin
            print,freq1_final(i),freq_sig1_final(i),order_number_final(i),angular_degree_final(i),azimuthal_number_final(i)
        endfor
        
        print,' ------------------------------------------'
        print,' '
    endif
endif

parameters = { catalog_id:     catalog_id,     $
               star_id:        star_id,        $
               run:            run,            $
               modality:       modality,       $
               numax:          numax,          $
               teff:           teff,           $
               best_dnu:       best_dnu,       $
               best_alpha:     best_alpha,     $
               best_epsi:      best_epsi,      $
               acf_dnu:        acf_dnu,        $
               local_epsi:     local_epsi,     $
               local_d02:      local_d02,      $
               local_dp:       local_dp,       $
               snr:            snr,            $
               fit_linewidth:      fit_linewidth,      $
               fwhm_radial_fit:    fwhm_radial_fit,    $
               avg_fwhm:       avg_fwhm,       $
               upper_height:   upper_height,   $
               threshold_asef: threshold_asef, $
               n_bins:         n_bins,         $
               tolerance:      0,              $
               n_freq:         n_freq,         $
               n_chunks:       n_chunks,       $
               n_radial_chunk: n_radial_chunk  $
             }

plot_summary,parameters,0

if info.save_eps eq 0 then begin
    read_jpeg,info.logo_filename,image
    TV, image, 0.855,0.865,TRUE = 1,/normal
endif else begin
    xyouts,0.987,0.2,sp.copyright_str+' FAMED',orientation=-90,charsize=lp.summary_charsize,charthick=lp.summary_charthick,/normal,color=250
endelse

if info.save_eps eq 1 then begin
    device,/close
    spawn,'open '+filename_star
    set_plot,'x'
endif
   
if info.save_png eq 1 then begin
    write_png, star_dir + info.figs_subdir + '/' + catalog_id + star_id + '_' + info.isla_subdir + '_' + run + '_' + modality + '.PNG', TVRD(/TRUE)
endif

; -------------------------------------------------------------------------------------------------------------------
; Save final outputs
; -------------------------------------------------------------------------------------------------------------------
if info.save_complete_lists eq 1 then begin
    ; Save total list of frequencies, without selection of significant peaks.
    
    get_lun,lun1
    openw,lun1,peakbagging_filename_chunk + run + '_' + modality + '.all.txt'
    printf,lun1,'# Col 1: n, Col 2: l, Col 3: Frequency (microHz), Col 4: 1-sigma Frequency (microHz), ',format='(A0)'
    printf,lun1,'# Col 5: Left frequency range (microHz), Col 6: Right frequency range (microHz), ',format='(A0)'
    printf,lun1,'# Col 7: Left frequency division (microHz), Col 8: Right frequency division (microHz), ',format='(A0)'
    printf,lun1,'# Col 9: ASEF maximum (nested iterations), Col 10: Sampling counts (counts), Col 11: SPSD maximum (ppm^2/microHz), ',format='(A0)'
    printf,lun1,'# Col 12: P (detection), Col 13: P (rotation), Col 14: P (duplicity), Col 15: Peak blending flag, Col 16: Sinc^2 profile flag.',format='(A0)'

    for i=0,n_freq-1 do begin
        printf,lun1,order_number(i),angular_degree(i),freq1(i),freq_sig1(i),range_maximum(0,i),range_maximum(1,i), $
               divisions_maximum(0,i),divisions_maximum(1,i),asef_maximum(i),sampling_counts(i),spsd_maximum(i),  $
               detection_probability(i),rotation_probability(i),duplicity_probability(i),blending_profile_flag(i),sinc_profile_flag(i),  $
               format='(I0,I5,F15.5,F12.5,F15.5,F15.5,F15.5,F15.5,F15.2,I15,E15.5,F12.3,F12.3,F12.3,I10,I10)'
    endfor

    free_lun,lun1
endif

if n_freq_final ne 0 then begin
    ; Save local asymptotic parameters and the final list of frequencies with peak significance and rotation tests applied.
    
    get_lun,lun1
    openw,lun1,peakbagging_filename_chunk + run + '_' + modality + '.txt'
    printf,lun1,'# Local epsilon, Local d02 (microHz), local DeltaP, FWHM radial mode (microHz), FWHM octupole mode (microHz), FWHM multi-modal fit (microHz)',format='(A0)'
    printf,lun1,local_epsi,local_d02,local_dp,fwhm_radial_fit,largest_octupole_fwhm,fit_linewidth,format='(F0.3,F10.3,F10.3,F10.4,F10.4,F10.4)'

    ; Save the individual oscillation frequencies from the global fit, their uncertainties, mode identification
    
    printf,lun1,'# Col 1: n, Col 2: l, Col 3: m, Col 4: Frequency (microHz), Col 5: 1-sigma Frequency (microHz), ',format='(A0)'
    printf,lun1,'# Col 6: Left frequency range (microHz), Col 7: Right frequency range (microHz), ',format='(A0)'
    printf,lun1,'# Col 8: Left frequency division (microHz), Col 9: Right frequency division (microHz), ',format='(A0)'
    printf,lun1,'# Col 10: ASEF maximum (nested iterations), Col 11: Sampling counts (counts), Col 12: SPSD maximum (ppm^2/microHz), Col 13: cosi',format='(A0)'
    printf,lun1,'# Col 14: P (detection), Col 15: P (rotation), Col 16: P (duplicity), Col 17: Peak blending profile flag, Col 18: Sinc^2 profile flag.',format='(A0)'

    for i=0,n_freq_final-1 do begin
        printf,lun1,order_number_final(i),angular_degree_final(i),azimuthal_number_final(i),freq1_final(i),freq_sig1_final(i),range_maximum_final(0,i),range_maximum_final(1,i), $
               divisions_maximum_final(0,i),divisions_maximum_final(1,i),asef_maximum_final(i),sampling_counts_final(i),spsd_maximum_final(i),cosi_final(i),  $
               detection_probability_final(i),rotation_probability_final(i),duplicity_probability_final(i),blending_profile_flag_final(i),sinc_profile_flag_final(i),  $
               format='(I2,I5,F8.1,F15.5,F12.5,F15.5,F15.5,F15.5,F15.5,F15.2,I15,E15.5,F15.4,F12.3,F12.3,F12.3,I10,I10)'
    endfor

    free_lun,lun1
endif
end