pro set_peakbagging,catalog_id,star_id,bgp,force=force
; -------------------------------------------------------------------------------------------------------
; Author:  Enrico Corsaro
; e-mail:  enrico.corsaro@inaf.it
; Date:    August 2019
; Place:   Catania, Italy
; Purpose: This procedure verifies that the peakbagging folder is properly set up for running DIAMONDS
;          PeakBagging fits. If the setup is missing, this function will prepare it for you, based on
;          the computation of a background model with DIAMONDS as a previous step. The background
;          parameters have to be provided as input.
;          
;          WARNING: If another fitting method other than DIAMONDS-Background has been used to fit the 
;          background in the input star, then change this routine for properly reading the necessary
;          information to setup the peak bagging folder.
; -------------------------------------------------------------------------------------------------------
COMMON STAR,info
COMMON CONFIG,cp

bg_par = bgp.parameters

if file_test(info.peakbagging_results_dir + catalog_id + star_id,/directory) eq 0 then begin
    file_mkdir, info.peakbagging_results_dir + catalog_id + star_id, /noexpand_path
endif

if (file_test(info.peakbagging_data_dir + catalog_id + star_id + '.txt') eq 0 or keyword_set(force)) then begin
    
    ; Check whether the background fit has to be read from an external source
    if info.external_background_results_dir eq '-99' then begin
        readcol,info.background_data_dir + catalog_id + star_id + '.txt',freq,psd,format='D,D',/silent 
    endif else begin
        readcol,info.external_background_results_dir + catalog_id + star_id + '.txt',freq,psd,format='D,D',/silent 
    endelse

    ; Trim the global PSD of the star, used for the background fit, to a narrower frequency region
    ; centered around nuMax

    numax = bg_par(n_elements(bg_par)-2)
    ;fwhm_env = 0.66*numax^0.88
    ;sig_env_pred = fwhm_env / (2 * sqrt(2*alog(2)))
    sig_env = bg_par(n_elements(bg_par)-1)
    dnu = compute_scaling_dnu(numax)
    
    if dnu le cp.dnu_cl then begin
        if dnu le cp.dnu_tip then begin
            width_factor = cp.n_sigma_envelope_tip
        endif else begin
            width_factor = cp.n_sigma_envelope_cl
        endelse
    endif else begin
        width_factor = cp.n_sigma_envelope
    endelse
    
    if cp.n_dnu_envelope*dnu lt sig_env*width_factor then begin
        env_width = sig_env*width_factor
    endif else begin
        env_width = cp.n_dnu_envelope*dnu
    endelse

    lower_bound = numax - env_width
    upper_bound = numax + env_width
    
    tmp_clipping = where(freq ge lower_bound and freq le upper_bound)
    maxpower = max(psd(tmp_clipping))
    freq_pb = freq(tmp_clipping)
    psd_pb = psd(tmp_clipping)
   
    get_lun,lun1
    openw,lun1,info.peakbagging_data_dir + catalog_id + star_id + '.txt'
    
    for i=0L,n_elements(freq_pb)-1 do begin
        printf,lun1,freq_pb(i),psd_pb(i),format='(F15.10,F20.10)'
    endfor
    
    free_lun,lun1
endif

if (file_test(info.peakbagging_results_dir + catalog_id + star_id + '/backgroundParameters.txt') eq 0 or keyword_set(force)) then begin
    get_lun,lun1
    openw,lun1,info.peakbagging_results_dir + catalog_id + star_id + '/backgroundParameters.txt'
    printf,lun1,'# Configuring parameters for background model of '+catalog_id+star_id,format='(A0)'
    printf,lun1,'# These parameters are the median values of the result derived',format='(A0)'
    printf,lun1,'# from the Background code based on DIAMONDS.',format='(A0)'
    
    for i=0,n_elements(bg_par)-4 do begin
        printf,lun1,bg_par(i),format='(F0.10)'
    endfor
    
    free_lun,lun1
endif

if (file_test(info.peakbagging_results_dir + catalog_id + star_id + '/gaussianEnvelopeParameters.txt') eq 0 or keyword_set(force)) then begin
    get_lun,lun1
    openw,lun1,info.peakbagging_results_dir + catalog_id + star_id + '/gaussianEnvelopeParameters.txt'
    printf,lun1,'# Parameters of the Gaussian envelope fit for background model of '+catalog_id+star_id,format='(A0)'
    printf,lun1,'# Row #1: height',format='(A0)'
    printf,lun1,'# Row #2: numax',format='(A0)'
    printf,lun1,'# Row #3: sigma_env',format='(A0)'
    
    for i=n_elements(bg_par)-3,n_elements(bg_par)-1 do begin
       printf,lun1,bg_par(i),format='(F0.10)'
    endfor
    
    free_lun,lun1
endif

if (file_test(info.peakbagging_results_dir + catalog_id + star_id + '/NyquistFrequency.txt') eq 0 or keyword_set(force)) then begin 
    if info.external_background_results_dir eq '-99' then begin
        spawn,'cp ' + info.background_results_dir + catalog_id + star_id + '/NyquistFrequency.txt ' + info.peakbagging_results_dir + catalog_id + star_id + '/'
    endif else begin
        spawn,'cp ' + info.external_background_results_dir + catalog_id + star_id + '_NyquistFrequency.txt ' + info.peakbagging_results_dir + catalog_id + star_id + '/NyquistFrequency.txt'
    endelse
endif
end
