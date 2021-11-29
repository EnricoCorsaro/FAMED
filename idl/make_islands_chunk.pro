pro make_islands_chunk,catalog_id,star_id,teff,chunk_id,external=external
; -------------------------------------------------------------------------------------------------------
; Auuthor:     Enrico Corsaro
; e-mail:      enrico.corsaro@inaf.it
; Date:        February 2019
; Place:       Catania, Italy
; Purpose:     Procedure to compute a chunk multi-modal fit in order to retrieve frequencies, their 
;              uncertainties, and mode identification. This procedure cannot be executed if the global
;              modality was not performed as a previous step.  
; Usage:       <catalog_id>: string specifying the Catalog name of the star (e.g. KIC, TIC, etc.).
;              <star_id>: string specifying the ID number of the star. 
;              <teff>: a value for the effective temperature of the star. Based on Teff, and nuMax
;              a proper l=0 linewidth estimate will be computed in order to run the fit.
;              <chunk_number>: an integer specifying the number of the chunk for which the sampling is
;              required. If < 0 or > the total number of chunks, by default it will compute the sampling 
;              for all the chunks identified in the GLOBAL module.
; -------------------------------------------------------------------------------------------------------
COMMON CONFIG,cp
COMMON STAR,info
COMMON DIAMONDS,dp

; Set up the computation only if the routine has not been called from an external procedure

bgp = get_background(catalog_id,star_id)

if ~keyword_set(external) then begin
    setup_computation
    set_peakbagging,catalog_id,star_id,bgp
endif

star_dir = info.peakbagging_results_dir + catalog_id + star_id + '/'
peakbagging_filename_global = info.peakbagging_results_dir + catalog_id + star_id + '/' + info.summary_subdir + '/' $
                              + catalog_id + star_id + info.peakbagging_filename_label + $
                                info.isla_subdir + '_' + info.global_subdir + '_' + 'GLOBAL.txt'

; Read input PSD and global asteroseismic parameters for the entire sample of LMLLRG

readcol,info.peakbagging_data_dir + catalog_id + star_id + '.txt',freq,psd,format='D,D',/silent,comment='#'
readcol,info.peakbagging_results_dir + catalog_id + star_id + '/gaussianEnvelopeParameters.txt',gauss_par,format='D',/silent 
numax = gauss_par(1)
freqbin = freq(1)-freq(0)
dnu = compute_scaling_dnu(numax)
maxpower = max(psd)
maxpower *= 1.1

; Run DIAMONDS to perform the multi-modal fit in the different chunks indicated by the global fit.
; This has to assume that a global multi-modal fit has already been done.

if file_test(peakbagging_filename_global) eq 0 then begin
    print,'Cannot perform multi-modal fit for individual chunks. Missing global fit peakbagging summary file.'
    return
endif else begin
    ; Load global parameters

    readcol,peakbagging_filename_global,acf_dnu,best_dnu,best_epsi,best_alpha,teff,n_chunks,format='x,D,D,D,D,D,I',numline=2,comment='#',/silent
    best_dnu = best_dnu(0)
    best_epsi = best_epsi(0)
    best_alpha = best_alpha(0)
    teff = teff(0)
    n_chunks = n_chunks(0)

    ; Load the frequency positions for each chunk

    readcol,peakbagging_filename_global,chunk_number,freq_left,freq_right,snr,format='I,D,D,F',numline=n_chunks,skipline=3,comment='#',/silent
    
    min_linewidths = fltarr(n_chunks)
    bg_names = strarr(n_chunks)

    if ((chunk_id lt 0) or (chunk_id ge n_chunks)) then begin
        first_it = 0
        last_it = n_chunks-1
        run_labels = indgen(n_chunks)
    endif else begin
        first_it = chunk_id
        last_it = chunk_id
        run_labels = [chunk_id]
    endelse
    
    for i=first_it,last_it do begin
        run_subdir = strcompress(string(i),/remove_all)
   
        ; Evaluate maximum PSD in the given chunk

        freq_index = where(freq le freq_right(i) and freq ge freq_left(i))
        psd_chunk = psd(freq_index)
        freq_chunk = freq(freq_index)
        max_psd_chunk = max(psd_chunk)
        min_freq = min(freq(freq_index))

        if best_dnu le cp.dnu_rg then begin
            ; In case the star is classified as a RG, distinguish among RGB, clump and RGB-tip
            
            if best_dnu le cp.dnu_tip then begin
                min_linewidths(i) = get_linewidth(min_freq,teff,numax)/cp.fwhm_chunk_scaling_tip
                avg_fwhm = freqbin
            endif else begin
                if best_dnu le cp.dnu_cl then begin
                    min_linewidths(i) = get_linewidth(min_freq,teff,numax)/cp.fwhm_chunk_scaling_cl
                endif else begin
                    min_linewidths(i) = get_linewidth(min_freq,teff,numax)/cp.fwhm_chunk_scaling_rg
                endelse
                
                avg_fwhm = min_linewidths(i)*cp.smoothing_fwhm_factor_rg
            endelse
        endif else begin
            if (best_dnu lt cp.dnu_sg and teff lt cp.teff_sg) then begin
                min_linewidths(i) = get_linewidth(mean(freq_chunk),teff,numax)/cp.fwhm_chunk_scaling_sg
            endif else begin
                min_linewidths(i) = get_linewidth(mean(freq_chunk),teff,numax)/cp.fwhm_chunk_scaling_ms
            endelse
            avg_fwhm = mean(get_linewidth([min(freq_chunk),max(freq_chunk)],teff,numax))
        endelse

        bg_names(i) = bgp.name
        
        ; Load the background level of the star

        readcol,star_dir + 'backgroundLevel.txt',bg_level,format='x,D',/silent,comment='#'
        
        ; Evaluate smoothed PSD to obtain more accurate prior height

        smth_bins = avg_fwhm/freqbin
        spsd = smooth(psd,smth_bins,/edge_truncate)
        spsd_chunk = spsd(freq_index)
        bg_level_chunk = bg_level(freq_index)
        
        get_lun,lun1
        openw,lun1,star_dir + info.isla_subdir + '/' + info.prior_filename + '_' + run_subdir + '.txt'
        printf,lun1,'#',format='(A0)'
        
        ; Apply overlapping chunks if the star is a subgiant or redgiant.
        ; This should solve the problem of having mixed modes that happen to fall very close to the
        ; radial mode of the previous chunk.

        if i ne 0 then begin
            if best_dnu lt cp.dnu_rg or (best_dnu lt cp.dnu_sg and teff lt cp.teff_sg) then begin
                printf,lun1,freq_left(i) - best_dnu*cp.dnu_overlap_fraction_rg,freq_right(i),format='(F0.5,F15.5)'
            endif else begin
                printf,lun1,freq_left(i) - best_dnu*cp.dnu_overlap_fraction_ms,freq_right(i),format='(F0.5,F15.5)'
            endelse
        endif else begin
            printf,lun1,freq_left(i),freq_right(i),format='(F0.5,F15.5)'
        endelse 
        
        printf,lun1,0,max(spsd_chunk)-mean(bg_level_chunk),format='(F0.5,F15.5)'
        free_lun,lun1
    endfor

    ; When running a single chunk, select only information related to that chunk

    if n_elements(run_labels) lt n_chunks then begin
        active_chunk_index = where(bg_names ne '')
        bg_names = bg_names(active_chunk_index)
        min_linewidths = min_linewidths(active_chunk_index)
    endif

    peakbagging_parameters = { subdir:     info.isla_subdir,     $
                               run:        run_labels,           $
                               background: bg_names,             $
                               fwhm:       min_linewidths,       $
                               duplet:     0                     $         
                             }

    flag_computation_completed = run_peakbagging(catalog_id,star_id,peakbagging_parameters,0,0,0,/merge)
endelse
end
