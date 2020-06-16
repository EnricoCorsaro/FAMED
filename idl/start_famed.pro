pro start_famed,catalog_id,star_id,teff,fit=fit,global=global,chunk=chunk,echelle=echelle,complete=complete,force=force
COMMON CONFIG,cp
COMMON STAR,info
COMMON DIAMONDS,dp

setup_computation
peakbagging_filename_global = info.peakbagging_results_dir + catalog_id + star_id + '/' + info.summary_subdir + '/' $
                                + catalog_id + star_id + info.peakbagging_filename_label + info.peakbagging_filename_global_label

; Perform GLOBAL modality

if keyword_set(global) then begin
    if info.print_on_screen eq 1 then begin
        print,'----------------------------'
        print,' Performing GLOBAL modality.'
        print,'----------------------------'
    endif
   
    bgp = get_background(catalog_id,star_id)
    if keyword_set(force) then begin
        set_peakbagging,catalog_id,star_id,bgp,/force
    endif else begin
        set_peakbagging,catalog_id,star_id,bgp
    endelse 

    readcol,info.peakbagging_results_dir + catalog_id + star_id + '/gaussianEnvelopeParameters.txt',gauss_par,format='D',/silent 
    numax = gauss_par(1)
    dnu = compute_scaling_dnu(numax)
    threshold_asef = cp.threshold_asef_global
    run_subdir = info.global_subdir
   
    if dnu lt cp.dnu_rg or (dnu lt cp.dnu_sg and teff lt cp.teff_sg) then begin
        if dnu gt cp.dnu_rg then begin
            tolerance = cp.skim_frequency_tolerance_sg
        endif else begin
            if dnu le cp.dnu_tip then begin
                tolerance = cp.skim_frequency_tolerance_tip
            endif else begin
                tolerance = cp.skim_frequency_tolerance_rg
            endelse
        endelse
    endif else begin
        tolerance = cp.skim_frequency_tolerance_ms
    endelse
   
    if keyword_set(fit) then begin
        ; Obtain the global modality fit
        make_islands_global,catalog_id,star_id,teff
    endif

    ; Find the actual frequencies
    if keyword_set(force) then begin
        find_islands_global,catalog_id,star_id,threshold_asef,tolerance,teff,/force
    endif else begin
        find_islands_global,catalog_id,star_id,threshold_asef,tolerance,teff
    endelse
endif

; Perform CHUNK modality

if keyword_set(chunk) then begin
    if info.print_on_screen eq 1 then begin
        print,'---------------------------'
        print,' Performing CHUNK modality.'
        print,'---------------------------'
    endif
    
    if keyword_set(fit) then begin
        ; Obtain the chunk modality fits
        make_islands_chunk,catalog_id,star_id,teff
    endif

    ; Load global parameters
    readcol,peakbagging_filename_global,best_dnu,n_chunks,format='x,x,D,x,x,x,I',numline=2,comment='#',/silent
    best_dnu = best_dnu(0)
    n_chunks = n_chunks(0)
    
    readcol,peakbagging_filename_global,snr,format='x,x,x,F',numline=n_chunks,skipline=3,comment='#',/silent
    
    ; Order chunks by decreasing SNR
    tmp_sort = reverse(sort(snr))
    
    ; Find the actual frequencies
    for j=0, n_elements(snr)-1 do begin
        chunk_number = tmp_sort(j)
        if best_dnu lt cp.dnu_rg then begin
            threshold_asef = cp.threshold_asef_chunk_rg
        endif else begin
            if (best_dnu lt cp.dnu_sg and teff lt cp.teff_sg) then begin
                threshold_asef = cp.threshold_asef_chunk_sg
            endif else begin
                threshold_asef = cp.threshold_asef_chunk_ms
            endelse
        endelse
    
        run_index = strcompress(string(chunk_number),/remove_all)
        
        if info.print_on_screen eq 1 then begin
            print,' Analyzing CHUNK '+run_index
            print,' ...'
            print,' ...'
            print,' '
        endif
    
        if keyword_set(force) then begin
            find_islands_chunk,catalog_id,star_id,run_index,threshold_asef,teff,/force
        endif else begin
            find_islands_chunk,catalog_id,star_id,run_index,threshold_asef,teff
        endelse
    endfor
endif

; Perform ECHELLE modality
if keyword_set(echelle) then begin
    if info.print_on_screen eq 1 then begin
        print,'-----------------------------'
        print,' Performing ECHELLE modality.'
        print,' (Currently not available)'
        print,'-----------------------------'
    endif
endif

; Perform COMPLETE modality
if keyword_set(complete) then begin
    if info.print_on_screen eq 1 then begin
        print,'-------------------------------'
        print,' Performing COMPLETE modality. '
        print,' (Currently not available). '
        print,'-------------------------------'
    endif
endif

end
