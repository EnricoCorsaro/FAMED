pro start_famed,catalog_id,star_id,teff,fit=fit,global=global,chunk=chunk,echelle=echelle,complete=complete,force=force
COMMON CONFIG,cp
COMMON STAR,info
COMMON DIAMONDS,dp

setup_computation
peakbagging_filename_global = info.peakbagging_results_dir + catalog_id + star_id + '/' + info.summary_subdir + '/' $
                              + catalog_id + star_id + info.peakbagging_filename_label + $
                                info.isla_subdir + '_' + info.global_subdir + '_GLOBAL.txt'

; Perform GLOBAL modality
if keyword_set(global) then begin
    if info.print_on_screen eq 1 then begin
        print,'-------------------------------------------------'
        print,' Performing GLOBAL modality for ' + catalog_id + star_id + '.'
        print,'-------------------------------------------------'
    endif
   
    if keyword_set(fit) then begin
        ; Obtain the global modality fit
        make_islands_global,catalog_id,star_id,teff,/force,/external
    endif

    ; Find the actual frequencies
    if keyword_set(force) then begin
        find_islands_global,catalog_id,star_id,teff,/force,/external
    endif else begin
        find_islands_global,catalog_id,star_id,teff,/external
    endelse
endif

; Perform CHUNK modality
if keyword_set(chunk) then begin
    if info.print_on_screen eq 1 then begin
        print,'-------------------------------------------------'
        print,' Performing CHUNK modality for ' + catalog_id + star_id + '.'
        print,'-------------------------------------------------'
    endif
    
    if keyword_set(fit) then begin
        ; Obtain the chunk modality fits (consider all the chunks)
        make_islands_chunk,catalog_id,star_id,-1,teff,/external
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
            print,' Analyzing CHUNK ' + run_index
            print,' ...'
            print,' ...'
            print,' '
        endif
    
        if keyword_set(force) then begin
            find_islands_chunk,catalog_id,star_id,run_index,threshold_asef,teff,/force,/external
        endif else begin
            find_islands_chunk,catalog_id,star_id,run_index,threshold_asef,teff,/external
        endelse
    endfor
endif

; Perform ECHELLE modality
if keyword_set(echelle) then begin
    if info.print_on_screen eq 1 then begin
        print,'-------------------------------------------------'
        print,' Performing ECHELLE modality for ' + catalog_id + star_id + '.'
        print,' (Currently not available)'
        print,'-------------------------------------------------'
    endif
endif

; Perform COMPLETE modality
if keyword_set(complete) then begin
    if info.print_on_screen eq 1 then begin
        print,'-------------------------------------------------'
        print,' Performing COMPLETE modality for ' + catalog_id + star_id + '.'
        print,' (Currently not available). '
        print,'-------------------------------------------------'
    endif
endif

end
