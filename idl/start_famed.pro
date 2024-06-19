pro start_famed,catalog_id,star_id,teff,background_run_number=background_run_number,fit=fit,global=global,chunk=chunk,echelle=echelle,complete=complete,force=force
; -------------------------------------------------------------------------------------------------------
; Auuthor:     Enrico Corsaro
; e-mail:      enrico.corsaro@inaf.it
; Date:        February 2019
; Place:       Catania, Italy
; Purpose:     Executes the FAMED pipeline in either of its modalities. The user can activate multiple 
;              modalities if desired, which will be automatically executed in a sequential order.
; Usage:       <catalog_id>: string specifying the Catalog name of the star (e.g. KIC, TIC, etc.).
;              <star_id>: string specifying the ID number of the star. 
;              <teff>: integer or float specifying the effective temperature of the star.
;              <background_run_number>: [OPTIONAL] string specifying the subfolder number containing
;              the results of the background fit obtained with the DIAMONDS+Background code. Omit this
;              input if the default value read by the input configuring parameters has to be adopted.
;              <fit>: [OPTIONAL] activate this keyword to force FAMED computing the multi-modal fit.
;              This keyword is unnecessary when the fitting part is already done and one wants to
;              process the multi-modal fit (e.g. by re-running a different set of input configuring
;              parameters of the pipeline)
;              <global>: activate this keyword to perform the analysis of the GLOBAL module.
;              <chunk>: activate this keyword to perform the analysis of the CHUNK module.
;              <force>: activate this keyword to force FAMED overwriting the results of the analysis 
;              over a previous run (these are not overwritten by default).
; -------------------------------------------------------------------------------------------------------
COMMON CONFIG,cp
COMMON STAR,info
COMMON DIAMONDS,dp

setup_computation
peakbagging_filename_global = info.peakbagging_results_dir + catalog_id + star_id + '/' + info.summary_subdir + '/' $
                              + catalog_id + star_id + info.peakbagging_filename_label + $
                                info.isla_subdir + '_' + info.global_subdir + '_GLOBAL.txt'

if keyword_set(background_run_number) then info.background_run_number = background_run_number

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
    readcol,peakbagging_filename_global,best_dnu,n_chunks,format='x,x,D,x,x,x,x,I',numline=2,comment='#',/silent
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
            print,' '
            print,' '
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

    if info.plot_total_solution eq 1 then begin
        plot_psd_total,catalog_id,star_id,teff
    endif
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
