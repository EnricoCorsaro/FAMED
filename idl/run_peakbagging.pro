function run_peakbagging,catalog_id,star_id,parameters,flag_peaktest,flag_asymptotic,flag_bglevel,merge=merge
; -------------------------------------------------------------------------------------------------------
; Author:  Enrico Corsaro
; e-mail:  enrico.corsaro@inaf.it
; Date:    August 2019
; Place:   Catania, Italy
; Purpose: This function executes the peakbagging code based on DIAMONDS for different kind of fits, 
;          as determined by the input flags. It requires a set of input parameters to set up the
;          computation, which also do vary depending on what kind of fit is required.
; -------------------------------------------------------------------------------------------------------
COMMON STAR,info
COMMON DIAMONDS,dp

; Set up DIAMONDS configuring parameters

star_dir = info.peakbagging_results_dir + catalog_id + star_id + '/'
nsmc_filename = star_dir + parameters.subdir + '/NSMC_configuringParameters.txt'
xmeans_filename = star_dir + parameters.subdir + '/Xmeans_configuringParameters.txt'
fwhm_minimum = min(parameters.fwhm)

; Uni-modal peak bagging or Peak testing for FWHM fit

if (fwhm_minimum lt 0 and flag_peaktest eq 0) or (flag_peaktest eq 2) then begin
    nsmc_parameters = [dp.pb.n_live, dp.pb.n_live_end, dp.pb.max_draw_attempts, dp.pb.n_initial_it, dp.pb.n_it_same_clust,    $
                       dp.pb.enlarg_fraction, dp.pb.shrinking_rate, dp.pb.termination_factor, dp.pb.max_nested_it]
    xmeans_parameters = [dp.pb.min_ncluster, dp.pb.max_ncluster]
endif

; Peak testing

if flag_peaktest eq 1 then begin
    nsmc_parameters = [dp.pb.n_live_test, dp.pb.n_live_end_test, dp.pb.max_draw_attempts, dp.pb.n_initial_it, dp.pb.n_it_same_clust,    $
                       dp.pb.enlarg_fraction, dp.pb.shrinking_rate, dp.pb.termination_factor_test, dp.pb.max_nested_it]
    xmeans_parameters = [dp.pb.min_ncluster, dp.pb.max_ncluster]
endif

; Multi-modal peak bagging

if fwhm_minimum ge 0 and flag_peaktest eq 0 and flag_asymptotic eq 0 then begin
    nsmc_parameters = [dp.isla.n_live, dp.isla.n_live_end, dp.isla.max_draw_attempts, dp.isla.n_initial_it, dp.isla.n_it_same_clust,    $
                       dp.isla.enlarg_fraction, dp.isla.shrinking_rate, dp.isla.termination_factor, dp.isla.max_nested_it]
    
    if parameters.duplet eq 0 then begin
        xmeans_parameters = [dp.isla.min_ncluster, dp.isla.max_ncluster]
    endif else begin
        xmeans_parameters = [dp.isla.min_ncluster, dp.isla.max_ncluster_duplet]
    endelse
endif

; Sliding pattern fit

if flag_asymptotic eq 1 then begin
    nsmc_parameters = [dp.slid.n_live, dp.slid.n_live_end, dp.slid.max_draw_attempts, dp.slid.n_initial_it, dp.slid.n_it_same_clust,    $
                      dp.slid.enlarg_fraction, dp.slid.shrinking_rate, dp.slid.termination_factor, dp.slid.max_nested_it]
    xmeans_parameters = [dp.slid.min_ncluster, dp.slid.max_ncluster]
endif

get_lun, lun1
openw, lun1, nsmc_filename
printf, lun1, nsmc_parameters, format = '(F0.2)'
free_lun, lun1

get_lun, lun1
openw, lun1, xmeans_filename
printf, lun1, xmeans_parameters, format = '(F0.2)'
free_lun, lun1

; Create subdirectories for each run, if not already present

n_runs = n_elements(parameters.run)

for i=0, n_runs-1 do begin
    if file_test(star_dir + parameters.subdir + '/' + strcompress(string(parameters.run(i)),/remove_all),/directory) eq 0 then $
       file_mkdir,star_dir + parameters.subdir + '/' + strcompress(string(parameters.run(i)),/remove_all),/noexpand_path
endfor

flag_peaktest = strcompress(string(flag_peaktest),/remove_all)
flag_asymptotic = strcompress(string(flag_asymptotic),/remove_all)
flag_bglevel = strcompress(string(flag_bglevel),/remove_all)

cd, info.peakbagging_path

if n_runs gt 1 then begin
    n_chunks = ceil(n_runs*1.0/info.n_threads)
    start_indices = intarr(n_chunks)
    end_indices = intarr(n_chunks)
    start_index = 0

    if n_runs ge info.n_threads then begin
        ; Estimate how many chunks to divide the total run set into
        
        end_index = info.n_threads-1
        start_indices(0) = start_index
        end_indices(0) = end_index

        for j=1, n_chunks-1 do begin
            if j eq n_chunks-1 then begin
                start_index = end_index + 1
                end_index = n_runs - 1
                start_indices(j) = start_index
                end_indices(j) = end_index
            endif else begin
                start_index = end_index + 1
                end_index = end_index + info.n_threads
                start_indices(j) = start_index
                end_indices(j) = end_index
            endelse
        endfor
    endif else begin
        start_indices(0) = start_index
        end_indices(0) = n_runs - 1
    endelse
   
    for j=0, n_chunks-1 do begin
        start_index = start_indices(j)
        end_index = end_indices(j)

        if keyword_set(merge) then begin
            get_lun,lun1
            filename_run = catalog_id + star_id + '_run_indices.txt'
            openw,lun1,filename_run
            printf,lun1,parameters.run(start_index:end_index),format='(I0)'
            free_lun,lun1
            
            get_lun,lun1
            filename_bg = catalog_id + star_id + '_bg_names.txt'
            openw,lun1,filename_bg
            printf,lun1,parameters.background(start_index:end_index),format='(A0)'
            free_lun,lun1
            
            get_lun,lun1
            filename_prior = catalog_id + star_id + '_prior_filenames.txt'
            openw,lun1,filename_prior

            for i=0, end_index-start_index do begin
                printf,lun1,info.prior_filename,format='(A0)'
            endfor
            
            free_lun,lun1
            
            get_lun,lun1
            filename_fwhm = catalog_id + star_id + '_linewidths.txt'
            openw,lun1,filename_fwhm
            printf,lun1,parameters.fwhm(start_index:end_index),format='(F0.3)'
            free_lun,lun1
            
            output_err_filename = catalog_id + star_id + '_' + parameters.subdir + '_parallel.out'
            spawn,'parallel ./peakbagging ::: ' + catalog_id + ' ::: ' + star_id + ' ::: ' + parameters.subdir + $
                   ' :::: ' + filename_run + ' ::::+ ' + filename_bg + ' ::::+ ' + filename_prior + ' ::::+ ' + filename_fwhm + $
                   ' ::: ' + flag_peaktest + ' ::: ' + flag_asymptotic + ' ::: ' + flag_bglevel,output,error
            
            file_delete,filename_run
            file_delete,filename_bg
            file_delete,filename_prior
            file_delete,filename_fwhm
        endif else begin
            get_lun,lun1
            openw,lun1,parameters.filename_run + '.txt'
            printf,lun1,parameters.run(start_index:end_index),format='(A0)'
            free_lun,lun1
            
            spawn,'parallel ./peakbagging ::: ' + catalog_id + ' ::: ' + star_id + ' ::: ' + parameters.subdir + $
                  ' :::: ' + parameters.filename_run + '.txt ::: ' + parameters.background + ' ::: ' + info.prior_filename +      $ 
                  ' ::: ' + strcompress(string(parameters.fwhm),/remove_all) + ' ::: ' + flag_peaktest + ' ::: ' + flag_asymptotic + ' ::: ' + flag_bglevel,output,error

            file_delete,parameters.filename_run + '.txt'
        endelse
    endfor
endif else begin
    spawn,'./peakbagging ' + catalog_id + ' ' + star_id + ' ' + $
          parameters.subdir + ' ' + strcompress(string(parameters.run),/remove_all) + ' ' + parameters.background + ' ' + $
          info.prior_filename + ' ' + strcompress(string(parameters.fwhm),/remove_all) + ' ' + flag_peaktest + ' ' + flag_asymptotic + ' ' + flag_bglevel,output,error
endelse

if (size(error))[0] ne 0 then begin
    if info.print_on_screen eq 1 then begin
        print,' Peakbagging fit could not be completed. '
    endif
endif

cd, info.famed_path

star_dir = info.peakbagging_results_dir + catalog_id + star_id + '/'
flag_computation_completed = intarr(n_elements(parameters.run))

for i=0, n_elements(parameters.run)-1 do begin
    if file_test(star_dir + parameters.subdir + '/' + strcompress(string(parameters.run(i)),/remove_all) + '/peakbagging_parameterSummary.txt') eq 1 then begin
        flag_computation_completed(i) = 1
    endif
endfor

return, flag_computation_completed
end
