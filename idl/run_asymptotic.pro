function run_asymptotic,catalog_id,star_id,parameters,numax,ell
; -------------------------------------------------------------------------------------------------------
; Author:  Enrico Corsaro
; e-mail:  enrico.corsaro@inaf.it
; Date:    August 2019
; Place:   Catania, Italy
; Purpose: This function executes the asymptotic code based on DIAMONDS for different asymptotic relations
;          as given by the input angular degree. It requires a set of input parameters to set up the
;          computation, which also do vary depending on what kind of fit is required.
; -------------------------------------------------------------------------------------------------------
COMMON STAR,info
COMMON DIAMONDS,dp

star_dir = info.peakbagging_results_dir + catalog_id + star_id + '/'
nsmc_filename = star_dir + parameters.subdir + '/NSMC_configuringParameters.txt'
xmeans_filename = star_dir + parameters.subdir + '/Xmeans_configuringParameters.txt'
nsmc_parameters = [dp.pb.n_live, dp.pb.n_live_end, dp.pb.max_draw_attempts, dp.pb.n_initial_it, $
                   dp.pb.n_it_same_clust, dp.pb.enlarg_fraction, dp.pb.shrinking_rate, dp.pb.termination_factor, dp.pb.max_nested_it]
xmeans_parameters = [dp.pb.min_ncluster, dp.pb.max_ncluster]

get_lun, lun1
openw, lun1, nsmc_filename
printf, lun1, nsmc_parameters, format = '(F0.2)'
free_lun, lun1

get_lun, lun1
openw, lun1, xmeans_filename
printf, lun1, xmeans_parameters, format = '(F0.2)'
free_lun, lun1

; Create subdirectory for the run, if not already present. Otherwise delete the folder content, if existing.
if file_test(star_dir + parameters.subdir + '/' + strcompress(string(parameters.run),/remove_all), /directory) eq 0 then begin
    file_mkdir,star_dir + parameters.subdir + '/' + strcompress(string(parameters.run),/remove_all), /noexpand_path
endif else begin
    if file_test(star_dir + parameters.subdir + '/' + strcompress(string(parameters.run),/remove_all) + '/*') eq 1 then begin
        spawn,'ls -1 '+star_dir + parameters.subdir + '/' + strcompress(string(parameters.run),/remove_all) + '/*',delete_filenames
        file_delete,delete_filenames,/recursive
    endif
endelse

cd, info.asymptotic_path

output_err_filename = catalog_id + star_id + '_' + parameters.subdir + '_' + strcompress(string(parameters.run),/remove_all) + '.out'
spawn,'./asymptotic ' + catalog_id + ' ' + star_id + ' ' + $
      parameters.subdir + ' ' + strcompress(string(parameters.run),/remove_all) + ' ' + info.prior_filename + $
      ' ' + strcompress(string(numax),/remove_all) + ' ' + strcompress(string(ell),/remove_all) + ' &> ' + output_err_filename,    $
      error,/stderr

file_delete, output_err_filename

cd, info.famed_path

star_dir = info.peakbagging_results_dir + catalog_id + star_id + '/'
flag_computation_completed = 0

if file_test(star_dir + parameters.subdir + '/' + strcompress(string(parameters.run),/remove_all) + '/peakbagging_parameterSummary.txt') eq 1 then begin
   flag_computation_completed = 1
endif

return, flag_computation_completed
end
