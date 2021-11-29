pro make_islands_global,catalog_id,star_id,teff,external=external
; -------------------------------------------------------------------------------------------------------
; Auuthor:     Enrico Corsaro
; e-mail:      enrico.corsaro@inaf.it
; Date:        February 2019
; Place:       Catania, Italy
; Purpose:     Procedure to compute a global multi-modal fit and divide the PSD of the star into chunks 
;              according to the number of radial modes identified. Each chunk will correspond to a 
;              sub-folder of the <isla_subdir> folder, with number labels in increasing frequency order. 
;              In <isla_subdir> will be stored the corresponding set of priors for the oscillation 
;              peaks falling in the given chunk, with same labeing as the individual run folders.
; Usage:       <catalog_id>: string specifying the Catalog name of the star (e.g. KIC, TIC, etc.).
;              <star_id>: string specifying the ID number of the star. 
;              <teff>: a value for the effective temperature of the star. Based on Teff, and nuMax
;              a proper l=0 linewidth estimate will be computed in order to run the fit.
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

; Create storing output directories if not already present

star_dir = info.peakbagging_results_dir + catalog_id + star_id + '/'

if file_test(star_dir + info.isla_subdir,/directory) eq 0 then file_mkdir, star_dir + info.isla_subdir, /noexpand_path
if file_test(star_dir + info.pb_subdir,/directory) eq 0 then file_mkdir, star_dir + info.pb_subdir, /noexpand_path
if file_test(star_dir + info.as_subdir,/directory) eq 0 then file_mkdir, star_dir + info.as_subdir, /noexpand_path
if file_test(star_dir + info.as_subdir + '/data', /directory) eq 0 then file_mkdir,star_dir + info.as_subdir + '/data', /noexpand_path
if file_test(star_dir + info.figs_subdir,/directory) eq 0 then file_mkdir, star_dir + info.figs_subdir,/noexpand_path
if file_test(star_dir + info.summary_subdir,/directory) eq 0 then file_mkdir, star_dir + info.summary_subdir,/noexpand_path

; Read input PSD and global asteroseismic parameters for the entire sample of LMLLRG

readcol,info.peakbagging_data_dir + catalog_id + star_id + '.txt',freq,psd,format='D,D',/silent,comment='#'
readcol,info.peakbagging_results_dir + catalog_id + star_id + '/gaussianEnvelopeParameters.txt',gauss_par,format='D',/silent 

numax = gauss_par(1)
freqbin = freq(1)-freq(0)
dnu = compute_scaling_dnu(numax)
maxpower = max(psd)
maxpower *= 1.1

; Run DIAMONDS in a global multi-modal fit in order to identify the different chunks and estimate global frequencies and the value of DeltaNu

avg_fwhm = mean(get_linewidth([min(freq),max(freq)],teff,numax))
smth_bins = avg_fwhm/freqbin
spsd = smooth(psd,smth_bins,/edge_truncate)

get_lun,lun1
openw,lun1,star_dir + info.isla_subdir + '/' + info.prior_filename + '_' + info.global_subdir + '.txt'
printf,lun1,'#',format='(A0)'
printf,lun1,min(freq),max(freq),format='(F0.5,F15.5)'
printf,lun1,0,mean(spsd),format='(F0.5,F15.5)'
free_lun,lun1

; Evaluate linewidth at nuMax and use it for the global multi-modal fit

if dnu le cp.dnu_threshold then begin
    if dnu le cp.dnu_tip then begin
        linewidth = get_linewidth(numax,teff,numax)/cp.fwhm_global_scaling_tip
    endif else begin
        linewidth = get_linewidth(numax,teff,numax)/cp.fwhm_global_scaling
    endelse
endif else begin
    linewidth = get_linewidth(numax,teff,numax)
endelse

linewidth = strcompress(string(linewidth),/remove_all)

; Perform the multi-modal fit. When performing the global fit, also save the background
; level as an output file

peakbagging_parameters = { subdir:     info.isla_subdir,     $
                           run:        info.global_subdir,   $
                           background: bgp.name,             $
                           fwhm:       linewidth,            $
                           duplet:     0                     $         
                         }

flag_computation_completed = run_peakbagging(catalog_id,star_id,peakbagging_parameters,0,0,1)
end