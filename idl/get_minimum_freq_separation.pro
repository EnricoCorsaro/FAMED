function get_minimum_freq_separation,dnu,flag_global
; -------------------------------------------------------------------------------------------------------
; Author:  Enrico Corsaro
; e-mail:  enrico.corsaro@inaf.it
; Date:    August 2019
; Place:   Catania, Italy
; Purpose: This function computes the minimum separation in frequency to define the number of bins in the
;          ASEF histogram. This is based on an estimate of the minimum width required to obtain a good
;          sampling of the actual frequency peaks obtained in the multi-modal fits.
; -------------------------------------------------------------------------------------------------------
COMMON CONFIG,cp

if flag_global eq 1 then begin
    min_separation = dnu/cp.min_sep_scaling_global
endif else begin
    if dnu lt cp.dnu_rg then begin
        ; !!!!!!!! to improve using a function of the radial FWHM such that when FWHM decreases
        ; the number of bins increases
        min_separation = dnu/cp.min_sep_scaling_chunk_rg        
    endif else begin
        d02 = cp.d02_unique_slope*dnu + cp.d02_unique_offset        ; Huber et al. 2010
        min_separation = d02/cp.min_sep_scaling_chunk_ms
    endelse
endelse

return, min_separation
end