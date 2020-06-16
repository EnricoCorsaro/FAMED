function compute_scaling_dnu,numax
; -------------------------------------------------------------------------------------------------------
; Author:  Enrico Corsaro
; e-mail:  enrico.corsaro@inaf.it
; Date:    March 2019
; Place:   Catania, Italy
; Purpose: This function computes a value of DeltaNu from an input nuMax as based on literature
;          scaling relations. It could be updated from time to time as soon as new, more accurate,
;          calibration of the relation nuMax-DeltaNu become available.
; -------------------------------------------------------------------------------------------------------
COMMON CONFIG,cp

; nuMax has to be in microHz. Following scaling relations calibrated by Huber et al. 2011

if numax lt cp.numax_threshold then begin
    dnu = cp.numax_coeff_low*numax^cp.numax_exponent_low
endif else begin
    dnu = cp.numax_coeff_high*numax^cp.numax_exponent_high
endelse

return, dnu
end