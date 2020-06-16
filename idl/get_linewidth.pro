function get_linewidth,freq,teff,numax
; -------------------------------------------------------------------------------------------------------
; Author:  Enrico Corsaro
; e-mail:  enrico.corsaro@inaf.it
; Date:    February 2019
; Place:   Catania, Italy
; Purpose: This function provides an empirical estimate of the linewidth of a radial mode, based on the
;          input value of the stellar Teff. It currently incorporates different relations depending on
;          which evolutionary stage the star belongs to. These empirical relations should be updated
;          with more accurate ones calibrated on a larger sample, and spanning a wider range of
;          nuMax and Teff.
; -------------------------------------------------------------------------------------------------------
COMMON CONFIG,cp

; !!!!!!! Need to improve the fit. Possibly replace with one common law, or improve the law for red giants. Corsaro et al. in prep.
; Does not work fine for Teff < 4400
if numax gt cp.numax_threshold then begin
    ; Bilinear fits coefficients from Ball et al. 2018

    offset = [-3.71*1.d0,-7.209*1.d1,-2.266*1.d-1,-2.190*1d3,-5.639*1d-1]
    teff_coeff = [1.073*1.d-3,1.543*1.d-2,5.083*1d-5,4.302*1d-1,1.138*1d-4]
    numax_coeff = [1.883*1.d-4,9.101*1d-4,2.715*1d-6,8.427*1d-1,1.312*1d-4]
    
    parameter = offset + teff_coeff*teff + numax_coeff*numax
    ln_fwhm = parameter(0)*alog(freq/numax) + alog(parameter(1)) + alog(parameter(2))/(1 + (2*alog(freq/parameter(3))/alog(parameter(4)/numax))^2)
endif else begin
    ; Polynomial fit from Corsaro et al. 2015b

    if teff lt 5500. then begin
        if teff gt 4900. then begin
            ln_gamma = alog(51.6371 -0.0295244*teff + 5.61255e-06*teff^2 -3.53773e-10*teff^3)
        endif else begin
            teff_local = 4900.
            ln_gamma = -2.12026         ; Reference constant value 0.12 muHz estimated by Lund+17 for the sample of red giant stars
            ;ln_gamma = alog(51.6371 -0.0295244*teff_local + 5.61255e-06*teff_local^2 -3.53773e-10*teff_local^3)
        endelse
    endif else begin
        ln_gamma = 1463.49 -1.03503*teff + 0.000271565*teff^2 -3.14139e-08*teff^3 + 1.35524e-12*teff^4        ; Total fit from RG to MS
    endelse
    
    ln_fwhm = fltarr(n_elements(freq))
    ln_fwhm += ln_gamma
endelse

return,exp(ln_fwhm)
end
