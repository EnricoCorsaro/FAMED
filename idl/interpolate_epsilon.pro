pro interpolate_epsilon,teff,dnu,output_epsi,output_epsi_array,output_dnu_array,output_teff_array,output_fit_array
; -------------------------------------------------------------------------------------------------------
; Author:  Enrico Corsaro
; e-mail:  enrico.corsaro@inaf.it
; Date:    August 2019
; Place:   Catania, Italy
; Purpose: This procedure computes the value of epsilon for either a main-sequence, subgiant star or a 
;          red giant star by exploiting empirical relations linking epsilon with the stellar parameters 
;          Teff, DeltaNu.
; -------------------------------------------------------------------------------------------------------
COMMON CONFIG,cp
COMMON GRAPHIC,pp,lp,sp,ppe,lpe

if dnu gt cp.dnu_threshold then begin
    ; Values from Lund+17 LEGACY
    
    epsi_array = [1.114, 0.911, 1.356,0.988,1.114,1.445,1.325,1.377,1.374,1.077,1.431,1.343,      $
    1.336,1.225,1.006,1.492,0.880,1.319,0.978,1.392,1.054,1.358,1.112,1.368,1.117,1.504,    $
    1.075,1.455,1.547,1.163,1.153,1.158,1.311,1.267,1.517,1.113,1.400,1.444,1.475,1.439,    $
    1.337,1.007,0.958,1.095,1.343,1.045,1.067,1.529,1.139,1.131,1.350,1.106,1.206,1.318,    $
    1.313,1.032,1.275,1.020,0.920,1.516,1.200,1.061,1.437,1.461,1.281,0.928]
  
    teff_array = [6326,6614,6045,6384,6193,5668,6107,5805,5846,6130,5853,6037,6033,6313,6331,    $
    5674,6479,5832,6344,6068,6305,5775,6171,5811,6248,5501,6235,5309,5488,6173,6343,6122,  $
    6067,6143,5719,6246,5873,5677,5270,5852,6302,6400,6538,6278,6047,6253,6321,5457,5860,  $
    6132,5949,6146,6177,5964,6045,6150,6140,6548,6642,5180,6179,6276,5825,5750,5964,6580]*1.d0

    
    tmp_sort = sort(teff_array)
    teff_array = teff_array(tmp_sort)
    epsi_array = epsi_array(tmp_sort)
    result = poly_fit(teff_array,epsi_array,4,yfit=fit_array,/double)
    good_epsi = interpol(fit_array,teff_array,teff,/spline)
    good_epsi = good_epsi(0)
    
    if ~(good_epsi gt 0) then begin
       good_epsi = interpol(fit_array,teff_array,teff)
       good_epsi = good_epsi(0)
    endif
    
    output_teff_array = teff_array
    output_epsi_array = epsi_array
endif else begin
    ; epsilon-DeltaNu relation by Corsaro et al. 2012b
    
    good_epsi = cp.epsilon_offset + cp.epsilon_slope*alog10(dnu)
    dnu_array = findgen(301)/300. * cp.dnu_threshold + 0.15
    fit_array = cp.epsilon_offset + cp.epsilon_slope*alog10(dnu_array)

    output_dnu_array = dnu_array
    output_fit_array = fit_array
endelse

output_epsi = good_epsi
output_fit_array = fit_array

end