function compute_asef,parameter,distribution,n_bins
; -------------------------------------------------------------------------------------------------------
; Author:  Enrico Corsaro
; e-mail:  enrico.corsaro@inaf.it
; Date:    August 2019
; Place:   Catania, Italy
; Purpose: This function computes an Averaged Shifted Envelope Function (ASEF) from an input sampling 
;          distribution. It requires that the number of bins of the distribution is provided as input.
; -------------------------------------------------------------------------------------------------------

; Sort by increasing parameter value

tmp = sort(parameter)
param = parameter(tmp)
distr = distribution(tmp)
min_par = min(param)
max_par = max(param)
binwidth = (max_par-min_par)/n_bins*1.d0

; Now perform average shifting for ASH

param_finalrebin = fltarr(n_bins)
distr_finalrebin = fltarr(n_bins)
n_sim = 20    ; Number of combined histograms
simsize = binwidth/n_sim

for k=0, n_sim-1 do begin
    distr_rebin  = fltarr(n_sim, n_bins)
    param_rebin = fltarr(n_sim, n_bins)
    
    for j=0L, n_bins-1 do begin
        upper_value = min_par + (j+1)*binwidth + k*simsize
        lower_value = min_par + j*binwidth + k*simsize
        tmp = where(param ge lower_value and param lt upper_value)
        
        if tmp(0) ne -1 then begin
            param_rebin(k,j) = min_par + j*binwidth + k*simsize + 0.5*binwidth
            distr_rebin(k,j) = max(distr(tmp))
        endif else begin
            param_rebin(k,j) = min_par + j*binwidth + k*simsize + 0.5*binwidth
            distr_rebin(k,j) = 0.0d
        endelse
    endfor
   
    param_finalrebin += param_rebin(k,*)
    distr_finalrebin += distr_rebin(k,*)
endfor
   
param_finalrebin /= n_sim*1.d0
distr_finalrebin /= n_sim*1.d0
param = param_finalrebin(0:n_bins-1)
distr = distr_finalrebin
distr = distr(0:n_bins-1)

ash = { x:     param,   $
        y:     distr    $
      }

return, ash
end