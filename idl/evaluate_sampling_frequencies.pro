function evaluate_sampling_frequencies,par0,par_hist,freq,spsd,maximum,range_maximum
; -------------------------------------------------------------------------------------------------------
; Author:  Enrico Corsaro
; e-mail:  enrico.corsaro@inaf.it
; Date:    August 2019
; Place:   Catania, Italy
; Purpose: This function computes the frequencies and their uncertainties from the multi-modal fit of
;          DIAMONDS. The extracted frequencies do not necessarily correspond to real oscillation peaks
;          but they correspond to the extracted local maxima from the ASEF histogram. The routine
;          also provides the sampling counts, SPSD local maximum for each frequency, as well as
;          an additional improvement to the definition of the frequency range for each ASEF peak.
; -------------------------------------------------------------------------------------------------------
COMMON CONFIG,cp

n_maxima = n_elements(maximum)
freq1 = dblarr(n_maxima)
freq_sig1 = dblarr(n_maxima)
sampling_counts = fltarr(n_maxima)
spsd_maximum = fltarr(n_maxima)
nest_iter = findgen(n_elements(par0))
range_maximum_old = range_maximum
freqbin = freq(1) - freq(0)

for i=0, n_maxima-1 do begin
    ; Consider the frequency range around the local maximum
    
    iterations = 0
    
    while iterations lt cp.max_iterations_frequency do begin
        upper_bound = range_maximum(1,i)
        lower_bound = range_maximum(0,i)
        tmp_freq_range = where(freq ge lower_bound and freq le upper_bound)
        
        ; Make sure that within each range, at least one frequency point is found.
       
        while tmp_freq_range(0) eq -1 do begin
            upper_bound += freqbin/2.
            lower_bound -= freqbin/2.
            tmp_freq_range = where(freq ge lower_bound and freq le upper_bound)
            range_maximum(1,i) = upper_bound
            range_maximum(0,i) = lower_bound
        endwhile
        
        tmp_range = where(par0 lt upper_bound and par0 ge lower_bound)
        par0_range = par0(tmp_range)
        spsd_range = spsd(tmp_freq_range)
        
        ; Count the total of the nested iteration values falling in this maximum bin. 
        ; Note this is not the number of nested sampling points, but it includes the actual nested iteration value from each point (which can be up to several thousands
        ; for an individual sampling point). In this way it is possible to better weight the local maxima by the actual level of likelihood that they have reached
        ; during the sampling.
        
        sampling_counts(i) = total(nest_iter(tmp_range))
        
        ; Save the maximum smoothed PSD corresponding to the region of this local maximum
        
        spsd_maximum(i) = max(spsd_range)
        
        ; Weighted by nested iteration value
        
        freq1(i) = total(par0_range*tmp_range^2)/total(tmp_range^2)
        freq_sig1(i) = sqrt(total((par0_range-freq1(i))^2*tmp_range^2)/total(tmp_range^2))
        
        if iterations eq cp.max_iterations_frequency-1 then break
        
        ; Improve frequency ranges for computation of uncertainties and evaluation of the number of sampling points
        ; Do not exceed ranges by more than the estimated sigma * X (usually 2).
       
        if upper_bound gt freq1(i)+freq_sig1(i)*cp.max_sigma_range then begin
            range_maximum(1,i) = freq1(i)+freq_sig1(i)*cp.max_sigma_range
        endif

        if lower_bound lt freq1(i)-freq_sig1(i)*cp.max_sigma_range then begin
            range_maximum(0,i) = freq1(i)-freq_sig1(i)*cp.max_sigma_range
        endif
        
        ; Try to make ranges extend by at least sigma * Y times (usually 1) on each side
        
        left_freq = freq1(i)-freq_sig1(i)*cp.min_sigma_range
        right_freq = freq1(i)+freq_sig1(i)*cp.min_sigma_range
        
        if lower_bound gt left_freq then begin
            if i eq 0 then begin
                if left_freq le min(par_hist) then range_maximum(0,i) = min(par_hist)
                if left_freq gt min(par_hist) then range_maximum(0,i) = left_freq
            endif else begin
                if left_freq le range_maximum(1,i-1) then range_maximum(0,i) = range_maximum(1,i-1)
                if left_freq gt range_maximum(1,i-1) then range_maximum(0,i) = left_freq
            endelse
        endif
        
        if upper_bound lt right_freq then begin
            if i eq n_maxima-1 then begin
                if right_freq ge max(par_hist) then range_maximum(1,i) = max(par_hist)
                if right_freq lt max(par_hist) then range_maximum(1,i) = right_freq
            endif else begin
                if right_freq ge range_maximum(0,i+1) then range_maximum(1,i) = range_maximum(0,i+1)
                if right_freq lt range_maximum(0,i+1) then range_maximum(1,i) = right_freq
            endelse
        endif
        
        ; Make sure to also include the original ASEF maximum frequency
        
        if range_maximum(0,i) gt maximum(i) then range_maximum(0,i) = maximum(i)
        if range_maximum(1,i) lt maximum(i) then range_maximum(1,i) = maximum(i)
        
        iterations++
    endwhile
endfor

sampled_estimates = {  freq1:             freq1,             $
                       freq_sig1:         freq_sig1,         $
                       sampling_counts:   sampling_counts,   $
                       spsd_maximum:      spsd_maximum       $
                    }

return, sampled_estimates
end
