function assess_freq_asymptotic,freq_list,enn,ell,dnu,epsilon,alpha,d01,numax,tolerance
; -------------------------------------------------------------------------------------------------------
; Author:  Enrico Corsaro
; e-mail:  enrico.corsaro@inaf.it
; Date:    February 2019
; Place:   Catania, Italy
; Purpose: This function verifies that the input list of frequencies for radial and dipole modes
;          follows the expected asymptotic pattern (as set up by the asymptotic parameters given as input)
;          and discard those frequencies that appear to exceed from the expected asymptotic pattern
;          by an input tolerance value. 
;          
;          WARNING: the input frequency list is expected to include only frequencies of the same angular
;          degree (ell). The tolerance is provided as a % of DeltaNu.
; -------------------------------------------------------------------------------------------------------
n_freq = n_elements(freq_list)
good_freq_index = [0]
tolerance_threshold = dnu*tolerance
good_freq_index = 0
residuals = dblarr(n_freq)

; When assessing the asymptotic position, also consider the curvature term alpha.
; Assume for simplicity that alpha is the same for both l=1 and l=0.

for i=0,n_freq-1 do begin
   if ell eq 0 then begin
       residuals(i) = abs(freq_list(i) - dnu*(enn(i) + epsilon + alpha/2.*(enn(i) - numax/dnu)^2))
   endif else begin
       residuals(i) = abs(freq_list(i) - dnu*(enn(i) + 0.5 + epsilon + alpha/2.*(enn(i) - numax/dnu)^2) + d01)
   endelse
   
   if (residuals(i) le tolerance_threshold) then begin
       good_freq_index = [good_freq_index,i]
   endif
endfor

good_freq_index = good_freq_index(1:*)

; Search for double l=0 peaks within the same radial order

if ell eq 0 then begin
   duplicate_freq_index = [-1]

   for i=0,n_freq-1 do begin
       actual_enn = enn(i)
       tmp_match = where(enn eq actual_enn)

       if n_elements(tmp_match) gt 1 then begin
           asymp_freq = dnu*(actual_enn + epsilon + alpha/2.*(actual_enn - numax/dnu)^2)
           closest_index = closest(asymp_freq,freq_list)
           duplicate_freq_index = [duplicate_freq_index,tmp_match(where(tmp_match ne closest_index))]
       endif
   endfor

   if n_elements(duplicate_freq_index) gt 1 then begin
       duplicate_freq_index = duplicate_freq_index(1:*)
       duplicate_freq_index = duplicate_freq_index(uniq(duplicate_freq_index))
       
       ; Remove these peaks from the list of good frequencies

       good_freq_index_final = [-1]

       for j=0, n_elements(good_freq_index)-1 do begin
           actual_good_freq_index = good_freq_index(j)
           tmp_match = where(duplicate_freq_index eq actual_good_freq_index)

           if tmp_match(0) eq -1 then begin
               good_freq_index_final = [good_freq_index_final, actual_good_freq_index]
           endif
       endfor
       if n_elements(good_freq_index_final) gt 1 then begin
           good_freq_index = good_freq_index_final(1:*)
       endif
   endif
endif

return,good_freq_index
end
