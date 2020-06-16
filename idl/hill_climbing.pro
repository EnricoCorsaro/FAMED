function hill_climbing,par_hist,asef_hist,threshold,minimum_bin_separation
; -------------------------------------------------------------------------------------------------------
; Author:  Enrico Corsaro
; e-mail:  enrico.corsaro@inaf.it
; Date:    August 2019
; Place:   Catania, Italy
; Purpose: This function locates the local maxima of a given input distribution using a simple 
;          hill-climbing algorithm.
; -------------------------------------------------------------------------------------------------------
n_bins = n_elements(par_hist)
index_maximum = [0]

new_i = 2
for i=2,n_bins-2 do begin
    if i ne new_i then continue
    
    distr_start = asef_hist(i)
    new_i = i + 1
    distr_next = asef_hist(new_i)
  
    ; Check if we are in a descending phase.
    ; If so, update the starting bin at each step, until a local minimum is reached.
    
    while (distr_next lt distr_start) and (new_i lt n_bins-2) do begin
        new_i = new_i + 1
        distr_start = distr_next
        distr_next = asef_hist(new_i)
    endwhile
   
    ; Now check if we are in an ascending phase.
    ; If so, update the final bin until a local maximum is reached.
    
    distr_act = distr_start
    while (distr_next gt distr_act) and (new_i lt n_bins-2) do begin
        new_i = new_i + 1
        distr_act = distr_next      ; Actual value of the distribution
        distr_next = asef_hist(new_i)   ; Next value of the distribution
    endwhile
    
    ; Store the local maximum only if the rising phase has an amount of distribution variation not less than a given threshold.
    
    if (distr_act - distr_start) ge threshold then begin
        index_maximum = [index_maximum,new_i-1]
    endif
endfor

index_maximum = index_maximum(1:*)
good_index_maximum = 0
index_bad = -1

; Remove local maxima if within minimum_bin_separation bins from one another

for i=0, n_elements(index_maximum)-2 do begin
    if i eq index_bad then continue

    index1 = index_maximum(i)      ;where(par_hist eq maximum(i))
    index2 = index_maximum(i+1)    ;where(par_hist eq maximum(i+1))
    
    if (index2 - index1) le minimum_bin_separation then begin
       distr1 = asef_hist(index1)
       distr2 = asef_hist(index2)
       
        ; Consider only the maximum with the best value of the distribution

        if distr1 ge distr2 then begin
            index_bad = i + 1      ; Right peak is bad, don't save it
            good_index_maximum = [good_index_maximum,index_maximum(i)]
        endif
    endif else begin
        good_index_maximum = [good_index_maximum,index_maximum(i)]
    endelse
endfor

; Add the last local maximum of the group if it is good

if index_bad ne n_elements(index_maximum)-1 then begin
    good_index_maximum = [good_index_maximum,index_maximum(n_elements(index_maximum)-1)]
endif

index_maximum = good_index_maximum(1:*)

return,index_maximum
end