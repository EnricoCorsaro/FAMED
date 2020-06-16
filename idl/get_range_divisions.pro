pro get_range_divisions,par_hist,asef_hist,index_maximum,range_maximum,divisions_maximum,chunk=chunk
; -------------------------------------------------------------------------------------------------------
; Author:  Enrico Corsaro
; e-mail:  enrico.corsaro@inaf.it
; Date:    August 2019
; Place:   Catania, Italy
; Purpose: This routine obtains both ranges and divisions for each local maximum peak found in the ASEF
;          histogram. The ranges are defined as those boundaries where the ASEF of the local maximum peak
;          starts to rise (left edge) and stops to decrease (right edge), while the divisions are defined
;          as the mid points between adjacent local maxima.
; -------------------------------------------------------------------------------------------------------
COMMON CONFIG,cp
COMMON DIAMONDS,dp

maximum = par_hist(index_maximum)
asef_maximum = asef_hist(index_maximum)
n_maxima = n_elements(maximum)
n_bins_tot = n_elements(par_hist)
n_bins_max = round(n_bins_tot/cp.n_bins_max_fraction)
asef_threshold = (dp.isla.max_nested_it + dp.isla.n_live) * 3./4.

; At least 2 maxima should be present to define a division array. Otherwise, take as first and last division only 
; the first and last frequency value of the given frequency range.

if n_maxima ge 2 then begin
    for i=0, n_maxima-1 do begin
        if i eq 0 then begin
            divisions_maximum(0,i) = min(par_hist)
        endif else begin
            divisions_maximum(0,i) = divisions_maximum(1,i-1)
        endelse

        if i eq n_maxima-1 then begin
            divisions_maximum(1,i) = max(par_hist)
        endif else begin
            divisions_maximum(1,i) = (maximum(i+1)+maximum(i))/2.d0
        endelse
    endfor

    avg_spacing = mean(maximum(1:*)-maximum(0:n_maxima-2))
    first_division = maximum(0) - avg_spacing
    last_division = maximum(n_maxima-1) + avg_spacing
    if first_division gt divisions_maximum(0,0) then divisions_maximum(0,0) = first_division
    if last_division lt divisions_maximum(1,n_maxima-1) then divisions_maximum(1,n_maxima-1) = last_division
endif else begin
    divisions_maximum(0,0) = min(par_hist)
    divisions_maximum(1,0) = max(par_hist)
endelse

; Find the indices for the ranges around each local maximum. Use a X% drop tolerance on the ASEF for those maxima corresponding
; to peaks that are close to the maximum value of the ASEF. This is done in order to optimize the identification
; of the ranges for possible plateau regions.

for i=0, n_maxima-1 do begin
    index_max_act = index_maximum(i)
   
    ; Left bound
    
    j = 0
    while (asef_hist(index_max_act-j) ge asef_hist(index_max_act-j-1)) and (index_max_act-j-1 gt 0) do begin
        j++
    endwhile
  
    if index_max_act-j lt 0 then begin
        range_maximum(0,i) = min(par_hist)
    endif else begin
        range_maximum(0,i) = par_hist(index_max_act-j)

        if i ge 1 then begin
            if range_maximum(0,i) le range_maximum(1,i-1) then begin
                range_maximum(0,i) = range_maximum(1,i-1)
            endif
        endif
    endelse

    ; Right bound
    
    if keyword_set(chunk) then begin
        ; Adopt a different drop tolerance for MS and SG+RG stars.

        if maximum(i) gt cp.dnu_sg then begin
            drop_tolerance = asef_hist(index_max_act)*cp.drop_tolerance_chunk_ms
        endif else begin
            drop_tolerance = asef_hist(index_max_act)*cp.drop_tolerance_chunk_rg
        endelse
    endif else begin
        drop_tolerance = asef_hist(index_max_act)*cp.drop_tolerance_global
    endelse

    j = 0
    for k=index_max_act, n_elements(par_hist)-2 do begin
        if (asef_hist(k) ge (asef_hist(k+1) - drop_tolerance)) then begin
            j++
        endif else begin
            break
        endelse
    endfor

    if index_max_act+j ge n_elements(par_hist)-1 then begin
        range_maximum(1,i) = max(par_hist)
    endif else begin
        range_maximum(1,i) = par_hist(index_max_act+j)
    endelse

    ; Check that the actual frequency range does not exceed the frequency division 

    if i le n_maxima-2 then begin
        if range_maximum(1,i) gt divisions_maximum(1,i) then begin
            range_maximum(1,i) = divisions_maximum(1,i)
        endif
    
        tmp_right_exceed = where(par_hist gt range_maximum(1,i) and par_hist le divisions_maximum(1,i))

        if tmp_right_exceed(0) ne -1 then begin
            par_hist_exceed = par_hist(tmp_right_exceed)
            index_right_maximum = closest(range_maximum(1,i),par_hist)
            asef_right_maximum = asef_hist(index_right_maximum)
            asef_hist_exceed = asef_hist(tmp_right_exceed)
     
            if min(asef_hist_exceed, min_asef_index) lt asef_right_maximum then begin
                if par_hist_exceed(min_asef_index) le divisions_maximum(1,i) then begin
                    range_maximum(1,i) = par_hist_exceed(min_asef_index)
                endif else begin
                    range_maximum(1,i) = divisions_maximum(1,i)
                endelse
            endif
        endif
    endif
  
    ; Extend right range beyond if the local maximum is the last one of the list and there is an ASEF value
    ; smaller than that at the right bound.

    if i eq n_maxima-1 then begin
        tmp_right_exceed = where(par_hist gt range_maximum(1,i))

        if tmp_right_exceed(0) ne -1 then begin
            par_hist_exceed = par_hist(tmp_right_exceed)
            index_right_maximum = closest(range_maximum(1,i),par_hist)
            asef_right_maximum = asef_hist(index_right_maximum)
            asef_hist_exceed = asef_hist(tmp_right_exceed)

            if min(asef_hist_exceed, min_asef_index) lt asef_right_maximum then begin
                range_maximum(1,i) = par_hist_exceed(min_asef_index)
            endif
        endif
    endif

    if keyword_set(chunk) then begin
        ; In the CHUNK modality verify that if the local maximum is not a prominent one, i.e. its ASEF is below 3/4,
        ; the ranges do not extend over a given maximum number of allowed bins.
        ; This will prevent from having very wide ranges if the peak is small, i.e. surrounded by a large flat ASEF region.

        tmp_bins_range = where(par_hist ge range_maximum(0,i) and par_hist le range_maximum(1,i))
        
        if tmp_bins_range(0) ne -1 then begin
            if n_elements(tmp_bins_range) ge n_bins_max and asef_maximum(i) lt asef_threshold then begin
                ; If here then find the range that is closest to the local maximum and make the other boundary not exceed the
                ; maximum number of bins allowed starting from the selected bound.
              
                index_range_closest = closest(range_maximum(*,i),maximum(i))
                if index_range_closest eq 0 then begin
                    index_left_range = closest(range_maximum(0,i),par_hist)
                    if index_left_range + n_bins_max le (n_elements(par_hist)-1) then begin
                        range_maximum(1,i) = par_hist(index_left_range + n_bins_max)
                    endif
                endif else begin
                    index_right_range = closest(range_maximum(1,i),par_hist)
                    if index_right_range - n_bins_max ge 0 then begin
                        range_maximum(0,i) = par_hist(index_right_range - n_bins_max)
                    endif
                endelse
            endif
        endif
    endif
endfor
end
