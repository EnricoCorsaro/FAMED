pro write_diamonds_data_range,filename,boundaries
; -------------------------------------------------------------------------------------------------------
; Author:  Enrico Corsaro
; e-mail:  enrico.corsaro@inaf.it
; Date:    August 2019
; Place:   Catania, Italy
; Purpose: This routine creates an input ASCII file in the format to be used by
;          DIAMONDS for its computation to provide the data boundaries for the fit. 
;          It only requires as inputs the full path (including the filename)
;          for the output prior file, and the boundaries values as an array.
; -------------------------------------------------------------------------------------------------------
get_lun, lun1
openw, lun1, filename
printf, lun1, boundaries, format = '(F0.5)'
free_lun, lun1
end