pro write_diamonds_prior,filename,boundaries
; -------------------------------------------------------------------------------------------------------
; Author:  Enrico Corsaro
; e-mail:  enrico.corsaro@inaf.it
; Date:    August 2019
; Place:   Catania, Italy
; Purpose: This routine creates a uniform prior hyper-parameter ASCII file in the format to be used by
;          DIAMONDS for its computation. It only requires as inputs the full path (including the filename)
;          for the output prior file, and the boundaries values as an array.
; -------------------------------------------------------------------------------------------------------
get_lun, lun1
openw, lun1, filename
printf, lun1, '# Hyper parameters used for setting up uniform priors.', format = '(A0)'
printf, lun1, '# Each line corresponds to a different free parameter (coordinate).', format = '(A0)'
printf, lun1, '# Column #1: Minima (lower boundaries)', format = '(A0)'
printf, lun1, '# Column #2: Maxima (upper boundaries)', format = '(A0)'
printf, lun1, boundaries, format = '(F0.3, F20.3)'
free_lun, lun1
end