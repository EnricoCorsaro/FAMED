function get_background,catalog_id,star_id
; -------------------------------------------------------------------------------------------------------
; Author:  Enrico Corsaro
; e-mail:  enrico.corsaro@inaf.it
; Date:    August 2019
; Place:   Catania, Italy
; Purpose: This function allows to retrieve the background name adopted during the background fit with 
;          DIAMONDS, as well as its resulting fitted parameters. 
;
;          WARNING: If another fitting method other than DIAMONDS-Background has been used to fit the 
;          background in the input star, then change this routine for properly returning the correct 
;          background name fitted in the given star.
; -------------------------------------------------------------------------------------------------------
COMMON STAR,info

bg_name_list = ['Flat','Original','OneHarvey','OneHarveyColor','TwoHarvey','TwoHarveyColor','ThreeHarvey','ThreeHarveyColor']

if info.external_background_results_dir eq '-99' then begin
    ; Read background model fitted parameters
    
    if info.print_on_screen eq 1 then begin
        print,' Reading background fit result from subfolder ', info.background_run_number
    endif

    readcol,info.background_results_dir + catalog_id + star_id + '/' + info.background_run_number + '/background_parameterSummary.txt',bg_par,format='D,x,x,x,x,x,x',   $
    comment='#',/silent

    ; Read background model name

    if file_test(info.background_results_dir + catalog_id + star_id + '/' + info.background_run_number + '/background_computationParameters.txt') eq 0 then begin
       readcol,info.background_results_dir + catalog_id + star_id + '/' + info.background_run_number + '/background_configuringParameters.txt',config,format='A',comment='#',/silent
    endif else begin
       readcol,info.background_results_dir + catalog_id + star_id + '/' + info.background_run_number + '/background_computationParameters.txt',config,format='A',comment='#',/silent
    endelse

    bg_name = config(n_elements(config)-2)

    ; Apply the following in case the background name is not listed among the computation parameters of the 
    ; Background fit (e.g. an older version of DIAMONDS was used for the Background fit).

    bg_name_npar = [1,3,3,5,5,7,7,9]
    flag_bg_name = [1,0,0,1,0,1,0,1]
    tmp_match = where(bg_name_list eq bg_name)

    if tmp_match(0) eq -1 then begin
        n_bg_par = n_elements(bg_par)-3
        tmp_match2 = where(bg_name_npar eq n_bg_par and flag_bg_name eq 1)
        bg_name = bg_name_list(tmp_match2(0))
    endif
endif else begin
    ; In this case the user has required to adopt an external background, not fitted using the Background code extension of Diamonds.

    ; Read background model fitted parameters and background name

    if info.print_on_screen eq 1 then begin
        print,' Using an external background fit provided by the user.'
    endif
    
    readcol,info.external_background_results_dir + catalog_id + star_id + info.external_background_filename_suffix,bg_variables,format='A',comment='#',/silent

    bg_name = bg_variables(0)
    bg_name = bg_name(0)
    tmp_match = where(bg_name_list eq bg_name)

    if tmp_match(0) eq -1 then begin
        if info.print_on_screen eq 1 then begin
            print,'There is no background model recognized by PeakBagging with the name ',bg_name
            print,'Quitting program.'
            stop
        endif
    endif else begin
        bg_par = float(bg_variables(1:*))
    endelse
endelse

bgp = { parameters:    bg_par,    $
        name:          bg_name    $
      }
return, bgp
end
