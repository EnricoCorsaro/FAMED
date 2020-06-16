function interpolate_epsilon,teff,dnu
; -------------------------------------------------------------------------------------------------------
; Author:  Enrico Corsaro
; e-mail:  enrico.corsaro@inaf.it
; Date:    August 2019
; Place:   Catania, Italy
; Purpose: This function computes the value of epsilon for either a main-sequence, subgiant star or a 
;          red giant star by exploiting empirical relations linking epsilon with the stellar parameters 
;          Teff, DeltaNu.
; -------------------------------------------------------------------------------------------------------
COMMON CONFIG,cp
COMMON GRAPHIC,pp,lp,sp,ppe,lpe
COMMON STAR,info

if (dnu gt cp.dnu_threshold) then begin
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
    result = poly_fit(teff_array,epsi_array,4,yfit=fit,/double)
    good_epsi = interpol(fit,teff_array,teff,/spline)
    good_epsi = good_epsi(0)
    
    if ~(good_epsi gt 0) then begin
       good_epsi = interpol(fit,teff_array,teff)
       good_epsi = good_epsi(0)
    endif
    
    if info.print_on_screen eq 1 then begin
        loadct,39,/silent
        plot,teff_array,epsi_array,yr=[1.7,0.7],xr=[7000.,4700.],/nodata,xtitle='!3T!Deff!n',ytitle=sp.epsi_str,    $
        charsize=pp.epsi.charsize,xcharsize=pp.epsi.xcharsize,ycharsize=pp.epsi.ycharsize,  $
        charthick=pp.epsi.charthick,xthick=pp.epsi.xthick,ythick=pp.epsi.ythick,position=pp.global.position_epsi,xticklen=0.05,yticklen=0.04
        
        for i=0,n_elements(teff_array)-1 do begin
            loadct,3,/silent
            plotsym,0,pp.symsize,/fill
            oplot,[teff_array(i)],[epsi_array(i)],psym=8,color=pp.sample_color
            loadct,39,/silent
            plotsym,0,pp.symsize+0.05,thick=pp.symthickness
            oplot,[teff_array(i)],[epsi_array(i)],psym=8,color=255
        endfor
        
        loadct,39,/silent
        oplot,teff_array,fit,color=pp.epsi.fit_color,thick=3
        plotsym,8,pp.epsi.symsize+0.2,/fill
        oplot,[teff],[good_epsi],psym=8,color=pp.epsi.cross_color
        plotsym,8,pp.epsi.symsize+0.25,thick=pp.epsi.symthickness
        oplot,[teff],[good_epsi],psym=8,color=255
        oplot,[teff,teff],[0,20],linestyle=2,color=pp.epsi.cross_color,thick=2
        oplot,[100,20d4],[good_epsi,good_epsi],linestyle=2,color=pp.epsi.cross_color,thick=2
    
        flag_plot_epsi = 0
        width=1.
    
        plotcolorfill,[5700,4800],[0.8,0.8],bottom=[0.95,0.95],color=pp.epsi.color_box,/noerase,/notrace,width=width
        xyouts,5600.,0.9,'!3'+sp.epsi_str+' = '+strcompress(string(good_epsi,format='(F0.3)'),/remove_all),    $
        charsize=pp.epsi.label_charsize,charthick=pp.epsi.charthick,font=-1,color=pp.epsi.cross_color,/data
    endif
endif else begin
    ; epsilon-DeltaNu relation by Corsaro et al. 2012b
    
    good_epsi = cp.epsilon_offset + cp.epsilon_slope*alog10(dnu)
    dnu_array = findgen(301)/300. * cp.dnu_threshold + 0.15
    epsi_array = cp.epsilon_offset + cp.epsilon_slope*alog10(dnu_array)

    if info.print_on_screen eq 1 then begin
        loadct,39,/silent
        plot,epsi_array,dnu_array,xr=[0.4,1.6],yr=[0.,cp.dnu_threshold],/nodata,ytitle=sp.dnu_str + '(' + sp.freq_unit_str + ')',xtitle=sp.epsi_str,    $
        charsize=pp.epsi.charsize,xcharsize=pp.epsi.xcharsize,ycharsize=pp.epsi.ycharsize,  $
        charthick=pp.epsi.charthick,xthick=pp.epsi.xthick,ythick=pp.epsi.ythick,position=pp.global.position_epsi,xticklen=0.05,yticklen=0.04
         
        oplot,epsi_array,dnu_array,color=pp.epsi.fit_color,thick=3

        plotsym,8,pp.epsi.symsize+0.2,/fill
        oplot,[good_epsi],[dnu],psym=8,color=pp.epsi.cross_color
        plotsym,8,pp.epsi.symsize+0.25,thick=pp.epsi.symthickness
        oplot,[good_epsi],[dnu],psym=8,color=255

        oplot,[good_epsi,good_epsi],[0,cp.dnu_threshold],linestyle=2,color=pp.epsi.cross_color,thick=2
        oplot,[0.2,1.8],[dnu,dnu],linestyle=2,color=pp.epsi.cross_color,thick=2
        flag_plot_epsi = 0

        width = 0.05
        plotcolorfill,[0.6,0.98],[cp.dnu_threshold*0.78,cp.dnu_threshold*0.78],bottom=[cp.dnu_threshold*0.67,cp.dnu_threshold*0.67],    $
        color=pp.epsi.color_box,/noerase,/notrace,width=width
        xyouts,0.6,cp.dnu_threshold*0.7,'!3'+sp.epsi_str+' = '+strcompress(string(good_epsi,format='(F0.3)'),/remove_all),    $
        charsize=pp.epsi.label_charsize,charthick=pp.epsi.charthick,font=-1,color=pp.epsi.cross_color,/data
    endif
endelse
 
return, good_epsi
end