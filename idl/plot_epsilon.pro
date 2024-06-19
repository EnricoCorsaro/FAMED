pro plot_epsilon,parameters
; -------------------------------------------------------------------------------------------------------
; Author:  Enrico Corsaro
; e-mail:  enrico.corsaro@inaf.it
; Date:    March 2023
; Place:   Catania, Italy
; Purpose: This routine plots the epsilon diagram (either epsilon-Dnu or epsilon-Teff) of the star by 
;          showing both the value predicted by the scaling relation and the value obtained from the 
; 		   asymptotic fit of the radial modes.
; -------------------------------------------------------------------------------------------------------
COMMON CONFIG,cp
COMMON GRAPHIC,pp,lp,sp,lpe
COMMON STAR,info

; Evaluate percentage variation with respect to the interpolated value
epsi_variation = (parameters.epsi_fit - parameters.epsi_int)/parameters.epsi_int*100

if parameters.dnu_fit gt cp.dnu_threshold then begin
    ; Plot the epsilon-Teff diagram

    loadct,39,/silent
    plot,parameters.teff_array,parameters.epsi_array,yr=[1.7,0.7],xr=[7000.,4700.],/nodata,xtitle='!3T!Deff!n',ytitle=sp.epsi_str,    $
    charsize=pp.epsi.charsize,xcharsize=pp.epsi.xcharsize,ycharsize=pp.epsi.ycharsize,  $
    charthick=pp.epsi.charthick,xthick=pp.epsi.xthick,ythick=pp.epsi.ythick,position=pp.global.position_epsi,xticklen=0.05,yticklen=0.04
    
    for i=0,n_elements(parameters.teff_array)-1 do begin
        loadct,0,/silent
        plotsym,0,pp.symsize-0.5,/fill
        oplot,[parameters.teff_array(i)],[parameters.epsi_array(i)],psym=8,color=pp.sample_color-120
        loadct,39,/silent
        plotsym,0,pp.symsize-0.45,thick=pp.symthickness-0.5
        oplot,[parameters.teff_array(i)],[parameters.epsi_array(i)],psym=8,color=255
    endfor
    
    loadct,39,/silent
    oplot,parameters.teff_array,parameters.fit_array,color=pp.epsi.fit_color,thick=3
    
    plotsym,8,pp.epsi.symsize+0.2,/fill
    oplot,[parameters.teff],[parameters.epsi_int],psym=8,color=pp.epsi.cross_color
    plotsym,8,pp.epsi.symsize+0.25,thick=pp.epsi.symthickness
    oplot,[parameters.teff],[parameters.epsi_int],psym=8,color=255
    oplot,[parameters.teff,parameters.teff],[0,20],linestyle=2,color=pp.epsi.cross_color,thick=2
    oplot,[100,20d4],[parameters.epsi_int,parameters.epsi_int],linestyle=2,color=pp.epsi.cross_color,thick=2
    
    loadct,18,/silent
	plotsym,0,pp.epsi.symsize+0.2,/fill
    oplot,[parameters.teff],[parameters.epsi_fit],psym=8,color=pp.epsi.cross_color-30
    plotsym,0,pp.epsi.symsize+0.25,thick=pp.epsi.symthickness
    oplot,[parameters.teff],[parameters.epsi_fit],psym=8,color=255
    oplot,[parameters.teff,parameters.teff],[0,20],linestyle=2,color=pp.epsi.cross_color-30,thick=2
    oplot,[100,20d4],[parameters.epsi_fit,parameters.epsi_fit],linestyle=2,color=pp.epsi.cross_color-30,thick=2

    width = 1.
   
    loadct,39,/silent 
    plotcolorfill,[5730,4890],[0.84,0.84],bottom=[0.93,0.93],color=pp.epsi.color_box,/noerase,/notrace,width=width
    xyouts,5650.,0.9,'!3'+sp.epsi_str+'!Dint!n = '+strcompress(string(parameters.epsi_int,format='(F0.3)'),/remove_all),    $
    charsize=pp.epsi.label_charsize,charthick=pp.epsi.charthick,font=-1,color=pp.epsi.cross_color,/data

    plotcolorfill,[5730,4890],[0.93,0.93],bottom=[1.04,1.04],color=pp.epsi.color_box,/noerase,/notrace,width=width
    loadct,18,/silent
    xyouts,5690.,1.00,'!3'+sp.epsi_str+'!Dech!n = '+strcompress(string(parameters.epsi_fit,format='(F0.3)'),/remove_all),    $
    charsize=pp.epsi.label_charsize,charthick=pp.epsi.charthick,font=-1,color=pp.epsi.cross_color-30,/data

    loadct,39,/silent
    plotcolorfill,[5570,4910],[1.04,1.04],bottom=[1.13,1.13],color=pp.epsi.color_box,/noerase,/notrace,width=width

	loadct,18,/silent
    xyouts,5550,1.10,'!3' + '(' + strcompress(string(epsi_variation,format='(F0.2)'),/remove_all) + ' %)',    $
    charsize=pp.epsi.label_charsize,charthick=pp.epsi.charthick,font=-1,color=pp.epsi.cross_color-30,/data
	loadct,39,/silent
endif else begin
    ; Plot the epsilon-DeltaNu relation by Corsaro et al. 2012b
    
    loadct,39,/silent
    plot,parameters.fit_array,parameters.dnu_array,xr=[0.4,1.6],yr=[0.,cp.dnu_threshold],/nodata,ytitle=sp.dnu_str + '(' + sp.freq_unit_str + ')',xtitle=sp.epsi_str,    $
    charsize=pp.epsi.charsize,xcharsize=pp.epsi.xcharsize,ycharsize=pp.epsi.ycharsize,  $
    charthick=pp.epsi.charthick,xthick=pp.epsi.xthick,ythick=pp.epsi.ythick,position=pp.global.position_epsi,xticklen=0.05,yticklen=0.04
     
    oplot,parameters.fit_array,parameters.dnu_array,color=pp.epsi.fit_color,thick=3

    plotsym,8,pp.epsi.symsize+0.2,/fill
    oplot,[parameters.epsi_int],[parameters.dnu_acf],psym=8,color=pp.epsi.cross_color
    plotsym,8,pp.epsi.symsize+0.25,thick=pp.epsi.symthickness
    oplot,[parameters.epsi_int],[parameters.dnu_acf],psym=8,color=255

    oplot,[parameters.epsi_int,parameters.epsi_int],[0,cp.dnu_threshold],linestyle=2,color=pp.epsi.cross_color,thick=2
    oplot,[0.2,1.8],[parameters.dnu_acf,parameters.dnu_acf],linestyle=2,color=pp.epsi.cross_color,thick=2

    loadct,18,/silent
    plotsym,0,pp.epsi.symsize+0.2,/fill
    oplot,[parameters.epsi_fit],[parameters.dnu_fit],psym=8,color=pp.epsi.cross_color-30
    plotsym,0,pp.epsi.symsize+0.25,thick=pp.epsi.symthickness
    oplot,[parameters.epsi_fit],[parameters.dnu_fit],psym=8,color=255

    oplot,[parameters.epsi_fit,parameters.epsi_fit],[0,cp.dnu_threshold],linestyle=2,color=pp.epsi.cross_color-30,thick=2
    oplot,[0.2,1.8],[parameters.dnu_fit,parameters.dnu_fit],linestyle=2,color=pp.epsi.cross_color-30,thick=2

    width = 0.05
    loadct,39,/silent
    plotcolorfill,[0.6,0.98],[cp.dnu_threshold*0.77,cp.dnu_threshold*0.77],bottom=[cp.dnu_threshold*0.66,cp.dnu_threshold*0.66],    $
    color=pp.epsi.color_box,/noerase,/notrace,width=width
    xyouts,0.6,cp.dnu_threshold*0.7,'!3'+sp.epsi_str+'!Dint!n = '+strcompress(string(parameters.epsi_int,format='(F0.3)'),/remove_all),    $
    charsize=pp.epsi.label_charsize,charthick=pp.epsi.charthick,font=-1,color=pp.epsi.cross_color,/data

    plotcolorfill,[0.6,0.98],[cp.dnu_threshold*0.66,cp.dnu_threshold*0.66],bottom=[cp.dnu_threshold*0.57,cp.dnu_threshold*0.57],    $
    color=pp.epsi.color_box,/noerase,/notrace,width=width

	loadct,18,/silent
    xyouts,0.58,cp.dnu_threshold*0.6,'!3' + sp.epsi_str + '!Dech!n = ' + strcompress(string(parameters.epsi_fit,format='(F0.3)'),/remove_all),     $
    charsize=pp.epsi.label_charsize,charthick=pp.epsi.charthick,font=-1,color=pp.epsi.cross_color-30,/data
    
    loadct,39,/silent
    plotcolorfill,[0.70,1.03],[cp.dnu_threshold*0.58,cp.dnu_threshold*0.58],bottom=[cp.dnu_threshold*0.47,cp.dnu_threshold*0.47],    $
    color=pp.epsi.color_box,/noerase,/notrace,width=width

	loadct,18,/silent
    xyouts,0.70,cp.dnu_threshold*0.5,'!3' + '(' + strcompress(string(epsi_variation,format='(F0.2)'),/remove_all) + ' %)',    $
    charsize=pp.epsi.label_charsize,charthick=pp.epsi.charthick,font=-1,color=pp.epsi.cross_color-30,/data

	loadct,39,/silent

endelse

end