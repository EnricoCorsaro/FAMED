pro plot_psd,freq,psd,spsd,bg_level,parameters,flag_global
; -------------------------------------------------------------------------------------------------------
; Author:  Enrico Corsaro
; e-mail:  enrico.corsaro@inaf.it
; Date:    August 2019
; Place:   Catania, Italy
; Purpose: This routine plots the PSD of the star by incorporating different features, depending on the
;          modality used by the FAMED pipeline. It overplots the extracted oscillation frequency peaks.
; -------------------------------------------------------------------------------------------------------
COMMON CONFIG,cp
COMMON GRAPHIC,pp,lp,sp,ppe,lpe
COMMON STAR,info

if flag_global eq 1 then begin
    h_fac = 17.
    position_psd = pp.global.position_psd
endif else begin
    h_fac = 6.
    position_psd = pp.chunk.position_psd
endelse

loadct,39,/silent
plot,freq,psd,ytitle='!3PSD (' + sp.psd_unit_str + ')',xr=[parameters.min_freq,parameters.max_freq],yr=[min(bg_level)*0.8,parameters.max_psd*h_fac],/nodata,   $
xcharsize=pp.xcharsize,ycharsize=pp.ycharsize,font=-1,thick=pp.thick,xthick=pp.xthick,ythick=pp.ythick,xticklen=0.04,yticklen=0.02,   $
charsize=pp.charsize,charthick=pp.charthick,position=position_psd,xtickn=replicate(' ',30),xstyle=1,ystyle=1,/ylog

loadct,0,/silent

if flag_global ne 1 then begin
    oplot,[parameters.lower_freq,parameters.lower_freq],[min(bg_level)*0.8,parameters.max_psd*h_fac],linestyle=3,color=130,thick=pp.thick
endif else begin
    plotcolorfill,[parameters.min_freq,parameters.max_freq],[parameters.max_psd*h_fac,parameters.max_psd*h_fac],color=pp.background_psd,/noerase,width=parameters.dnu/500.
endelse

colorbar = 24
colorFile = 'fsc_brewer.tbl'
LoadCT,colorbar,FILE=colorFile,/silent

for i=0,n_elements(parameters.separations)-2 do begin
    if flag_global eq 1 then begin
        if i mod 2 eq 0 then begin
            plotcolorfill,[parameters.separations(i),parameters.separations(i+1)],[parameters.max_psd*h_fac,parameters.max_psd*h_fac],  $
            color=pp.chunk_band_color1,/noerase,width=parameters.dnu/500.
        endif else begin
            plotcolorfill,[parameters.separations(i),parameters.separations(i+1)],[parameters.max_psd*h_fac,parameters.max_psd*h_fac],  $
            color=pp.chunk_band_color2,/noerase,width=parameters.dnu/500.
        endelse    
    endif else begin
        freq_new = findgen(n_elements(freq)*10)/(n_elements(freq)*10) * (max(freq) - min(freq)) + min(freq)
        
        if parameters.degree(i) eq 0 then begin
            tmpl = where(freq_new gt parameters.freq_list(i)-parameters.freq_sig_list(i)*3 and freq_new lt parameters.freq_list(i))
            tmpr = where(freq_new gt parameters.freq_list(i) and freq_new lt parameters.freq_list(i)+parameters.freq_sig_list(i)*3)
            freq_2l = freq_new(tmpl)
            freq_2r = freq_new(tmpr)
            nl = n_elements(freq_2l)
            nr = n_elements(freq_2r)
            y2l = fltarr(nl) + parameters.max_psd*h_fac
            y2r = fltarr(nr) + parameters.max_psd*h_fac

            if info.save_eps eq 1 then begin
                facl = nl/35.
                facr = nr/35.
                start_col_bin = 220.
            endif else begin
                facl = nl/110.
                facr = nr/110.
                start_col_bin = 110.
            endelse

            ncoll = floor(indgen(nl)/facl)
            ncolr = floor(indgen(nr)/facr)

            loadct,1,/silent

            if info.save_eps eq 1 then begin
                plotcolorfill,freq_2r,y2r,color=start_col_bin+ncolr,/noerase,/notrace
                plotcolorfill,freq_2l,y2l,color=start_col_bin+reverse(ncoll),/noerase,/notrace
            endif else begin
                plotcolorfill,freq_2r,y2r,color=start_col_bin-ncolr,/noerase,/notrace
                plotcolorfill,freq_2l,y2l,color=start_col_bin-reverse(ncoll),/noerase,/notrace
            endelse
        endif

        if parameters.degree(i) eq 1 then begin
            tmpl = where(freq_new gt parameters.freq_list(i)-parameters.freq_sig_list(i)*3 and freq_new lt parameters.freq_list(i))
            tmpr = where(freq_new gt parameters.freq_list(i) and freq_new lt parameters.freq_list(i)+parameters.freq_sig_list(i)*3)
            freq_2l = freq_new(tmpl)
            freq_2r = freq_new(tmpr)
            nl = n_elements(freq_2l)
            nr = n_elements(freq_2r)
            y2l = fltarr(nl) + parameters.max_psd*h_fac
            y2r = fltarr(nr) + parameters.max_psd*h_fac

            if info.save_eps eq 1 then begin
                facl = nl/20.
                facr = nr/20.
                start_col_bin = 235.
                loadct,7,/silent
            endif else begin
                facl = nl/50.
                facr = nr/50.
                start_col_bin = 50.
                loadct,3,/silent
            endelse

            ncoll = floor(indgen(nl)/facl)
            ncolr = floor(indgen(nr)/facr)

            if info.save_eps eq 1 then begin
                plotcolorfill,freq_2r,y2r,color=start_col_bin+ncolr,/noerase,/notrace
                plotcolorfill,freq_2l,y2l,color=start_col_bin+reverse(ncoll),/noerase,/notrace
            endif else begin
                plotcolorfill,freq_2r,y2r,color=start_col_bin-ncolr,/noerase,/notrace
                plotcolorfill,freq_2l,y2l,color=start_col_bin-reverse(ncoll),/noerase,/notrace
            endelse
        endif

        if parameters.degree(i) eq 2 then begin
            tmpl = where(freq_new gt parameters.freq_list(i)-parameters.freq_sig_list(i)*3 and freq_new lt parameters.freq_list(i))
            tmpr = where(freq_new gt parameters.freq_list(i) and freq_new lt parameters.freq_list(i)+parameters.freq_sig_list(i)*3)
            freq_2l = freq_new(tmpl)
            freq_2r = freq_new(tmpr)
            nl = n_elements(freq_2l)
            nr = n_elements(freq_2r)
            y2l = fltarr(nl) + parameters.max_psd*h_fac
            y2r = fltarr(nr) + parameters.max_psd*h_fac

            if info.save_eps eq 1 then begin
                facl = nl/35.
                facr = nr/35.
                start_col_bin = 220.
            endif else begin
                facl = nl/70.
                facr = nr/70.
                start_col_bin = 70.
            endelse

            ncoll = floor(indgen(nl)/facl)
            ncolr = floor(indgen(nr)/facr)

            loadct,8,/silent
            if info.save_eps eq 1 then begin
                plotcolorfill,freq_2r,y2r,color=start_col_bin+ncolr,/noerase,/notrace
                plotcolorfill,freq_2l,y2l,color=start_col_bin+reverse(ncoll),/noerase,/notrace
            endif else begin
                plotcolorfill,freq_2r,y2r,color=start_col_bin-ncolr,/noerase,/notrace
                plotcolorfill,freq_2l,y2l,color=start_col_bin-reverse(ncoll),/noerase,/notrace
            endelse            
        endif

        if parameters.degree(i) eq 3 then begin
            tmpl = where(freq_new gt parameters.freq_list(i)-parameters.freq_sig_list(i)*2 and freq_new lt parameters.freq_list(i))
            tmpr = where(freq_new gt parameters.freq_list(i) and freq_new lt parameters.freq_list(i)+parameters.freq_sig_list(i)*2)
            freq_2l = freq_new(tmpl)
            freq_2r = freq_new(tmpr)
            nl = n_elements(freq_2l)
            nr = n_elements(freq_2r)
            y2l = fltarr(nl) + parameters.max_psd*h_fac
            y2r = fltarr(nr) + parameters.max_psd*h_fac

            if info.save_eps eq 1 then begin
                facl = nl/35.
                facr = nr/35.
                start_col_bin = 220.
                loadct,3,/silent
            endif else begin
                facl = nl/60.
                facr = nr/60.
                start_col_bin = 60.
                loadct,0,/silent
            endelse

            ncoll = floor(indgen(nl)/facl)
            ncolr = floor(indgen(nr)/facr)

            if info.save_eps eq 1 then begin
                plotcolorfill,freq_2r,y2r,color=start_col_bin+ncolr,/noerase,/notrace
                plotcolorfill,freq_2l,y2l,color=start_col_bin+reverse(ncoll),/noerase,/notrace
            endif else begin
                plotcolorfill,freq_2r,y2r,color=start_col_bin-ncolr,/noerase,/notrace
                plotcolorfill,freq_2l,y2l,color=start_col_bin-reverse(ncoll),/noerase,/notrace
            endelse
        endif
      
    endelse
endfor

if (parameters.numax ge parameters.min_freq and parameters.numax le parameters.max_freq) then begin
    loadct,7,/silent
    arrow,parameters.numax,parameters.max_psd*h_fac,parameters.numax,parameters.max_psd*h_fac/2.,/data,thick=pp.numax_arrow_thick,/solid,hsize=pp.numax_arrow_hsize,color=pp.numax_arrow_color
endif

loadct,0,/silent

if flag_global eq 1 then begin
    oplot,freq,psd,color=pp.psd_color_global,thick=pp.psd_thick
    tickscolor = pp.tickscolor_global

    ; Overplot asymptotic frequencies for radial modes
    loadct,39,/silent

    for i=0, n_elements(parameters.freq_radial_asymp)-1 do begin
        arrow,parameters.freq_radial_asymp(i),parameters.max_psd*h_fac,parameters.freq_radial_asymp(i),parameters.max_psd*h_fac/2.2,/data,thick=pp.freq_asymp_arrow_thick,/solid, $
        hsize=pp.freq_asymp_arrow_hsize,color=pp.freq_asymp_arrow_color
    endfor
endif else begin
    oplot,freq,psd,color=pp.psd_color_chunk,thick=pp.psd_thick
    tickscolor = pp.tickscolor_chunk
endelse

loadct,39,/silent
if flag_global eq 1 then begin
    step = parameters.dnu*0.06
endif else begin
    step = parameters.dnu*0.01
endelse

for i=0,n_elements(parameters.freq_list)-1 do begin
    x_pos = parameters.freq_list(i)-step
    y_pos = parameters.max_psd*1.4

    enn = strcompress(string(parameters.order(i)),/remove_all)
    ell = strcompress(string(parameters.degree(i)),/remove_all)

    if info.save_eps eq 1 then begin
        if flag_global eq 1 then begin
            loadct,4,/silent
            xyouts,x_pos,y_pos,'!3('+enn+','+ell+')',charsize=pp.label_charsize-0.1,charthick=pp.label_charthick-0.2,/data,color=240,orientation=90
            oplot,[parameters.freq_list(i),parameters.freq_list(i)],[min(bg_level)*0.7,parameters.max_psd*h_fac],linestyle=1,color=220,thick=3
        endif else begin
            loadct,39,/silent
            xyouts,x_pos,y_pos,'!3('+enn+','+ell+')',charsize=pp.label_charsize,charthick=pp.label_charthick-0.2,/data,color=60,orientation=90
            oplot,[parameters.freq_list(i),parameters.freq_list(i)],[min(bg_level)*0.7,parameters.max_psd*h_fac],linestyle=1,color=40,thick=3
        endelse
    endif else begin
        loadct,4,/silent
        xyouts,x_pos,y_pos,'!3('+enn+','+ell+')',charsize=pp.label_charsize-0.4,charthick=pp.label_charthick-0.2,/data,color=240,orientation=90
        oplot,[parameters.freq_list(i),parameters.freq_list(i)],[min(bg_level)*0.7,parameters.max_psd*h_fac],linestyle=1,color=220,thick=3
    endelse
endfor

loadct,39,/silent
oplot,freq,spsd,color=pp.psd_smth_color,thick=pp.psd_smth_thick

if flag_global ne 1 then begin
    if parameters.n_radial ne 0 then begin
        if info.save_eps eq 1 then begin
            loadct,39,/silent
            oplot,[parameters.freq_radial-parameters.dnu/2.,parameters.freq_radial-parameters.dnu/2.],[min(bg_level)*0.8,parameters.max_psd*h_fac],linestyle=2,color=240,thick=3
            xyouts,parameters.freq_radial-parameters.dnu/2.+parameters.dnu*0.03,parameters.max_psd*0.1,sp.nu_str + '!D0!n - '+sp.dnu_str + '/2', $
            charsize=0.8,charthick=2.2,color=240,orientation=90,noclip=0
        endif else begin
            loadct,39,/silent
            oplot,[parameters.freq_radial-parameters.dnu/2.,parameters.freq_radial-parameters.dnu/2.],[min(bg_level)*0.8,parameters.max_psd*h_fac],linestyle=2,color=170,thick=3
            xyouts,parameters.freq_radial-parameters.dnu/2.+parameters.dnu*0.02,parameters.max_psd*0.1,sp.nu_str + '!D0!n - '+sp.dnu_str + '/2', $
            charsize=1.5,charthick=1.5,color=170,orientation=90,noclip=0
        endelse
    endif
endif

if flag_global eq 1 then begin
    freq_radial = parameters.freq_list(where(parameters.degree eq 0))
    freq_sig_radial = parameters.freq_sig_list(where(parameters.degree eq 0))

    for i=0,n_elements(parameters.freq_radial_org)-1 do begin
        tmp_radial = where(freq ge parameters.freq_radial_org(i)-2*freq_sig_radial(i) and freq le parameters.freq_radial_org(i)+2*freq_sig_radial(i))
        oplot,freq(tmp_radial),spsd(tmp_radial),color=190,thick=pp.psd_smth_thick
    endfor
endif

loadct,7,/silent
oplot,freq,bg_level,color=pp.bg_level_color,thick=pp.psd_smth_thick,linestyle=2

loadct,39,/silent

if flag_global eq 1 then begin
    ; Overplot chunk indices
    for jj=0, n_elements(parameters.separations)-2 do begin
        x_pos = (parameters.separations(jj)+parameters.separations(jj+1))/2.0
        y_pos = parameters.max_psd*h_fac/3.
        xyouts,x_pos,y_pos,strcompress(string(jj),/remove_all),charsize=pp.label_charsize,charthick=pp.label_charthick,/data,color=255
    endfor

    axis,xaxis=1,xticklen=0.04,xthick=pp.xthick,charsize=pp.charsize,charthick=pp.charthick,font=1,xtickn=replicate(' ',30),color=tickscolor
    axis,yaxis=0,yticklen=0.02,ythick=pp.ythick,charsize=pp.charsize,charthick=pp.charthick,font=1,ytickn=replicate(' ',30),color=tickscolor
    axis,yaxis=1,yticklen=0.02,ythick=pp.ythick,charsize=pp.charsize,charthick=pp.charthick,font=1,ytickn=replicate(' ',30),color=tickscolor
endif

axis,xaxis=0,xticklen=0.04,xthick=pp.xthick,charsize=pp.charsize,charthick=pp.charthick,font=1,xtickn=replicate(' ',30),color=pp.tickscolor_global
axis,xaxis=1,xticklen=0.04,xthick=pp.xthick,charsize=pp.charsize,charthick=pp.charthick,font=1,xtickn=replicate(' ',30),color=pp.tickscolor_global
!p.ticklen=0.0
axis,yaxis=0,ythick=pp.ythick,charsize=pp.charsize,charthick=pp.charthick,font=1,yticks=1,ytickn=replicate(' ',30)
axis,yaxis=1,ythick=pp.ythick,charsize=pp.charsize,charthick=pp.charthick,font=1,yticks=1,ytickn=replicate(' ',30)
axis,xaxis=0,xthick=pp.xthick,charsize=pp.charsize,charthick=pp.charthick,font=1,xtickn=replicate(' ',30)
axis,xaxis=1,xthick=pp.xthick,charsize=pp.charsize,charthick=pp.charthick,font=1,xtickn=replicate(' ',30)

if flag_global eq 1 then begin
    freq_range_inset = [parameters.numax - parameters.dnu*0.75,parameters.numax + parameters.dnu*0.75]
    upper_psd = max(spsd)*h_fac*1.5
   
    plot,freq,psd,xr=freq_range_inset, yr=[min(bg_level)*0.8,upper_psd],   $
    xcharsize=pp.xcharsize/1.5,ycharsize=pp.ycharsize,font=-1,thick=pp.thick,xthick=pp.xthick,ythick=pp.ythick,xticklen=0.09,yticklen=0.04,   $
    charsize=pp.charsize,charthick=pp.charthick+0.3,position=pp.global.position_inset,ytickn=replicate(' ',30),background=0,/nodata,color=255,/ylog, $
    xstyle=1,ystyle=1
   
    loadct,3,/silent
    colorbar = 22

    if info.save_eps eq 1 then LoadCT,colorbar,FILE=colorFile,/silent

    plotcolorfill,[min(freq_range_inset),max(freq_range_inset)],[upper_psd,upper_psd],color=pp.background_inset,/noerase,width=max(freq_range_inset)/10.

    loadct,0,/silent
    oplot,freq,psd,color=pp.psd_color_inset

    loadct,39,/silent
    oplot,freq,spsd,color=pp.psd_smth_inset_color,thick=pp.psd_smth_thick-0.1
    axis,xaxis=0,xticklen=0.04,xthick=pp.xthick,charsize=pp.charsize,charthick=pp.charthick,font=1,xtickn=replicate(' ',30)
 
    freq_step = mean(parameters.freq_list)*0.004 
    for i=0, n_elements(parameters.freq_list)-1 do begin
        if (parameters.freq_list(i) ge freq_range_inset(0)) and (parameters.freq_list(i) le freq_range_inset(1)) then begin
            xyouts,parameters.freq_list(i)-freq_step,max(spsd)*3,strcompress(string(parameters.degree(i)),/remove_all), $
            charsize=pp.label_charsize,charthick=pp.label_charthick,color=pp.label_color
        endif
    endfor

    axis,xaxis=0,xticklen=0.09,xthick=pp.xthick,charsize=pp.charsize,charthick=pp.charthick,font=1,xtickn=replicate(' ',30)
    axis,xaxis=1,xticklen=0.09,xthick=pp.xthick,charsize=pp.charsize,charthick=pp.charthick,font=1,xtickn=replicate(' ',30)
    axis,yaxis=0,yticklen=0.04,ythick=pp.ythick,charsize=pp.charsize,charthick=pp.charthick,font=1,ytickn=replicate(' ',30)
    axis,yaxis=1,yticklen=0.04,ythick=pp.ythick,charsize=pp.charsize,charthick=pp.charthick,font=1,ytickn=replicate(' ',30)
endif
end
