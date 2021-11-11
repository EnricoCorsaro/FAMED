pro plot_asef,par_hist,asef_hist,maximum,range_maximum,threshold,position_asef
; -------------------------------------------------------------------------------------------------------
; Author:  Enrico Corsaro
; e-mail:  enrico.corsaro@inaf.it
; Date:    August 2019
; Place:   Catania, Italy
; Purpose: This routine plots the ASEF histogram obtained from the multi-modal sampling of DIAMONDS and
;          highlights the frequency regions of each local maximum, as obtained from a hill-climbing
;          algorithm.
; -------------------------------------------------------------------------------------------------------
COMMON GRAPHIC,pp,lp,sp,lpe

bin_width = (max(par_hist)-min(par_hist))/n_elements(par_hist)
width = (par_hist(1)-par_hist(0))/10.d0
y_title = '!3ASEF (Nested iteration)'

plot,par_hist,asef_hist,xr=[min(par_hist),max(par_hist)],yr=[0,max(asef_hist)*1.30],       $
xcharsize=pp.xcharsize,ycharsize=pp.ycharsize,font=-1,thick=pp.thick,                      $
xthick=pp.xthick,ythick=pp.ythick,xticklen=0.03,yticklen=0.02,charsize=pp.charsize,        $
charthick=pp.charthick,ytitle=y_title,/nodata,xtitle='!3Frequency (' + sp.freq_unit_str + ')',     $
position=position_asef,xstyle=1,ystyle=1

loadct,39,/silent

for i=0,n_elements(maximum)-1 do begin
    oplot,[maximum(i),maximum(i)],[0,max(asef_hist)*1.4],linestyle=1,color=220,thick=3
endfor

loadct,0,/silent

for i=0,n_elements(par_hist)-1 do begin
    loadct,0,/silent
    plotcolorfill,[par_hist[i]-bin_width/2.,par_hist[i]+bin_width/2.],[asef_hist[i],asef_hist[i]],color=pp.color_hist,/noerase,width=width
    loadct,0,/silent
    oplot,[par_hist[i]-bin_width/2.,par_hist[i]-bin_width/2.],[0,asef_hist[i]],color=pp.color_line,thick=0.5
    oplot,[par_hist[i]+bin_width/2.,par_hist[i]+bin_width/2.],[0,asef_hist[i]],color=pp.color_line,thick=0.5
    oplot,[par_hist[i]-bin_width/2.,par_hist[i]+bin_width/2.],[asef_hist[i],asef_hist[i]],color=pp.color_line,thick=0.5
endfor

loadct,39,/silent
axis,xaxis=0,xticklen=0.03,xthick=pp.xthick,charsize=pp.charsize,charthick=pp.charthick,font=1,xtickn=replicate(' ',30)

n_maxima = n_elements(maximum)
color_step = floor(255./n_maxima)
color_start = fix(abs(255.-color_step*n_maxima))

; Plot the individual regions of each local maximum with a different color
for i=0, n_maxima-1 do begin
    upper_bound = range_maximum(1,i)
    lower_bound = range_maximum(0,i)

    upper_index = closest(upper_bound,par_hist)
    lower_index = closest(lower_bound,par_hist)
    par_hist_range = par_hist(lower_index:upper_index)
    asef_hist_range = asef_hist(lower_index:upper_index)
 
    colorFile = 'fsc_brewer.tbl'

    ; Plot histogram counts
    LoadCT,pp.brewer_colorbar_asef,FILE=colorFile,/silent
    for j=0,n_elements(par_hist_range)-1 do begin
        upper_bin = par_hist_range(j)+bin_width/2.
        lower_bin = par_hist_range(j)-bin_width/2.

        if lower_bound gt lower_bin then begin
            if upper_bound lt upper_bin then begin 
                plotcolorfill,[lower_bound,upper_bound],  $
                [asef_hist_range(j),asef_hist_range(j)],   $
                color=i*color_step+color_start,/noerase,width=width
                continue
            endif else begin
                plotcolorfill,[lower_bound,par_hist_range(j)+bin_width/2.],  $
                [asef_hist_range(j),asef_hist_range(j)],   $
                color=i*color_step+color_start,/noerase,width=width
                continue
            endelse
        endif

        if upper_bound lt upper_bin and lower_bound le lower_bin then begin
            plotcolorfill,[par_hist_range(j)-bin_width/2.,upper_bound],[asef_hist_range(j),asef_hist_range(j)],   $
            color=i*color_step+color_start,/noerase,width=width
            continue
        endif

        if upper_bound ge upper_bin and lower_bound le lower_bin then begin
           plotcolorfill,[par_hist_range(j)-bin_width/2.,par_hist_range(j)+bin_width/2.],  $
           [asef_hist_range(j),asef_hist_range(j)],   $
           color=i*color_step+color_start,/noerase,width=width
        endif

    endfor
    
    ; Plot histogram counts
    loadct,0,/silent
    for j=0,n_elements(par_hist_range)-1 do begin
        oplot,[par_hist_range(j)-bin_width/2.,par_hist_range(j)-bin_width/2.],[0,asef_hist_range(j)],color=pp.color_line,thick=0.5
        oplot,[par_hist_range(j)+bin_width/2.,par_hist_range(j)+bin_width/2.],[0,asef_hist_range(j)],color=pp.color_line,thick=0.5
        oplot,[par_hist_range(j)-bin_width/2.,par_hist_range(j)+bin_width/2.],  $
        [asef_hist_range(j),asef_hist_range(j)],color=pp.color_line,thick=0.5
    endfor
endfor

loadct,0,/silent
oplot,[par_hist[0]-bin_width/2.,par_hist[0]+bin_width/2.],[asef_hist[0],asef_hist[0]],color=pp.color_envelope,thick=pp.thick_border_asef

for i=1,n_elements(par_hist)-2 do begin
    loadct,0,/silent
    oplot,[par_hist[i]-bin_width/2.,par_hist[i]-bin_width/2.],[asef_hist[i-1],asef_hist[i]],color=pp.color_envelope,thick=pp.thick_border_asef
    oplot,[par_hist[i]+bin_width/2.,par_hist[i]+bin_width/2.],[asef_hist[i+1],asef_hist[i]],color=pp.color_envelope,thick=pp.thick_border_asef
    oplot,[par_hist[i]-bin_width/2.,par_hist[i]+bin_width/2.],[asef_hist[i],asef_hist[i]],color=pp.color_envelope,thick=pp.thick_border_asef
endfor

i = n_elements(par_hist)-1
oplot,[par_hist[i]-bin_width/2.,par_hist[i]+bin_width/2.],[asef_hist[i],asef_hist[i]],color=pp.color_envelope,thick=pp.thick_border_asef
loadct,39,/silent
axis,xaxis=0,xticklen=0.03,xthick=pp.xthick,charsize=pp.charsize,charthick=pp.charthick,font=1,xtickn=replicate(' ',30)

oplot,[min(par_hist),max(par_hist)],[threshold,threshold],linestyle=2,thick=pp.thick
axis,xaxis=0,xticklen=0.03,xthick=pp.xthick,charsize=pp.charsize,charthick=pp.charthick,font=1,xtickn=replicate(' ',30)
end
