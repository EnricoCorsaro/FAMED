function compute_acf_dnu,scaling_dnu,par_hist,asef_hist,plot=plot
; -------------------------------------------------------------------------------------------------------
; Author:  Enrico Corsaro
; e-mail:  enrico.corsaro@inaf.it
; Date:    August 2019
; Place:   Catania, Italy
; Purpose: This function computes the ACF^2 of the multi-modal sampling obtained with DIAMONDS to evaluate 
;          an estimate of the large frequency separation DeltaNu. It requires an estimate of DeltaNu
;          from scaling.
; -------------------------------------------------------------------------------------------------------
COMMON GRAPHIC,pp,lp,sp,ppe,lpe
COMMON CONFIG,cp 

dnu_range_side = scaling_dnu * cp.dnu_acf_range_side
top_dnu = scaling_dnu + dnu_range_side
bottom_dnu = scaling_dnu - dnu_range_side

freqbin = par_hist(1) - par_hist(0)
dnu_range_bins = round((top_dnu-bottom_dnu)/freqbin)
bottom_dnu_bins = round(bottom_dnu/freqbin)
lag = indgen(dnu_range_bins) + bottom_dnu_bins
result = a_correlate(asef_hist,lag,/double)
result = result - min(result)

n_interpol = 101
interpolated_dnu = findgen(n_interpol)/(n_interpol-1)*(top_dnu - bottom_dnu) + bottom_dnu
interpolated_acf = interpol(result,lag*freqbin,interpolated_dnu,/spline)
best_acf = max(interpolated_acf,j)
acf_dnu = interpolated_dnu(j)

; Perform a Gaussian fit to the ACF^2 from the smoothed PSD

x = interpolated_dnu
y = interpolated_acf^2
yfit = GAUSSFIT(x, y, coeff, NTERMS=6, ESTIMATES=[max(y),acf_dnu,0,min(y),0,0])
acf_dnu = coeff(1)

if keyword_set(plot) then begin
    plot,lag*freqbin,result^2,xr=[bottom_dnu,top_dnu],yr=[0,max(result^2)*1.2],xtitle=sp.dnu_str + ' (' + sp.freq_unit_str + ')',ytitle='!3ACF!u2!n',     $
    charsize=pp.acf.charsize,xcharsize=pp.acf.xcharsize,ycharsize=pp.acf.ycharsize,xthick=pp.acf.xthick,ythick=pp.acf.ythick, $
    charthick=pp.acf.charthick,position=pp.global.position_acf_dnu,xticklen=0.02,yticklen=0.05
    
    loadct,0,/silent
    oplot,lag*freqbin,result^2,color=80,thick=2
    
    loadct,39,/silent
    oplot,interpolated_dnu,interpolated_acf^2,color=pp.acf.sample_color,thick=2
    
    loadct,39,/silent
    oplot,[scaling_dnu,scaling_dnu],[0,max(result^2)*1.5],linestyle=1,color=240,thick=3
    
    axis,xaxis=0,charthick=pp.acf.charthick,xcharsize=pp.acf.xcharsize,xthick=pp.acf.xthick,xtickn=replicate(' ',30)
    oplot,[acf_dnu,acf_dnu],[0,max(result^2)*1.5],linestyle=2,color=pp.acf.cross_color,thick=3
endif

return, acf_dnu
end