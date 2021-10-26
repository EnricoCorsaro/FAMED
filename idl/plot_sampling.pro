pro plot_sampling,par0,position_sampling
; -------------------------------------------------------------------------------------------------------
; Author:  Enrico Corsaro
; e-mail:  enrico.corsaro@inaf.it
; Date:    August 2019
; Place:   Catania, Italy
; Purpose: This routine plots the sampling evolution obtained from the multi-modal fit with DIAMONDS.
; -------------------------------------------------------------------------------------------------------
COMMON GRAPHIC,pp,lp,sp,lpe

x_title = '!3Nested iteration'
loadct,39,/silent
plotsym,0,0.3,/fill,color=220
plot,par0,psym=8,charsize=pp.charsize,xcharsize=pp.xcharsize,ycharsize=pp.ycharsize,charthick=pp.charthick,                                $
xthick=pp.xthick,ythick=pp.ythick,thick=pp.thick,xtitle=x_title,ytitle='!3Frequency (' + sp.freq_unit_str + ')',yr=[min(par0),max(par0)],  $
position=position_sampling,xminor=5,xstyle=1,ystyle=1
axis,xaxis=0,xticklen=0.02,xthick=pp.xthick,charsize=pp.charsize,charthick=pp.charthick,font=1,xtickn=replicate(' ',30),xminor=5
axis,xaxis=1,xticklen=0.02,xthick=pp.xthick,charsize=pp.charsize,charthick=pp.charthick,font=1,xtickn=replicate(' ',30),xminor=5
axis,yaxis=0,yticklen=0.04,ythick=pp.xthick,charsize=pp.charsize,charthick=pp.charthick,font=1,ytickn=replicate(' ',30)
axis,yaxis=1,yticklen=0.04,ythick=pp.xthick,charsize=pp.charsize,charthick=pp.charthick,font=1,ytickn=replicate(' ',30)
end
