pro plot_summary,parameters,flag_global
; -------------------------------------------------------------------------------------------------------
; Author:  Enrico Corsaro
; e-mail:  enrico.corsaro@inaf.it
; Date:    August 2019
; Place:   Catania, Italy
; Purpose: This routine plots the summary information on the run (either global or chunk modality).
; -------------------------------------------------------------------------------------------------------
COMMON CONFIG,cp
COMMON STAR,info
COMMON GRAPHIC,pp,lp,sp,ppe,lpe

; General information on the run
xyouts,lp.block1_1[0],lp.block1_1[1],'!3' + parameters.catalog_id + ' ' + parameters.star_id + '!X', $
charsize=lp.summary_charsize_title,charthick=lp.summary_charthick_title,font=-1,color=lp.block1_color,/normal
xyouts,lp.block1_2[0],lp.block1_2[1],'!3Folder: '+strcompress(string(info.isla_subdir),/remove_all), $
charsize=lp.summary_charsize,charthick=lp.summary_charthick,font=-1,color=lp.block1_color,/normal
xyouts,lp.block1_3[0],lp.block1_3[1],'!3Run: '+strcompress(string(parameters.run),/remove_all),  $
charsize=lp.summary_charsize,charthick=lp.summary_charthick,font=-1,color=lp.block1_color,/normal
xyouts,lp.block1_4[0],lp.block1_4[1],'!3Modality: ' + parameters.modality,   $
charsize=lp.summary_charsize,charthick=lp.summary_charthick,font=-1,color=lp.block1_color,/normal

; Asymptotic
xyouts,lp.block2_1[0],lp.block2_1[1],'!3'+sp.nu_str+'!Dmax!n = '+strcompress(string(parameters.numax,format='(F0.3)'),/remove_all) + ' ' + sp.freq_unit_str,    $
charsize=lp.summary_charsize,charthick=lp.summary_charthick,font=-1,color=lp.block2_color,/normal
xyouts,lp.block2_2[0],lp.block2_2[1],'!3'+sp.dnu_str+'!Dfit!n = '+strcompress(string(parameters.best_dnu,format='(F0.3)'),/remove_all) + ' ' + sp.freq_unit_str,    $
charsize=lp.summary_charsize,charthick=lp.summary_charthick,font=-1,color=lp.block2_color,/normal
xyouts,lp.block2_7[0],lp.block2_7[1],'!3'+sp.alpha_str+' = '+strcompress(string(parameters.best_alpha,format='(F0.3)'),/remove_all),    $
charsize=lp.summary_charsize,charthick=lp.summary_charthick,font=-1,color=lp.block2_color,/normal

if flag_global eq 1 then begin
    xyouts,lp.block2_3[0],lp.block2_3[1],'!3'+sp.dnu_str+'!DACF!n = '+strcompress(string(parameters.acf_dnu,format='(F0.3)'),/remove_all) + ' ' + sp.freq_unit_str,   $
    charsize=lp.summary_charsize,charthick=lp.summary_charthick,font=-1,color=lp.block2_color,/normal
endif else begin
    xyouts,lp.block2_3[0],lp.block2_3[1],'!3'+sp.epsi_str+' = '+strcompress(string(parameters.local_epsi,format='(F0.2)'),/remove_all),    $
    charsize=lp.summary_charsize,charthick=lp.summary_charthick,font=-1,color=lp.block2_color,/normal
    xyouts,lp.block2_5[0],lp.block2_5[1],'!3'+sp.d02_str+' = '+strcompress(string(parameters.local_d02,format='(F0.2)'),/remove_all) + ' ' + sp.freq_unit_str,    $
    charsize=lp.summary_charsize,charthick=lp.summary_charthick,font=-1,color=lp.block2_color,/normal
endelse

xyouts,lp.block2_4[0],lp.block2_4[1],'!3SNR = '+strcompress(string(parameters.snr,format='(F0.1)'),/remove_all),    $
charsize=lp.summary_charsize,charthick=lp.summary_charthick,font=-1,color=lp.block2_color,/normal

if flag_global eq 1 then begin
    epsi_subscript = 'ech'
   
    if parameters.flag_bad_epsi eq 1 then begin
        xyouts,lp.block2_6[0],lp.block2_6[1],'!3'+sp.epsi_str+'!D'+epsi_subscript+'!n = '+strcompress(string(parameters.best_epsi,format='(F0.3)'),/remove_all),    $
        charsize=lp.summary_charsize,charthick=lp.summary_charthick,font=-1,color=245,/normal
    endif else begin
        xyouts,lp.block2_6[0],lp.block2_6[1],'!3'+sp.epsi_str+'!D'+epsi_subscript+'!n = '+strcompress(string(parameters.best_epsi,format='(F0.3)'),/remove_all),    $
        charsize=lp.summary_charsize,charthick=lp.summary_charthick,font=-1,color=lp.block2_color,/normal
    endelse
endif else begin
    xyouts,lp.block2_6[0],lp.block2_6[1],'!3'+sp.dp_str+' = '+strcompress(string(parameters.local_dp,format='(F0.1)'),/remove_all) + ' s',    $
    charsize=lp.summary_charsize,charthick=lp.summary_charthick,font=-1,color=lp.block2_color,/normal
endelse

; Priors and observables
xyouts,lp.block3_1[0],lp.block3_1[1],'!3'+sp.gamma_str+'!Dfit!n = '+strcompress(string(parameters.fit_linewidth,format='(F0.3)'),/remove_all) + ' ' + sp.freq_unit_str, $
charsize=lp.summary_charsize,charthick=lp.summary_charthick,font=-1,color=lp.block3_color,/normal

if flag_global eq 1 then begin
    xyouts,lp.block3_2[0],lp.block3_2[1],'!3'+sp.gamma_str+'!D'+sp.nu_str+'!Dmax!n = ' +  $
    strcompress(string(get_linewidth(parameters.numax,parameters.teff,parameters.numax),format='(F0.2)'),/remove_all) + ' ' + sp.freq_unit_str, $
    charsize=lp.summary_charsize,charthick=lp.summary_charthick,font=-1,color=lp.block3_color,/normal
endif else begin
    if parameters.n_radial_chunk ne 0 then begin
        xyouts,lp.block3_2[0],lp.block3_2[1],'!3'+sp.gamma_str+'!Dradial!n = '+strcompress(string(parameters.fwhm_radial_fit,format='(F0.2)'),/remove_all) + $
        ' ' + sp.freq_unit_str,charsize=lp.summary_charsize,charthick=lp.summary_charthick,font=-1,color=lp.block3_color,/normal
    endif else begin
        xyouts,lp.block3_2[0],lp.block3_2[1],'!3'+sp.gamma_str+'!Davg!n = '+strcompress(string(parameters.avg_fwhm,format='(F0.2)'),/remove_all) + ' ' + sp.freq_unit_str, $
        charsize=lp.summary_charsize,charthick=lp.summary_charthick,font=-1,color=lp.block3_color,/normal
    endelse
endelse

xyouts,lp.block3_3[0],lp.block3_3[1],'!3H!Dmax,prior!n = '+strcompress(string(parameters.upper_height,format='(E0.1)'),/remove_all) + ' ' + sp.psd_unit_str,  $
charsize=lp.summary_charsize,charthick=lp.summary_charthick,font=-1,color=lp.block3_color,/normal
xyouts,lp.block3_4[0],lp.block3_4[1],'!3T!Deff!n = '+strcompress(string(parameters.teff,format='(I0)'),/remove_all) + ' K',    $
charsize=lp.summary_charsize,charthick=lp.summary_charthick,font=-1,color=lp.block3_color,/normal

; Analysis configuration
xyouts,lp.block4_1[0],lp.block4_1[1],'!3ASEF!Dthreshold!n = '+strcompress(string(parameters.threshold_asef*100,format='(F0.2)'),/remove_all) + ' %',    $
charsize=lp.summary_charsize,charthick=lp.summary_charthick,font=-1,/normal
xyouts,lp.block4_2[0],lp.block4_2[1],'!3ASEF!Dbins!n = '+strcompress(string(parameters.n_bins),/remove_all),   $
charsize=lp.summary_charsize,charthick=lp.summary_charthick,font=-1,/normal

if flag_global eq 1 then begin
    xyouts,lp.block4_3[0],lp.block4_3[1],'!3'+sp.dnu_str+'!Dtolerance!n = '+strcompress(string(parameters.tolerance*100.,format='(F0.2)'),/remove_all) + ' %', $
    charsize=lp.summary_charsize,charthick=lp.summary_charthick,font=-1,/normal
    xyouts,lp.block4_4[0],lp.block4_4[1],'!3N!Dfreq!n = '+strcompress(string(parameters.n_freq),/remove_all),  $
    charsize=lp.summary_charsize,charthick=lp.summary_charthick,font=-1,/normal
    xyouts,lp.block4_5[0],lp.block4_5[1],'!3N!Dorders!n = '+strcompress(string(parameters.n_chunks),/remove_all),  $
    charsize=lp.summary_charsize,charthick=lp.summary_charthick,font=-1,/normal
endif else begin
    xyouts,lp.block4_4[0],lp.block4_3[1],'!3N!Dfreq!n = '+strcompress(string(parameters.n_freq),/remove_all),  $
    charsize=lp.summary_charsize,charthick=lp.summary_charthick,font=-1,/normal
    xyouts,lp.block4_5[0],lp.block4_3[1],'!3N!Dorders!n = '+strcompress(string(parameters.n_chunks),/remove_all),  $
    charsize=lp.summary_charsize,charthick=lp.summary_charthick,font=-1,/normal
endelse
end
