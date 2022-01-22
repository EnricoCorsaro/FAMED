import numpy as np
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import matplotlib.ticker

__all__ = ['global_plot',
           'nested_iterations_plot',
           'psd_plot',
           'psd_inset_plot',
           'epsilon_plot',
           'acf_plot',
           'asef_histogram',
           'text_panel']

plt.ion()

def global_plot(famed_obj):
    """
    Produce the multipanel GLOBAL output plot for FAMED.

    Parameters
    ----------
    famed_obj : FamedStar object
        A FamedStar class object that has been processed through the GLOBAL 
        module of FAMED.
    """
    fig = plt.figure(figsize=(18,8))
    plt.clf()
    ax_eps = fig.add_axes([0.05,0.7,0.15,0.27])
    epsilon_plot(famed_obj,ax_eps)
    ax_acf = fig.add_axes([0.25,0.7,0.1,0.27])
    acf_plot(famed_obj,ax_acf)
    ax_nested = fig.add_axes([0.05,0.07,0.3,0.55])
    nested_iterations_plot(famed_obj,ax_nested)
    ax_psd = fig.add_axes([0.42,0.435,0.56,0.42])
    psd_plot(famed_obj,ax_psd)
    ax_asef = fig.add_axes([0.42,0.07,0.56,0.35],sharex=ax_psd)
    asef_histogram(famed_obj,ax_asef)
    plt.setp(ax_psd.get_xticklabels(), visible=False)
    ax_inset = fig.add_axes([0.81,0.60,0.16,0.10],yticklabels=[])
    psd_inset_plot(famed_obj,ax_inset)
    ax_text = fig.add_axes([0.36,0.88,0.65,0.25],xticks=[],yticks=[])
    text_panel(famed_obj,ax_text)

def chunk_plot(famed_obj,chunk=None):
    """
    Produce the multipanel CHUNK output plot for FAMED.

    Parameters
    ----------
    famed_obj : FamedStar object
        A FamedStar class object that has been processed through the GLOBAL 
        module of FAMED.
    """
    fig = plt.figure(figsize=(18,8))
    plt.clf()

    ax_nested = fig.add_axes([0.05,0.07,0.3,0.9])
    nested_iterations_plot(famed_obj,ax_nested,chunk=chunk)
    ax_psd = fig.add_axes([0.42,0.435,0.56,0.42])
    psd_plot(famed_obj,ax_psd,chunk=chunk)
    ax_asef = fig.add_axes([0.42,0.07,0.56,0.35],sharex=ax_psd)
    asef_histogram(famed_obj,ax_asef,chunk=chunk)
    plt.setp(ax_psd.get_xticklabels(), visible=False)
    ax_text = fig.add_axes([0.36,0.88,0.65,0.25],xticks=[],yticks=[])
    text_panel(famed_obj,ax_text,chunk=chunk)

def complete_plot(famed_obj):
    print('This feature is not yet implemented')
    
def echelle_plot(famed_obj):
    print('This feature is not yet implemented')

def nested_iterations_plot(famed_obj,ax=None,chunk=None):
    """
    Plot the nested iterations from a multimodal sampling run of FAMED.

    Parameters
    ----------
    famed_obj : FamedStar object
        A FamedStar class object that has been processed through the GLOBAL 
        module of FAMED.
    ax : matplotlib.pyplot.axes object, default: None
        Axes to plot this figure to. If None, the current axes are cleared and
        the plot is created on the full figure.
    """
    if chunk is not None:
        nested_iters = famed_obj.nested_iters[chunk]
    else:
        nested_iters = famed_obj.nested_iters
    if ax is None:
        plt.clf()
        ax = plt.axes()
    ax.set_ylabel(r'Frequency ($\mu$Hz)')
    ax.set_xlabel('Nested iteration')
    ax.plot(nested_iters,'.',c=famed_obj.cp.nested_dots,zorder=1,ms=1)
    ax.set_xlim(0,len(nested_iters))
    plt.annotate('',xy=(.9,famed_obj.numax,),xycoords=transforms.blended_transform_factory(ax.transAxes,ax.transData),xytext=(1,famed_obj.numax),textcoords=transforms.blended_transform_factory(ax.transAxes,ax.transData),arrowprops=dict(color=famed_obj.cp.numax_arrow, lw=3, arrowstyle='-|>'),clip_on=True)
        
def psd_plot(famed_obj,ax=None,chunk=None):
    """
    Plot the PSD and identified modes from the GLOBAL module of FAMED.

    Parameters
    ----------
    famed_obj : FamedStar object
        A FamedStar class object that has been processed through the GLOBAL 
        module of FAMED.
    ax : matplotlib.pyplot.axes object, default: None
        Axes to plot this figure to. If None, the current axes are cleared and
        the plot is created on the full figure.
 """
    if chunk is not None:
        freq, psd, spsd, bg_level,low_cut = famed_obj.freq[chunk], famed_obj.psd[chunk], famed_obj.spsd[chunk], famed_obj.bg_level[chunk], famed_obj.low_cut_frequency[chunk]
        modes,orders,degrees,sigs = famed_obj.freqs[chunk], famed_obj.orders[chunk], famed_obj.degrees[chunk], famed_obj.freqs_sig[chunk]
    else:
        freq, psd, spsd, bg_level = famed_obj.freq, famed_obj.psd, famed_obj.spsd, famed_obj.bg_level
        modes,orders,degrees,sigs = famed_obj.freqs, famed_obj.orders, famed_obj.degrees, famed_obj.freqs_sig
    
    if ax is None:
        plt.clf()
        ax = plt.axes()
        ax.set_xlabel(r'Frequency ($\mu$Hz)')
    ax.set_ylabel(r'PSD (ppm$^2/\mu$Hz)')
    ax.semilogy(freq,psd,'-',c=famed_obj.cp.psd_psd)
    ax.semilogy(freq,spsd,'-',c=famed_obj.cp.psd_spsd,lw=2)
    ax.semilogy(freq,bg_level,'--',c=famed_obj.cp.psd_bg,lw=2)
    ax.set_xlim(min(freq),max(freq))
    ax.set_ylim(min(bg_level)*.85,max(spsd)*100)
    
    if chunk is None:
        # nu_max arrow
        ann = ax.annotate('',xy=(famed_obj.numax,.92),xycoords=transforms.blended_transform_factory(ax.transData,ax.transAxes),xytext=(famed_obj.numax,1),textcoords=transforms.blended_transform_factory(ax.transData,ax.transAxes),arrowprops=dict(color=famed_obj.cp.numax_arrow, lw=3, arrowstyle='-|>'))
        ann.arrow_patch.set_clip_box(ax.bbox)

        # Shading for global regions
        for i in range(0,len(famed_obj.separations)-1):
            if i%2 == 0:
                ax.axvspan(famed_obj.separations[i],famed_obj.separations[i+1],color=famed_obj.cp.psd_chunk1,alpha=1,ec='None')
                ax.text((famed_obj.separations[i]+famed_obj.separations[i+1])/2,.88,r'$%i$'%i,transform=transforms.blended_transform_factory(ax.transData,ax.transAxes),fontweight='heavy',fontsize='medium',color=famed_obj.cp.psd_chunkN,zorder=4,clip_on=True)
            else:
                ax.axvspan(famed_obj.separations[i],famed_obj.separations[i+1],color=famed_obj.cp.psd_chunk2,alpha=1,ec='None')
                ax.text((famed_obj.separations[i]+famed_obj.separations[i+1])/2,.88,r'$%i$'%i,transform=transforms.blended_transform_factory(ax.transData,ax.transAxes),fontweight='heavy',fontsize='medium',color=famed_obj.cp.psd_chunkN,zorder=4,clip_on=True)
    else:
        ## Lower cutoff for frequencies searched:
        if low_cut is not None:
            plt.axvline(low_cut,ls='--',color=famed_obj.cp.asef_bar)
        #### Fade-away shading of regions for l=0,1,2,3  at back
        colors = ['b','r','g','grey']
        intervals = np.linspace(.1,1,50)#**3
        if modes is not None:
            for mode,order,degree,sig in zip(modes,orders,degrees,sigs):
                for i in range(0,len(intervals)):
                    plt.axvspan(mode-intervals[i]*3*sig,mode+intervals[i]*3*sig,facecolor=colors[degree],edgecolor=colors[degree],alpha=.01,zorder=1)
        
        ### v0-dnu/2 line  on very top
        plt.axvline(famed_obj.freqs_radial_chunk[chunk]-famed_obj.best_dnu/2,ls='--',color=famed_obj.cp.text2,zorder=5)
        ax.text(famed_obj.freqs_radial_chunk[chunk]-famed_obj.best_dnu/2,.53,r'$\nu_0-\Delta\nu/2$',rotation='vertical',ha='right',fontsize='x-small',color=famed_obj.cp.text2,transform=transforms.blended_transform_factory(ax.transData,ax.transAxes),clip_on=True)
                
    # Mode lines, arrows, and id (n,l) 
    for mode,order,degree,sig in zip(modes,orders,degrees,sigs):
        ax.axvline(mode,ls=':',color=famed_obj.cp.psd_line,zorder=1)
        if degree== 0:
            if chunk is None:
                mask = np.where((freq>=mode-2*sig) & (freq<=mode+2*sig))[0]
                plt.plot(freq[mask],spsd[mask],lw=2,c=famed_obj.cp.psd_spsd0)
                ann = plt.annotate('',xy=(mode,.9),xycoords=transforms.blended_transform_factory(ax.transData,ax.transAxes),xytext=(mode,1),textcoords=transforms.blended_transform_factory(ax.transData,ax.transAxes),arrowprops=dict(color=famed_obj.cp.psd_line, lw=2, arrowstyle='-|>'))
                ann.arrow_patch.set_clip_box(ax.bbox)
        ax.text(mode,.73,'(%i,%i)'%(order,degree),rotation='vertical',ha='right',fontsize='x-small',color=famed_obj.cp.psd_modeid,transform=transforms.blended_transform_factory(ax.transData,ax.transAxes),clip_on=True)
        
def psd_inset_plot(famed_obj,ax=None,chunk=None):
    """
    Plot a small inset of the PSD centered around numax

    Parameters
    ----------
    famed_obj : FamedStar object
        A FamedStar class object that has been processed through the GLOBAL 
        module of FAMED.
    ax : matplotlib.pyplot.axes object, default: None
        Axes to plot this figure to. If None, the current axes are cleared and
        the plot is created on the full figure.
    """
    if ax is None:
        plt.clf()
        ax = plt.axes()
    mask =  (famed_obj.freq >= famed_obj.numax-famed_obj.acf_dnu/4.0*3) & (famed_obj.freq <= famed_obj.numax+famed_obj.acf_dnu/4.0*3)
    ax.semilogy(famed_obj.freq[mask],famed_obj.psd[mask],'-',lw=1,c=famed_obj.cp.psd_psd,alpha=.9)
    ax.semilogy(famed_obj.freq[mask],famed_obj.spsd[mask],'-',lw=1,c=famed_obj.cp.inset_spsd)
    ax.set_xlim(famed_obj.numax-famed_obj.acf_dnu/4.0*3,famed_obj.numax+famed_obj.acf_dnu/4.0*3)
    for mode,degree in zip(famed_obj.freqs,famed_obj.degrees):
        if (mode > min(famed_obj.freq[mask])) and (mode < max(famed_obj.freq[mask])):
            ax.text(mode,.8,'%i'%degree,color=famed_obj.cp.inset_id,transform=transforms.blended_transform_factory(ax.transData,ax.transAxes))
    plt.ylim(min(famed_obj.spsd[mask]),max(famed_obj.spsd[mask])*15)
    y_major = matplotlib.ticker.LogLocator(base = 10.0, numticks=10)
    ax.yaxis.set_major_locator(y_major)
    y_minor = matplotlib.ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
    ax.yaxis.set_minor_locator(y_minor)
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    plt.setp(ax.get_yticklabels(), visible=False)

def epsilon_plot(famed_obj,ax=None):
    """
    Plot either the Epsilon-Dnu or Epsilon-Teff relation.

    Parameters
    ----------
    famed_obj : FamedStar object
        A FamedStar class object that has been processed through the GLOBAL 
        module of FAMED.
    ax : matplotlib.pyplot.axes object, default: None
        Axes to plot this figure to. If None, the current axes are cleared and
        the plot is created on the full figure.
    """
    if ax is None:
        plt.clf()
        ax = plt.axes()
    if famed_obj.acf_dnu > famed_obj.cp.dnu_threshold:
        ax.set_ylabel(r'$\epsilon$')
        ax.set_xlabel(r'$T_{\mathrm{eff}}$')
        ax.plot(famed_obj.epsi_fit_teff_arr,famed_obj.epsi_fit_epsi_arr,'o',mfc=famed_obj.cp.epsi_dots, mec=famed_obj.cp.epsi_edge)
        ax.plot(famed_obj.epsi_fit_teff_arr,famed_obj.epsi_fit_arr,'-',c=famed_obj.cp.epsi_line)
        ax.axhline(famed_obj.epsi_fit,ls='--',color=famed_obj.cp.epsi_mark)
        ax.axvline(famed_obj.teff,ls='--',color=famed_obj.cp.epsi_mark)
        ax.plot([famed_obj.teff],[famed_obj.epsi_fit],'s',ms=8,mfc=famed_obj.cp.epsi_mark,mec=famed_obj.cp.epsi_edge)
        ax.invert_xaxis()
        ax.invert_yaxis()
        ax.text(.6,.8,r'$\epsilon = %.3f$'%famed_obj.epsi_fit,transform=ax.transAxes,color=famed_obj.cp.epsi_mark)        
    else:
        ax.set_xlabel(r'$\epsilon$')
        ax.set_ylabel(r'$\Delta\nu$ ($\mu$Hz)')
        ax.plot(famed_obj.epsi_fit_epsi_arr,famed_obj.epsi_fit_dnu_arr,'-',c=famed_obj.cp.epsi_line)
        ax.axvline(famed_obj.epsi_fit,ls='--',color=famed_obj.cp.epsi_mark)
        ax.axhline(famed_obj.acf_dnu,ls='--',color=famed_obj.cp.epsi_mark)
        ax.plot([famed_obj.epsi_fit],[famed_obj.acf_dnu],'s',ms=8,mfc=famed_obj.cp.epsi_mark,mec=famed_obj.cp.epsi_edge)
        ax.text(.2,.8,r'$\epsilon = %.3f$'%famed_obj.epsi_fit,transform=ax.transAxes,color=famed_obj.cp.epsi_mark)
        plt.xlim(min(famed_obj.epsi_fit_epsi_arr),max(famed_obj.epsi_fit_epsi_arr))
        plt.ylim(min(famed_obj.epsi_fit_dnu_arr),max(famed_obj.epsi_fit_dnu_arr))

def acf_plot(famed_obj,ax=None):
    """
    Plot the results of the autocorreltion to find `dnu`. 
    
    Parameters
    ----------
    famed_obj : FamedStar object
        A FamedStar class object that has been processed through the GLOBAL 
        module of FAMED.
    ax : matplotlib.pyplot.axes object, default: None
        Axes to plot this figure to. If None, the current axes are cleared and
        the plot is created on the full figure.
    """
    if ax is None:
        plt.clf()
        ax = plt.axes()
    ax.set_ylabel(r'ACF$^2$')
    ax.set_xlabel(r'$\Delta\nu$ ($\mu$Hz)')
    ax.plot(famed_obj.interpolated_dnu,famed_obj.interpolated_acf**2,'-',color=famed_obj.cp.acf_line,zorder=1)
    ax.axvline(famed_obj.scaling_dnu,ls=':',color=famed_obj.cp.acf_pos2,zorder=3)
    ax.axvline(famed_obj.acf_dnu,ls='--',color=famed_obj.cp.acf_pos1,zorder=2)
    plt.ylim(min(famed_obj.interpolated_acf),max(famed_obj.interpolated_acf**2)*1.15)

def asef_histogram(famed_obj,ax=None,chunk=None):
    """
    Plot the ASEF histogram and associated regions identified as modes.

    Parameters
    ----------
    famed_obj : FamedStar object
        A FamedStar class object that has been processed through the GLOBAL 
        module of FAMED.
    ax : matplotlib.pyplot.axes object, default: None
        Axes to plot this figure to. If None, the current axes are cleared and
        the plot is created on the full figure.
    """
    if chunk is not None:
        par_hist,asef_hist,maxima,ranges = famed_obj.par_hist[chunk],famed_obj.asef_hist[chunk],famed_obj.maxima[chunk],famed_obj.ranges[chunk]
        modes,errs,degrees = famed_obj.freqs[chunk],famed_obj.freqs_sig[chunk],famed_obj.degrees[chunk]
    else:
        par_hist,asef_hist,maxima,ranges = famed_obj.par_hist,famed_obj.asef_hist,famed_obj.maxima,famed_obj.ranges
        modes,errs,degrees = famed_obj.freqs,famed_obj.freqs_sig,famed_obj.degrees
                
    if ax is None:
        plt.clf()
        ax = plt.axes()
    ax.set_ylabel('ASEF (Nested iteration)')
    ax.set_xlabel(r'Frequency ($\mu$Hz)')

    # Plot histogram.
    patches = ax.bar(par_hist,asef_hist,width=np.diff(par_hist)[0],color=famed_obj.cp.asef_bar,zorder=1)
    upper_edge = plt.step(par_hist,asef_hist,color=famed_obj.cp.asef_edge,where='mid',zorder=4,lw=2)
    step_interp = interp1d(par_hist,asef_hist,'nearest')
    step_x = np.linspace(min(par_hist),max(par_hist),5000)
    step_y = step_interp(step_x)
    cm = plt.cm.get_cmap('RdYlBu')
    cm_min = min(maxima)
    cm_range = max(maxima)-cm_min
    for maximum, range_maximum in zip(maxima,ranges.T):
        lower,upper = range_maximum
        plt.axvline(maximum,ls=':',color=famed_obj.cp.asef_line,zorder=1)
        plt.fill_between(step_x,step_y,where=(step_x>lower) & (step_x<upper),color=cm((maximum-cm_min)/cm_range),zorder=2)

    patches = ax.bar(par_hist,asef_hist,width=np.diff(par_hist)[0],fill=False,edgecolor=famed_obj.cp.asef_sep,zorder=3)

    y1,y2=plt.ylim()
    plt.ylim(y1,y2*1.2)
    
    # Errorbars and modes.
    for mode,err,degree in zip(modes,errs,degrees):
        plt.errorbar(mode,.86,xerr=err,capthick=2,capsize=7,marker='o',ms=5,mfc=famed_obj.cp.asef_mark,mec=famed_obj.cp.asef_mec,ecolor=famed_obj.cp.asef_err,transform=transforms.blended_transform_factory(ax.transData,ax.transAxes))
        plt.text(mode,.92,'%i'%degree,fontweight='heavy',size='small',color=famed_obj.cp.asef_degree,transform=transforms.blended_transform_factory(ax.transData,ax.transAxes),zorder=4,clip_on=True)

    if chunk is not None:
        for mode in famed_obj.freqs_global[chunk]:
            plt.axvline(mode,ls='dashed',color=famed_obj.cp.psd_spsd)

        # Bracket showing location of octupole modes
        octu,octu_lower,octu_upper = famed_obj.octupole_freq_asymp[chunk],famed_obj.octupole_freq_lower[chunk],famed_obj.octupole_freq_upper[chunk]
        plt.axvline(octu_lower,.5,.55,zorder=5,color=famed_obj.cp.text4,lw=2)
        plt.axvline(octu_upper,.5,.55,zorder=5,color=famed_obj.cp.text4,lw=2)
        plt.plot([octu_lower,octu_upper],[.55,.55],transform=transforms.blended_transform_factory(ax.transData,ax.transAxes),zorder=5,color=famed_obj.cp.text4,lw=2)
        ann=plt.annotate('',xy=(octu,.58),xycoords=transforms.blended_transform_factory(ax.transData,ax.transAxes),xytext=(octu,.68),textcoords=transforms.blended_transform_factory(ax.transData,ax.transAxes),arrowprops=dict(color=famed_obj.cp.text4, lw=1, arrowstyle='-|>'),clip_on=True)
        #ann.arrow_patch.set_clip_box(ax.bbox
            
def text_panel(famed_obj,ax=None,chunk=None):
    """
    Print a text summary of GLOBAL results as text over and empty axes.

    Parameters
    ----------
    famed_obj : FamedStar object
        A FamedStar class object that has been processed through the GLOBAL 
        module of FAMED.
    ax : matplotlib.pyplot.axes object, default: None
        Axes to plot this figure to. If None, the current axes are cleared and
        the plot is created on the full figure.
    """
    if ax is None:
        plt.clf()
        ax = plt.axes()
  
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    # Logo
    logo = plt.imread(famed_obj.cp.famed_path/famed_obj.cp.logo_filename)
    ax.imshow(logo,extent=[0.78,.98,0,.45],aspect='auto')
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)

    if chunk is None:
        # Text 1
        ax.text(0.0,.05,r'$\mathbf{%s %s}$''\n Folder:%s\n Run:%s\n Modality:%s'%(famed_obj.catalog_id.replace('_', '\_'),famed_obj.star_id.replace('_', '\_'),famed_obj.cp.isla_subdir,famed_obj.cp.global_subdir,famed_obj.modality),fontsize='small',color=famed_obj.cp.text1)

        # Text 2 global
        if famed_obj.bad_epsi:
            ax.text(0.15,.05,r' $\nu_{\mathrm{max}}$ = %.3f $\mu$Hz''\n'r' $\Delta\nu_{\mathrm{fit}}$ = %.3f $\mu$Hz''\n'r' $\Delta\nu_{\mathrm{ACF}}$ = %.3f $\mu$Hz''\n'r' SNR = %.1f\, $\epsilon_{\mathrm{ech}}$ = %.3f'%(famed_obj.numax,famed_obj.dnu,famed_obj.acf_dnu,famed_obj.snr,famed_obj.epsilon),fontsize='small',color='red')
            ax.text(0.15,.05,r' $\nu_{\mathrm{max}}$ = %.3f $\mu$Hz''\n'r' $\Delta\nu_{\mathrm{fit}}$ = %.3f $\mu$Hz''\n'r' $\Delta\nu_{\mathrm{ACF}}$ = %.3f $\mu$Hz''\n'r' SNR = %.1f\,'%(famed_obj.numax,famed_obj.dnu,famed_obj.acf_dnu,famed_obj.snr),fontsize='small',color=famed_obj.cp.text2)
        else:
            ax.text(0.15,.05,r' $\nu_{\mathrm{max}}$ = %.3f $\mu$Hz''\n'r' $\Delta\nu_{\mathrm{fit}}$ = %.3f $\mu$Hz''\n'r' $\Delta\nu_{\mathrm{ACF}}$ = %.3f $\mu$Hz''\n'r' SNR = %.1f\, $\epsilon_{\mathrm{ech}}$ = %.3f'%(famed_obj.numax,famed_obj.dnu,famed_obj.acf_dnu,famed_obj.snr,famed_obj.epsilon),fontsize='small',color=famed_obj.cp.text2)

        # Text 3 global
        ax.text(0.36,.05,r' $\Gamma_{\mathrm{fit}}$ = %.3f $\mu$Hz''\n'r' $\Gamma_{\nu\mathrm{max}}$ = %.3f $\mu$Hz''\n'r' H$_{\mathrm{max,prior}}$ = %.1e ppm$^2/\mu$Hz''\n'r'  $\alpha$ = %.3f\, $T_{\mathrm{eff}}$ = %i K'%(famed_obj.linewidth_numax,famed_obj.linewidth_numax,famed_obj.hmax_prior,famed_obj.alpha,famed_obj.teff),fontsize='small',color=famed_obj.cp.text3)

        # Text 4
        ax.text(0.61,.05,r' ASEF$_{\mathrm{threshold}}$ = %.2f \%%''\n'r' ASEF$_{\mathrm{bins}}$ = %i''\n'r' $\Delta\nu_{\mathrm{tolerance}}$ = %.2f \%%''\n'r' N$_{\mathrm{freq}}$ = %i \,\,\,  N$_{\mathrm{orders}}$ = %i'%(100*famed_obj.cp.threshold_asef_global,famed_obj.asef_bins,100*famed_obj.dnu_tol,famed_obj.n_freq,famed_obj.n_chunks),fontsize='small',color=famed_obj.cp.text4)
    
    else:
        # Text 1
        ax.text(0.0,.05,r'$\mathbf{%s %s}$''\n Folder:%s\n Run:%s\n Modality:%s'%(famed_obj.catalog_id.replace('_', '\_'),famed_obj.star_id.replace('_', '\_'),famed_obj.cp.isla_subdir,famed_obj.chunk_number[chunk],famed_obj.modality),fontsize='small',color=famed_obj.cp.text1)

        # Text 2 chunk
        ax.text(0.15,.05,r' $\nu_{\mathrm{max}}$ = %.3f $\mu$Hz''\n'r' $\Delta\nu_{\mathrm{fit}}$ = %.3f $\mu$Hz''\n'r' $\epsilon$ = %.2f\, $\delta\nu_{02}$ = %.2f $\mu$Hz''\n'r' SNR = %.1f\, $\Delta P_1 = %.1f s$'%(famed_obj.numax,famed_obj.best_dnu,famed_obj.local_epsi[chunk],famed_obj.local_d02[chunk],famed_obj.snr[chunk],famed_obj.local_dp[chunk]),fontsize='small',color=famed_obj.cp.text2)

        # Text 3 chunk
        if famed_obj.n_radial_chunk[chunk] == 1:
            ax.text(0.36,.05,r' $\Gamma_{\mathrm{fit}}$ = %.3f $\mu$Hz''\n'r' $\Gamma_{radial}$ = %.3f $\mu$Hz''\n'r' H$_{\mathrm{max,prior}}$ = %.1e ppm$^2/\mu$Hz''\n'r'  $\alpha$ = %.3f\, $T_{\mathrm{eff}}$ = %i K'%(famed_obj.fit_linewidth[chunk],famed_obj.fwhm_radial_fit[chunk],famed_obj.upper_height[chunk],famed_obj.best_alpha,famed_obj.teff),fontsize='small',color=famed_obj.cp.text3)
        else:
            ax.text(0.36,.05,r' $\Gamma_{\mathrm{fit}}$ = %.3f $\mu$Hz''\n'r' $\Gamma_{radial}$ = %.3f $\mu$Hz''\n'r' H$_{\mathrm{max,prior}}$ = %.1e ppm$^2/\mu$Hz''\n'r'  $\alpha$ = %.3f\, $T_{\mathrm{eff}}$ = %i K'%(famed_obj.fit_linewidth[chunk],famed_obj.avg_fwhm[chunk],famed_obj.upper_height[chunk],famed_obj.best_alpha,famed_obj.teff),fontsize='small',color=famed_obj.cp.text3)
            
        # Text 4
        ax.text(0.61,.05,r' ASEF$_{\mathrm{threshold}}$ = %.2f \%%''\n'r' ASEF$_{\mathrm{bins}}$ = %i''\n'r' N$_{\mathrm{freq}}$ = %i \,\,\,  N$_{\mathrm{orders}}$ = %i'%(100*famed_obj.threshold_asef,famed_obj.asef_bins[chunk],famed_obj.n_freqs[chunk],famed_obj.n_chunks),fontsize='small',color=famed_obj.cp.text4)
