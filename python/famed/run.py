from .configuring_parameters import ConfiguringParameters
from .star_object import *
from .global_object import *
from .chunk_object import *
from .utils import *
from matplotlib.backends.backend_pdf import PdfPages


__all__ = ['GLOBAL','CHUNK','ECHELLE','COMPLETE']

def GLOBAL(catalog_id, star_id, teff, background_run_number=None, force=True, fit=True):
    """
    Helper function to run all steps, including plotting, of GLOBAL modality.

    Parameters
    ----------
    catalog_id : str
        Catalogue ID of the star (e.g. 'KIC' for Kepler).
    star_id : str
        ID of the star as a string (e.g. '0012008916' or '7037405').
    teff : float
        Effective temperature of the star in Kelvin.
    background_run_number : str or int
        Number of the background subfolder that contains the results of
        the background fit.
    force : bool
        Flag to force the computation of the sliding pattern fit.
    fit: bool
        Flag to force execution of the multi-modal sampling done by DIAMONDS.
    """
    famed_obj = Global(catalog_id,star_id,teff,background_run_number)

    famed_obj.make_islands(force=force,fit=fit)
    famed_obj.find_islands(force=force)
    famed_obj.make_global_plots()
    return famed_obj

def CHUNK(catalog_id, star_id, chunk_id=-1, background_run_number=None, force=True, fit=True):
    """
    Helper function to run all steps, including plotting, of CHUNK modality.

    Parameters
    ----------
    catalog_id : str
        Catalogue ID of the star (e.g. 'KIC' for Kepler).
    star_id : str
        ID of the star as a string (e.g. '0012008916' or '7037405').
    chunk_id : int
        The Chunk identifier number to analyze. If all the chunks must be processed, 
        leave this input set to default -1.
    background_run_number : str or int
        Number of the background subfolder that contains the results of
        the background fit.
    force : bool
        Flag to force the computation of the peak tests.
    fit: bool
        Flag to force execution of the multi-modal sampling done by DIAMONDS.
    """
    famed_obj = Chunk(catalog_id,star_id,background_run_number=background_run_number)
    result = famed_obj.make_islands(chunk_id,fit=fit)

    if result:
        if chunk_id != -1:
            print('\n\n NOW ANALYZING ONLY CHUNK: ',chunk_id,'\n\n')
            result = famed_obj.find_islands(chunk_id,force=force)
            if result:
                famed_obj.make_chunk_plots(chunk_id)
            else:
                print('CHUNK {} did not run.'.format(chunk_id))
        else:
            snr,chunks=result
            chunks = chunks[np.argsort(snr)]
            snr = snr[np.argsort(snr)]
            print(' Sorted chunks:')
            for i in range(0,len(snr))[::-1]:
                print(chunks[i],snr[i])
            pdf = PdfPages(famed_obj.star_dir/famed_obj.cp.figs_subdir/(famed_obj.catalog_id+famed_obj.star_id+'_'+famed_obj.cp.isla_subdir+'_all_CHUNK.pdf'))
            for chunk in chunks[::-1]:
                print('\n\n NOW ANALYZING CHUNK: ',chunk,'\n\n')

                result = famed_obj.find_islands(chunk,force=force)
                if result:
                    famed_obj.make_chunk_plots(chunk)
                    pdf.savefig()
                else:
                    print('CHUNK {} did not run.'.format(chunk))
            pdf.close()
    
    if famed_obj.cp.plot_total_solution and chunk_id == -1:
        famed_obj.make_chunk_plots()

    return famed_obj

def ECHELLE():
    print('This function is not yet implemented')

def COMPLETE():
    print('This function is not yet implemented')





