.. _file_system:

File System
===========
The file system of the pipeline follows a specific structure that is based on the adoption of DIAMONDS-related codes for the computation. The file system has a series of input files in order to be able to run each computation depending on which code is required. In addition the pipeline produces different outputs that are stored in the system according to different categories. An overview of the outputs is provided in the section :ref:`overview`.

The file system is structured into different levels. As shown in Figure 1, one file system is devoted to the Background folder. Here the star must be already analyzed, meaning that the background fit has to be performed and therefore the full PSD dataset of the star has to be available (see Corsaro et al. 2020, Sect. 3.1, for more details). In Figure 2 we can instead see the file system of the PeakBagging folder. The level 1 contains also the folders for Background and Asymptotic. Most of the input files and all of the outputs produced by FAMED are stored inside the PeakBagging folder. We note that the Asymptotic code has no input files inside its level 1 folder, and only contains the code that needs to be executed within the analysis.

At level 2 we find the results folder of PeakBagging and Background. This folder contains the folders of all the stars we intend to analyze.

At level 3 we find the folder of the star for both the Background and PeakBagging codes.

At level 4 we find the folder containing the fits from a single run of the Background. For the PeakBagging code instead, we find three different folders, each one related to a different modality of the PeakBagging code (see Corsaro et al. 2020, Sect. 3.3). 

At level 5 we have a series of sub-folders for each configuration of the PeakBagging code. In particular, for the folder ``isla``, we find sub-folders containing the multi-modal fits for the GLOBAL module and for the individual chunks within the CHUNK module. Then we have sub-folders containing the uni-modal peak fits for estimating the FWHM of radial and octupole modes, and for performing the model comparison for peak detection, blending, :math:`\mbox{sinc}^2` profile, rotation, and duplicity tests. Finally, we can find a series of sub-folders for fits related to the asymptotic pattern. In particular we have those containing the multi-modal fits of the sliding-pattern model, and those containing the datasets and results for the fitting of the asymptotic relation for radial modes by means of the Asymptotic code.

The yellow arrows mark the general workflow of the analysis, which shows that the first results to be produced are those of the multi-modal sampling in the ``isla`` folder. Subsequently, these results are used for the sliding-pattern fitting, the Asymptotic code fitting, and the peak testing phase.

Figure 1. The file system of the Background folder.
    .. image:: https://raw.githubusercontent.com/EnricoCorsaro/FAMED/master/docs/figures/background_file_system.png
        :width: 500 px

Figure 2. The file system of the PeakBagging folder.
    .. image:: https://raw.githubusercontent.com/EnricoCorsaro/FAMED/master/docs/figures/peakbagging_file_system.png
        :width: 600 px
