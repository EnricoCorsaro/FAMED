FAMED - Detailed asteroseismic analysis pipeline
================================================

.. image:: https://img.shields.io/badge/GitHub-FAMED-yellow
    :target: https://github.com/EnricoCorsaro/FAMED
.. image:: https://img.shields.io/badge/license-CC%20BY--SA-blue
    :target: https://github.com/EnricoCorsaro/FAMED/blob/master/LICENSE.txt
.. image:: https://img.shields.io/badge/arXiv-2006.08245-blueviolet
    :target: https://arxiv.org/abs/2006.08245
.. image:: https://readthedocs.org/projects/famed/badge/?version=latest
    :target: https://famed.readthedocs.io/en/latest/?badge=latest
.. image:: https://img.shields.io/github/issues-closed/EnricoCorsaro/FAMED
    :target: https://github.com/EnricoCorsaro/FAMED/issues
.. image:: https://raw.githubusercontent.com/EnricoCorsaro/FAMED/master/docs/figures/FAMED_LOGO_WHITE.jpg
    :width: 500 px

Author
^^^^^^
- `Enrico Corsaro <mailto:enrico.corsaro@inaf.it>`_

Python version development
^^^^^^^^^^^^^^^^^^^^^^^^^^
- `Jean McKeever <mailto:jean.mckeever@yale.edu>`_
- `James Kuszlewicz <mailto:kuszlewicz@mps.mpg.de>`_

Description
^^^^^^^^^^^
The **FAMED** (Fast and AutoMated pEak bagging with Diamonds) pipeline is a multi-platform parallelized software to perform and automated extraction and mode identification of oscillation frequencies for solar-like pulsators. This pipeline is based on the free code DIAMONDS for Bayesian parameter estimation and model comparison by means of the nested sampling Monte Carlo (NSMC) algorithm. The pipeline can be applied to a large variety of stars, ranging from hot F-type main sequence, subgiants, up to stars evolving along the red giant branch, settled into the core-Helium-burning main sequence, and even evolved beyond towards the early asymptotic giant branch.
The pipeline is organized in separate modules, each one performing different tasks in a different level of detail. The current version of FAMED includes two out of four modules, which are named GLOBAL and CHUNK. These two first modules are available in ``IDL`` version, and are being developed in ``Python`` too. The pipeline requires some system prerequisites and a dedicated installation procedure which can be found in the documentation below.

Navigation
^^^^^^^^^^
.. toctree::
   :maxdepth: 2
   
   package_content
   installation
   file_system
   configuring_parameters
   overview
   getting_started
   publications
   events
   logo