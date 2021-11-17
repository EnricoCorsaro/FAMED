.. _package_content:

Package Content
===============
The FAMED pipeline is available in a public GitHub repository. The repository, named FAMED, contains the package of the first two modules of the pipeline, namely GLOBAL and CHUNK. In addition, some tutorials are also provided in order to become more familiar with the pipeline running system. Currently, the pipeline is available in IDL only, but a Python version is under development too.


Downloading the Package
^^^^^^^^^^^^^^^^^^^^^^^
If you want to retrieve the code you have two options:

1. **[Recommended]** Open a Git Client (e.g. see `SourceTree <https://www.sourcetreeapp.com/>`_, available for free). In this case you will have to clone the public repository into your local machine. If you want to clone our public GitHub repository, you need to create a new repository in SourceTree by choosing the option "Clone from URL" and specifying as URL the one of the public repository of FAMED, namely https://github.com/EnricoCorsaro/FAMED. When you clone the repository, we suggest you provide as repository name **FAMED**. This option is the best way of retrieving the code because any update, even if small, will always be visible for you using the Git Client. You can then decide to pull the new change, so that the code will be updated without the need to re-download the entire package.

2. Download the code from the GitHub repository (by clicking on the download as a ZIP file button in the website, see image below). After your download is completed, you unzip the code package file, *FAMED-master.zip*, and you rename the folder FAMED-master into **FAMED**.

    .. image:: https://raw.githubusercontent.com/EnricoCorsaro/FAMED/master/docs/figures/download_zip.png
        :width: 300 px


FAMED Package
^^^^^^^^^^^^^
The content of the package is divided into different folders, which are described below:

* ``idl``: this folder contains the IDL routines that are necessary to run the pipeline under an IDL working system.

* ``python``: this folder contains the Python routines that are necessary to run the pipeline under a Python working system.

* ``tutorials``: this folder contains ready-to-use tutorials that can be promptly executed by the user. These tutorials can be useful as a starting point to run own application on different stars.

* ``docs``: this folder contains the readthedocs material for supporting documentation.
 
Additionally, three files are also included inside the FAMED package. These files are:

* **LICENSE.txt**. Contains the license of the code.
* **README.md**. Contains a short description of the code.
* **install.sh**. The installing shell script for an automatic installation of the DIAMONDS-related software, the GNUparallel tool, and for setting up the working paths used to run the FAMED pipeline.

Tutorials
^^^^^^^^^
Tutorials on how to run the pipeline can be found in the public GitHub repository `here <https://github.com/EnricoCorsaro/FAMED/tree/master/tutorials>`_.

