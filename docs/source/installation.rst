.. _installation:

Installation
============
The following installation guide is suited for Mac OS X users and may be subject to some variations for Linux-based OS, mostly concerning the compilation of the source files. If you encounter technical problems in installing the code, please send an `e-mail <mailto:enrico.corsaro@inaf.it>`_ specifying the error messages (possibly sending your compilation log file) and the commands you attempted to execute.

Prerequisites
^^^^^^^^^^^^^
The FAMED pipeline requires that other preliminary software is installed in your system. This software is mostly related to the `DIAMONDS Bayesian Inference code <https://github.com/EnricoCorsaro/DIAMONDS>_`, which FAMED is exploiting to perform all of its computations. Please follow the instructions below in order to get started with the installation.

Before continuiting we advise you to read the original code papers, namely the one related to DIAMONDS, `E. Corsaro & J. De Ridder 2014 Astronomy & Astrophysics, 571, 71 <https://www.aanda.org/articles/aa/abs/2014/11/aa24181-14/aa24181-14.html>`_, and the new ones related to the multi-modal sampling, `E. Corsaro 2019, Front. Astron. Space Sci, 6, 21 <https://www.frontiersin.org/articles/10.3389/fspas.2019.00021/full>`_ and the FAMED pipeline `E. Corsaro, J. M. McKeever, J. S. Kuszlewicz 2020, Astronomy & Astrophysics (in press) <https://arxiv.org/abs/2006.08245>`_.

1.
Before installing the DIAMONDS-related codes you need to install the `CMake <http://www.cmake.org/>`_ compiler, a compiler suited for C, C++ source files that is able to recognize the most suited compiler installed in your machine, depending on the platform you have. For Mac OS X it is clang, while for Linux-based OS it is gcc. For our purposes, we recommend you to install CMake 2.8 or later. You can find the dmg file of the version 2.8.12.2 (suggested for Mavericks OS) `here <http://www.cmake.org/files/v2.8/cmake-2.8.12.2-Darwin64-universal.dmg>`_, while more recent versions are required for El Captain OS or more recent OS X (we recommend CMake version 3.5.1 or later in this case). 

    .. warning:: 
        Make sure you install the CMake command line tool as well, since you need that to compile DIAMONDS via terminal. To do so, either open the CMake app and go to Tools/Install for Command Line Use if you have installed it already, or select the option during the installation phase. To avoid further compilation issues, we also recommend to update Xcode to its latest version if you are running under a Mac OSX.
	

You also have the possibility to install cmake directly from the terminal. If you are running on a Mac OS X system, then execute the command

.. code:: shell
    
    $ sudo brew install cmake

If you are running on a Unix system such as Ubuntu, then use the command

.. code:: shell

    $ sudo apt-get install cmake

Alternatively, CMake can be installed automatically with the pipeline by using the installing shell script, ``install_osx.sh`` for Mac OS X, or ``install_unix.sh`` for Unix OS, that is provided in the GitHub repository of FAMED.



2.
a)
Install git in your terminal system if this is not already installed. For this purpose you can visit `Git Download <https://git-scm.com/downloads>`_. This will allow you using the pre-configured installing shell script to run an automatic installation of all the required software as well as a setting up of the working paths to run the pipeline (see the section below). If you prefer to run a manual installation, then ignore this step and proceed with step 2b).

2.
b)
Retrieve the code package from the public GitHub repository. How to retrieve the package and a description of the content of the package are presented in the :ref:`package_content` section of this website. Then proceed with the manual installation (see the section below).


3.
In order to use the ``IDL`` version of the pipeline, you will need to have ``IDL`` installed in your OS. ``IDL`` is typically provided by your research institution under a license agreement.

   
4.
In oder to use the ``Python`` version of the pipeline, you will need to have a minimum working python distribution. The code has been built and tested on ``Python 3.8``. As an example, using ``conda``, you can create and use a working python environment with the commands

 .. code:: shell

    $ conda install -n famed_env scipy matplotlib pyaml
    $ conda activate famed_env

The directory ``FAMED/python/`` will need to be added to your ``$PYTHONPATH`` environment variable. You can type the following in your terminal or add to your shell config file, changing it to fit your local path:

For BASH:
 .. code:: shell

    $ export PYTHONPATH="$PYTHONPATH:/YOUR_LOCAL_PATH/FAMED/python"	  

For C-Shell:
 .. code:: shell

    $ setenv PYTHONPATH $PYTHONPATH':/YOUR_LOCAL_PATH/FAMED/python'
   

Shell script Installation (Mac OS X and Unix OS)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you decide to perform a shell script installation because you followed step 2a) then you need to execute the shell script ``install_osx.sh`` for Mac OS X, or ``install_unix.sh`` for Unix OS. The script is available in the GitHub repository of the pipeline, for `Mac OS X <https://github.com/EnricoCorsaro/FAMED/blob/master/install_osx.sh>`_ and for `Unix OS <https://github.com/EnricoCorsaro/FAMED/blob/master/install_unix.sh>`_.. Once you downloaded the script, place it under the main folder where you want all the software installed. Then we recommend to make it an executable by typing the terminal command (e.g. for the Mac OS version)

.. code:: shell
    
    $ chmod +x install_osx.sh

In order to start the installation from scratch, go to the directory where you want to place all the software and run the following command via terminal

.. code:: shell
    
    $ ./install_osx.sh -d -b -p -a -g

This will install the software DIAMONDS (-d), Background (-b), PeakBagging (-p), Asymptotic (-a), and the GNUparallel tool (-g) inside the folder where you ran the shell script. Additionally, the labels YOUR_LOCAL_ROOT_PATH_HERE inside the ``famed_configuring_parameters.txt`` file will be replaced with your local working path containing the DIAMONDS-related software.

The script is assuming that either curl or wget are available in your system as shell scripts to download the GNUparallel tool. We note that the ``install_osx.sh`` and ``install_unix.sh`` scripts can run using different options. If you happen to have any of the DIAMONDS, Background, and PeakBagging codes already installed, you can skip their installation by discarding the corresponding options when executing the installing shell script.

    .. warning:: 
        When installing FAMED without installing the DIAMONDS-related software, e.g. because already installed in your system, make sure that you have the latest versions of each software available in the corresponding GitHub repositories. If this is not the case, the FAMED pipeline will not be able to run.

Manual Installation (Mac OS X and Unix OS)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The manual installation requires a number of steps, which may take some time to be accomplished. We usually recommend using the shell script installation, but if you are already more familiar with the installing process of the DIAMONDS-related software, then this can also be accomplished quite easily.

1. Once the package is downloaded because you followed step 2b), you will have to install the codes DIAMONDS, Background, PeakBagging and Asymptotic. The requirement is that their corresponding folders will have to be placed under a single common directory. For installing instructions of these codes please visit their GitHub repositories at
`DIAMONDS <https://github.com/EnricoCorsaro/DIAMONDS>`_,
`Background <https://github.com/EnricoCorsaro/Background>`_,
`PeakBagging <https://github.com/EnricoCorsaro/PeakBagging>`_,
`Asymptotic <https://github.com/EnricoCorsaro/Asymptotic>`_. 

    .. warning:: 
        The FAMED package has to be placed inside the same main directory containing the codes DIAMONDS, Background, PeakBagging, and Asymptotic.


2. After Asymptotic is installed, make sure that its ``localPath.txt`` file, inside the ``Asymptotic/build/`` directory, contains the same path used for ``localPath.txt`` of the PeakBagging code. This is because the output files produced by Asymptotic will go into the PeakBagging file system. 

3. By completing the installation of the DIAMONDS-related software, you need to install the GNUparallel tool as a shell tool. For detailed instructions please visit `GNUparallel <https://www.gnu.org/software/parallel/>`_.

4. As a last step, you need to configure the working paths in your ``famed_configuring_parameters.txt`` file. For this purpose, open the file located under the ``FAMED/idl/`` directory of the FAMED package and replace the YOUR_LOCAL_ROOT_PATH_HERE labels with your actual local path containing the FAMED package. For more details please check the description of the configuring parameters presented in the :ref:`configuring_parameters` section of this website.

Windows OS 10
^^^^^^^^^^^^^
For Windows OS 10 we recommend using the free application for creating an Ubuntu virtual machine. For details on how to set up this environment, visit `Install Ubuntu on Windows 10 <https://ubuntu.com/tutorials/tutorial-ubuntu-on-windows#1-overview>`_. 

Once the Ubuntu VM is installed and running in Windows OS, simply follow the guidlines presented in the Linux OS section of this page. You can even decide to use the shell script installation with the ``install_unix.sh`` script inside the Ubuntu VM, making sure to have the basic ubuntu packages installed, which include the GCC compiler suite.
