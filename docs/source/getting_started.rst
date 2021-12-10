Getting Started
===============

In this section we provide a simple way to test that the FAMED pipeline is working properly. 

IDL Version
^^^^^^^^^^^
If you are using IDL in your system you can test that the pipeline is correctly installed by running the code and checking that the outputs (both text and plots) are created. This test is a run all at once mode, meaning that the pipeline will perform automatically each step of the analysis without your supervision. To run the test:

.. code :: shell

    $ cd /YOUR_LOCAL_PATH/FAMED/idl/
    $ idl

Once your IDL prompt is loaded, you can run the following commands:

.. code :: idl
    IDL> .rnew test_famed
    IDL> test_famed

This routine will execute the analysis of the MS star KIC 12069424 by performing both the GLOBAL and CHUNK modulesn in a step-by-step mode. All the output files will be produced inside the directory ``/YOUR_LOCAL_PATH/PeakBagging/results/KIC012069424/``.

Step-by-step
------------
You can run the GLOBAL module of the pipeline in a setp-by-step mode for the MS star KIC 12069424 (or for the RGB star KIC 6117517), making sure to follow each step in a sequential order as explained in this section. 

First, you need to perform some preliminary steps in order to correctly place your input data files. 

1. Copy the file ``/YOUR_LOCAL_PATH/FAMED/tutorials/data/Background/data/KIC012069424.txt`` inside the folder ``/YOUR_LOCAL_PATH/Background/data/``, i.e. the data folder of the Background code. If you don't have the Background code installed, make sure to create a folder named ``Background``, to be placed at the same level as FAMED (not inside). You can even decide to install the Background code as a preliminary step. See the `Background code <https://github.com/EnricoCorsaro/Background>`_ documentation for more information.

2. Copy the folder ``/YOUR_LOCAL_PATH/FAMED/tutorials/data/Background/results/KIC012069424`` inside the folder ``/YOUR_LOCAL_PATH/Background/results/``, i.e. the results folder of the Background code.

3. Go to the folder ``/YOUR_LOCAL_PATH/FAMED/idl`` and open an IDL prompt from there.

In order to execute the first part of the analyis using the GLOBAL module, type the following commands inside the IDL prompt:

.. code :: idl

    IDL> .rnew make_islands_global
    IDL> make_islands_global,'KIC','012069424',5825,/force
    IDL> .rnew find_islands_global
    IDL> find_islands_global,'KIC','012069424',5825,/force

The ``force`` option will force FAMED to repeat the entire peakbagging setting up (thus overriding any previous result) and it is convenient to adopt when repeating the analysis if one has been performed already.

Once the GLOBAL analysis is completed you can decide to move on by performing the CHUNK analysis. For this purpose type the following command inside the IDL prompt:

.. code :: idl

    IDL> start_famed,'KIC','012069424',5825,/fit,/chunk

The ``/fit`` option has to be activated when executing a module for the first time because this turns on the making of the multi-modal sampling with DIAMONDS. If one has to repeat the analysis starting from an existing multi-modal sampling (e.g. suppose you want to make a new test by changing any of the configuring parameters of the pipeline), then this option can be deactivated. The options ``/global`` and ``/chunk`` can also be used in combination so that the modules GLOBAL and CHUNK are executed one right after the other as soon as the former is completed, without requiring any intervention by the user. In this case, the entire analysis can be started from an IDL prompt by using the command:

.. code :: idl

    IDL> start_famed,'KIC','012069424',5825,/fit,/global,/chunk

Using an external background fit solution
-----------------------------------------
You have the possibility to run the pipeline by using a background fit solution that was obtained by means of a code different than the DIAMONDS-Background. This is because FAMED is capable of working without the DIAMONDS-Background code installed. If you want to test this functionality you can run the related test as follows:

.. code :: shell

    $ cd /YOUR_LOCAL_PATH/FAMED/idl/
    $ idl

Once your IDL prompt is loaded, then run:

.. code :: idl

    IDL> .rnew test_famed.pro
    IDL> test_external_background

The pipeline will make use of the background fit solution for the RGB star KIC 12008916, which is contained in the folder ``/YOUR_LOCAL_PATH/FAMED/tutorials/data/Background/results/KIC012008916/`` already in a format readable by FAMED. For more information about how to set up the input files when using an external background fit solution see the tutorial #2 provided `here <https://github.com/EnricoCorsaro/FAMED/tree/master/tutorials>`_.


Python Version
^^^^^^^^^^^^^^
To test that the code is installed and running properly, you can run ``test_famed.py`` and confirm that output, both text and plots have been created. To do this:

 .. code :: shell

    $ cd /YOUR_LOCAL_PATH/FAMED/python/test/
    $ python test_famed.py

It will ask confirmation to remove existing files for the two of the sample stars included. This is to provide a clean working directory for the test to run in. Once the test has completed and output has been verified, you are good to begin your own stars!

The configuring parameters for the Python version can be found in the files ``famed_config.yml`` and ``famed_configuring_parameters.txt``. The Python version reads ``famed_config.yml`` first and has a path to ``famed_configuring_parameters.txt`` within. You can customize the path to point to specific files for different stars if you need to change the parameters. However, the default values should work for most types of stars. Additionally, the Python version will look for a ``famed_config.yml`` file in your current working directory first, before using the default file under the source code directory.

There are two ways to perform the ``GLOBAL`` computations on your data: 1) Step-by-step and 2) All at once. 

Step-by-step
------------
In this method we load the data in, create the initial islands sampling, identifiy modes from the islands sampling, and plot the results as separate steps.  The required input are a catalog id, star id, and effective temperature. 

 .. code :: python

     >>> import famed as f
     >>> star = f.Global('KIC', '006117517', 4687)
     >>> star.make_islands()
     >>> star.find_islands()
     >>> star.make_global_plots()


To force ``make_islands`` to generate a data set with new configuring parameters we can use the ``force`` option:
 
 .. code :: python

     >>> star.make_islands(force=True)

To force ``find_islands`` to recompute the sliding pattern fit we can use the ``force`` option:
 
 .. code :: python

     >>> star.find_islands(force=True)

With this interactive method, you can change specific configuring parameters and recompute just the steps that you need to.

If the input configuring parameter ``save_progress_pickle`` is set to 1, a pickle of the star object is saved in the results directory of each star after both the ``make_islands()`` and ``find_islands()`` functions have been run. The keyword ``load_islands`` can be set to ``True`` when creating a ``Global`` object to load the pickled data. 

 .. code :: python

     >>> star = f.Global('KIC', '006117517', 4687, load_islands=True)

All at once
-----------
This method does everything in the step-by-step method with a single command. This is helpful if you do not need to examine individual steps of the process and just want to get the results and output created. By default this method has ``force=True`` for both ``make_islands`` and ``find_islands``. It will only produce plots if the ``save_png`` or ``save_eps`` flags are set in the configuring parameters.

 .. code :: python

     >>> import famed as f
     >>> f.run.GLOBAL('KIC', '012069424', 5825)

Using an external background fit solution
-----------------------------------------
An additional routine is provided in order to let you test the analysis of the RGB star KIC 12008916 by means of a background fit solution obtained by a code different than the DIAMONDS-Background. You can run this test in an interactive way following the steps below:

 .. code:: shell

    $ cd /YOUR_LOCAL_PATH/FAMED/python/test/

then open a Python prompt and from there:

 .. code :: python

     >>> from test_famed import *
     >>> test_external_background()

The test is also automatically executed when running the test_famed.py script directly through Python with

 .. code :: shell

    $ python test_famed.py

The input files contained in the folder ``/YOUR_LOCAL_PATH/FAMED/tutorials/data/Background/KIC012008916/`` provide an example of the format needed to execute an analysis with FAMED that is using an external background fit solution. For more information about how to set up the input files when using an external background fit solution see the tutorial #2 provided `here <https://github.com/EnricoCorsaro/FAMED/tree/master/tutorials>`_.