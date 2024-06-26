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

The ``/fit`` option has to be activated when executing a module for the first time because this turns on the making of the multi-modal sampling with DIAMONDS. If one has to repeat the analysis starting from an existing multi-modal sampling (e.g. suppose you want to make a new test by changing any of the configuring parameters of the pipeline), then this option can be deactivated. 

All at once
-----------
The options ``/global`` and ``/chunk`` of the ``start_famed`` procedure can thus be used in combination so that the modules GLOBAL and CHUNK are executed one right after the other as soon as the former is completed, without requiring any intervention by the user. In this case, the entire analysis can be started from an IDL prompt by using the command:

.. code :: idl

    IDL> start_famed,'KIC','012069424',5825,/fit,/global,/chunk

The all at once method offers an additional functionality. If desired you can specify an input integer or a string to change at runtime the reference subfolder containing the background fit solution obtained from the DIAMONDS+Background code. By default this subfolder is read from the input ``famed_configuring_parameters.txt`` file (where it is set to 00). This is possible once again by means of the ``start_famed`` procedure, using the following calling sequence:

.. code :: idl

    IDL> start_famed,'KIC','012069424',5825, background_run_number=10, /fit,/global,/chunk

where, in the example provided, the background fit solution will have to be contained within the subfolder 10 (instead of the default 00). This option can be useful if one wants to adopt and test the outcomes of the analysis for different background solutions for the same star, or in case a large number of targets (each one with a different subfolder) is being considered.

Using an external background fit solution
-----------------------------------------
You have the possibility to run the pipeline by using a background fit solution that was obtained by means of a code different than the DIAMONDS-Background. This is because FAMED is capable of working without the DIAMONDS-Background code installed. If you want to test this functionality you can run the related test as follows:

.. code :: shell

    $ cd /YOUR_LOCAL_PATH/FAMED/idl/
    $ idl

Once your IDL prompt is loaded, then run:

.. code :: idl

    IDL> .rnew test_famed
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

The ``force`` option overwrites any existing solution that was obtained in a previous run. With this interactive method, you can change specific configuring parameters and recompute just the steps that you need to.

If the input configuring parameter ``save_progress_pickle`` is set to 1, a pickle of the star object is saved in the results directory of each star after both the ``make_islands()`` and ``find_islands()`` functions have been run. The keyword ``load_islands`` can be set to ``True`` when creating a ``Global`` object to load the pickled data. 

.. code :: python

     >>> star = f.Global('KIC', '006117517', 4687, load_islands=True)

Similarly to the case of the IDL version, you can also decide to change at runtime the reference subfolder containing the background fit solution obtained with the DIAMONDS+Background code. For this purpose the general calling sequence of the GLOBAL module becomes:

.. code :: python

     >>> star = f.Global('KIC', '006117517', 4687, 10)

where in this example the last input, here set to 10, represents the subfolder name that will be considered when reading the background fit solution (instead of the default value 00 read from the input ``famed_configuring_parameters.txt`` file).

All at once
-----------
This method does everything in the step-by-step method with a single command. This is helpful if you do not need to examine individual steps of the process and just want to get the results and output created. By default this method has the option ``force=True`` for both ``make_islands`` and ``find_islands``. It will only produce plots if the ``save_png`` or ``save_eps`` flags are set in the configuring parameters.

.. code :: python

     >>> import famed as f
     >>> f.run.GLOBAL('KIC', '012069424', 5825)

In this calling sequence there is also the option ``fit=True`` set by default. The option ``fit`` executes the ``make_islands`` method to evaluate the multi-modal sampling. If repeating the multi-modal fit is not necessary, for example because one wants to re-execute the GLOBAL module only to attempt at improving the identification of the l=0,1 pairs, hence to repeat the sliding-pattern fit, we recommend using the following calling sequence:

.. code :: python

     >>> import famed as f
     >>> f.run.GLOBAL('KIC', '012069424', 5825, fit=False)

In addition, in this method one can decide to change at runtime the reference subfolder containing the background fit solution obtained with the DIAMONDS+Background code. The default subfolder is the one specified in the input configuring parameter file under the keyword ``background_run_number`` (normally set to ``00``). The general calling sequence to read at runtime the background fit solution contained in a different subfolder, e.g. the subfolder ``10``, thus becomes:

.. code :: python

     >>> f.run.GLOBAL('KIC', '012069424', 5825, background_run_number=10)

or simply

.. code :: python

     >>> f.run.GLOBAL('KIC', '012069424', 5825, 10)

For the analysis performed by the CHUNK module we recommend using the all at once method. Once the GLOBAL module has been completed, following the example above the user can activate the CHUNK with the following commands:

.. code ::

    >>> f.run.CHUNK('KIC', '012069424')

where we note that the input temperature is no longer required because this is obtained directly from the solution of the GLOBAL module. Similarly to the case of GLOBAL, the calling sequence for CHUNK also allows two optional flags, both activated by default, ``fit=True``, ``force=True``. The two options have an analogous meaning to that of GLOBAL. If repeating the multi-modal fit is not necessary, for example because one wants to re-execute the CHUNK module just for a better identification of the individual modes and/or to repeat the peak detection tests, then the following calling sequence should be used:

.. code ::

    >>> f.run.CHUNK('KIC', '012069424', fit=False)

There is also the possibility to recompute the analysis of an individual chunk instead of all the chunks together, for example because we only want to intervene in a specific region of the dataset. In this case, if for example one wants to re-analyze chunk #5, the following calling sequence should be used:

.. code ::

    >>> f.run.CHUNK('KIC', '012069424', chunk_id=5)

where ``chunk_id`` would be otherwise set to -1 by default, with -1 meaning that all the identified chunks should be processed. The number of the chunk to adopt matches the numbering provided by the GLOBAL module.
Similarly to GLOBAL, if a background fit solution has to be specified at runtime, the user can force its reading also in the CHUNK by means of the following calling sequence:

.. code ::

    >>> f.run.CHUNK('KIC', '012069424', background_run_number=10)  

Note that these are positional arguments, so in this case the argument name ``background_run_number`` has to be specified entirely when providing it as an explicit input, otherwise the code may confuse it with the other input ``chunk_id`` if the latter one is not specified instead. In summary, if one wants to re-execute the CHUNK analysis for the chunk number 5 by using the background fit solution contained in the sub-folder ``10`` but without repeating the multi-modal fit, then the user should adopt the following calling sequence:

.. code ::

    >>> f.run.CHUNK('KIC', '012069424', chunk_id=5, background_run_number=10, fit=False)  

Using an external background fit solution
-----------------------------------------
An additional routine is provided in order to let you test the analysis of the RGB star KIC 12008916 by means of a background fit solution obtained by a code different than the DIAMONDS-Background. You can run this test in an interactive way following the steps below:

.. code :: shell

    $ cd /YOUR_LOCAL_PATH/FAMED/python/test/

then open a Python prompt and from there:

.. code :: python

     >>> from test_famed import *
     >>> test_external_background()

The test is also automatically executed when running the test_famed.py script directly through Python with

.. code :: shell

    $ python test_famed.py

The input files contained in the folder ``/YOUR_LOCAL_PATH/FAMED/tutorials/data/Background/KIC012008916/`` provide an example of the format needed to execute an analysis with FAMED that is using an external background fit solution. For more information about how to set up the input files when using an external background fit solution see the tutorial #2 provided `here <https://github.com/EnricoCorsaro/FAMED/tree/master/tutorials>`_.