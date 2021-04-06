Getting Started
===============

To begin, you need data that has been fit using the ``Background`` module of ``DIAMONDS``, and a set of configuring parameters.

IDL Version
^^^^^^^^^^^

IDL instructions.



Python Version
^^^^^^^^^^^^^^
To test that the code is installed and running properly, you can run ``test_famed.py`` and confirm that output, both text and plots have been created. To do this:

 .. code:: shell

    $ cd /YOUR_LOCAL_PATH/FAMED/python/test/
    $ python test_famed.py

It will ask confirmation to remove existing files for the two sample stars included. This is to provide a clean working directory for the test to run in. Once the test has completed and output has been verified, you are good to begin your own stars!


The configuring parameters for the python version can be found in the files ``famed_config.yml`` and ``famed_configuring_parameters.txt``. The python version reads ``famed_config.yml`` first and has a path to ``famed_configuring_parameters.txt`` within. You can customize the path to point to specific files for different stars if you need to change the parameters, however, the default values should work for most types of stars. Additionally, the python version will look for a ``famed_config.yml`` file in your current working directory first, before using the default file under the source code directory.

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

If the configuring parameter flag ``save_progress_pickle`` is set, a pickle of the star object is saved in the results directory of each star after both the ``make_islands()`` and ``find_islands()`` functions have been run. The keyword ``load_islands`` can be set to ``True`` when creating a ``Global`` object to load the pickled data. 

 .. code :: python

     >>> star = f.Global('KIC', '006117517', 4687, load_islands=True)



     
All at once
-----------
This method does everything in the step-by-step method with a single command. This is helpful if you do not need to examine individual steps of the process and just want to get the results and output created. By default this method has ``force=True`` for both ``make_islands`` and ``find_islands``. It will only produce plots if the ``save_png`` or ``save_eps`` flags are set in the configuring parameters.

 .. code :: python
 
     >>> import famed as f
     >>> f.run.run_GLOBAL('KIC', '012069424', 5825)
   
     
