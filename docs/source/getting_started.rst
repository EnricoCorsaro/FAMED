Getting Started
===============


IDL Version
^^^^^^^^^^^

IDL instructions.



Python Version
^^^^^^^^^^^^^^

To begin, you need data that has been fit using the ``Background`` module of ``DIAMONDS``.

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

To force ``find_islands`` to recompute the sliding pattern fit we can use the ``force`` option:
 
 .. code :: python

     >>> star.find_islands(force=True)


 
     
All at once
-----------
This method does everything in the step-by-step method with a single command. Helpful if you do not need to examine individual steps of the process and just want to get a results and output created. By default this method has ``force=True``.

 .. code :: python
 
     >>> import famed as f
     >>> f.run.run_GLOBAL('KIC', '006117517', 4687)
   
     
