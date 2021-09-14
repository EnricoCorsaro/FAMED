# Tutorial to execute IDL-FAMED for the RGB star KIC 12008916
This tutorial shows how to run the IDL version of the FAMED pipeline using both the GLOBAL and CHUNK modules in order to extract the oscillation frequencies and mode identification of the RGB star KIC 12008916 observed by NASA Kepler for more than four years. 

Before running this tutorial make sure that the working paths have been setup properly according to your local directories. If you have installed the code using the install.sh script provided within the package this initial configuration is done automatically. Once these aspects are set you can run the tutorial using the steps listed below:

1. Copy the file `FAMED/tutorials/data/Background/data/KIC012069424.txt` inside the folder `../Background/data/`, i.e. the data folder of the Background code. If you don't have the Background code installed, make sure to create a folder named Background, to be placed at the same level as FAMED (not inside). You can even decide to install the Background code as a preliminary step. See [Background code](https://github.com/EnricoCorsaro/Background) for more information.
2. Copy the folder `FAMED/tutorials/data/Background/results/KIC012069424` inside the folder `../Background/results/`, i.e. the results folder of the Background code.
3. Go to the folder `FAMED/idl/` and open an IDL prompt from there.
4. In order to execute the first part of the analyis using the GLOBAL module, type the following commands inside the IDL prompt:

```idl
.rnew start_famed.pro
start_famed,'KIC','012008916',5825,/fit,/global
```	

5. Once the GLOBAL analysis is completed you can decide to move on by performing the CHUNK analysis. For this purpose type the following command inside the IDL prompt:

```idl
start_famed,'KIC','012008916',5825,/fit,/chunk
```	

The '/fit' option has to be activated when executing both modules for the first time because this activates the making of the multi-modal sampling with DIAMONDS. If one has to repeat the analysis starting from an existing multi-modal sampling (e.g. suppose you want to make a new test by changing any of the configuring parameters of the pipeline), then this option can be deactivated. The options `/global` and `/chunk` can also be used in combination so that the modules GLOBAL and CHUNK are executed one right after the other as soon as the former is completed, without requiring any intervention by the user. In this case, the entire analysis can be started from an IDL prompt by using the command:

```idl
start_famed,'KIC','012008916',5825,/fit,/global,/chunk
```	