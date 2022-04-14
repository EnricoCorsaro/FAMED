# Tutorial to execute IDL-FAMED for the MS star KIC 12069424

This tutorial shows how to run the IDL version of the FAMED pipeline using both the GLOBAL and CHUNK modules in order to extract the oscillation frequencies and mode identification of the MS star KIC 12069424 observed by NASA Kepler for more than four years. 

Before running this tutorial make sure that the working paths have been setup properly according to your local directories. If you have installed the code using the install.sh script provided within the package this initial configuration is done automatically. Once these aspects are set you can run the tutorial using the steps listed below:

1. Copy the file `FAMED/tutorials/data/Background/data/KIC012069424.txt` inside the folder `Background/data/`, i.e. the data folder of the Background code (note that the `Background` dir sits at the same level as the `FAMED` dir). If you don't have the Background code installed, make sure to create a folder named Background, to be placed at the same level as FAMED (not inside). You can even decide to install the Background code as a preliminary step. See [Background code](https://github.com/EnricoCorsaro/Background) for more information.
2. Copy the folder `FAMED/tutorials/data/Background/results/KIC012069424` inside the folder `Background/results/`, i.e. the results folder of the Background code.
3. Go to the folder `FAMED/idl/` and open an IDL prompt from there.
4. In order to execute the first part of the analyis using the GLOBAL module, type the following commands inside the IDL prompt:

```idl
.rnew start_famed.pro
start_famed,'KIC','012069424',5825,/fit,/global
```	

5. Once the GLOBAL analysis is completed you can decide to move on by performing the CHUNK analysis. For this purpose type the following command inside the IDL prompt:

```idl
start_famed,'KIC','012069424',5825,/fit,/chunk
```	

The '/fit' option has to be activated when executing both modules for the first time because this activates the making of the multi-modal sampling with DIAMONDS. If one has to repeat the analysis starting from an existing multi-modal sampling (e.g. suppose you want to make a new test by changing any of the configuring parameters of the pipeline), then this option can be deactivated. The options `/global` and `/chunk` can also be used in combination so that the modules GLOBAL and CHUNK are executed one right after the other as soon as the former is completed, without requiring any intervention by the user. In this case, the entire analysis can be started from an IDL prompt by using the command:

```idl
start_famed,'KIC','012069424',5825,/fit,/global,/chunk
```	

**NOTE**: you have the possibility to run this first tutorial without installing the DIAMONDS-Background code. If you plan to do so, we recommend to follow the guidelines described in the next tutorial.

# Tutorial to adopt a background fit obtained by a code different than the DIAMONDS-Background code

The file `famed_configuring_parameters.txt` contains two keywords that can be used to force FAMED to adopt a background fit solution that was not obtained using the DIAMONDS-based Background code available from [here](https://github.com/EnricoCorsaro/Background). When adopting a different software, the user has to readapt the output solution into a format that can be readable by FAMED. The main requirement in this case is that the adopted background model matches one of those implemented in the Background code. A list and description of the available background models can be found [here](https://famed.readthedocs.io/en/latest/background_models.html). Adopting a background solution from another software will not require that you install the DIAMONDS-Background code in your machine. FAMED can work properly independently of whether you have the DIAMONDS-Background code installed.

For this purpose the user has to supply three different files for each star that needs to be analyzed:
1. An ASCII file containing the dataset of the star.
2. An ASCII file containing the Nyquist frequency of the adopted dataset.
3. An ASCII file containing the solution of the background fit and the name of the model that was fitted (see below for more details). 

These three files must be placed under the same directory, which is specified by the keyword `external_background_results_dir`. This keyword, contained in the FAMED configuring parameter file `famed_configuring_parameters.txt`, has to be changed from its default value of -99 to the user specified directory (using the full path ending with `/`, see also [here](https://famed.readthedocs.io/en/latest/configuring_parameters.html)). 

Once this is done, the user can prepare the requested ASCII files in the proper format, which is described below: 
1. The file containing the dataset of the star has to be named with the label `<Catalog_ID + Star_ID>.txt` and be in a 2-column format with frequency (in muHz) and PSD (in ppm^2/muHz). The dataset contained in the first tutorial, `KIC012069424.txt` show already this file format for the star KIC 12069424.
2. The file containing the Nyquist frequency has to be named with the label `<Catalog_ID + Star_ID>_NyquistFrequency.txt` and specify the value in muHz. In the example of the tutorial, this file should be named `KIC012069424_NyquistFrequency.txt`. An example of the content of this file is provided in the first tutorial.
3. The file containing the fit solution (hereafter solution file) has to be named with the label `<Catalog_ID + Star_ID>` + `external_background_filename_suffix` and be in 1-column format. The keyword `external_background_filename_suffix` is set by default to `_backgroundParameters.txt` in the configuring parameter file of FAMED. For example, if we want to supply the solution file for the star KIC 12069424, then it has to be named `KIC012069424_backgroundParameters.txt`.

The content of the solution file is:
- Line 1. The name of the background model that was adopted, to match one of those implemented in the DIAMONDS-Background code (see [here](https://famed.readthedocs.io/en/latest/background_models.html) for a list of all the available model names).
- Lines from 2 to EOF. The estimates of each free parameter of the model as obtained from your own code but listed in the same order implemented in the DIAMONDS-Background code. See [here](https://famed.readthedocs.io/en/latest/background_models.html) for a complete description of the order of the free parameter estimates.

**NOTE**: it is mandatory that your background model incorporates the Gaussian envelope of the oscillations. This is encoded by three parameters which have to be placed as the last ones in the input file. FAMED always expect these parameters to be included because the region of the oscillations is fitted afterwards. The pipeline will not be usable for datasets that do not contain an oscillation envelope.