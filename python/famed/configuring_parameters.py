"""Defines a class to hold environment variables / configuring parameters determined from an input file"""

import numpy as np
from pathlib import Path
import yaml
import os

famed_path = Path(os.path.dirname(os.path.abspath(__file__)))


class ConfiguringParameters(object):
    """
    Class to hold fixed environment variables and configuring parameters.

    Contents of the class are all variables that have been set in the 
    `famed_configuing_parameters.txt` file (default name). Also included in 
    this file are the settings from the `config.yml` file. See the files and
    the online documentation for a list and description of all these quantities.

    Parameters
    ----------

    Returns
    -------
    """
    def __init__(self):

        # Using config file to setup plotting and paths
        # PATH FOR CONFIG FILE IS FIXED! THIS IS BAD AND MUST CHANGE!

        #Try to read local config file first, then default to famed_path.
        self.local_config = False
        try:
            with open('famed_config.yml', 'r') as f:
                params = yaml.safe_load(f)
                self.local_config = True
        except:
            with open(famed_path/'famed_config.yml', 'r') as f:
                params = yaml.safe_load(f)

        # Setup path to configuring parameters
        paths = params['PATHS']

        # Read in configuring parameters from a fixed filename
        self.configuring_parameters_file = paths['configuring_parameters_file']
        with open(self.configuring_parameters_file) as f:
            d = dict(x.rstrip().split(None,1) for x in f if not x.startswith('#'))

        # Set attributes from file
        for k,v in d.items():
            try:
                d[k] = int(v)
            except:
                try:
                    d[k] = float(v)
                except:
                    pass
            setattr(self, k, d[k])
        
       
        # Set up plotting parameters
        plotting = params['PLOTTING']
        self.plot_mode = plotting['plot_mode']
        for key in plotting[self.plot_mode]:
            setattr(self, key, plotting[self.plot_mode][key])
        
        # Set additional filepaths and filenames
        self.diamonds_path = Path(self.root_path)
        self.famed_path = famed_path
        
        self.peakbagging_filename_label = '_peakbagging_'
        self.peakbagging_filename_global_label = self.global_subdir+'.txt'
        self.peakbagging_filename_chunk_label = 'chunk_'
        
        # Set DIAMONDS parameters for different modalities
        diamonds = params['DIAMONDS']
        
        self.dp_isla = dict()
        for key in diamonds['isla']:
            self.dp_isla[key] = diamonds['isla'][key]
        self.dp_pb = dict()
        for key in diamonds['pb']:
            self.dp_pb[key] = diamonds['pb'][key]
        self.dp_slid = dict()
        for key in diamonds['slid']:
            self.dp_slid[key] = diamonds['slid'][key]

            

