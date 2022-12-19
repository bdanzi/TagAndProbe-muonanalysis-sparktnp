from __future__ import print_function
import os
import glob
import numpy as np
import pandas as pd
import itertools

#import uproot

from pyspark.ml.feature import Bucketizer
from pyspark.sql import SparkSession
from pyspark.sql import functions as F
from pyspark.sql.functions import pandas_udf, PandasUDFType

from iminuit import Minuit, describe
from scipy.stats import expon
from scipy.special import wofz, erfc

from muon_definitions import *

from registry import registry

import importlib.util
import sys
import math
import json
from array import array
import ctypes
import ROOT
import tdrstyle
import CMS_lumi

from dataset_allowed_definitions import get_allowed_sub_eras, get_data_mc_sub_eras

from compare_one_job import *

from config import Configuration


#
#
# CREATE JOB FILES AND SUBMIT TO CONDOR
#
#


def compare_multiple(particle, probe, resonance, era, config, **kwargs):
  
    # Input arguments

    _baseDir = kwargs.pop('baseDir', '')
    _subera1 = kwargs.pop('subera1', '')
    _subera2 = kwargs.pop('subera2', '')
    _era2    = kwargs.pop('era2', '')

    # If era2 is not specified, it is taken as era1

    if _era2=='':
        _era2 = era

    # Get subEras and check if running full era

    subEra_dic = get_allowed_sub_eras(resonance, era)
    subEra2_dic = get_allowed_sub_eras(resonance, _era2)

    subEra_dic.remove(era)
    subEra2_dic.remove(_era2)

    _FullEra = False

    if (_subera1 in subEra_dic) and (_subera2 in subEra2_dic):
        subEras = [_subera1, _subera2]
        print("\n")
        print("-------------------------------------")
        print("-------------------------------------")
        print("Procesing " + str(_subera1) + " vs " + str(_subera2) + " datasets")
        print("-------------------------------------")
        print("-------------------------------------")
        print("\n")
    elif _subera1 == '' and _subera2 == '':
        subEras = subEra_dic
        _FullEra = True
        print("\n")
        print("-------------------------------------")
        print("-------------------------------------")
        print("Procesing complete " + str(era) + " dataset")
        print("-------------------------------------")
        print("-------------------------------------")
        print("\n")
    else:
        print("\n")
        print("\n")
        print("The subEras selected are not availible")
        print("\n")
        print("\n")
        return 0

    
    # Get paths for parquet files

    fnames = []

    use_Data = False

    if _FullEra:
        for subEra in subEras:
            tmp = list(registry.parquet(particle, probe, resonance, era, subEra))
            fnames.append(tmp)
            
            if era.split("_")[0] in subEra:
                use_Data = True
    else:
        tmp = list(registry.parquet(particle, probe, resonance, era, _subera1))
        tmp2 = list(registry.parquet(particle, probe, resonance, _era2, _subera2))
        fnames.append(tmp)
        fnames.append(tmp2)


    res = [ele for ele in fnames if ele != []]  # Delete empty lists
    fnames = [ele[0] for ele in res]


    muon_IDs = []   

    runDir = os.getcwd()  # Current dir

    config_file = Configuration(config)  # Configuration file

    efficiencies = config_file.efficiencies()  # Efficiencies to plot

    if "./"in _baseDir:
        _baseDir = runDir + '/' + _baseDir.split("./")[1]
    elif _baseDir == '':
        _baseDir = runDir
    
    
    # Check efficiencies, one for each job

    if len(efficiencies) == 1:
        muon_IDs.append(efficiencies[0][0])
    else:
        for eff_pair in efficiencies:
            if len(eff_pair) == 1:
                muon_IDs.append(eff_pair[0])
            else:
                muon_IDs.append(eff_pair[0])
                muon_IDs.append(eff_pair[1])

    muon_IDs = list(dict.fromkeys(muon_IDs))
    
    if not os.path.exists('OUTPUT'):
        os.makedirs('OUTPUT')


    if _subera1 == '':
        _subera1 = 'subera1'
    if _subera2 == '':
        _subera2 = 'subera2'


    #### Compute lumi for both cases (full era or not)

    lumi = 0


    if _FullEra:
        
        if use_Data:
            for subEra in subEras:
                if (era.split('_')[0] not in subEra) and (era not in subEra):
                    continue
                else:
                    lumi = lumi + registry.luminosity(particle, probe, resonance, era, subEra)
        else:
            lumi = -1

    else:

        subera1_isMC = _subera1 in ['DY_madgraph', 'DY_powheg', 'JPsi_pythia8']
        subera2_isMC = _subera2 in ['DY_madgraph', 'DY_powheg', 'JPsi_pythia8']
        
        if subera1_isMC and subera2_isMC:
            lumi = -1
        elif subera1_isMC:
            lumi = registry.luminosity(particle, probe, resonance, _era2, _subera2)
        elif subera2_isMC:
            lumi = registry.luminosity(particle, probe, resonance, era, _subera1)
        else:
            lumi = registry.luminosity(particle, probe, resonance, _era2, _subera2) + registry.luminosity(particle, probe, resonance, era, _subera1)

    

    # For each muon ID or reconstruction, submit a job 
    
    for muon_ID in muon_IDs:
        
      
        args = [particle, probe, resonance, era, config, muon_ID, _baseDir, _subera1, _subera2, _era2, lumi] + fnames

        if era == _era2:
            files = ['env.sh', 'tdrstyle.py', 'CMS_lumi.py', 'dataset_allowed_definitions.py', 'registry.py', 'muon_definitions.py', 
                     'compare_one_job.py', 'config.py', config, 'pileup/data/'+era.split("_")[0]+'.root', 'pileup/mc/'+era+'.root']
        else:
            files = ['env.sh', 'tdrstyle.py', 'CMS_lumi.py', 'dataset_allowed_definitions.py', 'registry.py', 'muon_definitions.py',
                     'compare_one_job.py', 'config.py', config, 'pileup/data/'+era.split("_")[0]+'.root', 'pileup/mc/'+era+'.root', 
                     'pileup/data/'+_era2.split("_")[0]+'.root', 'pileup/mc/'+_era2+'.root']

        condorTag = ''

        submit_dir = ''

        JobFlavour = '"microcentury"'

        joblist = os.path.join(
            submit_dir,
            '{}joblist_{}_{}_{}_{}{}.txt'.format(
                particle,
                probe,
                resonance,
                era,
                muon_ID,
                '_'+condorTag if condorTag != '' else ''
            ))

        
        arguments = './compare_one_job.py'
        
        for i in args:
            arguments = arguments + ' ' + str(i)


        queue = 'queue 1'
        
        output = 'OUTPUT/job_'+ muon_ID + '.$(ClusterId).$(ProcId).out'
        error = 'OUTPUT/job_'+ muon_ID + '.$(ClusterId).$(ProcId).err' 
        log = 'OUTPUT/job_'+ muon_ID + '.$(ClusterId).$(ProcId).log' 
        
        
        config_condor = '''universe    = vanilla
executable  = condor_wrapper.sh
arguments   = {arguments}
transfer_input_files = {files}
output      = {output}
error       = {error}
log         = {log}
+JobFlavour = {JobFlavour}
{queue}'''.format(
            runDir=runDir,
            arguments=arguments,
            files=','.join(files),
            output=output,
            error=error,
            log=log,
            JobFlavour=JobFlavour,
            queue=queue,)
        
        configpath = 'condor_' + muon_ID + '.sub'
        
        with open(configpath, 'w') as f:
            f.write(config_condor)
            
        print('Condor submit script written to {}'.format(configpath))
        print('To submit:')
        print('    condor_submit {}'.format(configpath))
       

        os.system('condor_submit condor_'+muon_ID + '.sub') # Submit job

