#!/usr/bin/env python
from __future__ import print_function
import os
import glob
import numpy as np
import pandas as pd
import itertools

import uproot

from pyspark.ml.feature import Bucketizer
from pyspark.sql import SparkSession
from pyspark.sql import functions as F
from pyspark.sql.functions import pandas_udf, PandasUDFType

from iminuit import Minuit, describe
from scipy.stats import expon
from scipy.special import wofz, erfc

from muon_definitions import *

import importlib.util
import sys
import math
import json
from array import array
import ctypes
import ROOT
import tdrstyle
import CMS_lumi

from config import Configuration

from dataset_allowed_definitions import get_allowed_sub_eras, get_data_mc_sub_eras



# Get Pile Up. 
# Modified from muon_definitions.py to be able to run it from condor. Pile Up files are passed to condor as input files so paths are different.

def get_pileup_condor(resonance, era, subEra):
   '''
   Get the pileup distribution scalefactors to apply to simulation
   for a given era.
   '''
   # get the pileup
   dataPileup = {
       # Note: for now use ReReco version of pileup
       # TODO: need to redo splitting by 2016 B-F/F-H
       'Run2016_UL_HIPM': 'Run2016.root',
       'Run2016_UL': 'Run2016.root',
       'Run2017_UL': 'Run2017.root',
       'Run2018_UL': 'Run2018.root',
       'Run2016': 'Run2016.root',
       'Run2017': 'Run2017.root',
       'Run2018': 'Run2018.root',
       #'Run2022': 'Run2018.root'
       'Run2022': 'Run2022.root'
   }
   mcPileup = {
       # TODO: do the two eras have different profiles?
       'Run2016_UL_HIPM': 'Run2016_UL.root',
       'Run2016_UL': 'Run2016_UL.root',
       'Run2017_UL': 'Run2017_UL.root',
       'Run2018_UL': 'Run2018_UL.root',
       'Run2016': 'Run2016.root',
       'Run2017': 'Run2017.root',
       'Run2018': 'Run2018.root',
       #'Run2022': 'Run2018_UL.root'
       'Run2022': 'Run2022.root'
   }
   # get absolute path
   baseDir = os.path.dirname(__file__)
   dataPileup = {k: os.path.join(baseDir, dataPileup[k]) for k in dataPileup}
   mcPileup = {k: os.path.join(baseDir, mcPileup[k]) for k in mcPileup}
   with uproot.open(dataPileup[era]) as f:
       data_edges = f['pileup'].edges
       data_pileup = f['pileup'].values
       data_pileup /= sum(data_pileup)
   with uproot.open(mcPileup[era]) as f:
       mc_edges = f['pileup'].edges
       mc_pileup = f['pileup'].values
       mc_pileup /= sum(mc_pileup)
   pileup_edges = data_edges if len(data_edges) < len(mc_edges) else mc_edges
   pileup_ratio = [d/m if m else 1.0 for d, m in zip(
       data_pileup[:len(pileup_edges)-1], mc_pileup[:len(pileup_edges)-1])]

   return pileup_ratio, pileup_edges


# Get weighted dataframe
# Modified to use the pile up files from condor
# From muon_definitions.py

def get_weighted_dataframe_condor(df, doGen, resonance, era, subEra, shift=None):
   '''
   Produces a dataframe with a weight and weight2 column
   with weight corresponding to:
       1 for data
   or
       pileup for mc
   The optional shift parameter allows for a different
   systematic shift to the weights
   '''
   # TODO: implement systematic shifts in the weight such as PDF, pileup, etc.
   # get the pileup
   pileup_ratio, pileup_edges = get_pileup_condor(resonance, era, subEra)

   # build the weights (pileup for MC)
   # TODO: if there is a weight column (ie, gen weight) get that first
   if doGen:
       pileupMap = {e: r for e, r in zip(pileup_edges[:-1], pileup_ratio)}
       mapping_expr = F.create_map(
           [F.lit(x) for x in itertools.chain(*pileupMap.items())])
       # M.Oh: temporary solution for missing true PU branch in the new ntuples
       if 'pair_truePileUp' in df.columns:
           weightedDF = df.withColumn(
               'PUweight', mapping_expr.getItem(F.round('pair_truePileUp')))
       elif 'nTrueInteractions' in df.columns:
           weightedDF = df.withColumn(
               'PUweight', mapping_expr.getItem(F.round('nTrueInteractions')))
       elif 'nVertices' in df.columns:
           weightedDF = df.withColumn(
               'PUweight', mapping_expr.getItem(F.col('nVertices')))
       else:
           weightedDF = df.withColumn('PUweight', F.lit(1.0))
       # apply gen weights
       if 'genWeight' in weightedDF.columns:
           weightedDF = weightedDF.withColumn('genWeightSign', F.signum('genWeight'))
           weightedDF = weightedDF.withColumn('weight', F.col('genWeightSign') * F.col('PUweight'))
       elif 'pair_genWeight' in weightedDF.columns:
           weightedDF = weightedDF.withColumn('genWeightSign', F.signum('pair_genWeight'))
           weightedDF = weightedDF.withColumn('weight', F.col('genWeightSign') * F.col('PUweight'))
       else:
           weightedDF = weightedDF.withColumn('weight', F.col('PUweight'))
   else:
       weightedDF = df.withColumn('weight', F.lit(1.0))
   weightedDF = weightedDF.withColumn(
       'weight2', F.col('weight') * F.col('weight'))

   return weightedDF



# Get Pile Up ratio for Data/Data plots

def get_data_pileup(era, era2):
   '''                                                                                                                                                                                                      
   Get the pileup distribution scalefactors to apply to simulation                                                                                                                                          
   for a given era.                                                                                                                                                                                         
   '''
   # get the pileup                                                                                                                                                                                         
   dataPileup = {
       # Note: for now use ReReco version of pileup                                                                                                                                                         
       # TODO: need to redo splitting by 2016 B-F/F-H                                                                                                                                                       
       'Run2016_UL_HIPM': 'Run2016.root',
       'Run2016_UL': 'Run2016.root',
       'Run2017_UL': 'Run2017.root',
       'Run2018_UL': 'Run2018.root',
       'Run2016': 'Run2016.root',
       'Run2017': 'Run2017.root',
       'Run2018': 'Run2018.root'
       'Run2022': 'Run2022.root'
   }

   # get absolute path                                                                                                                                                                                      
   baseDir = os.path.dirname(__file__)
   dataPileup = {k: os.path.join(baseDir, dataPileup[k]) for k in dataPileup}
   with uproot.open(dataPileup[era]) as f:
       data1_edges = f['pileup'].edges
       data1_pileup = f['pileup'].values
       data1_pileup /= sum(data1_pileup)
   with uproot.open(dataPileup[era2]) as f:
       data2_edges = f['pileup'].edges
       data2_pileup = f['pileup'].values
       data2_pileup /= sum(data2_pileup)

   pileup_edges = data1_edges if len(data1_edges) < len(data2_edges) else data2_edges
   pileup_ratio = [d1/d2 if d2 else 1.0 for d1, d2 in zip(
       data1_pileup[:len(pileup_edges)-1], data2_pileup[:len(pileup_edges)-1])]


   return pileup_ratio, pileup_edges


# Get weighted dataframe
# Modified to use the pile up files from condor
# From muon_definitions.py
# By this moment, nTrueInteractions and nPUInteractions not implemented in the ntuples (-1 value for all events), nVertices is used instead.

def get_weighted_data(df, era, era2, shift=None):
                                                                                                                            
   # get the pileup                                                                                                                                                                                         
   pileup_ratio, pileup_edges = get_data_pileup(era, era2)

   # build the weights (pileup for Data2)                                                                                                                                               

   pileupMap = {e: r for e, r in zip(pileup_edges[:-1], pileup_ratio)}
   mapping_expr = F.create_map(
       [F.lit(x) for x in itertools.chain(*pileupMap.items())])
   
   # M.Oh: temporary solution for missing true PU branch in the new ntuples                                                                                                                                 
   if 'pair_truePileUp' in df.columns:
      weightedDF = df.withColumn('PUweight', mapping_expr.getItem(F.round('pair_truePileUp')))
   #elif 'nTrueInteractions' in df.columns:                                                                                                                                                                 
   #   weightedDF = df.withColumn('PUweight', mapping_expr.getItem(F.round('nTrueInteractions')))                                                                                                           
   #   test = weightedDF.withColumn('PU', F.round('nTrueInteractions'))                                                                                                                                     
   #   test = test.select("PU")                                                                                                                                                                             
   #   test.show()                                                                                                                                                                                          
   elif 'nVertices' in df.columns:
      weightedDF = df.withColumn('PUweight', mapping_expr.getItem(F.col('nVertices')))
   else:
      weightedDF = df.withColumn('PUweight', F.lit(1.0))


   weightedDF = weightedDF.withColumn('weight', F.col('PUweight'))
   weightedDF = weightedDF.withColumn('weight2', F.col('weight') * F.col('weight'))

   return weightedDF


#
# Run_files. It gets the parquet files, process it and generates pandas dataframes from binned ones using spark
#

def run_files(particle, probe, resonance, era, subEra, config, spark, muon_ID, doDataRew, era1, fnames):
    

    # Load parquet files (or root)
    print('Loading parquet files:', fnames)
    if isinstance(fnames, list):
        baseDF = spark.read.parquet(*fnames)
    else:
        baseDF = spark.read.parquet(fnames)
    
    # Load definitions and filter events
    
    doGen = subEra in ['DY_madgraph', 'DY_powheg', 'JPsi_pythia8']
    
    #definitions = config['definitions']
    definitions = config.definitions()

    defDF = baseDF
    for d in definitions:
        defDF = defDF.withColumn(d, F.expr(definitions[d]))

    
    #tagsDF = defDF.filter(config['selection'])
    tagsDF = defDF.filter(config.selection() +' and ' + muon_ID)
    
    if doGen:
        if 'mc_selection' in config:
            tagsDF = tagsDF.filter(config.mc_selection())
    else:
        if 'data_selection' in config:
            tagsDF = tagsDF.filter(config.data_selection())

            
    # Weight data and MC with the PileUp
    if doDataRew:
       weightedDF = get_weighted_data(tagsDF, era1, era, shift='Nominal')
    else:
       weightedDF = get_weighted_dataframe_condor(tagsDF, doGen, resonance, era, subEra, shift='Nominal')
        
    binning = config.binning()
    variables = config.variables()
    binVariables = config.binVariables()

    binningSet = set()

    for bvs in binVariables:
        binningSet = binningSet.union(set(bvs))

    binnedDF = weightedDF
        
    for bName in binningSet:
        binnedDF = get_binned_dataframe(
            binnedDF, bName+"Bin",
            variables[bName]['variable'],
            binning[bName])                
        
    
    yields = {}
    for binVars in binningSet:
        key = binVars
        yields[key] = binnedDF.groupBy(key+'Bin', *[key+'Bin']).agg({'weight': 'sum'})                                            


    realized = {}
    for binVars in yields:
        realized[binVars] = yields[binVars].toPandas()
        
    return realized


#
# Compare_one. Initialized spark, call run_files and draw and save the histograms.
#

def compare_one(particle, probe, resonance, era, config_name, muon_ID, _baseDir, _subera1, _subera2, _era2, lumi, fnames):
    

   # Start spark

    _useLocalSpark = False
    useParquet = True

    print("\n")
    print("\n")
    print("*************************************")
    print("******* Initializing Spark **********")
    print("*************************************")
    print("\n")
    print("\n")


    local_jars = ','.join([
        './laurelin-1.0.0.jar',
        './log4j-api-2.13.0.jar',
        './log4j-core-2.13.0.jar',
    ])
    
    spark = SparkSession\
        .builder\
        .appName("TnP")
    
    if useParquet == False:
        spark = spark\
        .config("spark.jars", local_jars)\
        .config("spark.driver.extraClassPath", local_jars)\
        .config("spark.executor.extraClassPath", local_jars)\
        .config("spark.dynamicAllocation.maxExecutors", "100")\
        .config("spark.driver.memory", "6g")\
        .config("spark.executor.memory", "4g")\
        .config("spark.executor.cores", "2")

    if _useLocalSpark == True:
        spark = spark.master("local")

    spark = spark.getOrCreate()

    sc = spark.sparkContext
    print(sc.getConf().toDebugString())


    config_real = config_name.split('/')[1]

    config = Configuration(config_real)

    if _era2 == 'era2':
        _era2 = era

    if _subera1 == 'subera1':
        _subera1 = ''
    if _subera2 == 'subera2':
        _subera2 = ''


    _FullEra = False

    subEra_dic = get_allowed_sub_eras(resonance, era)
    subEra2_dic = get_allowed_sub_eras(resonance, _era2)

    subEra_dic.remove(era)
    subEra2_dic.remove(_era2)

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
    
        

    efficiencies = config.efficiencies()

    realized = {}
    

    ### luminosity 
    lumi = float(lumi)


    if _FullEra:
        

        use_Data = False
        use_MC = False
        
        Z_peak = False
        JPsi_peak = False

        for subEra in subEras:
            if (era.split('_')[0] in subEra) or (era in subEra):
                use_Data = True
            elif 'DY_madgraph' in subEra:
                use_MC = True
                Z_peak = True
            elif 'DY_powheg' in subEra:
                use_MC = True
                Z_peak = True
            elif 'JPsi_pythia8' in subEra:
                use_MC = True
                JPsi_peak = True
       
            
        for subEra in subEras:
            key_name = ''
            for fname in fnames:
                if subEra in fname:
                    key_name = fname

                    realized[subEra] = run_files(particle, probe, resonance, era, subEra, config, spark, muon_ID, False, '', key_name) # Dictionary with pandas dataframes for each subEra
        
   
        binning = config.binning()
        variables = config.variables()
        binVariables = config.binVariables()
        
        ROOT.gROOT.SetBatch()
        ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")
        tdrstyle.setTDRStyle()
        
        
        # For each variable, draw and save a plot

        for binVar in binVariables:
            
            if len(binVar) == 1:
                
                # Data array
                fill_array = np.zeros(len(binning[binVar[0]])+1) 
                
                # MC array
                fill_array_mc = np.zeros(len(binning[binVar[0]])+1) 
                
                var_name = binVar[0]
                
                # Fill the arrays for each subEra
                for subEra in subEras:
                    
                    df = realized[subEra][var_name]
                    df = df.T.drop_duplicates().T
                    df = df.sort_values(by=[var_name+'Bin'])
                    df = df.reset_index()
                    
                    values = pd.Series(np.zeros(len(binning[var_name])+1))
                    values[df.index] = df['sum(weight)']
                    
                    if (subEra in ['DY_madgraph', 'DY_powheg', 'JPsi_pythia8']):
                        for i in df.index:
                            fill_array_mc[i] = fill_array_mc[i] + float(values[i])
                    else:
                        for i in df.index:
                            fill_array[i] = fill_array[i] + float(values[i])
                    
                           
                #
                # Initialize histogram
                #
                
                bins = np.array(binning[var_name])
                bins = bins.astype(np.float32)
                
                print("\n")
                print("Binning for " + var_name + " is ", bins)
                print("\n")
                
                hist = ROOT.TH1F(var_name, var_name, len(bins)-1, bins) #ROOT histogram
                hist_mc = ROOT.TH1F(var_name+'_mc', var_name+'_mc', len(bins)-1, bins) #ROOT histogram
                
                
                # Axis titles and options
                
                hist_mc.GetYaxis().SetTitle('Events')
                hist_mc.GetYaxis().SetTitleOffset(1.33)
                hist_mc.GetXaxis().SetTitle(variables[var_name]['pretty'])
                
                
                # Fill histogram
                
                for i in range(0, len(fill_array)):
                    hist.SetBinContent(i, fill_array[i])
                    hist_mc.SetBinContent(i, fill_array_mc[i])
                

                # Normalize histograms
                
                hist_mc.Scale(hist.Integral()/hist_mc.Integral())
                
                # Canvas initialization    
                
                cName = var_name        
                canvas = ROOT.TCanvas(cName, cName, 900, 800)
                ROOT.gStyle.SetPaintTextFormat("5.3f")
                #ROOT.gStyle.SetPaintTextFormat("4.1f")
                canvas.SetRightMargin(0.24)
            
            
                # Make up 
                
                hist.SetMarkerStyle(ROOT.kFullCircle)
                
                hist_mc.SetLineWidth(3)
                hist_mc.SetLineColor(ROOT.kAzure-4)
                hist_mc.SetFillColor(ROOT.kAzure-4)
                
                
                
                if use_MC:
                    hist_mc.Draw("hist")
                    
                if use_Data:
                    hist.Draw("e1 same")
                else:
                    hist.Draw("e1")


            
                #plotPath = os.path.join(plotDir, h)
                #canvas.SetLogy()
                #canvas.SetGrid()
                
                max1 = hist.GetMaximum()
                max2 = hist_mc.GetMaximum()
                
                if max1 > max2:
                    hist.SetMaximum(max1*0.2 + max1)
                    hist_mc.SetMaximum(max1*0.2 + max1)
                else:
                    hist.SetMaximum(max2*0.2 + max2)
                    hist_mc.SetMaximum(max2*0.2 + max2)
        
                
                legend = ROOT.TLegend(0.9, 0.9, 0.7, 0.78)
                #legend = ROOT.TLegend(0.5, 0.70, 0.92, 0.92)
                
                if use_Data:
                    legend.AddEntry(hist, "Data")
                if use_MC:
                    legend.AddEntry(hist_mc, "Simulation", "l")
                    
                legend.SetTextFont(42)
                legend.SetBorderSize(0)
                legend.SetFillColor(0)
                legend.Draw()
                
                canvas.Modified()
                canvas.Update()
                
                # CMS title and lumi
                
                CMS_lumi.cmsText = 'CMS'
                CMS_lumi.writeExtraText = True
                CMS_lumi.extraText = 'Preliminary'
                if lumi!=-1:
                    CMS_lumi.lumi_13TeV = "%0.1f fb^{-1}" % (lumi)
                else:
                    CMS_lumi.lumi_13TeV = ""
                CMS_lumi.CMS_lumi(canvas, 4, 11)
            

                #
                # Draw and save no ratio plots    
                #
                
                directory = _baseDir + "/plots/" + particle + "/" + probe + "/" + resonance + "/" + era + "/" + muon_ID + "/"
                
                if not os.path.exists(directory):
                    os.makedirs(directory)
                    
                canvas.Draw()        
                canvas.SaveAs(_baseDir + "/plots/" + particle + "/" + probe + "/" + resonance + "/" +
                              era + "/" + muon_ID + "/c_" + var_name + "_muon_val.png")  #Save .png file 
            
            
                print("\n")
                print(str(var_name) + " distribution saved at: " + _baseDir + "/plots/" + particle + "/" + probe + "/" + resonance + "/" +
                      era + "/" + muon_ID + "/c_" + var_name + "_muon_val.png")
                
                
                canvas.SetLogy()
                canvas.Update()
                canvas.Draw()
                canvas.SaveAs(_baseDir + "/plots/" + particle + "/" + probe + "/" + resonance + "/" +
                              era + "/" + muon_ID + "/log_c_" + var_name + "_muon_val.png")  #Save .png file 

                print(str(var_name) + " distribution saved at: " + _baseDir + "/plots/" + particle + "/" + probe + "/" + resonance + "/" +
                      era + "/" + muon_ID + "/log_c_" + var_name + "_muon_val.png")
                
                print("\n")
                
                #
                # Produce ratio plots
                #
                
                
                rcanvas = ROOT.TCanvas("r"+cName, "r"+cName, 900, 800)
                rcanvas.SetRightMargin(0.24)
                
                rcanvas.Divide(1,2)
                rcanvas.cd(1)
                
                plotPad = rcanvas.GetPad(1)
                plotPad.SetPad(0.,0.2,1.,1.)
                
                hist.GetXaxis().SetLabelSize(0.)
                hist.GetXaxis().SetTitleSize(0.)
                
                hist_mc.GetXaxis().SetLabelSize(0.)
                hist_mc.GetXaxis().SetTitleSize(0.)
            
                if use_MC:
                    hist_mc.Draw("HIST")
                    
                if use_Data:
                    hist.Draw("E1 SAME")
                else:
                    hist.Draw("E1")

                    
                #plotPath = os.path.join(plotDir, h)                                                                                                                                               
                #rcanvas.SetLogy()                                                                                                                                                                
                #rcanvas.SetGrid()                                                                                                                                                                
                    
                max1 = hist.GetMaximum()
                max2 = hist_mc.GetMaximum()
                
                if max1 > max2:
                    hist.SetMaximum(max1*0.2 + max1)
                    hist_mc.SetMaximum(max1*0.2 + max1)
                    hist_mc.GetYaxis().SetRangeUser(1, max1*0.2 + max1)
                    hist.GetYaxis().SetRangeUser(1, max1*0.2 + max1)
                else:
                    hist.SetMaximum(max2*0.2 + max2)
                    hist_mc.SetMaximum(max2*0.2 + max2)
                    hist_mc.GetYaxis().SetRangeUser(1, max2*0.2 + max2)
                    hist.GetYaxis().SetRangeUser(1, max2*0.2 + max2)


                legend = ROOT.TLegend(0.95, 0.9, 0.75, 0.78)
                #legend = ROOT.TLegend(0.5, 0.70, 0.92, 0.92)                                                                                                                                          
                
                if use_Data:
                    legend.AddEntry(hist, "Data")
                if use_MC:
                    legend.AddEntry(hist_mc, "Simulation", "l")
                    
                legend.SetTextFont(42)
                legend.SetBorderSize(0)
                legend.SetFillColor(0)
                legend.Draw()
                
                CMS_lumi.cmsText = 'CMS'
                CMS_lumi.writeExtraText = True
                CMS_lumi.extraText = 'Preliminary'
                if lumi!=-1:
                    CMS_lumi.lumi_13TeV = "%0.1f fb^{-1}" % (lumi)
                else:
                    CMS_lumi.lumi_13TeV = ""
                CMS_lumi.CMS_lumi(plotPad, 4, 11)


                #
                # Ratio pad
                #
                
                rcanvas.cd(2)
                
                ratioPad = rcanvas.GetPad(2)
                ratioPad.SetPad(0.,0.,1.,0.31)
                
                ratioPad.SetFillStyle(4000)
                ratioPad.SetBottomMargin(0.2)
                
                hRatio = hist.Clone()
                hRatio.SetTitle(" ")
                
                hRatio.GetXaxis().SetLabelSize(0.1)
                hRatio.GetXaxis().SetTitleSize(0.1)
                hRatio.GetXaxis().SetTitleOffset(.85)
                hRatio.GetXaxis().SetTitle(variables[var_name]['pretty'])
                
                hRatio.GetYaxis().SetLabelSize(0.07)
                hRatio.GetYaxis().SetTitleSize(0.1)
                hRatio.GetYaxis().SetTitleOffset(0.5)
                hRatio.GetYaxis().SetTitle("Data/MC")
                hRatio.GetYaxis().SetRangeUser(0.5,1.5)
                
                hRatio.Divide(hist_mc)
                hRatio.Draw()
                
                Xmax = hRatio.GetXaxis().GetXmax()
                Xmin = hRatio.GetXaxis().GetXmin()
                
                l = ROOT.TLine(Xmin, 1, Xmax, 1)
                l.SetLineColor(1) 
                l.Draw("same") 
                
                
                rcanvas.Draw()
                rcanvas.SaveAs(_baseDir + "/plots/" + particle + "/" + probe + "/" + resonance + "/" + era + "/" + muon_ID + "/c_ratio_" + var_name + "_muon_val.png")  #Save .png file           
                
                print("\n")
                print(str(var_name) + " distribution saved at: " + _baseDir + "/plots/" + particle + "/" + probe + "/" + resonance + "/" +
                      era + "/" + muon_ID + "/c_ratio" + var_name + "_muon_val.png")
                

                rcanvas.cd(1)
                
                plotPad.SetLogy()
                rcanvas.Update()
                rcanvas.Draw()
                
                rcanvas.SaveAs(_baseDir + "/plots/" + particle + "/" + probe + "/" + resonance + "/" + era + "/" + muon_ID +  "/log_c_ratio_" + var_name + "_muon_val.png")  #Save .png file           

                print(str(var_name) + " distribution saved at: " + _baseDir + "/plots/" + particle + "/" + probe + "/" + resonance + "/" +
                      era + "/" + muon_ID + "/log_c_ratio" + var_name + "_muon_val.png")

                print("\n")
                    
            if len(binVar) == 2:
                continue
                
    else:
       
        print(fnames)
        if _subera1 in fnames[0]:
           key1 = fnames[0]
           key2 = fnames[1]
        elif _subera1 in fnames[1]:
           key1 = fnames[1]
           key2 = fnames[0]
        else:
           print("Error: suberas not in path")
            
        subera1_isMC = _subera1 in ['DY_madgraph', 'DY_powheg', 'JPsi_pythia8']
        subera2_isMC = _subera2 in ['DY_madgraph', 'DY_powheg', 'JPsi_pythia8']
           
        realized[_subera1] = run_files(particle, probe, resonance, era, _subera1, config, spark, muon_ID, False, '', key1)
      
        if not subera1_isMC and not subera2_isMC:
            realized[_subera2] = run_files(particle, probe, resonance, _era2, _subera2, config, spark, muon_ID, True, era, key2) 
        else:
            realized[_subera2] = run_files(particle, probe, resonance, _era2, _subera2, config, spark, muon_ID, False, '', key2)

      
        binning = config.binning()
        variables = config.variables()
        binVariables = config.binVariables()
        
        ROOT.gROOT.SetBatch()
        ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")
        tdrstyle.setTDRStyle()
        
        
        for binVar in binVariables:
            
            if len(binVar) == 1:
                
                # subera1 array
                fill_array = np.zeros(len(binning[binVar[0]])+1) 
                
                # subera2 array
                fill_array_2 = np.zeros(len(binning[binVar[0]])+1) 
                
                var_name = binVar[0]
                
                # Fill the arrays for each subEra
                for subEra in subEras:
                    
                    df = realized[subEra][var_name]
                    df = df.T.drop_duplicates().T
                    df = df.sort_values(by=[var_name+'Bin'])
                    df = df.reset_index()
                    
                    values = pd.Series(np.zeros(len(binning[var_name])+1))
                    values[df.index] = df['sum(weight)']
                    
                    if (_subera1 == subEra):
                        for i in df.index:
                            fill_array[i] = fill_array[i] + float(values[i])
                    else:
                        for i in df.index:
                            fill_array_2[i] = fill_array_2[i] + float(values[i])
                            
                    
                #
                # Initialize histogram
                #
                
                bins = np.array(binning[var_name])
                bins = bins.astype(np.float32)
                
                hist = ROOT.TH1F(var_name, var_name, len(bins)-1, bins) #ROOT histogram
                hist_2 = ROOT.TH1F(var_name+'_2', var_name+'_2', len(bins)-1, bins) #ROOT histogram
                
                
                # Axis titles and options
                
                hist.GetYaxis().SetTitle('Events')
                hist.GetYaxis().SetTitleOffset(1.33)
                hist.GetXaxis().SetTitle(variables[var_name]['pretty'])
                
                hist_2.GetYaxis().SetTitle('Events')
                hist_2.GetYaxis().SetTitleOffset(1.33)
                hist_2.GetXaxis().SetTitle(variables[var_name]['pretty'])
                
                # Fill histogram
                
                for i in range(0, len(fill_array)):
                    hist.SetBinContent(i, fill_array[i])
                    hist_2.SetBinContent(i, fill_array_2[i])

                    
                # Canvas initialization
                
                
                cName = var_name
                canvas = ROOT.TCanvas(cName, cName, 900, 800)
                ROOT.gStyle.SetPaintTextFormat("5.3f")
                canvas.SetRightMargin(0.24)
                
                
                # Normalize histograms
                
                if subera1_isMC and not subera2_isMC:
                    hist_2.SetMarkerStyle(ROOT.kFullCircle)
                    
                    hist.SetLineWidth(3)
                    hist.SetLineColor(ROOT.kAzure-2)
                    hist.SetFillColor(ROOT.kAzure-2)
                    
                    his.Scale(hist_2.Integral()/hist.Integral())
                    hist.Draw("HIST")
                    hist_2.Draw("E1 SAME")
                elif subera2_isMC and not subera1_isMC:
                    hist.SetMarkerStyle(ROOT.kFullCircle)

                    hist_2.SetLineWidth(3)
                    hist_2.SetLineColor(ROOT.kAzure-4)
                    hist_2.SetFillColor(ROOT.kAzure-4)

                    hist_2.Scale(hist.Integral()/hist_2.Integral())
                    hist_2.Draw("HIST")
                    hist.Draw("E1 SAME")
                elif (not subera2_isMC) and (not subera1_isMC):
                    hist.SetMarkerStyle(ROOT.kFullCircle)

                    hist_2.SetMarkerStyle(ROOT.kFullCircle)
                    hist_2.SetMarkerColor(ROOT.kRed)

                    hist.Scale(1/hist.Integral())
                    hist_2.Scale(1/hist_2.Integral())
                    hist_2.Draw("E1")
                    hist.Draw("E1 SAME")

                else:
                    hist.SetMarkerStyle(ROOT.kFullCircle)

                    hist_2.SetMarkerStyle(ROOT.kFullCircle)
                    hist_2.SetMarkerColor(ROOT.kRed)

                    hist_2.Scale(hist.Integral()/hist_2.Integral())
                    hist_2.Draw("E1")
                    hist.Draw("E1 SAME")
                    
                    
                #plotPath = os.path.join(plotDir, h)
                #canvas.SetLogy()
                #canvas.SetGrid()
                
                max1 = hist.GetMaximum()
                max2 = hist_2.GetMaximum()
                
                if max1 > max2:
                    hist.SetMaximum(max1*0.2 + max1)
                    hist_2.SetMaximum(max1*0.2 + max1)
                else:
                    hist.SetMaximum(max2*0.2 + max2)
                    hist_2.SetMaximum(max2*0.2 + max2)
                    
                    
                legend = ROOT.TLegend(0.9, 0.9, 0.7, 0.78)
                
                
                if subera1_isMC and not subera2_isMC:
                    legend.AddEntry(hist, _subera1, "l")
                    legend.AddEntry(hist_2, _subera2)
                elif (not subera2_isMC) and (not subera1_isMC):
                    legend.AddEntry(hist, _subera1)
                    legend.AddEntry(hist_2, _subera2)
                else:
                    legend.AddEntry(hist, _subera1)
                    legend.AddEntry(hist_2, _subera2, "l")
                    
                legend.SetTextFont(42)
                legend.SetBorderSize(0)
                legend.SetFillColor(0)
                legend.Draw()
                
                canvas.Modified()
                canvas.Update()
                
                # CMS title and lumi
                
                CMS_lumi.cmsText = 'CMS'
                CMS_lumi.writeExtraText = True
                CMS_lumi.extraText = 'Preliminary'
                if lumi!=-1:
                    CMS_lumi.lumi_13TeV = "%0.1f fb^{-1}" % (lumi)
                else:
                    CMS_lumi.lumi_13TeV = ""
                CMS_lumi.CMS_lumi(canvas, 4, 11)
                
                # Draw    
                # Saved as file: ./baseDir/plots/muon/generalTracks/Z/Run2018A_vs_Run2018B/TightID/muon_pt_Run2018A_vs_Run2018B.png 
                
                directory = _baseDir + "/plots/" + particle + "/" + probe + "/" + resonance + "/" + _subera1 + "_vs_" + _subera2 + "/" + muon_ID + "/"
                
                if not os.path.exists(directory):
                    os.makedirs(directory)
                    
                canvas.Draw()        
                canvas.SaveAs(_baseDir + "/plots/" + particle + "/" + probe + "/" + resonance + "/" 
                              + _subera1 + "_vs_" + _subera2 + "/" + muon_ID + "/c_" + var_name + "_" + _subera1 + "_vs_" + _subera2 + "_muon_val.png")  #Save .png file 
            
                print("\n")
                print(str(var_name) + " distribution saved at: " + _baseDir + "/plots/" + particle + "/" + probe + "/" + resonance
                      + "/" + _subera1 + "_vs_" + _subera2 + "/" + muon_ID  + "/c_" + var_name + "_" + _subera1 + "_vs_" + _subera2 + "_muon_val.png")

                
                canvas.SetLogy()
                canvas.Update()
                canvas.Draw()
                
                canvas.SaveAs(_baseDir + "/plots/" + particle + "/" + probe + "/" + resonance + "/" +
                              _subera1 + "_vs_" + _subera2 + "/" + muon_ID  + "/log_c_" + var_name + "_" + _subera1 + "_vs_" + _subera2 + "_muon_val.png")
                
                print(str(var_name) + " distribution saved at: " + _baseDir + "/plots/" + particle + "/" + probe + "/" + resonance + "/" +
                      _subera1 + "_vs_" + _subera2 + "/" + muon_ID + "/log_c_" + var_name + "_" + _subera1 + "_vs_" + _subera2 + "_muon_val.png")
                
                print("\n")

                
                #                                                                                                                                                                                      
                # Produce ratio plots                                                                                                                                                                  
                #                                                                                                                                                                                     
                
                
                rcanvas = ROOT.TCanvas("r"+cName, "r"+cName, 900, 800)
                rcanvas.SetRightMargin(0.24)
                
                rcanvas.Divide(1,2)
                rcanvas.cd(1)
                
                plotPad = rcanvas.GetPad(1)
                plotPad.SetPad(0.,0.2,1.,1.)
                
                hist.GetXaxis().SetLabelSize(0.)
                hist.GetXaxis().SetTitleSize(0.)
                
                hist_2.GetXaxis().SetLabelSize(0.)
                hist_2.GetXaxis().SetTitleSize(0.)
            
                
                if subera1_isMC and not subera2_isMC:
                    hist_2.SetMarkerStyle(ROOT.kFullCircle)
                    
                    hist.SetLineWidth(3)
                    hist.SetLineColor(ROOT.kAzure-2)
                    hist.SetFillColor(ROOT.kAzure-2)
                    
                    his.Scale(hist_2.Integral()/hist.Integral())
                    hist.Draw("HIST")
                    hist_2.Draw("E1 SAME")

                elif subera2_isMC and not subera1_isMC:
                    hist.SetMarkerStyle(ROOT.kFullCircle)

                    hist_2.SetLineWidth(3)
                    hist_2.SetLineColor(ROOT.kAzure-4)
                    hist_2.SetFillColor(ROOT.kAzure-4)

                    hist_2.Scale(hist.Integral()/hist_2.Integral())
                    hist_2.Draw("HIST")
                    hist.Draw("E1 SAME")
                elif (not subera2_isMC) and (not subera1_isMC):
                    hist.SetMarkerStyle(ROOT.kFullCircle)

                    hist_2.SetMarkerStyle(ROOT.kFullCircle)
                    hist_2.SetMarkerColor(ROOT.kRed)

                    hist.Scale(1/hist.Integral())
                    hist_2.Scale(1/hist_2.Integral())
                    hist_2.Draw("E1")
                    hist.Draw("E1 SAME")
                else:
                    hist.SetMarkerStyle(ROOT.kFullCircle)

                    hist_2.SetMarkerStyle(ROOT.kFullCircle)
                    hist_2.SetMarkerColor(ROOT.kRed)

                    #hist_2.SetLineWidth(3)
                    #hist_2.SetLineColor(ROOT.kAzure-4)
                    #hist_2.SetFillColor(ROOT.kAzure-4)
                    
                    hist_2.Scale(hist.Integral()/hist_2.Integral())
                    hist_2.Draw("E1")
                    hist.Draw("E1 SAME")
                    
                #plotPath = os.path.join(plotDir, h)                                                                                                                                                  
                #rcanvas.SetLogy()                                                                                                                                                                     
                #rcanvas.SetGrid()                                                                                                                                                                     
                
                max1 = hist.GetMaximum()
                max2 = hist_2.GetMaximum()
                
                if max1 > max2:
                    hist.SetMaximum(max1*0.2 + max1)
                    hist_2.SetMaximum(max1*0.2 + max1)
                    hist_2.GetYaxis().SetRangeUser(1, max1*0.2 + max1)
                    hist.GetYaxis().SetRangeUser(1, max1*0.2 + max1)
                else:
                    hist.SetMaximum(max2*0.2 + max2)
                    hist_2.SetMaximum(max2*0.2 + max2)
                    hist_2.GetYaxis().SetRangeUser(1, max2*0.2 + max2)
                    hist.GetYaxis().SetRangeUser(1, max2*0.2 + max2)
                    
                    
                legend.Draw()
                
                
                CMS_lumi.cmsText = 'CMS'
                CMS_lumi.writeExtraText = True
                CMS_lumi.extraText = 'Preliminary'
                if lumi!=-1:
                    CMS_lumi.lumi_13TeV = "%0.1f fb^{-1}" % (lumi)
                else:
                    CMS_lumi.lumi_13TeV = ""
                CMS_lumi.CMS_lumi(plotPad, 4, 11)


                #                                                                                                                                                                                     
                # Ratio pad                                                                                                                                                                            
                #                                                                                                                                                                                  


                rcanvas.cd(2)
                
                ratioPad = rcanvas.GetPad(2)
                ratioPad.SetPad(0.,0.,1.,0.31)
                
                ratioPad.SetFillStyle(4000)
                ratioPad.SetBottomMargin(0.2)
                
                if subera1_isMC and not subera2_isMC:
                    hRatio = hist_2.Clone()
                    hRatio.Divide(hist)
                else:
                    hRatio = hist.Clone()
                    hRatio.Divide(hist_2)


                hRatio.SetTitle(" ")
                
                hRatio.GetXaxis().SetLabelSize(0.1)
                hRatio.GetXaxis().SetTitleSize(0.1)
                hRatio.GetXaxis().SetTitleOffset(.85)
                hRatio.GetXaxis().SetTitle(variables[var_name]['pretty'])

                hRatio.GetYaxis().SetTitle('Events')
                hRatio.GetYaxis().SetTitleOffset(1.33)
                
                hRatio.GetYaxis().SetLabelSize(0.07)
                hRatio.GetYaxis().SetTitleSize(0.1)
                hRatio.GetYaxis().SetTitleOffset(0.5)
            
                if subera1_isMC or subera2_isMC:
                    hRatio.GetYaxis().SetTitle("Data/MC")
                else:
                    hRatio.GetYaxis().SetTitle("Data 1/Data 2")
                hRatio.GetYaxis().SetRangeUser(0.5,1.5)
                    
                hRatio.Draw()
                
                Xmax = hRatio.GetXaxis().GetXmax()
                Xmin = hRatio.GetXaxis().GetXmin()
                
                l = ROOT.TLine(Xmin, 1, Xmax, 1)
                l.SetLineColor(1)
                l.Draw("same")
            
                
                rcanvas.Draw()
                rcanvas.SaveAs(_baseDir + "/plots/" + particle + "/" + probe + "/" + resonance + "/" + _subera1 + "_vs_" + _subera2 + "/" + muon_ID  
                               + "/c_ratio_" + var_name + "_" + _subera1 + "_vs_" + _subera2 + "_muon_val.png")  #Save .png file
                
                print("\n")
                print(str(var_name) + " distribution saved at: " + _baseDir + "/plots/" + particle + "/" + probe + "/" + resonance + "/" +
                      _subera1 + "_vs_" + _subera2 + "/" + muon_ID  + "/c_ratio_" + var_name + "_" + _subera1 + "_vs_" + _subera2 + "_muon_val.png")
                
                
                rcanvas.cd(1)
                
                plotPad.SetLogy()
                rcanvas.Update()
                rcanvas.Draw()
                
                rcanvas.SaveAs(_baseDir + "/plots/" + particle + "/" + probe + "/" + resonance + "/" + _subera1 + "_vs_" + _subera2 + "/" +  muon_ID  + "/log_c_ratio_" 
                               + var_name + "_" + _subera1 + "_vs_" + _subera2 + "_muon_val.png")
                
                print(str(var_name) + " distribution saved at: " + _baseDir + "/plots/" + particle + "/" + probe + "/" + resonance + "/" +
                      _subera1 + "_vs_" + _subera2 + "/" + muon_ID  + "/log_c_ratio_" + var_name + "_" + _subera1 + "_vs_" + _subera2 + "_muon_val.png")
                
                print("\n")

                
            if len(binVar) == 2:
                pass
                

                
    spark.stop()
    
    
if __name__ == "__main__":
    argv = sys.argv[1:]
    fnames = argv[11:]
    print(argv)
    compare_one(argv[0], argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9], argv[10], fnames)
    
