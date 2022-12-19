from __future__ import print_function
import os
import importlib.util
import sys
import math
import itertools
import json
from array import array
import ctypes
import ROOT
import tdrstyle
import CMS_lumi

from dataset_allowed_definitions import get_data_mc_sub_eras
from muon_definitions import (get_full_name, get_eff_name,
                              get_bin_name,
                              get_extended_eff_name,
                              get_variables_name)
from registry import registry

#ROOT.gROOT.LoadMacro('make_ratioplots.C+')
ROOT.gROOT.SetBatch()
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")
tdrstyle.setTDRStyle()

def generateClopperPearsonInterval(num,den):
    confidenceLevel = 0.68
    alpha = 1 - confidenceLevel
    
    lowerLimit = round(ROOT.Math.beta_quantile(alpha/2,num,den-num + 1),4)
    if num==den:
        upperLimit=1
    else:
        upperLimit = round(ROOT.Math.beta_quantile(1-alpha/2,num + 1,den-num),4)
    return lowerLimit,upperLimit


def computeEff(n1, n2, e1, e2):
    eff = n1 / (n1 + n2)
    err = 1 / (n1 + n2) * math.sqrt(e1 * e1 * n2 * n2 + e2 * e2 * n1 * n1) / (n1 + n2)
    effD, effU = generateClopperPearsonInterval(n1,n1+n2)
    effD = err
    effU = err
    return eff, effD, effU



def getEff(binName, fname, massRanges, shift=None, cutAndCount=False, resonance='Z'):
    try:
        # MC Eff is always cut and count now
        tfile = ROOT.TFile(fname, 'read')
        if cutAndCount:
            hP = tfile.Get('{}_Pass'.format(binName))
            hF = tfile.Get('{}_Fail'.format(binName))
        else:
            hP = tfile.Get('{}_GenPass'.format(binName))
            hF = tfile.Get('{}_GenFail'.format(binName))
        # hard code Z for now (same as in run_single_fit.py)
        # AF: TODO make the Z vs JPsi treatment less hacky perhaps
        if resonance == 'JPsi':
            if shift == 'massRangeUp':
                blow, bhigh = 2.96, 3.36
            elif shift == 'massRangeDown':
                blow, bhigh = 2.84, 3.24
            else:
                blow, bhigh = 2.90, 3.30
        else:
            if shift == 'massRangeUp':
                blow, bhigh = 75, 135
            elif shift == 'massRangeDown':
                blow, bhigh = 65, 125
            else:
                blow, bhigh = 70, 115
        #blow, bhigh = massRanges[resonance].get(shift, massRanges[resonance]["nominal"])

        bin1 = hP.GetXaxis().FindBin(blow)
        bin2 = hP.GetXaxis().FindBin(bhigh)
        eP = ctypes.c_double(-1.0)
        eF = ctypes.c_double(-1.0)
        nP = hP.IntegralAndError(bin1, bin2, eP)
        nF = hF.IntegralAndError(bin1, bin2, eF)
        eff, effD, effU = computeEff(nP, nF, eP, eF)
        #errD = abs(eff-effD)
        #errU = abs(eff-effU)
        errD = effD
        errU = effU
        tfile.Close()
        return eff, errD, errU
    except Exception as e:
        print('Exception for getEff', binName)
        print(e)
        # raise e
        return 1., 0., 0.


def getMCEff(binName, fname, massRanges, shift=None, cutAndCount=False, resonance='Z'):
    # hard code Z for now (same as in run_single_fit.py)
        # AF: TODO make the Z vs JPsi treatment less hacky perhaps
    if resonance == 'JPsi':
        if shift == 'massRangeUp':
            blow, bhigh = 2.96, 3.36
        elif shift == 'massRangeDown':
            blow, bhigh = 2.84, 3.24
        else:
            blow, bhigh = 2.90, 3.30
    else:
        if shift == 'massRangeUp':
            blow, bhigh = 75, 135
        elif shift == 'massRangeDown':
            blow, bhigh = 65, 125
        else:
            blow, bhigh = 70, 115
    #blow, bhigh = massRanges[resonance].get(shift, massRanges[resonance]["nominal"])
    try:
        tfile = ROOT.TFile(fname, 'read')
        if cutAndCount:
            hP = tfile.Get('{}_GenPass'.format(binName))
            hF = tfile.Get('{}_GenFail'.format(binName))

            bin1 = hP.GetXaxis().FindBin(blow)
            bin2 = hP.GetXaxis().FindBin(bhigh)
            eP = ctypes.c_double(-1.0)
            eF = ctypes.c_double(-1.0)
            nP = hP.IntegralAndError(bin1, bin2, eP)
            nF = hF.IntegralAndError(bin1, bin2, eF)
            eff, effD, effU = computeEff(nP, nF, eP, eF)
        else:
            fitresP = tfile.Get('{}_resP'.format(binName))
            fitresF = tfile.Get('{}_resF'.format(binName))

            fitP = fitresP.floatParsFinal().find('nSigP')
            fitF = fitresF.floatParsFinal().find('nSigF')

            nP = fitP.getVal()
            nF = fitF.getVal()
            eP = fitP.getError()
            eF = fitF.getError()

            hP = tfile.Get('{}_Pass'.format(binName))
            hF = tfile.Get('{}_Fail'.format(binName))

            bin1 = hP.GetXaxis().FindBin(blow)
            bin2 = hP.GetXaxis().FindBin(bhigh)
            ePalt = ctypes.c_double(-1.0)
            eFalt = ctypes.c_double(-1.0)
            hP.IntegralAndError(bin1, bin2, ePalt)
            hF.IntegralAndError(bin1, bin2, eFalt)

            #eP = min(eP, ePalt.value)
            #eF = min(eF, eFalt.value)
            eff, effD, effU = computeEff(nP, nF, eP, eF)
        #errD = abs(eff-effD)
        #errU = abs(eff-effU)
        errD = effD
        errU = effU
        tfile.Close()
        return eff, errD, errU

    except Exception as e:
        print('Exception for getMCEff', binName)
        print(e)
        # raise e
        return 1., 1., 1.



def getDataEff(binName, fname, massRanges, shift=None, cutAndCount=False, resonance='Z'):
    # hard code Z for now (same as in run_single_fit.py)
        # AF: TODO make the Z vs JPsi treatment less hacky perhaps
    if resonance == 'JPsi':
        if shift == 'massRangeUp':
            blow, bhigh = 2.96, 3.36
        elif shift == 'massRangeDown':
            blow, bhigh = 2.84, 3.24
        else:
            blow, bhigh = 2.90, 3.30
    else:
        if shift == 'massRangeUp':
            blow, bhigh = 75, 135
        elif shift == 'massRangeDown':
            blow, bhigh = 65, 125
        else:
            blow, bhigh = 70, 115
    #blow, bhigh = massRanges[resonance].get(shift, massRanges[resonance]["nominal"])
    try:
        tfile = ROOT.TFile(fname, 'read')
        if cutAndCount:
            hP = tfile.Get('{}_Pass'.format(binName))
            hF = tfile.Get('{}_Fail'.format(binName))

            bin1 = hP.GetXaxis().FindBin(blow)
            bin2 = hP.GetXaxis().FindBin(bhigh)
            eP = ctypes.c_double(-1.0)
            eF = ctypes.c_double(-1.0)
            nP = hP.IntegralAndError(bin1, bin2, eP)
            nF = hF.IntegralAndError(bin1, bin2, eF)
            eff, effD, effU = computeEff(nP, nF, eP, eF)
        else:
            fitresP = tfile.Get('{}_resP'.format(binName))
            fitresF = tfile.Get('{}_resF'.format(binName))

            fitP = fitresP.floatParsFinal().find('nSigP')
            fitF = fitresF.floatParsFinal().find('nSigF')

            nP = fitP.getVal()
            nF = fitF.getVal()
            eP = fitP.getError()
            eF = fitF.getError()

            hP = tfile.Get('{}_Pass'.format(binName))
            hF = tfile.Get('{}_Fail'.format(binName))

            bin1 = hP.GetXaxis().FindBin(blow)
            bin2 = hP.GetXaxis().FindBin(bhigh)
            ePalt = ctypes.c_double(-1.0)
            eFalt = ctypes.c_double(-1.0)
            hP.IntegralAndError(bin1, bin2, ePalt)
            hF.IntegralAndError(bin1, bin2, eFalt)

            #eP = min(eP, ePalt.value)
            #eF = min(eF, eFalt.value)
            eff, effD, effU = computeEff(nP, nF, eP, eF)
        #errD = abs(eff-effD)
        #errU = abs(eff-effU)
        errD = effD
        errU = effU
        tfile.Close()
        return eff, errD, errU

    except Exception as e:
        print('Exception for getDataEff', binName)
        print(e)
        # raise e
        return 1., 1., 1.


def getSF(binName, fname, massRanges, shift=None, resonance='Z'):
    mcEff, mcErrD, mcErrU = getMCEff(binName, fname, massRanges, shift, False, resonance)
    dataEff, dataErrD, dataErrU = getDataEff(binName, fname, massRanges, shift, False, resonance)
    sf = dataEff / mcEff if mcEff else 0.0
    sf_err = 0.0
    if dataEff and mcEff:
        dataErr = max(dataErrD, dataErrU)
        mcErr = max(mcErrD, mcErrU)
        sf_err = sf * ((dataErr / dataEff)**2 + (mcErr / mcEff)**2)**0.5
    return sf, sf_err, dataEff, dataErrD, dataErrU, mcEff, mcErrD, mcErrU


def getSF_fitdataMC(binName, fname_data, fname_mc,massRanges, shift=None, resonance='Z'):
    mcEff, mcErrD, mcErrU = getMCEff(binName, fname_mc, massRanges, shift, False, resonance)
    dataEff, dataErrD, dataErrU = getDataEff(binName, fname_data, massRanges, shift, False, resonance)
    if not 'Matched_fakerate_' in fname_data:
        fname_mc_fake = '%s' % fname_mc
        fname_data_fake = '%s' % fname_data
        binName_fake = '%s' % binName
        fname_mc_fake = fname_mc_fake.replace('Matched_', 'Matched_fakerate_')
        fname_data_fake = fname_data_fake.replace('Matched_', 'Matched_fakerate_')
        binName_fake = binName_fake.replace('Matched_', 'Matched_fakerate_')
        mcEff_fake, mcErr_fakeD, mcErr_fakeU  = getMCEff(binName_fake, fname_mc_fake, massRanges, shift, False, resonance)
        dataEff_fake, dataErr_fakeD, dataErr_fakeU = getDataEff(binName_fake, fname_data_fake, massRanges, shift, False, resonance)
        if mcEff_fake != 1 and dataEff_fake != 1 :
            mcEff = (mcEff - mcEff_fake)/(1 - mcEff_fake) #(   y   -    yf   )/(1-    yf   )
            mcErrU =   ((mcEff+mcErrU)-(mcEff_fake-mcErr_fakeD))/(1-(mcEff_fake-mcErr_fakeD)) - mcEff
            mcErrD =   mcEff - ((mcEff-mcErrD)-(mcEff_fake+mcErr_fakeU))/(1-(mcEff_fake+mcErr_fakeU)) #double ycorr_lo = ((y-eyl)-(yf+eyhf))/(1-(yf+eyhf));
            dataEff = (dataEff - dataEff_fake)/(1 - dataEff_fake) 
            dataErrU =   ((dataEff+dataErrU)-(dataEff_fake-dataErr_fakeD))/(1-(dataEff_fake-dataErr_fakeD)) - dataEff
            dataErrD =   dataEff - ((dataEff-dataErrD)-(dataEff_fake+dataErr_fakeU))/(1-(dataEff_fake+dataErr_fakeU))
    sf = dataEff / mcEff if mcEff else 0.0
    sf_err = 0.0
    if dataEff and mcEff:
        dataErr = max(dataErrD, dataErrU)
        mcErr = max(mcErrD, mcErrU)
        sf_err = sf * ((dataErr / dataEff)**2 + (mcErr / mcEff)**2)**0.5
    return sf, sf_err, dataEff, dataErrD, dataErrU, mcEff, mcErrD, mcErrU


def getSF_cutAndCount(binName, fnameData, fnameMC, massRanges, shift=None, resonance='Z'):
    mcEff, mcErrD, mcErrU = getEff(binName, fnameMC, massRanges, shift, True, resonance)
    dataEff, dataErrD, dataErrU = getDataEff(binName, fnameData, massRanges, shift, True, resonance)
    sf = dataEff / mcEff if mcEff else 0.0
    sf_err = 0.0
    if dataEff and mcEff:
        dataErr = max(dataErrD, dataErrU)
        mcErr = max(mcErrD, mcErrU)
        sf_err = sf * ((dataErr / dataEff)**2 + (mcErr / mcEff)**2)**0.5
    return sf, sf_err, dataEff, dataErrD, dataErrU, mcEff, mcErrD, mcErrU


#def getSyst(binName, fname, fitTypes, shiftTypes, massRanges, resonance='Z'):
#    sf, sf_err, dataEff, dataErr, mcEff, mcErr = getSF(binName, fname, massRanges, resonance=resonance)
def getSyst(binName, fname_data, fname_mc,fitTypes, shiftTypes, massRanges, resonance='Z'):
    sf, sf_err, dataEff, dataErrD, dataErrU, mcEff, mcErrD, mcErrU = getSF_fitdataMC(binName, fname_data, fname_mc, massRanges, resonance=resonance)
    syst = {}
    for isyst in fitTypes:
        systfname = fname.replace('Nominal', isyst)
        # sf, sf_err, dataEff, dataErrD, dataErrU, mcEff, mcErrD, mcErrU
        tmp = getSF(binName, systfname, massRanges, isyst, resonance=resonance)
        syst[isyst] = {
            'sf': tmp[0],
            'err': abs(tmp[0]-sf),
            'dataEff': tmp[2],
            'dataErr': abs(tmp[2]-dataEff),
            'mcEff': tmp[5],
            'mcErr': abs(tmp[5]-mcEff),
        }

    for isyst in shiftTypes:
        systUpfname = fname.replace('Nominal', isyst+'Up')
        systDnfname = fname.replace('Nominal', isyst+'Down')
        # sf, sf_err, dataEff, dataErrD, dataErrU, mcEff, mcErrD, mcErrU
        tmpUp = getSF(binName, systUpfname, massRanges, isyst+'Up'  , resonance=resonance)
        tmpDn = getSF(binName, systDnfname, massRanges, isyst+'Down', resonance=resonance)
        tmp = [
            (tmpUp[0]+tmpDn[0])/2,
            (abs(tmpUp[0]-sf)+abs(tmpDn[0]-sf))/2,
            (tmpUp[2]+tmpDn[2])/2,
            (abs(tmpUp[2]-dataEff)+abs(tmpDn[2]-dataEff))/2,
            (tmpUp[5]+tmpDn[5])/2,
            (abs(tmpUp[5]-mcEff)+abs(tmpDn[5]-mcEff))/2,
        ]
        syst[isyst] = {
            'sf': tmp[0],
            'err': tmp[1],
            'dataEff': tmp[2],
            'dataErr': tmp[3],
            'mcEff': tmp[4],
            'mcErr': tmp[5],
        }
        syst[isyst+'Up'] = {
            'sf': tmpUp[0],
            'err': abs(tmpUp[0]-sf),
            'dataEff': tmpUp[2],
            'dataErr': abs(tmpUp[2]-dataEff),
            'mcEff': tmpUp[5],
            'mcErr': abs(tmpUp[5]-mcEff),
        }
        syst[isyst+'Down'] = {
            'sf': tmpDn[0],
            'err': abs(tmpDn[0]-sf),
            'dataEff': tmpDn[2],
            'dataErr': abs(tmpDn[2]-dataEff),
            'mcEff': tmpDn[5],
            'mcErr': abs(tmpDn[5]-mcEff),
        }

    return syst


def getSyst_cutAndCount(binName, fnameData, fnameMC, fitTypes, shiftTypes, massRanges, resonance='Z'):
    sf, sf_err, dataEff, dataErrD, dataErrU, mcEff, mcErrD, mcErrU = getSF_cutAndCount(
        binName, fnameData, fnameMC, massRanges, resonance=resonance)

    syst = {}
    for isyst in fitTypes:
        systfnameData = fnameData.replace('Nominal', isyst)
        systfnameMC = fnameMC.replace('Nominal', isyst)
        # sf, sf_err, dataEff, dataErrD, dataErrU, mcEff, mcErrD, mcErrU
        tmp = getSF_cutAndCount(binName, systfnameData, systfnameMC, massRanges, isyst, resonance=resonance)
        syst[isyst] = {
            'sf': tmp[0],
            'err': abs(tmp[0]-sf),
            'dataEff': tmp[2],
            'dataErr': abs(tmp[2]-dataEff),
            'mcEff': tmp[5],
            'mcErr': abs(tmp[5]-mcEff),
        }

    for isyst in shiftTypes:
        systUpfnameData = fnameData.replace('Nominal', isyst+'Up')
        systDnfnameData = fnameData.replace('Nominal', isyst+'Down')
        systUpfnameMC = fnameMC.replace('Nominal', isyst+'Up')
        systDnfnameMC = fnameMC.replace('Nominal', isyst+'Down')
         # sf, sf_err, dataEff, dataErrD, dataErrU, mcEff, mcErrD, mcErrU
        tmpUp = getSF_cutAndCount(binName, systUpfnameData,
                                  systUpfnameMC, massRanges, isyst+'Up'  , resonance=resonance)
        tmpDn = getSF_cutAndCount(binName, systDnfnameData,
                                  systDnfnameMC, massRanges, isyst+'Down', resonance=resonance)
        tmp = [
            (tmpUp[0]+tmpDn[0])/2,
            (abs(tmpUp[0]-sf)+abs(tmpDn[0]-sf))/2,
            (tmpUp[2]+tmpDn[2])/2,
            (abs(tmpUp[2]-dataEff)+abs(tmpDn[2]-dataEff))/2,
            (tmpUp[5]+tmpDn[5])/2,
            (abs(tmpUp[5]-mcEff)+abs(tmpDn[5]-mcEff))/2,
        ]
        syst[isyst] = {
            'sf': tmp[0],
            'err': tmp[1],
            'dataEff': tmp[2],
            'dataErr': tmp[3],
            'mcEff': tmp[4],
            'mcErr': tmp[5],
        }
        syst[isyst+'Up'] = {
            'sf': tmpUp[0],
            'err': abs(tmpUp[0]-sf),
            'dataEff': tmpUp[2],
            'dataErr': abs(tmpUp[2]-dataEff),
            'mcEff': tmpUp[5],
            'mcErr': abs(tmpUp[5]-mcEff),
        }
        syst[isyst+'Down'] = {
            'sf': tmpDn[0],
            'err': abs(tmpDn[0]-sf),
            'dataEff': tmpDn[2],
            'dataErr': abs(tmpDn[2]-dataEff),
            'mcEff': tmpDn[5],
            'mcErr': abs(tmpDn[5]-mcEff),
        }

    return syst


def prepare_corrected(baseDir, particle, probe, resonance, era,
            config, num, denom, variableLabels, lumi,
            skipPlots=False, cutAndCount=False):
    
    hists = {}

    effType = config.type() if 'type' in config else ''
    effName = get_eff_name(num, denom)
    extEffName = get_extended_eff_name(num, denom, variableLabels)
    binning = config.binning()
    dataSubEra, mcSubEra = get_data_mc_sub_eras(resonance, era)

    systList = config.get('systematics',
                          {x: {'fitTypes': [],
                               'shiftTypes': []}
                           for x in ['SF', 'dataEff', 'mcEff']})

    massRanges = config.get("massRanges") # Start from a sensible default, but allow customization

    def get_variable_name_pretty(variableLabel):
        variables = config.variables()
        return variables.get(variableLabel, {}).get('pretty', variableLabel)

    # create output histograms
    nVars = len(variableLabels)
    if nVars == 1:
        THX = ROOT.TH1F
    elif nVars == 2:
        THX = ROOT.TH2F
    elif nVars == 3:
        THX = ROOT.TH3F
    else:
        raise NotImplementedError(
            'More than 3 dimensions are not supported for scale factors'
        )

    hargs = [extEffName, extEffName]
    for variableLabel in variableLabels:
        hargs += [len(binning[variableLabel]) - 1,
                  array('d', binning[variableLabel])]
    hist = THX(*hargs)
    axes = [hist.GetXaxis(), hist.GetYaxis(), hist.GetZaxis()]
    for vi, variableLabel in enumerate(variableLabels):
        axes[vi].SetTitle(get_variable_name_pretty(variableLabel))
    if nVars == 1:
        hist.GetYaxis().SetTitle('Scalefactor')
    if nVars == 2:
        hist.SetOption('colz')
        hist.GetZaxis().SetTitle('Scalefactor')
    hist_stat = hist.Clone(extEffName+'_stat')
    hist_syst = hist.Clone(extEffName+'_syst')
    histList_syst = {
        'combined_syst': hist.Clone(extEffName+'_combined_syst'),
    }
    if nVars == 2:
        histList_syst['combined_syst'].GetZaxis().SetTitle('Uncertainty')

    hist_dataEff = hist.Clone(extEffName+'_efficiencyData')
    if nVars == 1:
        hist_dataEff.GetYaxis().SetTitle('Efficiency')
    if nVars == 2:
        hist_dataEff.GetZaxis().SetTitle('Efficiency')
    hist_dataEff_errD = hist_dataEff.Clone(extEffName+'_efficiencyData_errD')
    hist_dataEff_errU = hist_dataEff.Clone(extEffName+'_efficiencyData_errU')
    hist_dataEff_stat = hist_dataEff.Clone(extEffName+'_efficiencyData_stat')
    hist_dataEff_syst = hist_dataEff.Clone(extEffName+'_efficiencyData_syst')
    histList_dataEff_syst = {
        'combined_syst': hist_dataEff.Clone(
            extEffName+'_efficiencyData_combined_syst'),
    }
    if nVars == 2:
        histList_dataEff_syst['combined_syst'].GetZaxis().SetTitle('Uncertainty')
    hist_mcEff = hist_dataEff.Clone(extEffName+'_efficiencyMC')
    hist_mcEff_errD = hist_dataEff.Clone(extEffName+'_efficiencyMC_errD')
    hist_mcEff_errU = hist_dataEff.Clone(extEffName+'_efficiencyMC_errU')
    hist_mcEff_stat = hist_dataEff.Clone(extEffName+'_efficiencyMC_stat')
    hist_mcEff_syst = hist_dataEff.Clone(extEffName+'_efficiencyMC_syst')
    histList_mcEff_syst = {
        'combined_syst': hist_dataEff.Clone(
            extEffName+'_efficiencyMC_combined_syst'),
    }
    if nVars == 2:
        histList_mcEff_syst['combined_syst'].GetZaxis().SetTitle('Uncertainty')

    # the individual systematics
    for iSyst in itertools.chain(systList['SF']['fitTypes'],
                                 systList['SF']['shiftTypes']):
        histList_syst[iSyst] = hist.Clone(extEffName+'_'+iSyst)
        histList_syst[iSyst+'_syst'] = hist.Clone(extEffName+'_'+iSyst+'_syst')
        if nVars == 2:
            histList_syst[iSyst+'_syst'].GetZaxis().SetTitle('Uncertainty')
    for iSyst in itertools.chain(systList['dataEff']['fitTypes'],
                                 systList['dataEff']['shiftTypes']):
        histList_dataEff_syst[iSyst] = hist_dataEff.Clone(extEffName+'_'+iSyst)
        histList_dataEff_syst[iSyst+'_syst'] = hist_dataEff.Clone(
            extEffName+'_'+iSyst+'_syst')
        if nVars == 2:
            histList_dataEff_syst[iSyst+'_syst'].GetZaxis().SetTitle('Uncertainty')
    for iSyst in itertools.chain(systList['mcEff']['fitTypes'],
                                 systList['mcEff']['shiftTypes']):
        histList_mcEff_syst[iSyst] = hist_mcEff.Clone(extEffName+'_'+iSyst)
        histList_mcEff_syst[iSyst+'_syst'] = hist_mcEff.Clone(
            extEffName+'_'+iSyst+'_syst')
        if nVars == 2:
            histList_mcEff_syst[iSyst+'_syst'].GetZaxis().SetTitle('Uncertainty')

    varName = get_variables_name(variableLabels)

    # iterate through the bin indices
    # this does nested for loops of the N-D binning (e.g. pt, eta)
    # binning starts at 1 (0 is underflow), same as ROOT
    indices = [list(range(1, len(binning[variableLabel])))
               for variableLabel in variableLabels]
    output = {effName: {varName: {}}}
    all_systematics = {}
    for index in itertools.product(*indices):
        binName = get_full_name(num, denom, variableLabels, index)
        subVarKeys = [
            '{}:[{},{}]'.format(
                variableLabels[i],
                binning[variableLabels[i]][ind-1],
                binning[variableLabels[i]][ind]
            ) for i, ind in enumerate(index)
        ]
        _out = output[effName][varName]

        # add binning definitions
        _out['binning'] = [
            {
                'variable': vl,
                'binning': binning[vl].tolist(),
            }
            for vl in variableLabels
        ]

        for subVarKey in subVarKeys:
            if subVarKey not in _out:
                _out[subVarKey] = {}
            _out = _out[subVarKey]

        # the fitted distributions
        fitType = 'Nominal'
        dataFNameFit = os.path.join(baseDir, 'fits_data',
                                    particle, probe,
                                    resonance, era,
                                    fitType, effName,
                                    binName + '.root')
        mcFNameFit = os.path.join(baseDir, 'fits_mc',
                                    particle, probe,
                                    resonance, era,
                                    fitType, effName,
                                    binName + '.root')
        dataFNameCNC = os.path.join(baseDir, 'flat',
                                    particle, probe,
                                    resonance, era,
                                    dataSubEra, 'Nominal',
                                    extEffName + '.root')
        mcFNameCNC = os.path.join(baseDir, 'flat',
                                  particle, probe,
                                  resonance, era,
                                  mcSubEra, 'Nominal',
                                  extEffName + '.root')
        if cutAndCount:
            sf, sf_stat, dataEff, dataStatD, dataStatU, mcEff, mcStatD, mcStatU = getSF_cutAndCount(
                binName, dataFNameCNC, mcFNameCNC, massRanges, resonance=resonance)
        else:
            #sf, sf_stat, dataEff, dataStatD, dataStatU, mcEff, mcStatD, mcStatU = getSF(
            #    binName, dataFNameFit, massRanges, resonance=resonance)
            sf, sf_stat, dataEff, dataStatD, dataStatU, mcEff, mcStatD, mcStatU = getSF_fitdataMC(
                binName, dataFNameFit, mcFNameFit, massRanges, shift=None,resonance=resonance)
        #As statistical error, we take the maximum between up and down asym. errors
        dataStat = max(dataStatD, dataStatU)
        mcStat = max(mcStatD, mcStatU)
        fitTypes = set(systList['SF']['fitTypes']
                       + systList['dataEff']['fitTypes']
                       + systList['mcEff']['fitTypes'])
        shiftTypes = set(systList['SF']['shiftTypes']
                         + systList['dataEff']['shiftTypes']
                         + systList['mcEff']['shiftTypes'])
        if cutAndCount:
            sf_syst = getSyst_cutAndCount(binName, dataFNameCNC, mcFNameCNC,
                                          fitTypes, shiftTypes, massRanges, resonance=resonance)
        else:
            #sf_syst = getSyst(binName, dataFNameFit,
            #                  fitTypes, shiftTypes, massRanges, resonance=resonance)
            sf_syst = getSyst(binName, dataFNameFit, mcFNameFit,
                                fitTypes,shiftTypes,massRanges,resonance=resonance)
        combined_syst = {}
        for kind in ['SF', 'dataEff', 'mcEff']:
            combined_syst[kind] = 0
            errKey = 'err'
            if kind == 'dataEff':
                errKey = 'dataErr'
            if kind == 'mcEff':
                errKey = 'mcErr'
            for t in itertools.chain(systList[kind]['fitTypes'],
                                     systList[kind]['shiftTypes']):
                combined_syst[kind] += sf_syst[t][errKey]**2
            combined_syst[kind] = combined_syst[kind]**0.5

        sf_err = (sf_stat**2 + combined_syst['SF']**2)**0.5
        dataErr = (dataStat**2 + combined_syst['dataEff']**2)**0.5
        dataErrD = dataStatD
        if (dataEff + dataStatU) >= 1.0:
            dataErrU = 1.0 - dataEff
        else:
            dataErrU = dataStatU
        #dataErrD = (dataStatD**2 + combined_syst['dataEff']**2)**0.5
        #dataErrU = (dataStatU**2 + combined_syst['dataEff']**2)**0.5
        mcErr = (mcStat**2 + combined_syst['mcEff']**2)**0.5
        #mcErrD = (mcStatD**2 + combined_syst['mcEff']**2)**0.5
        #mcErrU = (mcStatU**2 + combined_syst['mcEff']**2)**0.5
        mcErrD = dataStatD
        #print("dataEff:")
        #print(dataEff)
        #print("dataEffU:")
        #print(dataErrU)
        #print("dataEffD:")
        #print(dataErrD)
        if (mcEff + mcStatU) >= 1.0:
            mcErrU = 1.0 - mcEff
        else:
            mcErrU = mcStatU
        _out['value'] = sf
        _out['stat'] = sf_stat
        _out['syst'] = combined_syst['SF']
        for s in itertools.chain(systList['SF']['fitTypes'],
                                 systList['SF']['shiftTypes']):
            _out[s] = sf_syst[s]['err']

        # copy systs for later schema
        all_systematics[index] = _out.copy()

        def set_bin(hist, index, val, err):
            index = list(index)
            val_args = index + [val]
            err_args = index + [err]
            hist.SetBinContent(*val_args)
            if err >= 0:
                hist.SetBinError(*err_args)

        set_bin(hist, index, sf, sf_err)
        set_bin(hist_stat, index, sf, sf_stat)
        set_bin(hist_syst, index, sf, combined_syst['SF'])
        set_bin(histList_syst['combined_syst'], index,
                combined_syst['SF'], -1)
        set_bin(hist_dataEff, index, dataEff, dataErr)
        set_bin(hist_dataEff_errD, index, dataEff, dataErrD)
        set_bin(hist_dataEff_errU, index, dataEff, dataErrU)
        set_bin(hist_dataEff_stat, index, dataEff, dataStat)
        set_bin(hist_dataEff_syst, index, dataEff, combined_syst['dataEff'])
        set_bin(histList_dataEff_syst['combined_syst'], index,
                combined_syst['dataEff'], -1)
        set_bin(hist_mcEff, index, mcEff, mcErr)
        set_bin(hist_mcEff_errD, index, mcEff, mcErrD)
        set_bin(hist_mcEff_errU, index, mcEff, mcErrU)
        set_bin(hist_mcEff_stat, index, mcEff, mcStat)
        set_bin(hist_mcEff_syst, index, mcEff, combined_syst['mcEff'])
        set_bin(histList_mcEff_syst['combined_syst'], index,
                combined_syst['mcEff'], -1)
        for iKey in sf_syst.keys():
            if iKey in histList_syst:
                set_bin(histList_syst[iKey], index,
                        sf_syst[iKey]['sf'], sf_syst[iKey]['err'])
                set_bin(histList_syst[iKey+'_syst'], index,
                        sf_syst[iKey]['err'], -1)

            if iKey in histList_dataEff_syst:
                set_bin(histList_dataEff_syst[iKey], index,
                        sf_syst[iKey]['dataEff'], sf_syst[iKey]['dataErr'])
                set_bin(histList_dataEff_syst[iKey+'_syst'], index,
                        sf_syst[iKey]['dataErr'], -1)

            if iKey in histList_mcEff_syst:
                set_bin(histList_mcEff_syst[iKey], index,
                        sf_syst[iKey]['mcEff'], sf_syst[iKey]['mcErr'])
                set_bin(histList_mcEff_syst[iKey+'_syst'], index,
                        sf_syst[iKey]['mcErr'], -1)

    hists[extEffName] = hist
    hists[extEffName+'_stat'] = hist_stat
    hists[extEffName+'_syst'] = hist_syst
    hists[extEffName+'_efficiencyData'] = hist_dataEff
    hists[extEffName+'_efficiencyData_errD'] = hist_dataEff_errD
    hists[extEffName+'_efficiencyData_errU'] = hist_dataEff_errU
    hists[extEffName+'_efficiencyData_stat'] = hist_dataEff_stat
    hists[extEffName+'_efficiencyData_syst'] = hist_dataEff_syst
    hists[extEffName+'_efficiencyMC'] = hist_mcEff
    hists[extEffName+'_efficiencyMC_errD'] = hist_mcEff_errD
    hists[extEffName+'_efficiencyMC_errU'] = hist_mcEff_errU
    hists[extEffName+'_efficiencyMC_stat'] = hist_mcEff_stat
    hists[extEffName+'_efficiencyMC_syst'] = hist_mcEff_syst
    for iKey in histList_syst.keys():
        hname = extEffName+'_'+iKey
        hists[hname] = histList_syst[iKey]
    for iKey in histList_dataEff_syst.keys():
        hname = extEffName+'_efficiencyData_'+iKey
        hists[hname] = histList_dataEff_syst[iKey]
    for iKey in histList_mcEff_syst.keys():
        hname = extEffName+'_efficiencyMC_'+iKey
        hists[hname] = histList_mcEff_syst[iKey]

    # save the efficiency
    plotDir = os.path.join(baseDir, 'plots',
                           particle, probe,
                           resonance, era,
                           effName, 'efficiency')
    os.makedirs(plotDir, exist_ok=True)

    effDir = os.path.join(baseDir, 'efficiencies',
                          particle, probe,
                          resonance, era,
                          effName)
    os.makedirs(effDir, exist_ok=True)
    effPath = os.path.join(effDir, extEffName)

    # JSON format
    with open('{}.json'.format(effPath), 'w') as f:
        f.write(json.dumps(output, indent=4, sort_keys=True))

    # Now build the new xPOG schema v1 if correctionlib and pydantic installed

    schemav1 = None
    libname = 'correctionlib.schemav1'
    if libname in sys.modules:
        schemav1 = sys.modules[libname]
    elif importlib.util.find_spec('correctionlib.schemav1') is not None:
        spec = importlib.util.find_spec('correctionlib.schemav1')
        schemav1 = importlib.util.module_from_spec(spec)
        sys.modules[libname] = schemav1
        spec.loader.exec_module(schemav1)

    if schemav1 is not None:

        def build_schema(dim, index):
            # If we reach recursion bottom, build and return the systematics node
            if dim == len(variableLabels) + 1:
                keys, content = [], []
                for syst, value in all_systematics[index].items():
                    keys.append(syst)
                    content.append(value) 
                return schemav1.Category.parse_obj({
                    "nodetype": "category",
                    "keys": keys,
                    "content": content
                })
            # If not, build a binning node
            edges = list(map(float, binning[variableLabels[dim-1]]))
            content = [build_schema(dim+1, tuple(list(index)[0:dim-1]+[i]+list(index)[dim:])) for i in indices[dim-1]]
            return schemav1.Binning.parse_obj({
                "nodetype": "binning",
                "edges": edges,
                "content": content
            })

        inputs = [{"name": vl, "type": "real"} for vl in variableLabels]
        inputs += [{"name": "uncertainties", "type": "string"}]

        corr = schemav1.Correction.parse_obj({
                "version": 1,
                "name": effName,
                "description": effName,
                "inputs": inputs,
                "output": {"name": "weight", "type": "real"},
                "data": build_schema(1, tuple([1]*len(variableLabels)))
        })
        cset = schemav1.CorrectionSet.parse_obj({
            "schema_version": 1,
            "corrections": [corr]
        })

        # Write out schema json
        with open('{}_schemaV1.json'.format(effPath), "w") as fout:
            fout.write(cset.json(exclude_unset=True, indent=4))

    else:
        print("Warning: correctionlib not installed. Not producing schema jsons.")

    # ROOT histogram format
    tfile = ROOT.TFile.Open('{}.root'.format(effPath), 'recreate')
    for h in sorted(hists):
        hists[h].Write(h)

        if skipPlots:
            continue

        def setLog(canvas, hist, thr = 110.):
            if hist.GetXaxis().GetBinLowEdge(hist.GetXaxis().GetNbins()) > thr:
                hist.GetXaxis().SetMoreLogLabels()
                canvas.SetLogx()
            if hist.GetYaxis().GetBinLowEdge(hist.GetYaxis().GetNbins()) > thr:
                hist.GetYaxis().SetMoreLogLabels()
                canvas.SetLogy()
            if hist.GetZaxis().GetBinLowEdge(hist.GetZaxis().GetNbins()) > thr:
                canvas.SetLogz()

        if nVars == 2:
            cName = 'c' + h
            canvas = ROOT.TCanvas(cName, cName, 1000, 800)
            ROOT.gStyle.SetPaintTextFormat("5.3f")
            canvas.SetRightMargin(0.24)
            hists[h].Draw('colz text')
            plotPath = os.path.join(plotDir, h)
            canvas.Modified()
            canvas.Update()

            CMS_lumi.cmsText = 'CMS'
            CMS_lumi.writeExtraText = True
            CMS_lumi.extraText = 'Preliminary'
            #CMS_lumi.extraText = 'Work in progress'
            #CMS_lumi.lumi_13TeV = "%0.1f fb^{-1}" % (lumi)
            CMS_lumi.lumi_13p6TeV = "%0.1f fb^{-1}" % (lumi)
            CMS_lumi.CMS_lumi(canvas, 4, 11)

            if effType == 'trig':
                setLog(canvas, hists[h])

            canvas.Print('{}.png'.format(plotPath))
            canvas.Print('{}.pdf'.format(plotPath))
            canvas.Print('{}.root'.format(plotPath))

        elif nVars == 3:
            axes = [hists[h].GetXaxis(),
                    hists[h].GetYaxis(),
                    hists[h].GetZaxis()]
            axislabels = ['x', 'y', 'z']

            def zAxisTitle(effName):
                for iSyst in itertools.chain(systList['SF']['fitTypes'],
                                             systList['SF']['shiftTypes'],
                                             systList['dataEff']['fitTypes'],
                                             systList['dataEff']['shiftTypes'],
                                             systList['mcEff']['fitTypes'],
                                             systList['mcEff']['shiftTypes']):
                    if effName.endswith(iSyst+'_syst'):
                        return 'Uncertainty'

                if effName.endswith('combined_syst'):
                    return 'Uncertainty'
                elif 'efficiency' in effName:
                    return 'Efficiency'
                else:
                    return 'Scalefactor'

            for vi, variableLabel in enumerate(variableLabels):
                # if there are more than 2 bins, skip this variable
                if len(binning[variableLabel]) > 3:
                    continue

                projOpt = 'zyxe'.replace(axislabels[vi], '')
                for ibin in range(1, len(binning[variableLabel])):
                    axes[vi].SetRange(ibin, ibin)
                    projEffName = h.replace(variableLabel, variableLabel+'_{}'.format(ibin))
                    hist_proj = hists[h].Project3D(projOpt).Clone(projEffName)

                    hist_proj.GetZaxis().SetTitle(zAxisTitle(projEffName))

                    cName = 'c' + projEffName
                    canvas = ROOT.TCanvas(cName, cName, 1000, 800)
                    ROOT.gStyle.SetPaintTextFormat("5.3f")
                    canvas.SetRightMargin(0.24)
                    hist_proj.Draw('colz text')
                    plotPath = os.path.join(plotDir, projEffName)
                    canvas.Modified()
                    canvas.Update()

                    CMS_lumi.cmsText = 'CMS'
                    CMS_lumi.writeExtraText = True
                    CMS_lumi.extraText = 'Preliminary'
                    #CMS_lumi.extraText = 'Work in progress'
                    #CMS_lumi.lumi_13TeV = "%0.1f fb^{-1}" % (lumi)
                    CMS_lumi.lumi_13p6TeV = "%0.1f fb^{-1}" % (lumi)
                    CMS_lumi.CMS_lumi(canvas, 4, 11)

                    if effType == 'trig':
                        setLog(canvas, hist_proj)

                    canvas.Print('{}.png'.format(plotPath))
                    canvas.Print('{}.pdf'.format(plotPath))
                    canvas.Print('{}.root'.format(plotPath))

    tfile.Close()

    if skipPlots:
        return

    # gets a graph projection of an ND histogram for a given axis
    # with axis index (ie x,y,z = 0,1,2) and other dimensions ind
    def get_graph(histD, histU, axis, axis_ind, *ind):
        ind = list(ind)
        ni = axis.GetNbins()
        xvals = [axis.GetBinCenter(i+1) for i in range(ni)]
        xvals_errLow = [xvals[i]-axis.GetBinLowEdge(i+1) for i in range(ni)]
        xvals_errHigh = [axis.GetBinUpEdge(i+1)-xvals[i] for i in range(ni)]
        yvals = [
            histD.GetBinContent(
                *ind[:axis_ind]
                + [i+1]
                + ind[axis_ind:]
            ) for i in range(ni)]
        yvals_errLow = [
            histD.GetBinError(
                *ind[:axis_ind]
                + [i+1]
                + ind[axis_ind:]
            ) for i in range(ni)]
        yvals_errHigh = [
            histU.GetBinError(
                *ind[:axis_ind]
                + [i+1]
                + ind[axis_ind:]
            ) for i in range(ni)]
        #print("y vals")
        #print(yvals)
        #print("y vals low")
        #print(array('d', yvals_errLow))
        #print("y vals up")
        #print(array('d', yvals_errHigh))
        graph = ROOT.TGraphAsymmErrors(
            ni,
            array('d', xvals),
            array('d', yvals),
            array('d', xvals_errLow),
            array('d', xvals_errHigh),
            array('d', yvals_errLow),
            array('d', yvals_errHigh),
        )
        return graph

    # plot the efficiencies
    # some default colors for plots
    #colors = [ROOT.kBlack, ROOT.kViolet-2,ROOT.kViolet-2, ROOT.kRed, ROOT.kGreen+2, #for trackeronlyseed
    colors = [ROOT.kBlack, 65, 62, ROOT.kRed, ROOT.kGreen+2, #for alltracks
              ROOT.kMagenta+1, ROOT.kOrange+1, ROOT.kTeal-1,
              ROOT.kRed-3, ROOT.kBlue]

    def plot_1d_eff(savename, graphs,
                    labels=['Data', 'Drell-Yan (MC)'],
                    colors=colors,
                    xlabel='', ylabel='Efficiency',
                    xRange=[], additional_text=[]):
        ng = len(graphs)
        mg = ROOT.TMultiGraph()
        MCfilled = 0
        ratiofilled = 0
        meaneff_MC = 0
        meaneff_data = 0
        for gi in range(ng):
            if gi == 1:
                #graphs[gi].SetMarkerStyle(21)
                graphs[gi].SetFillColor(colors[gi])
                mg.Add(graphs[gi],'AE2')
                graphs[gi].SetLineColor(colors[gi+1])
                graphs[gi].SetMarkerColor(colors[gi])
                graphs[gi].SetLineWidth(1)
                MCfilled = 1
                meaneff_MC = graphs[gi].GetMean(2)
                meaneff_MC = meaneff_MC*100
                meaneff_MC = "%.2f" % meaneff_MC
            if  MCfilled: #gi == 0:
                mg.Add(graphs[0],'AP0')
                graphs[0].SetLineColor(colors[0])
                graphs[0].SetMarkerColor(colors[0])
                graphs[0].SetMarkerSize(1.6)
                graphs[0].SetLineWidth(2)
                meaneff_data = graphs[0].GetMean(2)
                meaneff_data = meaneff_data*100
                meaneff_data = "%.2f" % meaneff_data
                #if ratiofilled == 0:
                #    ratio = ROOT.DividTGraphs(graphs[0], graphs[1])
                #    ratiofilled = 1
        #graphs[0].SetLineColor(colors[0])
        #graphs[0].SetMarkerColor(colors[0])
        #graphs[0].SetLineWidth(2)
        #mg.Add(graphs[0],'AP0')
        #graphs[1].SetLineColor(colors[1])
        #graphs[1].SetFillColor(colors[1])
        #graphs[1].SetLineWidth(2)
        #mg.Add(graphs[1],'AE2')
        
        canvas = ROOT.TCanvas(savename, savename, 800, 800)
        # Added 2
        #pad1 = ROOT.TPad("pad1", "pad1", 0, 0.2, 1, 1.0)
        #pad1.SetBottomMargin(0.05)
        #pad1.SetTopMargin(0.06)
        #pad1.Draw()
        #pad1.cd()
        #ratio = ROOT.TH1F()
        mg.Draw('A')
        mg.GetXaxis().SetTitle(xlabel)
        if xRange:
            #mg.GetXaxis().SetLimits(*xRange)
            #mg.GetXaxis().SetLabelSize(0.) # to make ratio plots
            mg.GetXaxis().SetLabelSize(0.04)
            mg.GetXaxis().SetLabelOffset(0.01)
            mg.GetXaxis().SetTitleSize(0.05)
            mg.GetXaxis().SetTitleOffset(1)
            mg.GetXaxis().SetRangeUser(*xRange)
            #mg.GetXaxis().CenterLabels(ROOT.kTRUE)
        if not 'fake' in savename:
            mg.GetYaxis().SetTitle(ylabel)
        else:
            mg.GetYaxis().SetTitle('Fake matching rate')
        mg.GetYaxis().SetTitleSize(0.05)
        mg.GetYaxis().SetLabelSize(0.04)
        mg.GetYaxis().SetLabelOffset(0.01)
        mg.GetYaxis().SetTitleOffset(1.5)
        if effType == 'trig':
            mg.GetYaxis().SetRangeUser(0.6, 1.20)
        else:
            if not 'fake' in savename:
                mg.GetYaxis().SetRangeUser(0.85, 1.10)
                #mg.GetYaxis().SetRangeUser(0.97, 1.02 )
            else:
                mg.GetYaxis().SetRangeUser(0.000001, 0.4)
                ROOT.gPad.SetLogy()
        if not 'fake' in savename:
            legend = ROOT.TLegend(0.20, 0.3, 0.5, 0.45)
        else:
            legend = ROOT.TLegend(0.20, 0.6, 0.5, 0.75)
        legend.SetTextFont(42)
        legend.SetTextSize(0.03)
        legend.SetBorderSize(0)
        legend.SetFillColor(0)
        meanvalue = ROOT.TPaveText(0.65, 0.82, 0.9, 0.92,"NDC")
        meanvalue.SetTextAlign(13)
        meanvalue.SetFillColor(0)
        meanvalue.SetTextFont(52)
        meanvalue.SetBorderSize(0)
        #meanvalue.AddText("< {}_{{ Data}} > = {} %}".format('\epsilon',meaneff_data))
        #meanvalue.AddText("< \epsilon_{{ MC  }} > = {} %}".format(meaneff_MC))
        meanvalue.SetTextFont(42)
        meanvalue.SetTextSize(0.03)
        meanvalue.Draw()
        
        for gi in range(ng):
            #legend.AddEntry(graphs[gi], labels[gi], 'LP')
            if gi == 0:
                legend.AddEntry(graphs[gi], labels[gi], 'LP')
            if gi == 1:
                legend.AddEntry(graphs[gi], labels[gi], 'F')
        #legend.SetHeader('Tracker-only Seeded Tracks','C')# {} / {}'.format(num, denom))
        #legend.SetHeader('All Tracks','h')# {} / {}'.format(num, denom))
        legend.Draw()

        if additional_text:
            nother = len(additional_text)+1
            #dims = [0.55, 0.52, 0.93, 0.28]
            #dims = [0.55, 0.65-nother*0.04-0.02, 0.75+nother*0.055+0.02, 0.3]
            if not 'fake' in savename:
                dims = [0.57, 0.3, 0.85, 0.5]#alltracks
                #dims = [0.56, 0.3, 0.85, 0.5] # trackeronly
            else:
                #dims = [0.56, 0.6, 0.85, 0.8] # trackeronly
                dims = [0.56, 0.6, 0.85, 0.8] #alltracks
            text = ROOT.TPaveText(*dims+['NB NDC'])
            text.AddText("All Tracks") #alltracks
            #text.AddText("Tracker-only Tracks") #trackeronly
            text.SetTextSize(0.04)
            text.SetTextAlign(12)
            text.SetTextFont(61)
            for rtext in additional_text:
                text.AddText(rtext)
                text.GetListOfLines().Last().SetTextFont(42)
                text.GetListOfLines().Last().SetTextAlign(11)
                text.GetListOfLines().Last().SetTextSize(0.03)
            text.SetBorderSize(0)
            text.SetFillColor(0)
            text.Draw()

        CMS_lumi.cmsText = 'CMS'
        CMS_lumi.writeExtraText = True
        CMS_lumi.extraText = 'Preliminary'
        #CMS_lumi.extraText = 'Work in progress'
        CMS_lumi.lumi_13p6TeV = "%0.1f fb^{-1}" % (lumi)
        #CMS_lumi.lumi_13TeV = "%0.1f fb^{-1}" % (lumi)
        CMS_lumi.CMS_lumi(canvas, 4, 11)

        # Added 2 for Make ratioplots
        #pad2 = ROOT.TPad("pad2", "pad2", 0, 0., 1, 0.21)
        #pad2.SetTopMargin(0.01) 
        #pad2.SetBottomMargin(0.35) 
        #pad2.SetGridy() 
        #pad2.Draw()
        #pad2.cd()
        #ratio.SetStats(0)
        #ratio.SetTitle("")
        #ratio.SetLineWidth(2)
        #ratio.SetLineColor(1)
        #ratio.SetMarkerStyle(20)
        #ratio.SetMarkerColor(1)
        #ratio.GetXaxis().SetTitle(xlabel)
        #if xRange:
        #    mg.GetXaxis().SetLimits(*xRange)
        #    ratio.GetXaxis().SetLabelSize(0.2)
        #    ratio.GetXaxis().SetLabelOffset(0.01)
        #    ratio.GetXaxis().SetTitleSize(0.2)
        #    ratio.GetXaxis().SetTitleOffset(1)
        #    ratio.GetXaxis().SetRangeUser(*xRange)
        #    #mg.GetXaxis().CenterLabels(ROOT.kTRUE)
        #ratio.GetYaxis().SetRangeUser(*xRange)
        #ratio.GetYaxis().SetRangeUser(0.96,1.04)
        #ratio.GetYaxis().SetTitle("Data/MC")
        #ratio.GetYaxis().SetNdivisions(505)
        #ratio.GetYaxis().SetTitleSize(0.1)
        #ratio.GetYaxis().SetLabelSize(0.1)
        #ratio.GetYaxis().SetLabelOffset(0.05)
        #ratio.GetYaxis().SetTitleOffset(0.01)
        #ratio.GetXaxis().SetTitle(xlabel)
        #ratio.GetXaxis().SetTitleFont(42)
        #ratio.Draw("same")
        #canvas.Update()
        canvas.Modified()
        canvas.Update()
        canvas.Print('{}.png'.format(savename))
        canvas.Print('{}.pdf'.format(savename))
        canvas.Print('{}.root'.format(savename))

        # save each graph
        tfile = ROOT.TFile('{}.root'.format(savename), 'update')
        for gi in range(ng):
            graphs[gi].SetTitle(labels[gi])
            graphs[gi].Write('g_{}_{}'.format(gi, labels[gi]))
        tfile.Close()

    # enumerate over the axis/variable to plot
    axes = [hists[extEffName].GetXaxis(),
            hists[extEffName].GetYaxis(),
            hists[extEffName].GetZaxis()]
    for vi, variableLabel in enumerate(variableLabels):

        # iterate over the other axis indices
        otherVariableLabels = [ovl for ovl in variableLabels
                               if ovl != variableLabel]
        otherVariableIndices = [ovi for ovi, ovl in enumerate(variableLabels)
                                if ovl != variableLabel]
        indices = [list(range(1, len(binning[vl])))
                   for vl in otherVariableLabels]
        if indices:
            for index in itertools.product(*indices):
                #print("get Graph Data")
                graph_data = get_graph(hists[extEffName+'_efficiencyData_errD'],
                                       hists[extEffName+'_efficiencyData_errU'],
                                       axes[vi], vi, *index)
                #print("get MC Data")
                graph_mc = get_graph(hists[extEffName+'_efficiencyMC_errD'],
                                     hists[extEffName+'_efficiencyMC_errU'],
                                     axes[vi], vi, *index)
                xlabel = get_variable_name_pretty(variableLabel)
                ylabel = 'Efficiency'
                xRange = [axes[vi].GetBinLowEdge(1),
                          axes[vi].GetBinUpEdge(axes[vi].GetNbins())]
                additional_text = []
                for novi, (ovi, ovl) in enumerate(zip(otherVariableIndices,
                                                      otherVariableLabels)):
                    xlow = axes[ovi].GetBinLowEdge(index[novi])
                    xhigh = axes[ovi].GetBinUpEdge(index[novi])
                    rtext = '{} < {} < {}'.format(
                        xlow, get_variable_name_pretty(ovl), xhigh)
                    additional_text += [rtext]
                plotDir = os.path.join(baseDir, 'plots',
                                       particle, probe,
                                       resonance, era,
                                       effName, 'efficiency')
                os.makedirs(plotDir, exist_ok=True)
                otherVariableLabel = get_bin_name(otherVariableLabels, index)
                plotName = '{}_{}_vs_{}'.format(effName,
                                                otherVariableLabel,
                                                variableLabel)
                plotPath = os.path.join(plotDir, plotName)
                plot_1d_eff(plotPath, [graph_data, graph_mc],
                            xlabel=xlabel, ylabel=ylabel,
                            xRange=xRange, additional_text=additional_text)

                # dataEfficiency systs
                graphs = [get_graph(hists[extEffName+'_efficiencyData_errD'],
                                    hists[extEffName+'_efficiencyData_errU'],
                                    axes[vi], vi, *index)]
                labels = ['Nominal']
                for iSyst in itertools.chain(
                        systList['dataEff']['fitTypes'],
                        systList['dataEff']['shiftTypes']):
                    graphs += [get_graph(
                        hists[extEffName+'_efficiencyData_'+iSyst],
                        hists[extEffName+'_efficiencyData_'+iSyst],
                        axes[vi], vi, *index)]
                    labels += [iSyst]
                plotName = '{}_{}_vs_{}_efficiencyData_syst'.format(
                    effName,
                    otherVariableLabel,
                    variableLabel,
                )
                plotPath = os.path.join(plotDir, plotName)
                plot_1d_eff(plotPath, graphs,
                            labels=labels,
                            xlabel=xlabel, ylabel=ylabel,
                            xRange=xRange, additional_text=additional_text)

                # mcEfficiency systs
                graphs = [get_graph(hists[extEffName+'_efficiencyMC_errD'],
                                    hists[extEffName+'_efficiencyMC_errU'],
                                    axes[vi], vi, *index)]
                labels = ['Nominal']
                for iSyst in itertools.chain(
                        systList['mcEff']['fitTypes'],
                        systList['mcEff']['shiftTypes']):
                    graphs += [get_graph(
                        hists[extEffName+'_efficiencyMC_'+iSyst],
                        hists[extEffName+'_efficiencyMC_'+iSyst],
                        axes[vi], vi, *index)]
                    labels += [iSyst]
                plotName = '{}_{}_vs_{}_efficiencyMC_syst'.format(
                    effName,
                    otherVariableLabel,
                    variableLabel,
                )
                plotPath = os.path.join(plotDir, plotName)
                plot_1d_eff(plotPath, graphs,
                            labels=labels,
                            xlabel=xlabel, ylabel=ylabel,
                            xRange=xRange, additional_text=additional_text)

        # if no indices, easier, just itself
        else:
            
            graph_data = get_graph(hists[extEffName+'_efficiencyData_errD'],
                                   hists[extEffName+'_efficiencyData_errU'],
                                   axes[vi], vi)
         
            graph_mc = get_graph(hists[extEffName+'_efficiencyMC_errD'],
                                 hists[extEffName+'_efficiencyMC_errU'],
                                 axes[vi], vi)

            xlabel = get_variable_name_pretty(variableLabel)
            ylabel = 'Efficiency'
            xRange = [axes[0].GetBinLowEdge(1),
                      axes[0].GetBinUpEdge(axes[0].GetNbins())]
            plotDir = os.path.join(baseDir, 'plots',
                                   particle, probe,
                                   resonance, era,
                                   effName, 'efficiency')
            os.makedirs(plotDir, exist_ok=True)
            plotName = '{}_vs_{}'.format(effName, variableLabel)
            plotPath = os.path.join(plotDir, plotName)
            plot_1d_eff(plotPath, [graph_data, graph_mc],
                        xlabel=xlabel, ylabel=ylabel,
                        xRange=xRange)


def build_prepare_corrected_jobs(particle, probe, resonance, era,
                       config, **kwargs):
    _baseDir = kwargs.pop('baseDir', '')
    _numerator = kwargs.pop('numerator', [])
    _denominator = kwargs.pop('denominator', [])
    _registry = kwargs.pop('registry', None)

    if _registry is not None:
        registry.reset()
        registry.load_json(_registry)

    subEra = era.split('_')[0]  # data subera is beginning of era
    lumi = registry.luminosity(particle, probe, resonance, era, subEra, **kwargs)

    jobs = []
    # iterate through the efficiencies
    efficiencies = config.efficiencies()
    for num, denom in efficiencies:
        if _numerator and num not in _numerator:
            continue
        if _denominator and denom not in _denominator:
            continue

        # iterate through the output binning structure
        for variableLabels in config.binVariables():

            jobs += [[_baseDir, particle, probe, resonance, era,
                     config, num, denom, tuple(variableLabels), lumi]]

    return jobs
