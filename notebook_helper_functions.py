import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

from flattener import *
from config import Configuration
from notebook_helper_functions import *

def create_binning_structure(config, subEra):
    doGen = subEra in ['DY_madgraph', 'JPsi_pythia8']
    
    fitVariable = config.fitVariable()
    binningSet = set([fitVariable])
    if doGen:
        fitVariableGen = config.fitVariableGen()
        binningSet = binningSet.union(set([fitVariableGen]))
    binVariables = config.binVariables()
    for bvs in binVariables:
        binningSet = binningSet.union(set(bvs))
    return binningSet

def plot(dfs, plotVar, edges, subEras, xlabel='', additional_text=''):
    data_subera, mc_subera = subEras[0], subEras[1]
    
    def reduce_df(df):
        df = get_binned_dataframe(df, 'tmpBin', plotVar, edges)
        df = df.groupBy('tmpBin').agg({'weight': 'sum', 'weight2': 'sum'})
        df = df.withColumnRenamed('sum(weight)','sumw')
        df = df.withColumnRenamed('sum(weight2)','sumw2')
        df = df.toPandas().set_index('tmpBin').sort_index()
        return df
    
    df_data = reduce_df(dfs[data_subera])
    df_mc = reduce_df(dfs[mc_subera])
    
    x_center = [(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)]
    x_width = [abs(x_center[i]-edges[i]) for i in range(len(edges)-1)]
    
    def make_plot_content(df):
        y = pd.Series(np.zeros(len(edges)+1))
        yerr = pd.Series(np.zeros(len(edges)+1))
        y[df.index] = df['sumw']
        integral = sum(y)
        y /= integral
        y = y[1:len(edges)]
        yerr[df.index] = np.sqrt(df['sumw2']) 
        yerr /= integral
        yerr = yerr[1:len(edges)]
        return y, yerr
    
    y_data, yerr_data = make_plot_content(df_data)
    y_mc, yerr_mc = make_plot_content(df_mc)
    
    plt.errorbar(x_center, y_data, xerr=x_width, yerr=yerr_data, fmt='.k')
    plt.errorbar(x_center, y_mc, xerr=x_width, yerr=yerr_mc, fmt='.b')

    plt.legend(['Data', 'Simulation'], frameon=False)        
    plt.ylabel('Unit normalized')
    plt.xlabel(xlabel)
    
    if additional_text != '':
        plt.text(edges[int(len(edges)/2)], 0.90*max(*y_data, *y_mc), additional_text)
        
    plt.show()
    

def plot_eff(dfs, num, den, config, plotVar, subEras, xlabel='', additional_text=''):
    
    data_subera, mc_subera = subEras[0], subEras[1]
    edges = config.binning()[plotVar]
    
    df_data = dfs[data_subera]
    df_mc = dfs[mc_subera]
    
    x_center = [(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)]
    x_width = [abs(x_center[i]-edges[i]) for i in range(len(edges)-1)]
    
    # remove under and overflow bins 
    # careful: this is a df index slicing, not array slicing, so *does have* inclusive range
    y_data = df_data['eff'][1:len(edges)-1].copy()
    yerr_data = df_data['eff_err'][1:len(edges)-1].copy()
    y_mc = df_mc['eff'][1:len(edges)-1].copy()
    yerr_mc = df_mc['eff_err'][1:len(edges)-1].copy()
    
    fig, (ax0, ax1) = plt.subplots(2, sharex=True, gridspec_kw={'hspace': 0, 'height_ratios':[3, 1]})
    
    ax0.errorbar(x_center, y_data, xerr=x_width, yerr=yerr_data, fmt='.k', label='Data')
    ax0.errorbar(x_center, y_mc, xerr=x_width, yerr=yerr_mc, fmt='.b', label='Simulation')
    ax0.legend(frameon=False)
    ax0.set_ylabel('Efficiency')
    ax0.set_ylim(0.88,1.1)
    ax0.yaxis.set_minor_locator(AutoMinorLocator())
    ax0.set_xlabel(xlabel)
    ax0.text(0.05, 0.90, f'{num} / {den}', transform=ax0.transAxes)
    ax0.text(0.05, 0.82, additional_text, transform=ax0.transAxes)
    ax0.label_outer()
    
    ratio = y_data/y_mc
    err_ratio = ratio * ((yerr_data / y_data)**2 + (yerr_mc / y_mc)**2)**0.5
    ax1.errorbar(x_center, y_data/y_mc, xerr=x_width, yerr=err_ratio, fmt='.r')
    ax1.set_ylabel('Data/Simulation')
    ax1.set_ylim(0.95,1.02)
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax1.set_xlabel(xlabel)
    ax1.label_outer()
    
    fig.show()