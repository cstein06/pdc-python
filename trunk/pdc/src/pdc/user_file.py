# -*- coding:utf-8 -*-
"""
This is a file you user could use for your code.

Mini instructions:

Use an_.set_params() to set parameters. 
See available parameters in globals.py->Params.
Call an_.measure_full(data)

"""

import pdc.analysis as an_
import pdc.ar_data as ar_
import pdc.ar_fit as fit_
import pdc.plotting as pl_
import pdc.asymp as ass
from pdc.globals import *
import pdc.states as sta_

from numpy import *

import matplotlib.pyplot as pp

def simple_simulation_analysis():
    '''Simple test of connectivity routines
    
    Here we simulate some data using the ar_data module, from a given
    MVAR matrix A and covariance matrix er.
    
    '''
    
    #Definition of the MVAR model
    A = array([[[0.2, 0],[0.3,-0.2],[0.3,-0.2]], 
               [[0, 0],[0.8,-0.1],[0.4,-0.1]],
               [[0, 0],[0.3,0.2],[0.4,0.1]]], dtype = float) 
    er = identity(3)
    
    #number of data points generated
    nd = 2000
    
    #number of frequency points analyzed
    nf = 64
    
    #error of the confidence interval and threshold
    alpha = 0.05
    
    #model order parameters
    n = A.shape[0]
    maxp = A.shape[2]
    
    #type of PDC used (refer to manual to see what it means)
    metric = 'diag'
    
    #Generate data from AR
    data = ar_.ar_data(A, er, nd)
    
    # Set parameters for analysis
    an_.set_params(nf = nf, ss = True, metric = metric,
                   detrend = True)
    
    #Call any connectivity routine routine. 
    #Here are some calling examples, uncomment your preferred one for use:
    
    # Set which connectivity method to use
    an_.set_params(alg = 'coh')
    # Call the method
    an_.measure_full(data)
    
    #Instead of setting the 'alg' parameter, you can call directly:
    #an_.coh_full(data), an_.pdc_full(data), etc.
        
    pp.show()
    
def states_simulation_analysis():
    
    root = 'G:\\stein\\dados\\teste edu\\gotas psd\\'
        
    input = 'ES59_13_07_09_melhores3.txt'
    
    instates = 'ES59_13_07_09_estagiamentojanela10s_limpo.txt'
    
    #algoritmo = 'pdc'
    algoritmo = 'coh'
    
    window_size = 10
    n_frequencies = 250
    sampling_rate = 500
    
    plot_labels = ['Ca1e', 'Ca3e', 'Ca1d', 'Ca2d', 'Ca3d']
    plot_states = array([1,2,3,4,5,6])
    plot_freq = 150
    plota = True
    do_window_log = True
    
    #pr_.ss = True
    #pr_.logss = False
    #pr_.plot_diag = True
    valid_states = [1,2,3,4,5,6]
    ordem_max = 25
    ordem_fixa = True
    detrend = True
    espectro_em_potencia = True
    metrica_pdc = 'diag'
    
    ####################################################
    
    #nao mexer daqui pra frente
    
    set_params(alg = algoritmo, window_size = window_size, 
               nf = n_frequencies, sample_f = sampling_rate,
               maxp = ordem_max, fixp = ordem_fixa, detrend = detrend, 
               do_states_log = do_window_log,
               power = espectro_em_potencia, metric = metrica_pdc, 
               do_plot = plota,  plot_labels = plot_labels, plotf = plot_freq,
               root_dir = root, stinput = input, plot_states = plot_states,
               valid_states = valid_states)
    
    data = loadtxt(root+input).T
    states = loadtxt(root+instates)
    
    print 'Data loaded:', input
    
    res, mea, stds, nstates = sta_.states_analysis(data, states)
    
    return res, mea, stds, nstates

def window_simulation_analysis():
    
    root = 'G:\\stein\\dados\\teste edu\\gotas psd\\'
        
    input = 'ES59_13_07_09_melhores3.txt'
    
    #algoritmo = 'pdc'
    algoritmo = 'coh'
    
    window_size = 10
    n_frequencies = 250
    sampling_rate = 500
    
    plot_labels = ['Ca1e', 'Ca3e', 'Ca1d', 'Ca2d', 'Ca3d']
    plot_states = array([1,2,3,4,5,6])
    plot_freq = 150
    plota = True
    do_window_log = True
    
    #pr_.ss = True
    #pr_.logss = False
    #pr_.plot_diag = True
    valid_states = [1,2,3,4,5,6]
    ordem_max = 25
    ordem_fixa = True
    detrend = True
    espectro_em_potencia = True
    metrica_pdc = 'diag'
    
    ####################################################
    
    #nao mexer daqui pra frente
    
    set_params(alg = algoritmo, window_size = window_size, 
               nf = n_frequencies, sample_f = sampling_rate,
               maxp = ordem_max, fixp = ordem_fixa, detrend = detrend, 
               do_states_log = do_window_log,
               power = espectro_em_potencia, metric = metrica_pdc, 
               do_plot = plota,  plot_labels = plot_labels, plotf = plot_freq,
               root_dir = root, stinput = input, plot_states = plot_states,
               valid_states = valid_states)
    
    data = loadtxt(root+input).T
    
    print 'Data loaded:', input
    
    res = sta_.window_analysis(data)
    
    return res

if __name__ == '__main__':
    
    simple_simulation_analysis()
    #states_simulation_analysis()
    #window_simulation_analysis()
    
    
