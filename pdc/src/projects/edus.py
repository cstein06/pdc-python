from numpy import *

import time
import os

from scipy.io import savemat

import pdc.states as sta_

from pdc.globals import *

def main_analysis():
    
    #root = 'G:\\stein\\dados\\teste edu\\'
    root = 'G:\\stein\\dados\\teste edu\\gotas psd\\'
        
    #input = root + 'ES60_21_07_09_melhores4.txt'
    input = 'ES57_13_02_09_melhores3_test.txt'
    #input = root + 'test.txt'
    
    inestag = 'ES57_13_02_09_estagiamentojanela10s_limpo.txt'
    
    #algoritmo = 'pdc'
    algoritmo = 'ar'
    
    window_size = 10
    n_frequencies = 250
    sampling_rate = 500
    
    plot_labels = ['Ca1e', 'Ca2e', 'Ca1d', 'Ca2d', 'Ca3d']
    #plot_labels = ['Ca1e', 'Ca2e', 'Ca1d']
    plot_states = array([1,2,3,4,5,6])
    #plot_states = array([2])
    plot_freq = 150
    plota = True
    do_window_log = True
    
    valid_states = [1,2,3,4,5,6]
    ordem_max = 25
    ordem_fixa = True
    detrend = True
    espectro_em_potencia = True
    metrica_pdc = 'diag'
    
    ####################################################
    
    #nao mexer daqui pra frente
    
    data = loadtxt(root+input).T
    estag = loadtxt(root+inestag)
    
    print 'Data loaded:', input
    

    set_params(alg = algoritmo, window_size = window_size, 
               nf = n_frequencies, sample_f = sampling_rate,
               maxp = ordem_max, fixp = ordem_fixa, detrend = detrend, 
               do_window_log = do_window_log,
               power = espectro_em_potencia, metric = metrica_pdc, 
               do_plot = plota, plotf = plot_freq, plot_labels = plot_labels,
               root_dir = root, stinput = input, plot_states = plot_states,
               valid_states = valid_states)
    
    res, meds, stds, nstates = sta_.states_analysis(data, estag)

        #read with: medias2 = permute(reshape(medias', shape(4), shape(3), 
        # shape(2), shape(1)), [4,3,2,1]);

    return res, meds, stds

def batch_analysis():
    
    #root = 'C:\\Documents and Settings\\Stein\\My Documents\\dados\\teste edu\\'
    #root = 'G:\\stein\\dados\\teste edu\\'
    root = 'G:\\stein\\dados\\teste edu\\gotas psd\\'
        
    #input = 'C:/Documents and Settings/Stein/Desktop/teste edu/baccala2001a_ex4_sim_01.txt'
    #input = root + 'ES57_09_02_09_medias.txt'
    #input = root + 'ES57_09_02_09_medias_curto.txt'

    #inputs = ['ES57_09_02_09_medias_test.txt', 'ES57_09_02_09_medias_test2.txt']
    inputs = ['ES60_21_07_09_melhores4.txt', 'ES59_16_07_09_melhores3.txt']
    inestags = ['ES57_09_02_09_estagiamentojanela10s_limpo.txt','ES57_09_02_09_estagiamentojanela10s_limpo.txt']
    
    #for i in arange(10):
    #    inestags.append('ES57_09_02_09_estagiamentojanela10s_limpo.txt')
    
    algoritmo = 'pdc'
    #algoritmo = 'coh'
    
    window_size = 10
    n_frequencies = 250
    sampling_rate = 500
    
    plot_states = array([1,2,3,4,5,6])
    #plot_states = array([2])
    plot_freq = 150
    plota = False
    do_window_log = True
    
    ordem_max = 25
    ordem_fixa = True
    detrend = True
    espectro_em_potencia = True
    metrica_pdc = 'diag'
    
    #########################################
    
    #nao mexer daqui pra frente
    
    for i in arange(size(inputs)):
        
        data = loadtxt(root+inputs[i]).T
        estag = loadtxt(root+inestags[i])
        
        
        
        print 'Data loaded:', inputs[i]
    
        set_params(alg = algoritmo, window_size = window_size, 
                   nf = n_frequencies, sample_f = sampling_rate,
                   maxp = ordem_max, fixp = ordem_fixa, detrend = detrend, 
                   do_window_log = do_window_log,
                   power = espectro_em_potencia, metric = metrica_pdc, 
                   do_plot = plota, plotf = plot_freq, plot_states = plot_states,
                   stinput = inputs[i])
        
        res, meds, stds = sta_.states_analysis(data, estag)
    

    #return res, meds, stds
  


    
    
if __name__ == '__main__':
    
    main_analysis()
    
    #batch_analysis()
