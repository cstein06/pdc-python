from numpy import *

import time
import os

from scipy.io import savemat

import projects.edu_estagios as edu

from pdc.globals import *

def main_analysis():
    
    #root = 'C:\\Documents and Settings\\Stein\\My Documents\\dados\\teste edu\\'
    root = 'G:\\stein\\dados\\teste edu\\'
        
    #input = 'C:/Documents and Settings/Stein/Desktop/teste edu/baccala2001a_ex4_sim_01.txt'
    #input = root + 'ES57_09_02_09_medias.txt'
    #input = root + 'ES57_09_02_09_medias_curto.txt'
    #input = root + 'ES57_09_02_09_medias_test.txt'
    #input = root + 'ES57_09_02_09_melhores.txt'
    #input = root + 'ES60_21_07_09_melhores4.txt'
    input = root + 'ES60_21_07_09_melhores4.txt'
    input = root + 'test.txt'
    
    inestag = root + 'ES57_09_02_09_estagiamentojanela10s_limpo.txt'
    
    outputres = root + 'ES57_09_02_09_medias_res'
    outputmeds = root + 'ES57_09_02_09_medias_meds'

    algoritmo = 'pdc'
    #algoritmo = 'coh'
    
    window_size = 10
    n_frequencies = 250
    sampling_rate = 500
    
    #plot_labels = ['Ca1e', 'Ca2e', 'Ca1d', 'Ca2d', 'Ca3d']
    plot_labels = ['Ca1e', 'Ca2e', 'Ca1d']
    plot_states = array([1,2,3,4,5,6])
    #plot_states = array([2])
    plot_freq = 150
    plota = True
    do_window_log = True
    
    ordem_max = 25
    ordem_fixa = True
    detrend = True
    espectro_em_potencia = True
    metrica_pdc = 'diag'
    
    ####################################################
    
    #nao mexer daqui pra frente
    
    data = loadtxt(input).T
    estag = loadtxt(inestag)
    
    print 'Data loaded:', input

    set_params(alg = algoritmo, window_size = window_size, 
               nf = n_frequencies, sample_f = sampling_rate,
               maxp = ordem_max, fixp = ordem_fixa, detrend = detrend, 
               do_window_log = do_window_log,
               power = espectro_em_potencia, metric = metrica_pdc, 
               do_plot = plota, plotf = plot_freq, plot_labels = plot_labels)
    
    res, meds, stds = edu.states_analysis(data, estag, plot_states = plot_states)

    if pr_.do_window_log:
            
        if os.path.isfile(outputres + '.mat'):
            print '\nOverwriting result .mat file!'
            
        if os.path.isfile(outputmeds + '.mat'):
            print '\nOverwriting mean .mat file!'
        
        savemat(outputres, {'result':res, 'shape':res.shape, 'time':time.ctime()})
        
        savemat(outputmeds, {'medias':meds, 'stds':stds, 
                             'shape':meds.shape, 'time':time.ctime()})
        #read with: medias2 = permute(reshape(medias', shape(4), shape(3), 
        # shape(2), shape(1)), [4,3,2,1]);

    return res, meds, stds

def batch_analysis():
    
    #root = 'C:\\Documents and Settings\\Stein\\My Documents\\dados\\teste edu\\'
    root = 'G:\\stein\\dados\\teste edu\\'
        
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
        inputs[i] = root + inputs[i]
    #input = root + 'ES57_09_02_09_melhores.txt'
    
    for i in arange(size(inestags)):
        inestags[i] = root + inestags[i]
        
        
    outputres = []
    outputmeds = []
    for i in arange(size(inputs)):
        if inputs[i][-4:] == '.txt':
            auxin = inputs[i][:-4]
        else:
            auxin = inputs[i]
        outputres.append(auxin + '_' + algoritmo + '_res')
        outputmeds.append(auxin + '_' + algoritmo + '_meds')
    
    
    for i in arange(size(inputs)):
        
        data = loadtxt(inputs[i]).T
        estag = loadtxt(inestags[i])
        
        print 'Data loaded:', inputs[i]
    
        set_params(alg = algoritmo, window_size = window_size, 
                   nf = n_frequencies, sample_f = sampling_rate,
                   maxp = ordem_max, fixp = ordem_fixa, detrend = detrend, 
                   do_window_log = do_window_log,
                   power = espectro_em_potencia, metric = metrica_pdc, 
                   do_plot = plota, plotf = plot_freq)
        
        res, meds, stds = edu.states_analysis(data, estag, plot_states = plot_states)
    
        if pr_.do_window_log:
            
            import os
            if os.path.isfile(outputres[i] + '.mat'):
                print '\nOverwriting result .mat file!'
                
            if os.path.isfile(outputmeds[i] + '.mat'):
                print '\nOverwriting mean .mat file!'
            
            savemat(outputres[i], {'result':res, 'shape':res.shape, 'time':time.ctime()})
            
            savemat(outputmeds[i], {'mean':meds, 'stds':stds, 
                                    'shape':meds.shape, 'time':time.ctime()})
            #read with: medias2 = permute(reshape(medias', shape(4), shape(3), 
            # shape(2), shape(1)), [4,3,2,1]);

    #return res, meds, stds
  


    
    
if __name__ == '__main__':
    
    main_analysis()
    
    #batch_analysis()
