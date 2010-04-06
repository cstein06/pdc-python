from numpy import *

import pdc.states as sta_
from pdc.globals import *

import matplotlib.pyplot as pp
from scipy.io.matlab.mio import loadmat

def main_analysis():
    
    #root = 'G:\\stein\\dados\\teste edu\\'
    root = 'G:\\stein\\dados\\teste edu\\gotas psd\\'
        
    #input = root + 'ES60_21_07_09_melhores4.txt'
    #input = 'ES57_13_02_09_melhores3_test.txt'
    input = 'ES59_13_07_09_melhores3.txt'
    #input = root + 'test.txt'
    
    instates = 'ES59_13_07_09_estagiamentojanela10s_limpo.txt'
    
    #algoritmo = 'pdc'
    algoritmo = 'coh'
    
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
    
    pr_.ss = True
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
    
    res, mea, stds, nstates = sta_.states_analysis_bind(data, states)

        #read with: medias2 = permute(reshape(medias', shape(4), shape(3), 
        # shape(2), shape(1)), [4,3,2,1]);

    return res, mea, stds, nstates

def batch_analysis():
    
    #root = 'C:\\Documents and Settings\\Stein\\My Documents\\dados\\teste edu\\'
    root = 'G:\\stein\\dados\\teste edu\\'
    #root = 'G:\\stein\\dados\\teste edu\\gotas psd\\'
        
    #input = 'C:/Documents and Settings/Stein/Desktop/teste edu/baccala2001a_ex4_sim_01.txt'
    #input = root + 'ES57_09_02_09_medias.txt'
    #input = root + 'ES57_09_02_09_medias_curto.txt'

    #inputs = ['ES57_09_02_09_medias_test.txt', 'ES57_09_02_09_medias_test2.txt']
    #inputs = ['ES60_21_07_09_melhores4.txt', 'ES59_16_07_09_melhores3.txt']
    inestags = ['ES57_09_02_09_estagiamentojanela10s_limpo.txt','ES57_09_02_09_estagiamentojanela10s_limpo.txt']
    
#    d1 = loadmat(root+'ES57_AD23_dia_09_02_09.mat')
#    d2 = loadmat(root+'ES57_AD23combed_dia_09_02_09.mat')
#    d3 = loadmat(root+'ES57_09_02_09_AD23_combed_runicado.mat')
#    d4 = loadmat(root+'ES57_09_02_09_AD23_limpo_sem_filt.mat')
#    d5 = loadmat(root+'ES57_09_02_09_AD23_runicado.mat')
    
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
    
        set_params(alg = algoritmo, window_size = window_size, 
                   nf = n_frequencies, sample_f = sampling_rate,
                   maxp = ordem_max, fixp = ordem_fixa, detrend = detrend, 
                   do_window_log = do_window_log,
                   power = espectro_em_potencia, metric = metrica_pdc, 
                   do_plot = plota, plotf = plot_freq, plot_states = plot_states,
                   stinput = inputs[i], root_dir = root)    
        
        data = loadtxt(root+inputs[i]).T
        estag = loadtxt(root+inestags[i])
                
        print 'Data loaded:', inputs[i]
        res, meds, stds, nstates = sta_.states_analysis(data, estag)
    
    return res, meds, stds, nstates
  


def check_hist():
    
    root = "G:\\stein\\dados\\edu_comp\\"
    
    d1 = loadmat(root+'ES57_AD23_dia_09_02_09.mat')
    d2 = loadmat(root+'ES57_AD23combed_dia_09_02_09.mat')
    d3 = loadmat(root+'ES57_09_02_09_AD23_combed_runicado.mat')
    d4 = loadmat(root+'ES57_09_02_09_AD23_limpo_sem_filt.mat')
    d5 = loadmat(root+'ES57_09_02_09_AD23_runicado.mat')
    
    d5 = d5['ADrunicado']
    d4 = d4['ADlimpo']
    d3 = d3['ADcombed_runicado']
    d2 = d2['ad_combed']
    d1 = d1['ad']
    
    pp.subplot(2,3,1)
    pp.specgram(d1[:1000000:2,0], NFFT = 1000)
    pp.subplot(2,3,2)
    pp.specgram(d2[0,:500000], NFFT = 1000)
    pp.subplot(2,3,3)
    pp.specgram(d3[0,:500000], NFFT = 1000)
    pp.subplot(2,3,4)
    pp.specgram(d4[0,:500000], NFFT = 1000)
    pp.subplot(2,3,5)
    pp.specgram(d5[0,:500000], NFFT = 1000)
    
    d1 = loadmat(root+'ES59_AD22_dia_13_07_09.mat')
    d2 = loadmat(root+'ES59_AD22combed_dia_13_07_09.mat')
    d3 = loadmat(root+'ES59_13_07_09_AD22_limpo_sem_filt.mat')
    d3 = d3['ADlimpo']
    d2 = d2['ad_combed']
    d1 = d1['ad_downsampled']

    pp.figure()
    pp.subplot(2,3,1)
    pp.specgram(d1[:1000000:2,0], NFFT = 1000)
    pp.subplot(2,3,2)
    pp.specgram(d2[0,:500000], NFFT = 1000)
    pp.subplot(2,3,3)
    pp.specgram(d3[0,:500000], NFFT = 1000)
    
if __name__ == '__main__':
    
    main_analysis()
    
    #batch_analysis()
