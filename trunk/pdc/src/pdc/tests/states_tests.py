
from numpy import *

import time

import matplotlib.pyplot as pp

import scipy.stats as st
from scipy.linalg.basic import det
from scipy.io import savemat

from pdc import *


def windows_test():
    
    data = ar_data(m = 20000, model = 0)
    
    #algoritmo = 'pdc'
    algoritmo = 'coh'
    
    window_size = 100
    n_frequencies = 64
    sampling_rate = 1
    
    #plot_labels = ['Ca1e', 'Ca3e', 'Ca1d', 'Ca2d', 'Ca3d']
    #plot_labels = ['Ca1e', 'Ca2e', 'Ca1d']
    #plot_states = array([1,2,3,4,5,6])
    #plot_states = array([2])
    #plot_freq = 150
    plota = False
    do_window_log = False
    
    pr_.ss = True
    #pr_.logss = False
    #pr_.plot_diag = True
    valid_states = [1,2,3,4,5,6]
    ordem_max = 3
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
               do_plot = plota,
               stinput = input,
               valid_states = valid_states)
    
    #data = loadtxt(root+input).T
    
    print 'Data loaded:', input
    
    res = window_analysis(data)
    
    plot_coherogram(res)

    
    
    pp.show()

    

        #read with: medias2 = permute(reshape(medias', shape(4), shape(3), 
        # shape(2), shape(1)), [4,3,2,1]);

    #return res, mea, stds, nstates

def main_analysis():
    
    #root = 'G:\\stein\\dados\\teste edu\\'
    #root = 'G:\\stein\\dados\\teste edu\\gotas psd\\'
    root = '/home/stein/dados/teste edu/gotas psd/'
        
    #input = root + 'ES60_21_07_09_melhores4.txt'
    input = 'ES57_13_02_09_melhores3_test.txt'
    #input = 'ES59_13_07_09_melhores3.txt'
    #input = root + 'test.txt'
    
    instates = 'ES59_13_07_09_estagiamentojanela10s_limpo.txt'
    
    #algoritmo = 'pdc'
    algoritmo = 'coh'
    
    window_size = 10
    n_frequencies = 250
    sampling_rate = 500
    
    plot_labels = ['Ca1e', 'Ca3e', 'Ca1d', 'Ca2d', 'Ca3d']
    #plot_labels = ['Ca1e', 'Ca2e', 'Ca1d']
    plot_states = array([1,2,3,4,5,6])
    #plot_states = array([2])
    plot_freq = 150
    plota = False
    do_window_log = False
    
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
    
    res, mea, stds, nstates = states_analysis(data, states)
    
    plot_coherogram(res, states)

    
    
    pp.show()

    

        #read with: medias2 = permute(reshape(medias', shape(4), shape(3), 
        # shape(2), shape(1)), [4,3,2,1]);

    return res, mea, stds, nstates

def group_analysis():
    
    root = 'C:\\Documents and Settings\\Stein\\My Documents\\dados\\teste edu\\'
        
    #input = 'C:/Documents and Settings/Stein/Desktop/teste edu/baccala2001a_ex4_sim_01.txt'
    input = root + 'ES57_09_02_09_medias.txt'
    #input = root + 'ES57_09_02_09_medias_curto.txt'
    #input = root + 'ES57_09_02_09_medias_test.txt'
    #input = root + 'ES57_09_02_09_melhores.txt'
    
    inestag = root + 'ES57_09_02_09_estagiamentojanela10s_limpo.txt'
    
    plot_states = array([1,2,3,4,5,6])
    plot_states = array([6])

    algoritmo = 'pdc'
    #algoritmo = 'coh'
    window_size = 10
    n_frequencies = 250
    plot_freq = 150
    sampling_rate = 500
    ordem_max = 25
    ordem_fixa = True
    detrend = True
    espectro_em_potencia = True
    metrica_pdc = 'diag'
    plota = True
    
    
    #nao mexer daqui pra frente
    
    data = loadtxt(input).T[:2]
    estag = loadtxt(inestag)
    
    print 'data loaded'

    set_params(alg = algoritmo, window_size = window_size, 
               nf = n_frequencies, sample_f = sampling_rate,
               maxp = ordem_max, fixp = ordem_fixa, detrend = detrend, 
               power = espectro_em_potencia, metric = metrica_pdc, do_plot = plota)
    
    n,T = data.shape
    
    win = int(window_size*sampling_rate)
    nwin = T/win
    dm = win*(nwin/2)
    data1 = data[:,:dm]
    data2 = data[:,dm:]
    estag1 = estag[:nwin/2]
    estag2 = estag[nwin/2:]
    
    r1 = states_analysis(data1, estag1, plot_states = plot_states, plot_freq = plot_freq)
    r2 = states_analysis(data2, estag2, plot_states = plot_states, plot_freq = plot_freq)

    return r1, r2

def testa_std_asymp():
    
    root = 'C:\\Documents and Settings\\Stein\\My Documents\\dados\\teste edu\\'
        
    #input = 'C:/Documents and Settings/Stein/Desktop/teste edu/baccala2001a_ex4_sim_01.txt'
    input = root + 'ES57_09_02_09_medias.txt'
    #input = root + 'ES57_09_02_09_medias_curto.txt'
    #input = root + 'ES57_09_02_09_medias_test.txt'
    #input = root + 'ES57_09_02_09_melhores.txt'
    
    inestag = root + 'ES57_09_02_09_estagiamentojanela10s_limpo.txt'
    
    plot_states = array([1,2,3,4,5,6])
    #plot_states = array([])

    algoritmo = 'pdc'
    #algoritmo = 'coh'
    window_size = 10
    n_frequencies = 250
    plot_freq = 150
    sampling_rate = 500
    ordem_max = 25
    ordem_fixa = True
    detrend = True
    espectro_em_potencia = True
    metrica_pdc = 'diag'
    plota = False
    alpha = 0.05
    
    est = 4
    
    #nao mexer daqui pra frente
    
    data = loadtxt(input).T[:2,:]
    estag = loadtxt(inestag)
    
    print 'data loaded'

    set_params(alg = algoritmo, window_size = window_size, 
               nf = n_frequencies, sample_f = sampling_rate,
               maxp = ordem_max, fixp = ordem_fixa, detrend = detrend, 
               power = espectro_em_potencia, metric = metrica_pdc, do_plot = plota)
    
    res, mer, str = states_analysis(data, estag, plot_states = plot_states, plot_freq = plot_freq)

    win = int(window_size*sampling_rate)

    res, merA, strA, merE, strE = states_analysis(data, estag, plot_states = plot_states, plot_freq = plot_freq)

    #Aest, er = fit_.ar_estim(data[:,:win], ordem_max, fixp = True) 

    #simd = ar_data(Aest, er, win)
    simd = ar_data(merA[est], merE[est], win)
    
    
    asy = vars()[algoritmo + '_full'](simd)
     
    plot_all()
    
    set_results(mes = mer[est])
    for i in arange(size(mer[est])):
        res_.ic1.flat[i] = st.norm.ppf(alpha, mer[est].flat[i], str[est].flat[i])
        res_.ic2.flat[i] = st.norm.ppf(1-alpha, mer[est].flat[i], str[est].flat[i])
    set_results(th = zeros(mer[est].shape))
    
    plot_all()
    
    #pp.show()
    
    return res, mer, str, asy
        
    
def testa_aic():
    
#    input = 'C:/Documents and Settings/Stein/Desktop/teste edu/ES57_09_02_09_medias.txt'
    input = 'C:/Documents and Settings/Stein/Desktop/teste edu/ES57_09_02_09_medias_curto.txt'
           
    
    data = loadtxt(input).T
    
    print 'Data loaded'
    
    n, T = data.shape
    
    aicn = 50
    
    vaic = zeros(aicn)
    
    win = 10*1000/2
    
    pp.ion()
    
    for j in arange(10):
        
        print 'j', j
        
        daux = data[:,j*win:(j+1)*win]
        
        for i in arange(1,aicn+1):
    
            [npf, na, npb, nb, nef, neb, ISTAT]=nstrand(daux,i,False)
        
            vaic[i-1]=win*log(det(npf))+2*n*n*i
            
            if i > 1: 
                if vaic[i-1] > vaic[i-2]:
                    print 'min:', i-1
                    
        print 'min global:', argmin(vaic)
        
        pp.plot(vaic)
        
    pp.show()
    
    
def testa_ordens():
    
#    input = 'C:/Documents and Settings/Stein/Desktop/teste edu/ES57_09_02_09_medias.txt'
    input = 'C:/Documents and Settings/Stein/Desktop/teste edu/ES57_09_02_09_medias_curto.txt'
    
    data = loadtxt(input).T
    
    print 'leu dados'
    
    n, T = data.shape
    
    #ords = array([10, 15, 20, 25, 30])
    #ords = array([10, 20, 30])
    ords = array([20])
    
    win = 10*1000/2
    sf = 500
    
    
    
    #data2 = data[:,::2]
    #win2 = win/2
    #ords2 = ords/2
    #sf2 = sf/2
    
    wins = 5
    
    for ord in ords:
        for t in arange(wins):
            daux = data[:,100000+t*win:100000+t*win+win]
            #daux2 = data[:,100000:100000+win/2]
            #daux3 = data[:,100000+win/2:100000+win]
            #daux4 = data[:,100000:100000+win/4]
            #daux5 = data[:,100000:100000+win/10]
            #print an_.igct(daux, maxp = 20, fixp = True)
            #print an_.white_test(daux, maxp = 20, h = 50, fixp = True)
        
            #pdcr = an_.pdc_full(daux, maxp = ord, metric = 'gen',
            #                         fixp = True, ss = True, sample_f = sf)
            pdcr = pdc_full(daux, maxp = ord, metric = 'diag',
                                     fixp = True, ss = True, sample_f = sf)
            #pdcr = an_.coh_and_plot(daux2, maxp = ord, 
            #                         fixp = True, ss = False, sample_f = sf)
            #pdcr = an_.coh_and_plot(daux4, maxp = ord, 
            #                         fixp = True, ss = False, sample_f = sf)
            #pdcr = an_.coh_and_plot(daux5, maxp = ord, 
            #                         fixp = True, ss = False, sample_f = sf)
            #pdcr = an_.coh_and_plot(data2[:win2], maxp = ord/2, 
            #                         fixp = True, ss = False, sample_f = sf/2)
        
    pp.show()
    
if __name__ == '__main__':
    pass
    #testa_std_asymp()
    #testa_aic()
    #testa_ordens()
    #main_analysis()
    windows_test()
