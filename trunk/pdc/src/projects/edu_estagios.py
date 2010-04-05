import pdc.analysis as an_

from numpy import *

import time

import matplotlib.pyplot as pp

import scipy.stats as st
from scipy.linalg.basic import det
from scipy.io import savemat

import pdc.plotting as pl
import pdc.ar_fit as fit_
from pdc.globals import *
from pdc.ar_data import ar_data 

def window_analysis(data, **args):
    
    pr_.do_log = False 
    
    read_args(args)
    
    print '\nTotal data shape:', data.shape
    n, T = data.shape
    
    win = int(pr_.window_size*pr_.sample_f)
    
    nwins = T/win
    
    print '\nCalculating first window verbosely (others will be silent):\n'
    aux_v = pr_.v
    
    aux = an_.measure(data[:,0*win:0*win+win], ss = False)
    
    resp = []
    if type(aux) is 'tuple':
        for i in arange(len(aux)):
            resp.append(zeros([nwins]+list(aux[i].shape), dtype = aux[i].dtype))
            resp[i][0] = aux
    else:
        resp = [zeros([nwins]+list(aux.shape), dtype = aux.dtype)]
            
    for i in arange(0,nwins):
        
        #print i*win
        
        aux = an_.measure(data[:,i*win:i*win+win], ss = False)
        for j in arange(len(resp)):
            resp[j][i] = aux[j]
        
        print '\nProcessed', i+1, 'of', nwins, 'windows:', 100*(i+1.0)/nwins, '%'
        
        pr_.v = False
    
    pr_.v = aux_v
    
    return resp

def mean_states(mes, states):
        
    nwins,n,n,nf = mes.shape
    
    if mes.dtype == 'complex':
        print 'Taking mean of complex measure!`'
    
    mpdc = zeros([len(pr_.valid_states), n, n, nf], dtype = mes.dtype)
    spdc = zeros([len(pr_.valid_states), n, n, nf], dtype = mes.dtype)
    
    states = states[:nwins]
    
    if (nwins > size(states)):
        print 'More windows than states in the states file!'
    
    for i in pr_.valid_states:
        if sum(states == i) > 0:
            auxi = pr_.st_dict[i]
            mpdc[auxi] = mean(mes[states == i], 0)
            spdc[auxi] = std(mes[states == i], 0)
            
    return mpdc, spdc

def mean_states_list(mes, estag, maxe = 6):
        
    mpdc = []
    spdc = []
    for i in arange(len(mes)):
        aux = mean_states(mes[i], estag)
        mpdc.append(aux[0])
        spdc.append(aux[1])
            
    return mpdc, spdc
        
def states_analysis(data, states, **args):
        
    read_args(args)

    tim = time.clock()
    
    
    result = window_analysis(data)
    
    nwins = result[0].shape[0]
    
    states = states[:nwins]
    
    if pr_.valid_states is None:
        set_params(valid_states = list(sorted(set(states[states >= 0]))))
    
    nstates = zeros(len(pr_.valid_states))
    
    pr_.st_dict = {}
    for i in arange(len(pr_.valid_states)):
        pr_.st_dict[pr_.valid_states[i]] = i
        nstates[i] = sum(states == pr_.valid_states[i])
    
    #nstates = histogram(states, bins=arange(1,8))[0]
    
    print '\nNumber of windows used:', states.shape[0]
    
    mpdc, spdc = mean_states_list(result, states)
    
    #nf = mpdc.shape[3]
    
    print '\nNumber of windows per state:', nstates
        
    if pr_.do_plot:
        if pr_.plot_states is None:
            pr_.plot_states = pr_.valid_states
        for st in pr_.plot_states:
            i = pr_.st_dict[st]
            if nstates[i] > 0:
                pr_.plot_color = pr_.state_colors[i]
                res_.mes = mpdc[0][i]
                pl.plot_all()
        pp.show()
        
    if len(mpdc) == 1:
        result = result[0]
        mpdc = mpdc[0]
        spdc = spdc[0]
                
    print '\nTotal time in secs:', time.clock() - tim
    
    if pr_.do_states_log:
        log_windows_results(result, mpdc, spdc, nstates)
    
    return result, mpdc, spdc, nstates



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

    res, merA, strA, merE, strE = states_analysis_A(data, estag, plot_states = plot_states, plot_freq = plot_freq)

    #Aest, er = fit_.ar_fit(data[:,:win], ordem_max, fixp = True) 

    #simd = ar_data(Aest, er, win)
    simd = ar_data(merA[est], merE[est], win)
    
    asy = vars(an_)[algoritmo + '_full'](simd)
     
    pl.plot_all()
    
    set_results(mes = mer[est])
    for i in arange(size(mer[est])):
        res_.ic1.flat[i] = st.norm.ppf(alpha, mer[est].flat[i], str[est].flat[i])
        res_.ic2.flat[i] = st.norm.ppf(1-alpha, mer[est].flat[i], str[est].flat[i])
    set_results(th = zeros(mer[est].shape))
    
    pl.plot_all()
    
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
    
            [npf, na, npb, nb, nef, neb, ISTAT]=fit_.nstrand(daux,i,False)
        
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
            pdcr = an_.pdc_full(daux, maxp = ord, metric = 'diag',
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
