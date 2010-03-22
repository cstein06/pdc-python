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
    
    aux = an_.measure(data[:,0*win:0*win+win], ss = False)
    
    
    resp = zeros([nwins, n, n, pr_.nf], dtype = aux.dtype)
    resp[0] = aux
    
    for i in arange(1,nwins):
        set_params(v = False)
        
        aux = an_.measure(data[:,i*win:i*win+win], ss = False)
        resp[i] = aux
        
        print '\nProcessed', i+1, 'of', nwins, 'windows:', 100*(i+1.0)/nwins, '%'
    
    return resp

def mean_estag(mes, estag, maxe = 6, nulle = -1):
        
    nwins,n,n,nf = mes.shape
    
    if mes.dtype == 'complex':
        print 'Taking mean of complex measure!`'
    
    mpdc = zeros([maxe, n, n, nf], dtype = mes.dtype)
    spdc = zeros([maxe, n, n, nf], dtype = mes.dtype)
    m2 = zeros([maxe, n, n, nf], dtype = mes.dtype)
    s2 = zeros([maxe, n, n, nf], dtype = mes.dtype)
    
    estag = estag[:nwins]
    
    for i in arange(maxe):
        if sum(estag-1 == i) > 0:
            mpdc[i] = mean(mes[estag-1 == i], 0)
            spdc[i] = std(mes[estag-1 == i], 0)
            
            
    return mpdc, spdc
        
def states_analysis(data, states, plot_states = None, plot_freq = None, **args):
        
    read_args(args)

    tim = time.clock()
    
    
    result = window_analysis(data)
    
    nwins = result.shape[0]
    
    states = states[:nwins]
    
    nstates = histogram(states, new = True, bins=arange(1,8))[0]
    
    print '\nNumber of windows used:', states.shape
    
    mpdc, spdc = mean_estag(result, states)
    
    nf = mpdc.shape[3]
    
    print '\nNumber of windows per state:', nstates
    
    if pr_.do_plot:
        for i in plot_states:
            if nstates[i-1] > 0:
                pr_.plot_color = state_colors[i-1]
                pl.pdc_plot(mpdc[i-1,:,:,:plot_freq])
        
    print '\nTotal time in secs:', time.clock() - tim
    
    pp.show()
    
    return result, mpdc, spdc


def mean_estag_A(mes, estag, maxe = 6, nulle = -1):
        
    nwins,n,n,dum = mes[0].shape
    
    mpdcA = zeros([maxe, n, n, pr_.maxp])
    spdcA = zeros([maxe, n, n, pr_.maxp])
    mpdcE = zeros([maxe, n, n])
    spdcE = zeros([maxe, n, n])
    
    estag = estag[:nwins]
    
    for i in arange(maxe):
        mpdcA[i] = mean(mes[0][estag-1 == i], 0)
        spdcA[i] = std(mes[0][estag-1 == i], 0)
        
    for i in arange(maxe):
        mpdcE[i] = mean(mes[1][estag-1 == i], 0)
        spdcE[i] = std(mes[1][estag-1 == i], 0)
        
    return mpdcA, spdcA, mpdcE, spdcE

def window_analysis_A(data, **args):
    
    read_args(args)
    
    print data.shape
    n, T = data.shape
    
    win = int(pr_.window_size*pr_.sample_f)
    
    nwins = T/win
    
    aux = fit_.ar_fit(data[:,0*win:0*win+win], pr_.maxp, criterion = 1)
    
    respA = zeros([nwins, n, n, pr_.maxp])
    respE = zeros([nwins, n, n])
    respA[0] = aux[0]
    respE[0] = aux[1]
    
    for i in arange(1,nwins):
        
        aux = fit_.ar_fit(data[:,i*win:i*win+win], pr_.maxp, criterion = 1)
        #aux = an_.measure(data[:,i*win:i*win+win], ss = False)
        respA[i] = aux[0]
        respE[i] = aux[1]
        
        print 'Processed', i+1, 'of', nwins, 'windows:', 100*(i+1.0)/nwins, '%'
        set_params(v = False)
    
    return respA,respE

def states_analysis_A(data, states, plot_states = None, plot_freq = None, **args):
        
    read_args(args)

    tim = time.clock()
    
    result = window_analysis_A(data)
    
    nwins = result[0].shape[0]
    
    states = states[:nwins]
    
    nstates = histogram(states, new = True, bins=arange(1,8))[0]
    
    print states.shape
    
    mpdcA, spdcA, mpdcE, spdcE = mean_estag_A(result, states)
    
#    dum,a,b,c = mpdcA.shape
    
#    mpdc = zeros([6, a, b, c])
#    for i in arange(6):
#        mpdc[i] = an_.pdc_alg(mpdcA[i], mpdcE[i], ss = False)
#        mpdc[i] = abs(mpdc[i])**2
#        
#        
#    print nstates
#    if pr_.do_plot:
#        for i in plot_states:
#            if nstates[i-1] > 0:
#                pl.pdc_plot(mpdc[i-1,:,:,:plot_freq])
#        
#    pp.show()
    
    print 'tempo gasto:', time.clock() - tim
    
    #savetxt('C:/Documents and Settings/Stein/Desktop/teste edu/baccala2001pdc.txt', pdcan)
    #savetxt('C:/Documents and Settings/Stein/Desktop/teste edu/baccala2001coh.txt', cohan)
    
    return result, mpdcA, spdcA, mpdcE, spdcE


def main_analysis():
    
    root = 'C:\\Documents and Settings\\Stein\\My Documents\\dados\\teste edu\\'
        
    #input = 'C:/Documents and Settings/Stein/Desktop/teste edu/baccala2001a_ex4_sim_01.txt'
    #input = root + 'ES57_09_02_09_medias.txt'
    #input = root + 'ES57_09_02_09_medias_curto.txt'
    input = root + 'ES57_09_02_09_medias_test.txt'
    #input = root + 'ES57_09_02_09_melhores.txt'
    
    inestag = root + 'ES57_09_02_09_estagiamentojanela10s_limpo.txt'
    
    
    outputres = root + 'ES57_09_02_09_medias_res'
    outputmeds = root + 'ES57_09_02_09_medias_meds'
    
    plot_states = array([1,2,3,4,5,6])
    #plot_states = array([2])

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
    
    
    #nao mexer daqui pra frente
    
    data = loadtxt(input).T
    estag = loadtxt(inestag)
    
    print 'Data loaded'

    set_params(alg = algoritmo, window_size = window_size, 
               nf = n_frequencies, sample_f = sampling_rate,
               maxp = ordem_max, fixp = ordem_fixa, detrend = detrend, 
               power = espectro_em_potencia, metric = metrica_pdc, do_plot = plota)
    
    res, meds, stds = states_analysis(data, estag, plot_states = plot_states, plot_freq = plot_freq)

    if pr_.do_window_log:
        savemat(outputres, {'result':res})
        
        savemat(outputmeds, {'medias':meds, 'stds':stds, 'shape':meds.shape})
        #read with: medias2 = permute(reshape(medias', shape(4), shape(3), 
        # shape(2), shape(1)), [4,3,2,1]);

    return res, meds, stds

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

    #Aest, er = fit_.ar_fit(data[:,:win], ordem_max, criterion=1) 

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
    
    #testa_std_asymp()
    main_analysis()
    #testa_aic()
    #testa_ordens()
    
    