import pdc.analysis as an_

from numpy import *

import time

import matplotlib.pyplot as pp

import pdc.ar_fit as fit_

from scipy.linalg.basic import det

import pdc.plotting as pl

from pdc.globals import *

def window_analysis(data, **args):
    
    read_args(args)
    
    print data.shape
    n, T = data.shape
    
    win = pr_.window_size*pr_.sample_f
    
    nwins = T/win
    
    resp = zeros([nwins, n, n, pr_.nf])
    
    for i in arange(nwins):
        resp[i] = an_.measure(data[:,i*win:i*win+win], ss = False)
        
        print 'Processed', i+1, 'of', nwins, 'windows:', 100*(i+1.0)/nwins, '%'
    
    return resp

def mean_estag(mes, estag, maxe = 6, nulle = -1):
    
    nwins,n,n,nf = mes.shape
    
    mpdc = zeros([maxe, n, n, nf])
    npdc = zeros(maxe)
    
    for i in arange(nwins):
        if (estag[i] == nulle):
            continue
        mpdc[estag[i]-1] += mes[i]
        npdc[estag[i]-1] += 1
    
    print npdc
    
    for i in arange(maxe):
        if npdc[i] > 0:            
            mpdc[i] = mpdc[i]/npdc[i]
            
    return mpdc
        
def main_analysis(data = None):
    
    root = 'C:\\Documents and Settings\\Stein\\My Documents\\dados\\teste edu\\'
        
    #input = 'C:/Documents and Settings/Stein/Desktop/teste edu/baccala2001a_ex4_sim_01.txt'
    #input = root + 'ES57_09_02_09_medias.txt'
    #input = root + 'ES57_09_02_09_medias_curto.txt'
    input = root + 'ES57_09_02_09_medias_test.txt'
    #input = root + 'ES57_09_02_09_melhores.txt'
    
    inestag = root + 'ES57_09_02_09_estagiamentojanela10s_limpo.txt'
    
    plot_estag = array([1,2,3,4,5,6])

    algoritmo = 'pdc'
    #algoritmo = 'coh'
    window_size = 10
    n_frequencies = 250
    plot_freq = 150
    sampling_rate = 500
    ordem_max = 25
    ordem_fixa = True
    detrend = True
    espectro_em_potencia = False
    metrica_pdc = 'diag'


    #nao mexer a partir daqui

    tim = time.clock()
    
    set_params(alg = algoritmo, window_size = window_size, 
               nf = n_frequencies, sample_f = sampling_rate,
               maxp = ordem_max, fixp = ordem_fixa, detrend = detrend, 
               power = espectro_em_potencia, metric = metrica_pdc)
      
    if (data == None):  
        data = loadtxt(input).T
    
    #data  = data[:,:5000]
    
    print 'data loaded'
    
    result = window_analysis(data)
    
    estag = loadtxt(inestag)
    
    print estag.shape
    mpdc = mean_estag(result, estag)
    
    nf = mpdc.shape[3]
    
    for i in plot_estag-1:
        pl.pdc_plot(mpdc[i,:,:,:plot_freq])
        
    pp.show()
    
    print 'tempo gasto:', time.clock() - tim
    
    #savetxt('C:/Documents and Settings/Stein/Desktop/teste edu/baccala2001pdc.txt', pdcan)
    #savetxt('C:/Documents and Settings/Stein/Desktop/teste edu/baccala2001coh.txt', cohan)
    
    return result, mpdc
    
def testa_aic():
    
#    input = 'C:/Documents and Settings/Stein/Desktop/teste edu/ES57_09_02_09_medias.txt'
    input = 'C:/Documents and Settings/Stein/Desktop/teste edu/ES57_09_02_09_medias_curto.txt'
           
    
    data = loadtxt(input).T
    
    print 'leu dados'
    
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
    
    main_analysis()
    #testa_aic()
    #testa_ordens()
    
    