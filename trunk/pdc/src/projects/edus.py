from numpy import *

import pdc.states as sta_
from pdc.globals import *

import matplotlib.pyplot as pp
from scipy.io.matlab.mio import loadmat
import pdc.plotting as pl_
import pdc.analysis as an

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
    
    #root = 'G:\\stein\\dados\\teste edu\\'
    #root = 'G:\\stein\\dados\\teste edu\\gotas psd\\'
    root = "G:\\stein\\dados\\edu_comp\\"
        
    #input = 'C:/Documents and Settings/Stein/Desktop/teste edu/baccala2001a_ex4_sim_01.txt'
    #input = root + 'ES57_09_02_09_medias.txt'
    #input = root + 'ES57_09_02_09_medias_curto.txt'

    inputs = ['ES57_09_02_09_medias_test.txt', 'ES57_09_02_09_medias_test2.txt']
    #inputs = ['ES60_21_07_09_melhores4.txt', 'ES59_16_07_09_melhores3.txt']
    instates = ['ES57_09_02_09_estagiamentojanela10s_limpo.txt','ES57_09_02_09_estagiamentojanela10s_limpo.txt']
    
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
    
    set_params(alg = algoritmo, window_size = window_size, 
               nf = n_frequencies, sample_f = sampling_rate,
               maxp = ordem_max, fixp = ordem_fixa, detrend = detrend, 
               do_window_log = do_window_log,
               power = espectro_em_potencia, metric = metrica_pdc, 
               do_plot = plota, plotf = plot_freq, plot_states = plot_states,
               root_dir = root)   
    
    for i in arange(size(inputs)): 
        
        pr_.stinput = inputs[i]
        
        data = loadtxt(root+inputs[i]).T
        state = loadtxt(root+instates[i])
                
        print 'Data loaded:', inputs[i]
        res, meds, stds, nstates = sta_.states_analysis(data, state)
    
    return res, meds, stds, nstates
  
ds = [['09_02_09', '10_02_09', '11_02_09', '12_02_09', '13_02_09'],
      ['13_07_09', '14_07_09', '15_07_09', '16_07_09', '17_07_09'],
      ['19_07_09', '20_07_09', '21_07_09', '22_07_09']]

files = [['ES%(r)s_%(d)s_AD%(c)d_limpo_sem_filt.mat',
         'ES%(r)s_AD%(c)d_dia_%(d)s.mat',
         'ES%(r)s_AD%(c)dcombed_dia_%(d)s.mat',
         'ES%(r)s_%(d)s_AD%(c)d_combed_runicado.mat',
         'ES%(r)s_%(d)s_AD%(c)d_runicado.mat'],
         ['ES%(r)s_%(d)s_AD%(c)d_limpo_sem_filt.mat',
          'ES%(r)s_AD%(c)d_dia_%(d)s.mat',
          'ES%(r)s_AD%(c)dcombed_dia_%(d)s.mat'],
         ['ES%(r)s_%(d)s_AD%(c)d_limpo_sem_filt.mat',
          'ES%(r)s_AD%(c)d_dia_%(d)s.mat',
          'ES%(r)s_AD%(c)dcombed_dia_%(d)s.mat']]

rs = ['57', '59', '60']

chs = [12, 21, 46, 38, 53]
         
cdata = 'ES%(r)s_%(d)s_AD%(c)d_limpo_sem_filt.mat'
         
         
#    d1 = loadmat(root+'ES59_AD22_dia_13_07_09.mat')
#    d2 = loadmat(root+'ES59_AD22combed_dia_13_07_09.mat')
#    d3 = loadmat(root+'ES59_13_07_09_AD22_limpo_sem_filt.mat')
#    d3 = d3['ADlimpo']
#    d2 = d2['ad_combed']
#    d1 = d1['ad_downsampled']

fields = [['ADlimpo',
          'ad',
          'ad_combed',
          'ADcombed_runicado',
          'ADrunicado'],
          ['ADlimpo',
           'ad_downsampled',
          'ad_combed'],
          ['ADlimpo',
          'ad_downsampled',
          'ad_combed']]

chs = [[12,21,46,38,53],
       [14,27,22,47,38,54],
       [14,12,22,47,38,54]]

ls = [['Ca1e', 'Ca3e', 'Ca1d', 'Ca2d', 'Ca3d'],
      ['Ca1e', 'Ca2e', 'Ca3e','Ca1d', 'Ca2d', 'Ca3d'],
      ['Ca1e', 'Ca2e', 'Ca3e','Ca1d', 'Ca2d', 'Ca3d']]

def get_data(t = 0, r = 0, d = 0, cs = None):
    root = "G:\\stein\\dados\\edu_comp\\"

    if cs is None:
        cs = arange(len(chs[r]))
        
    #pr_.plot_labels = ls[r][cs]
    #print ls[r], ls[r][cs], cs
    
    
    
    field = fields[r][t]
    fi = files[r][t]
    
    da1 = loadmat(root+fi % {'d':ds[r][d], 'r':rs[r], 'c':chs[r][cs[0]]})
    d1s = size(da1[field])
    
    d1 = zeros([len(cs), d1s])

    d1[0] = da1[field].ravel()
    
    for c in arange(len(cs)):
        fst = root+fi % {'d':ds[r][d], 'r':rs[r], 'c':chs[r][cs[c]]}
        da1 = loadmat(fst)
        d1[c] = da1[field].ravel()
        
    if t == 1:
        #downsample raw data that comes in 1000HZ to 500HZ like the others
        return d1[:,::2]
        
    return d1

def get_res(r = 0, d = 0, alg = 'coh'):
    #root = "G:\\stein\\dados\\edu_comp\\results\\"
    
    fi = 'R%s_D%s_%s_res.pic' % (rs[r], ds[r][d], alg) 
    
    return load_win_pic(fi)


def get_state(r = 0, d = 0):
    root = "G:\\stein\\dados\\edu_comp\\"
    
    fst = 'ES%s_%s_estagiamentojanela10s_limpo.txt' % (rs[r], ds[r][d])
    
    return loadtxt(root+fst)

def set_def():
    #algoritmo = 'pdc'
    pr_.alg = 'coh'
    
    pr_.window_size = 10
    pr_.nf = 250
    pr_.sample_f = 500
    
    pr_.plot_labels = ['Ca1e', 'Ca2e', 'Ca1d', 'Ca2d', 'Ca3d']
    #plot_labels = ['Ca1e', 'Ca2e', 'Ca1d']
    #plot_states = array([1,2,3,4,5,6])
    #plot_states = array([2])
    pr_.plotf = 120
    pr_.plota = True
    pr_.plot_ic = True
    pr_.do_window_log = False
    
    pr_.ss = True
    pr_.logss = True
    #pr_.plot_diag = True
    #valid_states = [1,2,3,4,5,6]
    pr_.maxp = 25
    pr_.fixp = True
    pr_.detrend = True
    pr_.power = True
    pr_.metric = 'diag'
    
    
def check_filt(data):
    
    set_def()
    pr_.alg = 'coh'
    
    res = sta_.window_analysis(data[:,:100000])[0]
    
    pl_.plot_coherogram(res)

    return res
  


def check_hist():
    
    si = 100000
    nf = 1024
    
    for i in arange(5):
    
        pp.subplot(2,3,i+1)
        pp.specgram(get_data(t = i)[:si,0], NFFT = nf)
    
    pp.show()
    
def all_cohero():
    
    for i in arange(5):
        pp.figure(i+1)
        check_filt(get_data(i))
        
    pp.show()
    

def all_coh():
       
    set_def()
       
    for i in arange(5):
        pp.figure(i+1)
        an.pdc_full(get_data(i)[:,:30000])
        
    pp.show()
    
def all_psd():
       
    set_def()
       
    for i in arange(5):
        pp.figure(i+1)
        pp.psd(get_data(i)[1,:40000])
        
    pp.show()
    
def final_analysis():
    
    root = "G:\\stein\\dados\\edu_comp\\"
    pr_.output_dir = "G:\\stein\\dados\\edu_comp\\results\\"
    
    algoritmo = 'pdc'
    #algoritmo = 'coh'
    
    window_size = 10
    n_frequencies = 250
    sampling_rate = 500
    
    plot_states = array([1,2,3,4,5,6])
    #plot_states = array([2])
    plot_freq = 150
    plota = False
    pr_.do_states_log = True
    
    ordem_max = 25
    ordem_fixa = True
    detrend = True
    espectro_em_potencia = True
    metrica_pdc = 'diag'
    
    pr_.ar_fit = 'yw'
    
    set_params(alg = algoritmo, window_size = window_size, 
               nf = n_frequencies, sample_f = sampling_rate,
               maxp = ordem_max, fixp = ordem_fixa, detrend = detrend, 
               power = espectro_em_potencia, metric = metrica_pdc, 
               do_plot = plota, plotf = plot_freq, plot_states = plot_states,
               root_dir = root)    
        
    for r in [0,1,2]:
        for d in arange(len(ds[r])):
            
            data = get_data(r = r, d = d)
            
            #data = data[:,:10000]
            
            print 'Data loaded:', 'r=', r, ' d=', d
        
            pr_.plot_labels = ls[r]
            pr_.stinput = 'R%s_D%s' % (rs[r], ds[r][d]) 
                
            state = get_state(r = r, d = d)
            
            pr_.alg = 'pdc'
            res, meds, stds, nstates = sta_.states_analysis(data, state)
            
            return
            
            pr_.alg = 'coh'
            res, meds, stds, nstates = sta_.states_analysis(data, state)
    
    return res, meds, stds, nstates
    
    
if __name__ == '__main__':
    
    #main_analysis()
    
    #check_hist()
    
    #batch_analysis()

    #all_filters()
    
    #all_coh()
    
    #all_psd()
    
    final_analysis()
    