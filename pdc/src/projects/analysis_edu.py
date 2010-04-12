# -*- coding:utf-8 -*-
"""
Created on Apr 7, 2010

@author: Stein
"""

from numpy import *
import matplotlib.pyplot as pp

from projects.edus import get_data, get_res, get_state
import pdc.plotting as pl_
import pdc.states as sta_
import pdc.analysis as an_

from pdc.globals import *

from projects.edus import ds, rs, ls, chs

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

def check_specgram():
    
    si = 100000
    nf = 1024
    
    for i in arange(5):
    
        pp.subplot(2,3,i+1)
        pp.specgram(get_data(t = i)[:si,0], NFFT = nf, origin='lower')
    
    pp.show()
    
def all_cohero():
    
    rootim = 'G:\\stein\\dados\\edu_comp\\results\\images\\'
    fileim = 'R%s_D%s_%s.png'
    
    alg = 'coh'
    
    set_def()
    
    #pr_.state_colors = ['k', 'lightblue', 'darkblue', 
    #                    'pink', 'red', 'green']
    pr_.state_colors = ['lightblue', 'lightblue', 'lightblue', 
                        'lightblue', 'red', 'green']
    k = 1
    
    for i in [0,1,2]:
        for j in arange(len(ds[i])):
            pr_.plot_labels = ls[i]
            fig = pp.figure(k)
            k += 1
            sta = get_state(r = i, d = j)
            res = get_res(r = i, d = j, alg = alg)[0]
            pl_.plot_coherogram(res, sta)
            #check_filt(get_data(i))
            
            fig.savefig(rootim + fileim % (rs[i], ds[i][j], alg))
                        
            pp.close()
        
    pp.show()
    

def all_coh():
       
    set_def()
    
    k = 1
    for i in arange(len(rs)):
        for j in arange(len(ds[i])):
            pp.figure(k)
            k += 1
            #an_.pdc_full(get_data(i)[:,:30000])
        
    pp.show()
    
def all_psd():
       
    set_def()
       
    for i in arange(5):
        pp.figure(i+1)
        pp.psd(get_data(i)[1,:40000])
        
    pp.show()
    
    
if __name__ == '__main__':
    
    pass
    #all_cohero()
    check_specgram()



