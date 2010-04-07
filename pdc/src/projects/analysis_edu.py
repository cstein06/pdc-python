# -*- coding:utf-8 -*-
"""
Created on Apr 7, 2010

@author: Stein
"""



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