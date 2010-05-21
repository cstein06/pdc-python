'''
Created on May 20, 2010

@author: stein
'''

from numpy import *

def sort_spike_ch(chs, ts, chans = None):
    
    if chans == None:
        chans = r_[33,65]
    
    res = []
    for i in chans:
        aux = where(chs == i)
        auxv = ts[aux]
        
        #res append(auxv)
        
    return res
    
    