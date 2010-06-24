'''
Created on May 20, 2010

@author: stein
'''

from numpy import *

def sort_spike_ch(chs, ts, chans = None):
    
    if chans == None:
        chans = r_[33:65]
    
    res = []
    for i in chans:
        aux = where(chs == i)
        auxv = ts[aux]
        
        res.append(auxv)
        
        print i, size(auxv)
        
    return res
    
def spiketotime(spi, maxtime = 2000000, subs = 20):
    
    spi = spi[spi > 0]
    
    #mt = min(spi.max()/20, maxtime)
    mt = maxtime
    res = zeros(mt)
    for i in arange(size(spi)):
        if spi[i]/20 >= mt or spi[i] == 0:
            break
        
        res[spi[i]/20] += 1
        
    return res
