# -*- coding:utf-8 -*-
"""
Created on 27/10/2009

@author: Carlos Stein
"""

__all__ = ['subsample']

from numpy import *

def subsample(data, k, ave=True):
    if data.ndim == 1:
        if ave:
            return convolve(data,ones(k)/float(k),'same')[k-1::k]
        else:
            data[::k]
    
    res = empty((data.shape[0], data.shape[1]/k))
    for i in arange(data.shape[0]):
        if ave:
            res[i] = convolve(data[i],ones(k)/float(k),'same')[k-1::k]
        else:
            res[i] = data[i,::k]
            
    return res
