'''
Created on Jul 8, 2010

@author: stein
'''
import cPickle
import matplotlib.pyplot as pp

from pdc import *
from numpy.lib.io import loadtxt

root = '/media/8c8a676c-a8cd-4a18-ae81-0ad35333149b/dados/mestre/'

def cogram():
    ind = [194,254]
    ins = [163,223]
    
    f = open(root + 'coh.pic')
    co = cPickle.load(f)[194:254]
    f.close()
    
    f = open(root + 'pdc.pic')
    pd = cPickle.load(f)[194:254]
    f.close()
    
    #f = open(root + 'ES60_19_07_09_est_10s_limpo.txt')
    est = loadtxt(root + 'ES60_19_07_09_est_10s_limpo.txt')[163:223]
    f.close()
    
    pr_.window_size = 10
    pr_.nf = 250
    pr_.sample_f = 500
    pr_.plotf = 50
    pr_.plota = True
    pr_.plot_ic = True
    pr_.ss = True
    pr_.logss = True
    pr_.plot_labels = ['Ca1e', 'Ca3e', 'Ca1d', 'Ca3d']
    pr_.alg = 'coh'
    
    pp.figure()
    plot_coherogram(co, est)
    
    pr_.alg = 'pdc'
    pp.figure()
    plot_coherogram(pd, est)
    
    pp.show()
    
if __name__ == '__main__':
    cogram()