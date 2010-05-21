# -*- coding:utf-8 -*-
"""
Created on 18/05/2010

Simulações para qualificação

Plots simples do modelo

@author: Carlos Stein
"""

from pdc import *

import projects.quali.test_asymp_distr_alpha as qua

from numpy import *

import matplotlib.pyplot as pp

def plots(alg):
    
    pr_.alg = alg
    
    pr_.detrend = False
    
    pr_.plot_ic = True
    
    pr_.alpha = 0.05
    
    A, e = qua.model1()
    
    data = ar_data(A, e, nd = 10000) 
    
    measure_full(data)

def plot_simples():

    plots('pdc')
    
    pp.show()
    
def pergunta_daniel():
    
    plots('dtf')
    
    pp.show()
    

if __name__ == "__main__":
    
    #plot_simples()
    
    pergunta_daniel()
    
    pass
