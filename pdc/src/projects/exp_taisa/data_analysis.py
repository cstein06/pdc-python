'''
Created on May 5, 2010

@author: Marcelo Gomes Mattar
'''

from numpy import loadtxt
import scipy.stats as st
import pdc.analysis as an
from numpy.matlib import randn
import matplotlib.pyplot as pp
import cProfile
import load_plexon

def import_data(filename):
    filepath = './' + filename
    data = loadtxt(filepath)
    #print(data)
    output = load_plexon.readPLXspikes('D:/Documents/Mestrado/dados_taisa/R01_15_09_09_7POLM_ESTEIRA.plx')
    channels = output[1]
    spiketimes = output[3]

if __name__ == '__main__':
   
    import_data('data.txt')