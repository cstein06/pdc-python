'''
Created on 18/03/2010

@author: Marcelo
'''

from numpy import *
import scipy.stats as st
import pdc.analysis as an
from numpy.matlib import randn
import matplotlib.pyplot as pp
import cProfile

def fp_pdc(numint, timewindow, spikerate_A, spikerate_B, alpha_pdc):
    fp = 0
    for i in range(0,numint):
        print i
        neuronA = st.bernoulli.rvs(spikerate_A,size=timewindow)
        neuronB = st.bernoulli.rvs(spikerate_B,size=timewindow)
        res = an.pdc_full([neuronA,neuronB],alpha=alpha_pdc)
        
        mes = res[0]
        th = res[1]
    #    ic1 = res[2]
    #    ic2 = res[3]
        
        for freq in range(0,mes.shape[2]):
            if mes[0,1,freq] > th[0,1,freq] :
                fp = fp + 1
                break
        
    return fp



def fp_gct(numint, timewindow, spikerate_A, spikerate_B):
    fp = 0
    for i in range(0,numint):
        print i
        neuronA = st.bernoulli.rvs(spikerate_A,size=timewindow)
        neuronB = st.bernoulli.rvs(spikerate_B,size=timewindow)
        res = an.gct([neuronA,neuronB])
        
        mes = res[0]
        
        if (mes[0,1] < 0.05) | (mes[1,0] < 0.05):
            fp = fp + 1   
    
    print fp
    print 'O Granger Causality Test apresenta ', fp*100/numint, '% de falsos positivos'     
    return fp




#def fp_gct_convolved(numint, timewindow, spikerate_A, spikerate_B):
#    fp = 0
#    t = arange(-0.25,0.25,0.005)
#    h = (sin(pi*t/0.05)

def fp_gct_convolved(numint, timewindow, timestep, spikerate_A, spikerate_B, threshold):
    
    #numint = number of iterations
    #timewindow = lenght of the data (in seconds)
    #timestep = time of each step, and duration of a spike (in seconds)
    #spikerate_A = spike rate of the first neuron (in Hertz)
    #spikerate_B = spike rate of the second neuron (in Hertz)
    #threshold (p value) used to indicate connectivity in the Granger Causality Test
    
    
    fp = 0
    
    #parameters of the kernel
    t = arange(-0.075,0.075,timestep)
    h = (sin(pi*t/0.05))/(pi*t/0.05)
                          
                                                    
    for i in range(0,numint):
        print i
        neuronA = st.bernoulli.rvs(spikerate_A*timestep,size=timewindow/timestep)
        neuronB = st.bernoulli.rvs(spikerate_B*timestep,size=timewindow/timestep)
        convolvedA = convolve(neuronA, h)
        convolvedB = convolve(neuronB, h)
#        pp.plot(convolvedA)
#        pp.show()
#        return
        
#        res = an.gct([convolvedA,convolvedB], maxp = 10, fixp = True, test_allp = False, v = False)
        res = an.gct([neuronA,neuronB], maxp = 5, fixp = False, test_allp = False, v = False)
#        res = an.gct(randn(2,timewindow), maxp = 10, fixp = True, v = False)
        
        mes = res[0]

        if (mes[0,1] < threshold) | (mes[1,0] < threshold):
            fp = fp + 1   
    
#    print fp
    print 'O Granger Causality Test apresenta ', fp*100/numint, '% de falsos positivos'     
    return fp



if __name__ == '__main__':
   
#    cProfile.run('fp_gct_convolved(10,24000,0.005,0.005)')
#    fp_gct(1000, 24000, 0.005, 0.005)
    fp_gct_convolved(1000,120,0.005,10,10,0.01)