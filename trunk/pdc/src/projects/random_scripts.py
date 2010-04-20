from numpy import *
import matplotlib.pyplot as pp

from pdc.sim_data import ar_data
import pdc.analysis as pdc_

def test_subsample():
    
    A = array([[0.6, 0.2], [-0.2, 0.5]])
    er = array([[2,0],[0,1]])
    m = 300

    data = ar_data(A, er, m, 40)
    
    sdata = data[0::3]
    
    pdc_.pdc_and_plot(data, 2, 30, metric='gen')
    pdc_.pdc_and_plot(sdata, 2, 30, metric='gen')
    
def test_plot_rc():
    pp.rcParams['ps.usedistiller'] = 'xpdf'
    pp.plot(arange(10), 'k+')
    pp.title('test ps')
    pp.show()
    
    
if __name__ == "__main__":
    
    #test_subsample()
    
    test_plot_rc()