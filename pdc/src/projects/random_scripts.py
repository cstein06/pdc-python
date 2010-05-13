from numpy import *
import matplotlib.pyplot as pp

from pdc.sim_data import ar_data
import pdc.analysis as pdc_
import cPickle

def test_subsample():
    
    A = array([[0.6, 0.2], [-0.2, 0.5]])
    er = array([[2,0],[0,1]])
    m = 300

    data = ar_data(A, er, m, 40)
    
    sdata = data[0::3]
    
    pdc_.pdc_and_plot(data, 2, 30, metric='info')
    pdc_.pdc_and_plot(sdata, 2, 30, metric='info')
    
def test_plot_rc():
    pp.rcParams['ps.usedistiller'] = 'xpdf'
    pp.plot(arange(10), 'k+')
    pp.title('test ps')
    pp.show()
    
    

root = '/home/stein/producao/quali/simulados/'
def test_pickle():
    a = ones(10)
    b = 5
    d = ones([3,3])
    e = [3, arange(2)]
    
    a.dump(root +'test.pic')
    d.dump(root +'test.pic')
    f = open(root +'test.pic', 'ab')
    cPickle.dump(b, f)
    cPickle.dump(d, f)
    #f.close()
    d.dump(f)
    f.close()
    f = open(root +'test.pic', 'ab')
    cPickle.dump(e, f)
    f.close()
    
    f = open(root +'test.pic', 'rb')
    print cPickle.load(f)
    print cPickle.load(f)
    print cPickle.load(f)
    print cPickle.load(f)
    print cPickle.load(f)
    print cPickle.load(f)
    
if __name__ == "__main__":
    
    #test_subsample()
    
    #test_plot_rc()
    
    test_pickle()