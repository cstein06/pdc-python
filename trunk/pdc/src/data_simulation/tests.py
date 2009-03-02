
from numpy import *
from numpy.random import randn

from data_simulation.ar_data import ar_data

def test_ar_data():
    
    A = array([[2,3],[0,5]])/10.0
    A.resize([2,2,1])
    data = ar_data(A)
    print data.shape
    print data[:,-10:]

if __name__ == '__main__':
    test_ar_data()