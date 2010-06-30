'''
Created on Jun 29, 2010

@author: stein
'''

from pdc import *
from numpy import *

def model(b = 0):
    A = zeros([3,3,1])
    A[0,0,0] = 0.5
    A[0,1,0] = 1.0
    A[2,1,0] = b
    er = identity(3)
    return A, er

#def pdttestar():
    
