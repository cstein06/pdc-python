'''
Created on May 20, 2010

@author: stein
'''

from numpy import *

from pdc import *

import projects.data_io.plxio as plx

root = "/media/8c8a676c-a8cd-4a18-ae81-0ad35333149b/dados/dados taisa/"

def test_plxio():
    
    filename = root + 'R16_17_03_10_42POLM_ESTEIRA'
    unit = r_[0:5]
    channel = r_[33:65]
    maxspikes = 30000
    
    return plx.readPLXspikes(filename, unit, channel, maxspikes = maxspikes,
                             maxtime = 2000000)


if __name__ == '__main__':
    
    test_plxio()