from scipy.io import read_array
from numpy import *

def read_hor_eeg(file_s):
    fil = open(file_s, 'r')
    data = array([])
    header = []
    for line in fil:
        header.append(line[:4].strip())
        d = fromstring(line[4:], sep=' ')
        data = concatenate((data, d))
    return data.reshape(32,-1)

if __name__ == "__main__":
    
    s = 'D:\work\dados\\fmri\einstein\CarlosVIS\Carlos_VIS_MRI Cardiobalistogram Correction_OK.dat'.replace('\\', '/')
    print read_hor_eeg(s).shape