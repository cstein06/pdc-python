
from numpy import *
from numpy.random import randn

from data_simulation.ar_data import ar_data

from scipy import weave
from scipy.weave import converters

#===============================================================================
#   4         """Takes a time step using inlined C code -- this version uses
#   5         blitz arrays."""
#   6         g = self.grid
#   7         nx, ny = g.u.shape
#   8         dx2, dy2 = g.dx**2, g.dy**2
#   9         dnr_inv = 0.5/(dx2 + dy2)
#  10         u = g.u
#  11 
#  12         code = """
#  13                #line 120 "laplace.py" (This is only useful for debugging)
#  14                double tmp, err, diff;
#  15                err = 0.0;
#  16                for (int i=1; i<nx-1; ++i) {
#  17                    for (int j=1; j<ny-1; ++j) {
#  18                        tmp = u(i,j);
#  19                        u(i,j) = ((u(i-1,j) + u(i+1,j))*dy2 +
#  20                                  (u(i,j-1) + u(i,j+1))*dx2)*dnr_inv;
#  21                        diff = u(i,j) - tmp;
#  22                        err += diff*diff;
#  23                    }
#  24                }
#  25                return_val = sqrt(err);
#  26                """
#  27         # compiler keyword only needed on windows with MSVC installed
#  28         err = weave.inline(code,
#  29                            ['u', 'dx2', 'dy2', 'dnr_inv', 'nx', 'ny'],
#  30                            type_converters=converters.blitz,
#  31                            compiler = 'gcc')
#  32         return err
#===============================================================================


def test_ar_data():
    
    A = array([[2,3],[0,5]])/10.0
    A.resize([2,2,1])
    data = ar_data(A)
    print data.shape
    print data[:,-10:]

def test_weave():
    
    a = arange(10).reshape(2,5)
    code3 = '''
            double b = 0.0;
           for (int i = 0; i < 5; i++) {
                b = b + a(1,i);
               // b = b + 1;
               a(1,i) = 7;
           }
           return_val = b; '''
    code = '''
        return_val = 8; '''
    code2 = '''
        return_val = 9;'''
    print code
    ret = weave.inline(code3, ['a'], #verbose = 1,
                       type_converters=converters.blitz,
                       compiler = 'gcc')
    #weave.inline('printf("%d\\n",a);',['a'])

    print a
    print ret

if __name__ == '__main__':
    test_ar_data()
    
    