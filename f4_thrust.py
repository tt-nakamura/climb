# referece: A. E. Bryson, M. N. Desai and W. C. Hoffman
#  "The Energy-State Approximation in Performance Optimizaiton
#   of Supersonic Aircraft" Journal of Aircraft 6 (1969) 481

import numpy as np
from scipy.interpolate import SmoothBivariateSpline
from scipy.constants import foot,pound,g

# F4 jet engine data
M = [0.0, # Mach number
     0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
     0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
     0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
     0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
     1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2,
          1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4,
                    1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6,
                              1.8, 1.8, 1.8, 1.8, 1.8]
z = [0, # altitude / kft
     0, 5, 10, 15, 20, 25, 30,
     0, 5, 10, 15, 20, 25, 30, 40, 50,
     0, 5, 10, 15, 20, 25, 30, 40, 50,
     0, 5, 10, 15, 20, 25, 30, 40, 50, 70,
     0, 5, 10, 15, 20, 25, 30, 40, 50, 70,
     0, 5, 10, 15, 20, 25, 30, 40, 50, 70,
        5, 10, 15, 20, 25, 30, 40, 50, 70,
               15, 20, 25, 30, 40, 50, 70,
                       25, 30, 40, 50, 70]
T = [24.2, # thrust / klb
     28.0, 24.6, 21.1, 18.1, 15.2, 12.8, 10.7,
     28.3, 25.2, 21.9, 18.7, 15.9, 13.4, 11.2, 7.3,  4.4,
     30.8, 27.2, 23.8, 20.5, 17.3, 14.7, 12.3, 8.1,  4.9,
     34.5, 30.3, 26.6, 23.2, 19.8, 16.8, 14.1, 9.4,  5.6,  1.1,
     37.9, 34.3, 30.4, 26.8, 23.3, 19.8, 16.8, 11.2, 6.8,  1.4,
     36.1, 38.0, 34.9, 31.3, 27.3, 23.6, 20.1, 13.4, 8.3,  1.7,
           36.6, 38.5, 36.1, 31.6, 28.1, 24.2, 16.2, 10.0, 2.2,
                       38.7, 35.7, 32.0, 28.1, 19.3, 11.9, 2.9,
                                   34.6, 31.1, 21.7, 13.3, 3.1]
z = np.asarray(z)*1000*foot # altitude / m
T = np.asarray(T)*1000*pound*g # thrust / newton
T_interp = SmoothBivariateSpline(M, z, T, kx=4, ky=4)

def thrust(M, z, dx=0, dy=0, grid=False):
    """ spline interpolation of F4 propulsion data
    M = Mach number
    z = altitude / m
    dx = order of derivative with respect to M
    dy = order of derivative with respect to z
    grid = whether to output meshed grid in (M,z)
    return maximun thrust / newton
    """
    return T_interp(M, z, dx, dy, grid)
