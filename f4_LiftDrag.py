# referece: A. E. Bryson, M. N. Desai and W. C. Hoffman
#  "The Energy-State Approximation in Performance Optimizaiton
#   of Supersonic Aircraft" Journal of Aircraft 6 (1969) 481

from scipy.interpolate import PchipInterpolator

# F4 aerodynamic data
# Lift coefficient C_L = cLa * alpha
# Drag coefficient C_D = cD0 + kappa*cLa*alpha**2
# alpha = angle of attack / radian
M = [0, 0.4, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8]# Mach number
cLa = [3.44,  3.44,  3.44,  3.58,  4.44,  3.44,  3.01,  2.86,  2.44]
cD0 = [0.013, 0.013, 0.013, 0.014, 0.031, 0.041, 0.039, 0.036, 0.035]
kappa = [0.54,  0.54,  0.54,  0.75,  0.79,  0.78,  0.89,  0.93,  0.93]

LiftDragCoeffs = PchipInterpolator(M, [cLa,cD0,kappa], axis=-1)
