import numpy as np
import matplotlib.pyplot as plt
from f4_thrust import thrust

M = np.linspace(0, 1.8, 37)
z = np.linspace(0, 20000, 41)
T = thrust(M, z, grid=True)
T[M[:,np.newaxis] - z/5e3 > 1] = np.nan
levels = np.linspace(10, 160, 16)
cs = plt.contour(M, z/1e3, T.T/1e3, levels=levels)
plt.clabel(cs, fmt='%d kN')
plt.xlabel(r'$M$ = Mach number', fontsize=14)
plt.ylabel(r'$z$ = altitude / km', fontsize=14)
plt.tight_layout()
plt.savefig('fig3.eps')
plt.show()
