import numpy as np
import matplotlib.pyplot as plt
from f4_climb import f4_edot, l_unit, f_unit, t_unit
from scipy.constants import pound,g
from os.path import exists

E_unit = f_unit*l_unit
Edot_unit = E_unit/t_unit

m = 37400*pound
v = np.linspace(130, 550, 50)
z = np.linspace(0, 20020, 50)
v,z = np.meshgrid(v,z)
E = m*(v**2/2 + g*z)/E_unit
Edot = f4_edot(v,z,m)/Edot_unit

plt.contour(v, z/1000, Edot,
            levels=np.linspace(-0.5, 2.5, 40))
plt.contour(v, z/1000, E, colors='k',
            linestyles='dotted',
            levels=np.linspace(1, 50, 20))

if not exists('f4_climb.npy'):
    exec(open('f4_climb.py').read())

_,z,v,_,_,_,_ = np.load('f4_climb.npy')

plt.plot(v, z/1000, 'r')
plt.xlabel(r'$v$ = velocity / m/s', fontsize=14)
plt.ylabel(r'$z$ = altitude / km', fontsize=14)
plt.tight_layout()
plt.savefig('fig7.eps')
plt.show()
