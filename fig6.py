import numpy as np
import matplotlib.pyplot as plt
from os.path import exists

if not exists('f4_climb.npy'):
    exec(open('f4_climb.py').read())

_,z,v,_,_,_,_ = np.load('f4_climb.npy')

plt.plot(v, z/1000)
plt.xlabel(r'$v$ = velocity / m/s', fontsize=14)
plt.ylabel(r'$z$ = altitude / km', fontsize=14)
plt.tight_layout()
plt.savefig('fig6.eps')
plt.show()
