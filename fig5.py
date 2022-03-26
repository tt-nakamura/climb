import numpy as np
import matplotlib.pyplot as plt
from os.path import exists

if not exists('f4_climb.npy'):
    exec(open('f4_climb.py').read())

_,z,_,_,_,x,_ = np.load('f4_climb.npy')

plt.plot(x/1000, z/1000)
plt.xlabel(r'$x$ = horizontal distance / km', fontsize=14)
plt.ylabel(r'$z$ = altitude / km', fontsize=14)
plt.tight_layout()
plt.savefig('fig5.eps')
plt.show()
