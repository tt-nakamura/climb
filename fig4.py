import numpy as np
import matplotlib.pyplot as plt
from os.path import exists

if not exists('f4_climb.npy'):
    exec(open('f4_climb.py').read())

t,z,v,gamma,alpha,x,m = np.load('f4_climb.npy')

plt.figure(figsize=(5,10))

plt.subplot(5,1,1)
plt.plot(t,v)
plt.ylabel(r'$v$ = velocity / m/s')

plt.subplot(5,1,2)
plt.plot(t,z/1000)
plt.ylabel(r'$z$ = altitude / km')

plt.subplot(5,1,3)
plt.plot(t,alpha)
plt.ylabel(r'$\alpha$ = attack angle / rad')

plt.subplot(5,1,4)
plt.plot(t,gamma)
plt.ylabel(r'$\gamma$ = path angle / rad')

plt.subplot(5,1,5)
plt.plot(t,m/1000)
plt.ylabel(r'$m$ = mass / $10^3$kg')
plt.xlabel(r'$t$ = time / sec')

plt.tight_layout()
plt.savefig('fig4.eps')
plt.show()
