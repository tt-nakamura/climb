import numpy as np
import matplotlib.pyplot as plt
from atmosphere import density, sound_speed

z = np.linspace(0, 30000, 61)
rho = density(z)
a = sound_speed(z)

plt.figure(figsize=(5,5))
plt.subplot(2,1,1)
plt.plot(z/1e3, rho)
plt.ylabel(r'$\rho$ = air density / kg/m$^3$')
plt.subplot(2,1,2)
plt.plot(z/1e3, a)
plt.ylabel(r'$a$ = sound speed / m/s')
plt.xlabel(r'$z$ = altitude / km')

plt.tight_layout()
plt.savefig('fig2.eps')
plt.show()
