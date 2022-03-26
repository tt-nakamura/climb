import numpy as np
import matplotlib.pyplot as plt
from f4_LiftDrag import LiftDragCoeffs, M, cLa, cD0, kappa

plt.figure(figsize=(5,6))

plt.subplot(3,1,1); plt.plot(M, cLa, 'x')
plt.subplot(3,1,2); plt.plot(M, cD0, 'x')
plt.subplot(3,1,3); plt.plot(M, kappa, 'x')

M = np.linspace(0, 1.8, 65)
cLa, cD0, kappa = LiftDragCoeffs(M)

plt.subplot(3,1,1)
plt.plot(M, cLa)
plt.ylabel(r'$C_{\rm L\alpha}$')

plt.subplot(3,1,2)
plt.plot(M, cD0)
plt.ylabel(r'$C_{\rm D0}$')

plt.subplot(3,1,3)
plt.plot(M, kappa)
plt.ylabel(r'$\kappa$')
plt.xlabel(r'$M$ = Mach number')

plt.tight_layout()
plt.savefig('fig1.eps')
plt.show()

