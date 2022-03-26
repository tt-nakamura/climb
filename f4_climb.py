# F4 min time to climb using inverse optimization
# reference: A. E. Bryson
#  "Dynamic Optimization" problem 9.3.15

import numpy as np
from scipy.integrate import cumtrapz
from scipy.optimize import fmin_slsqp
from scipy.constants import foot,pound,g,degree
from atmosphere import density, sound_speed
from f4_thrust import thrust
from f4_LiftDrag import LiftDragCoeffs

# Aircraft parameters (F4 aircraft)
m0 = 41998*pound # initial weight / kg
S = 530*foot**2 # wing area / m^2
epsilon = 3*degree # angle of zero lift / radian
c = 1600 # fuel comsumption time / sec
rho0 = density(0) # air density at ground / kg/m^3

l_unit = 2*m0/rho0/S # length unit / m
v_unit = np.sqrt(g*l_unit) # velocity unit / m/s
t_unit = l_unit/v_unit # time unit / sec
f_unit = m0*g # force unit / newton

def f4_climb(init, zf, z0=0, N=128, verbose=2, maxiter=200):
    """ F4 min time to climb using inverse optimization
    z0,zf = initial and final altitude / m (specified)
    init = either (v0,vf) or (v,gamma) or
               (tf,v0,vf) or (tf,v,gamma) where
      tf = initial guess of min time / sec
      v0,vf = initial and final velocity / m/s (specified)
      v = velocity history (shape(N,)) / m/s
      gamma = path angle history (shape(N,)) / radian
    N = number of mesh points
    verbose = passed to fmin_slsqp as iprint
    maxiter = passed to fmin_slsqp as iter
    if (v,gamma) are given, N is ignored
    if tf is not given, tf is set to 2*(zf-z0)/(v0+vf)
    return t,z,v,gamma,alpha,x,m (each shape(N,)) where
      t = time / sec, z = altitude / m,
      v = velocity / m/s
      gamma = path angle / radian
      alpha = angle of attack / radian,
      x = horizontal distance / m,
      m = mass of aircraft / kg
      for optimized trajectory such that
      t[0],t[-1] = 0,tf; z[0],z[-1] = z0,zf
      v[0],v[-1] = v0,vf; gamma[0],gamma[-1] = 0,0;
      x[0] = 0; m[0] = m0
    reference: A. E. Bryson
     "Dynamic Optimization" problem 9.3.15
    """
    if len(init)==3: tf,a,b = init
    elif len(init)==2:  a,b = init
    else: raise RuntimeError('bad init')

    if np.isscalar(a):
        v0,vf = a,b
        v = np.linspace(v0,vf,N)
        gamma = np.zeros(N)
    else:
        v,gamma = a,b
        v0,vf = v[0],v[-1]
        N = len(v)

    if len(init)==2: tf = 2*(zf-z0)/(v0+vf)
    N1 = N-1

    def f_eqcons(p):# equality constraints
        dt = p[0]*t_unit/N1
        v = np.r_[v0, p[1:N1]*v_unit, vf]
        gamma = np.r_[0, p[N1:], 0]
        z = cumtrapz(v*np.sin(gamma), dx=dt, initial=0) + z0
        con0 = z[-1] - zf # final altitude specified

        vd = np.diff(v)/dt
        gd = np.diff(gamma)/dt
        z = (z[:-1] + z[1:])/2 # mid points
        v = (v[:-1] + v[1:])/2
        gamma = (gamma[:-1] + gamma[1:])/2

        M = v/sound_speed(z) # Mach number
        T = thrust(M,z)
        m = m0 - np.r_[0, np.cumsum(T/c/g)*dt] # mass
        m = (m[:-1] + m[1:])/2
        cLa, cD0, kappa = LiftDragCoeffs(M)
        rho = density(z)
        e = rho*v**2*S/2
        # alpha = angle of attack, determined from
        # equation of motion perpendicular to
        # velocity (small angle approximation)
        alpha = m*(v*gd + g*np.cos(gamma)) - T*epsilon
        alpha /= T + cLa*e
        D = (cD0 + kappa*cLa*alpha**2)*e # Drag
        # equation of motion parallel to velocity
        con1 = m*(vd + g*np.sin(gamma)) + D
        con1 -= T*np.cos(alpha + epsilon)
        return np.r_[con0/l_unit, con1/f_unit]

    def f_ieqcons(p):# inequality constraints
        dt = p[0]*t_unit/N1
        v = np.r_[v0, p[1:N1]*v_unit, vf]
        gamma = np.r_[0, p[N1:], 0]
        z = cumtrapz(v*np.sin(gamma), dx=dt) + z0
        return z[:-1]/l_unit # altitude >= 0

    p0 = np.r_[tf/t_unit, v[1:-1]/v_unit, gamma[1:-1]]
    p = fmin_slsqp(lambda p: p[0], p0,
                   f_eqcons=f_eqcons,
                   f_ieqcons=f_ieqcons,
                   iprint=verbose, iter=maxiter)

    t = np.linspace(0, p[0]*t_unit, N)
    v = np.r_[v0, p[1:N1]*v_unit, vf]
    gamma = np.r_[0, p[N1:], 0]
    cg,sg = np.cos(gamma), np.sin(gamma)
    z = cumtrapz(v*sg, t, initial=0) + z0
    x = cumtrapz(v*cg, t, initial=0)
    M = v/sound_speed(z)
    cLa,_,_ = LiftDragCoeffs(M)
    T = thrust(M,z)
    m = m0 - cumtrapz(T/c/g, t, initial=0)
    rho = density(z)
    e = rho*v**2*S/2
    gd = np.gradient(gamma, t)
    # angle of attack (small angle approximation)
    alpha = (m*(v*gd + g*cg) - T*epsilon)/(T + cLa*e)
    return t,z,v,gamma,alpha,x,m


def f4_edot(v,z,m):
    """ energy-state approximation for dE/dt
    v = speed / m/s, z = altitude / km, m = mass (const)
    return dE/dt = (d/dt)((1/2)mv^2 + mgz)
    reference: A. E. Bryson
     "Dynamic Optimization" example 4.5.4
    """
    M = v/sound_speed(z)
    rho = density(z)
    T = thrust(M,z)
    cLa, cD0, kappa = LiftDragCoeffs(M)
    e = rho*v**2*S/2
    alpha = (m*g - T*epsilon)/(T + cLa*e)
    D = (cD0 + kappa*cLa*alpha**2)*e
    return (T-D)*v


if __name__ == '__main__':
    zf = 65673*foot # final altitude / m
    v0 = 440*foot # initial velocity / m/s
    vf = 968*foot # final velocity / m/s

    t,z,v,gamma,alpha,x,m = f4_climb((v0,vf),zf)
    np.save('f4_climb.npy', (t,z,v,gamma,alpha,x,m))
    print('minimun time to climb: %f sec' % t[-1])
