##
# @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
# @Id $Id: classBlaser.py 135 2016-09-12 11:55:18Z mmbosschaert $
##

import numpy as np
from pydelay import dde23
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import colors

# number of time units
tfinal = 3000

# define DDE
eqns = {
	'x' : 'omega-x-x*(y+gamma*y(t-tau))',
  'y' : 'k*(x-1)*y'
}

# set parameters
params = {
        'omega':1.5,
        'k':50,
        'gamma':0.763036858,
        'tau':1.60017288
}

# solve DDE
dde = dde23(eqns=eqns, params=params)
dde.set_sim_params(tfinal=tfinal, dtmax=0.1, AbsTol=1e-8, RelTol=1e-6)
dde.set_sim_params(tfinal=tfinal)
histfunc = {'x': lambda t: 1.1, 'y': lambda t: 0.4} 
dde.hist_from_funcs(histfunc, 51)
dde.run()

# subtract solution components for time series
tau=1
sol0 = dde.sample(tfinal-tfinal/3, tfinal, 0.01)
sol1 = dde.sample(tfinal-tfinal/3-tau, tfinal-tau, 0.01)
t = sol0['t']
x  = sol0['x']
y  = sol0['y']
z  = sol1['x']

# plot time series
fig = plt.figure()
plt.figure(1)
plt.subplot(121)
plt.xlabel('$t$')
plt.ylabel('$x(t)$')
plt.plot(t, x, c='royalblue')

plt.subplot(122)
plt.xlabel('$t$')
plt.ylabel('$y(t)$')
plt.plot(t, y, c='royalblue')
fig.set_size_inches(10, 4)
plt.show()

# subtract solution components for phase-space plot
tau=1
sol0 = dde.sample(tfinal-tfinal/6, tfinal, 0.01)
sol1 = dde.sample(tfinal-tfinal/6-tau, tfinal-tau, 0.01)
x  = sol0['x']
y  = sol0['y']
z  = sol1['x']

# plot the solution in phase-space
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x, y, z, c='royalblue')

ax.set_xlabel('$x(t)$')
ax.set_ylabel('$y(t)$')
ax.set_zlabel('$x(t-\\tau)$')

ax.view_init(22, 46)
fig.set_size_inches(18.5, 10.5)
plt.show()
