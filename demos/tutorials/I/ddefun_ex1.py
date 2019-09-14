##
# @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
# @Id $Id: ddefun_ex1.py 135 2016-09-12 11:55:18Z mmbosschaert $
##

import numpy as np
from pydelay import dde23
import matplotlib.pyplot as plt
from matplotlib import colors

# number of time units
tfinal = 4

# define DDE
eqns = {
	'x' : 'x(t-tau)*(x-1)'
}

# set parameters
params = {
        'tau':1
}

# solve DDE
dde = dde23(eqns=eqns, params=params)
dde.set_sim_params(tfinal=tfinal)
histfunc = {'x': lambda t: np.cos(t)} 
dde.hist_from_funcs(histfunc, 51)
dde.run()

# subtract solution components for time series
sol = dde.sample(0, tfinal, 0.01)
t = sol['t']
x1 = sol['x']

# solve with second history function
histfunc = {'x': lambda t: np.exp(t)} 
dde.hist_from_funcs(histfunc, 51)
dde.run()
sol = dde.sample(0, tfinal, 0.01)
x2 = sol['x']

# solve with thrid history function
histfunc = {'x': lambda t: 1-t} 
dde.hist_from_funcs(histfunc, 51)
dde.run()
sol = dde.sample(0, tfinal, 0.01)
x3 = sol['x']

# plot time series
fig = plt.figure()
plt.figure(1)
plt.plot(t, x1, c='royalblue')
plt.plot(t, x2, ls=':', c='red')
plt.plot(t, x3, ls='--', c='gold')

# add history functions to the plot
t=np.linspace(-1,0, num=51)
plt.plot(t, np.cos(t), c='royalblue',label='$\cos(t)$')
plt.plot(t, np.exp(t), ls=':', c='red', label='$\exp(t)$')
plt.plot(t, 1-t, ls='--', c='gold', label='$1-t$')

plt.legend();
plt.xlabel('$t$')
plt.ylabel('$x$')

fig.set_size_inches(10, 4)
plt.show()

