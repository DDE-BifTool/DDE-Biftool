# Simulation of the neural DDE analyzed in
# Visser, S. (2013). From spiking neurons to brain waves.
#
# @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
# @Id $Id: neural.py 135 2016-09-12 11:55:18Z mmbosschaert $
#
import matplotlib.pyplot as plt
from matplotlib import colors
from pydelay import dde23
import time

# number of time units
t = time.time()
tfinal = 1000000

# define DDE
eqns = {
	'x1' : '-x1-a*(tanh(b*x1(t-tau1)-1)+tanh(1))*pow(cosh(1),2)'
		'+c*(tanh(d*x2(t-tau2)-1)+tanh(1))*pow(cosh(1),2)',
	'x2' : '-x2-a*(tanh(b*x2(t-tau1)-1)+tanh(1))*pow(cosh(1),2)'
		'+c*(tanh(d*x1(t-tau2)-1)+tanh(1))*pow(cosh(1),2)'
}

# set parameters
# point with stable trivial equalibrium
params1 = {
        'a':0.279899956451376,
        'b':2.0,
        'c':0.573986945974371,
        'd':1.2,
        'tau1':12.99,
        'tau2':20.15
}

# point where Sebastiaan found a torus
params2 = {
	'a':0.5598/2,
	'b':2.0,
	'c':0.6890/1.2,
	'd':1.2,
	'tau1':12.99,
	'tau2':20.15
}

# solve DDE
dde = dde23(eqns=eqns, params=params2)
dde.set_sim_params(tfinal=tfinal, dtmax=0.1, AbsTol=1e-8, RelTol=1e-6)
histfunc = {'x1': lambda t: 0.05, 'x2': lambda t: 0.075 } 
dde.hist_from_funcs(histfunc, 51)
dde.run()
elapsed = time.time() - t
elapsed

# subtract solution components
tfinal=30000
sol1 = dde.sample(tfinal-2500, tfinal, 0.1)
sol2 = dde.sample(tfinal-2500, tfinal, 0.1)
t = sol1['t']
x1 = sol1['x1']
x2 = sol2['x2']

# plot the solution in phase-space
plt.plot(x1, x2,c='royalblue')
plt.xlabel('$x_1(t)$')
plt.ylabel('$x_2(t)$')
filename='/home/maikel/Dropbox/Documents/Uni/Master/2015/' + \
  'Reading course DDE 2/VI/images/tutorial_VI_torus_30000.png'
plt.savefig(filename, bbox_inches='tight')

# subtract solution components
tfinal=300000
sol1 = dde.sample(tfinal-2500, tfinal, 0.1)
sol2 = dde.sample(tfinal-2500, tfinal, 0.1)
t = sol1['t']
x1 = sol1['x1']
x2 = sol2['x2']

# plot the solution in phase-space
plt.plot(x1, x2, c='forestgreen')
plt.xlabel('$x_1(t)$')
plt.ylabel('$x_2(t)$')
filename='/home/maikel/Dropbox/Documents/Uni/Master/2015/' + \
  'Reading course DDE 2/VI/images/tutorial_VI_torus_300000.png'
plt.savefig(filename, bbox_inches='tight')

# subtract solution components
tfinal=1000000
sol1 = dde.sample(tfinal-2500, tfinal, 0.1)
sol2 = dde.sample(tfinal-2500, tfinal, 0.1)
t = sol1['t']
x1 = sol1['x1']
x2 = sol2['x2']

# plot the solution in phase-space
plt.plot(x1, x2)
plt.xlabel('$x_1(t)$')
plt.ylabel('$x_2(t)$')
filename='/home/maikel/Dropbox/Documents/Uni/Master/2015/' + \
  'Reading course DDE 2/VI/images/tutorial_VI_torus_1000000.png'
plt.savefig(filename, bbox_inches='tight')
plt.show()

# subtract solution components
tfinal=1000000
sol1 = dde.sample(tfinal-1000, tfinal, 0.1)
sol2 = dde.sample(tfinal-1000, tfinal, 0.1)
t = sol1['t']
x1 = sol1['x1']
x2 = sol2['x2']

# plot time series x1
plt.figure(figsize=(10,3))
plt.plot(t,x1,c='royalblue')
plt.xlabel('$t$')
plt.ylabel('$x_1(t)$')
filename='/home/maikel/Dropbox/Documents/Uni/Master/2015/' + \
  'Reading course DDE 2/VI/images/tutorial_VI_torus_1000000_time_series_x1.png'
plt.savefig(filename, bbox_inches='tight')
plt.show()

# plot time series x2
plt.figure(figsize=(10,3))
plt.plot(t,x2,c='royalblue')
plt.xlabel('$t$')
plt.ylabel('$x_2(t)$')
filename='/home/maikel/Dropbox/Documents/Uni/Master/2015/' + \
  'Reading course DDE 2/VI/images/tutorial_VI_torus_1000000_time_series_x2.png'
plt.savefig(filename, bbox_inches='tight')
plt.show()

