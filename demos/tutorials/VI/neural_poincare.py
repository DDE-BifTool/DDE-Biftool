# Simulation of the neural DDE analyzed in
# Visser, S. (2013). From spiking neurons to brain waves.
#
# @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
# @Id $Id: neural_poincare.py 135 2016-09-12 11:55:18Z mmbosschaert $
#
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from pydelay import dde23

# number of time units
tfinal = 1000000

# define DDE
eqns = {
	'x1' : '-x1-a*(tanh(b*x1(t-tau1)-1)+tanh(1))*pow(cosh(1),2)'
		'+c*(tanh(d*x2(t-tau2)-1)+tanh(1))*pow(cosh(1),2)',
	'x2' : '-x2-a*(tanh(b*x2(t-tau1)-1)+tanh(1))*pow(cosh(1),2)'
		'+c*(tanh(d*x1(t-tau2)-1)+tanh(1))*pow(cosh(1),2)'
}

# set parameters
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

# subtract solution components
tfinal=30000
tau1=1
sol1 = dde.sample(tfinal-20000, tfinal, 0.001)
sol2 = dde.sample(tfinal-20000, tfinal, 0.001)
sol3 = dde.sample(tfinal-20000-tau1, tfinal-tau1, 0.001)
t = sol1['t']
x1 = sol1['x1']
x2 = sol2['x2']
x3 = sol3['x1']

# Poincare section
def poincaresection(x1, x2, x3,x1_label, x2_label, val, color):
  zero_cross = np.where(np.diff(np.sign(x3-val)))[0]
  plt.xlabel(x1_label)
  plt.ylabel(x2_label)
  plt.plot(x1[zero_cross], x2[zero_cross],'.', c=color)
  return;

x1_label='$x_1(t)$'
x2_label='$x_2(t)$'
poincaresection(x1,x2,x3,x1_label,x2_label,0,'royalblue')

# subtract solution components
tfinal=300000
sol1 = dde.sample(tfinal-20000, tfinal, 0.001)
sol2 = dde.sample(tfinal-20000, tfinal, 0.001)
sol3 = dde.sample(tfinal-20000-tau1, tfinal-tau1, 0.001)
t = sol1['t']
x1 = sol1['x1']
x2 = sol2['x2']
x3 = sol3['x1']

poincaresection(x1,x2,x3,x1_label,x2_label,0,'forestgreen')

# subtract solution components
tfinal=1000000
sol1 = dde.sample(tfinal-2000, tfinal, 0.001)
sol2 = dde.sample(tfinal-2000, tfinal, 0.001)
sol3 = dde.sample(tfinal-2000-tau1, tfinal-tau1, 0.001)
t = sol1['t']
x1 = sol1['x1']
x2 = sol2['x2']
x3 = sol3['x1']

poincaresection(x1,x2,x3,x1_label,x2_label,0,'red')
plt.show()
