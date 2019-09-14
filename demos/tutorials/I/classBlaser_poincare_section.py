##
# @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
# @Id $Id: classBlaser_poincare_section.py 135 2016-09-12 11:55:18Z mmbosschaert $
##

import numpy as np
from pydelay import dde23
import matplotlib.pyplot as plt
from matplotlib import colors

# number of time units
tfinal = 9000

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
dde.set_sim_params(tfinal=tfinal)
histfunc = {'x': lambda t: 1.1, 'y': lambda t: 0.4} 
dde.hist_from_funcs(histfunc, 51)
dde.run()

# subtract solution components
tau=1.60017288
sol0 = dde.sample(tfinal-tfinal/6, tfinal, 0.0001)
sol1 = dde.sample(tfinal-tfinal/6-tau, tfinal-tau, 0.0001)
x  = sol0['x']
y  = sol0['y']
z  = sol1['x']

# Poincare section
def poincaresection(x1, x2, x3,x1_label, x2_label, val, fname):
  zero_cross = np.where(np.diff(np.sign(x3-val)))[0]
  fig = plt.figure()
  plt.figure(1)
  plt.xlabel(x1_label)
  plt.ylabel(x2_label)
  plt.plot(x1[zero_cross], x2[zero_cross],'.', c='royalblue')
  plt.show()
  return;

x1_label='$x(t)$'
x2_label='$y(t)$'
x3_label='$x(t-\\tau)$'
poincaresection(x,y,z,x1_label,x2_label,1,'torus_poincare_section_1.eps')
poincaresection(y,z,x,x2_label,x3_label,1.05,'torus_poincare_section_2.eps')
