#
# @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
# @Id
#

import numpy as np
import pylab as pl
from pydelay import dde23
from matplotlib import colors

# number of time units
tfinal = 90000

# define DDE
eqns = {
    'x1' : 'tau*x2',
    'x2' : 'tau*(-x1-0.1*x1(t-1)-2*zeta*x2-0.52*x2(t-1)+0.1*pow(x1(t-1),3))'
}

# set parameters
tau=5.901783308978358
zeta1=-0.015485728828307 # periodic orbit
zeta2=zeta1-0.0002 # torus
zeta3=zeta1-0.0004 # torus near bifurcation to 3d torus
zeta4=zeta1-0.000445 # 3d torus
zeta5=zeta1-0.0004496 # 3d torus near destruction

# solve DDE
dde = dde23(eqns=eqns, params={'zeta':zeta1, 'tau':tau})
dde.set_sim_params(tfinal=tfinal, dtmax=0.1, AbsTol=1e-08, RelTol=1e-06)
histfunc = {'x1': lambda t: 0.01, 'x2': lambda t: 0 } 
dde.hist_from_funcs(histfunc, 51)
dde.run()

# subtract solution components
a1=1;
dt=1e-04
sol0 = dde.sample(tfinal-tfinal/10, tfinal, dt)
sol1 = dde.sample(tfinal-tfinal/10-a1, tfinal-a1, dt)
t = sol0['t']
x1 = sol0['x1']
x2 = sol0['x2']
x3 = sol1['x1']

zero_crossings = np.where(np.diff(np.sign(x3)))[0]
# scatter plot
x1_cros=x1[zero_crossings]
x2_cros=x2[zero_crossings]

params = {'figure.figsize': (15, 15),
         'axes.labelsize': 'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pl.rcParams.update(params)

fig = pl.figure()
pl.figure(1)
pl.scatter(x1_cros,x2_cros,s=0.8,color='royalblue')
pl.xlabel('$x_1(t)$')
pl.ylabel('$x_2(t)$')
pl.show()
