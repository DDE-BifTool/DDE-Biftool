#
# @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
# @Id $Id: torus.py 135 2016-09-12 11:55:18Z mmbosschaert $
#

from pydelay import dde23
from matplotlib import colors
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# number of time units
tfinal = 3000

# define DDE
eqns = {
'u' : '-gamma*u-kappa1*u(t-tau1)-kappa2*u(t-tau2)'
 '-kappa1*c*u*(gamma*u(t-tau1)+kappa1*u(t-tau3)+kappa2*u(t-tau4))'
 '-kappa2*c*u*(gamma*u(t-tau2)+kappa1*u(t-tau4)+kappa2*u(t-tau5))'
 '-kappa1*kappa1*pow(c,2)*u*u(t-tau1)*(gamma*u(t-tau3)'
 '+kappa1*u(t-tau6)+kappa2*u(t-tau7))'
 '-kappa1*kappa2*pow(c,2)*u*u(t-tau1)*(gamma*u(t-tau4)'
 '+kappa1*u(t-tau7)+kappa2*u(t-tau8))'
 '-kappa2*kappa1*pow(c,2)*u*u(t-tau2)*(gamma*u(t-tau4)'
 '+kappa1*u(t-tau7)+kappa2*u(t-tau8))'
 '-kappa2*kappa2*pow(c,2)*u*u(t-tau2)*(gamma*u(t-tau5)'
 '+kappa1*u(t-tau8)+kappa2*u(t-tau9))'
 '-1/2*pow(c*u,2)*kappa1*(pow(gamma,2)*u(t-tau1)'
 '+2*gamma*(kappa1*u(t-tau3)+kappa2*u(t-tau4))'
 '+kappa1*kappa1*u(t-tau6)+2*kappa1*kappa2*u(t-tau7)'
 '+kappa2*kappa2*u(t-tau8))'
 '-1/2*pow(c*u,2)*kappa2*(pow(gamma,2)*u(t-tau2)'
 '+2*gamma*(kappa1*u(t-tau4)+kappa2*u(t-tau5))'
 '+kappa1*kappa1*u(t-tau7)+2*kappa1*kappa2*u(t-tau8)'
 '+kappa2*kappa2*u(t-tau9))'
}

# set parameters
# point with stable trivial equalibrium
a1=1.3
a2=6
params = {
        'kappa1':2.757858545579159,
        'kappa2':3.383471633934356,
        'c':1,
        'gamma':4.75,
        'tau1':a1,
        'tau2':a2,
        'tau3':2*a1,
        'tau4':a1+a2,
        'tau5':2*a2,
        'tau6':3*a1,
        'tau7':a2+2*a1,
        'tau8':2*a2+a1,
        'tau9':3*a2
}

# solve DDE
dde = dde23(eqns=eqns, params=params)
dde.set_sim_params(tfinal=tfinal, dtmax=0.1, AbsTol=1e-8, RelTol=1e-6)
histfunc = {'u': lambda t: 0.036964714041287} 
dde.hist_from_funcs(histfunc, 51)
dde.run()

# subtract solution components
sol0 = dde.sample(tfinal-400, tfinal, 0.01)
sol1 = dde.sample(tfinal-400-a1, tfinal-a1, 0.01)
sol2 = dde.sample(tfinal-400-a2, tfinal-a2, 0.01)
t = sol0['t']
u0 = sol0['u']
u1 = sol1['u']
u2 = sol2['u']

# plot the solution in phase-space
mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(u0, u1, u2,c='royalblue')
ax.legend()

ax.set_xlabel('$u(t)$')
ax.set_ylabel('$u(t-a_1)$')
ax.set_zlabel('$u(t-a_2)$')

ax.view_init(16, 31)
fig.set_size_inches(18.5, 18.5)

filename='/home/maikel/Dropbox/Documents/Uni/Master/2015/Reading course DDE 2' + \
  '/VI/images/tutorial_VI_stable_torus_30000.png'
plt.savefig(filename, bbox_inches='tight')
plt.show()
