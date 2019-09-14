import numpy as np
from pydelay import dde23
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import colors

# number of time units
tfinal = 40000

# define DDE
eqns = {
        'x1' : '(tau0+mu2)*x2',
        'x2' : '(tau0+mu2)*((1+mu1)*x1(t-tau)-0.2*x2(t-tau)-0.2*pow(x1(t-tau),2)'
                '-0.2*x1(t-tau)*x2(t-tau)-0.2*pow(x2(t-tau),2)'
		'-epsilon*(pow(x1,2)-1)*x2-x1+0.5*pow(x1,3))'
}

# set parameters
# period cycle 
params1 = {
        'mu1':-0.0080,
        'mu2':0.0024,
	'tau0':1.757290761249588,
        'epsilon':0.3,
        'tau':1
}

# torus (ns2_br.point(15))
params2 = {
        'mu1':-0.0086,
        'mu2':0.0047,
	'tau0':1.757290761249588,
        'epsilon':0.3,
        'tau':1
}

# torus (ns1_br.point(15))
params3 = {
        'mu1':-0.0086,
        'mu2':0.0047+0.0000032,
	'tau0':1.757290761249588,
        'epsilon':0.3,
        'tau':1
}

# solve DDE
dde = dde23(eqns=eqns, params=params3)
dde.set_sim_params(tfinal=tfinal, dtmax=0.1, AbsTol=1e-8, RelTol=1e-6)
histfunc1 = {'x1': lambda t: -0.0366, 'x2': lambda t: -0.0028 }
# hist function for second torus (param3)
histfunc2 = {'x1': lambda t: 0.0049, 'x2': lambda t: 0.0034 }
dde.hist_from_funcs(histfunc2, 51)
dde.run()

# subtract solution components
sol1 = dde.sample(tfinal-tfinal/2, tfinal, 0.01)
t = sol1['t']
x1 = sol1['x1']
x2 = sol1['x2']
sol3 = dde.sample(tfinal-tfinal/2-1.757290761249588, tfinal-1.757290761249588, 0.01)
t1 = sol3['t']
x3 = sol3['x1']

# plot time series
fig = plt.figure()
plt.figure(1)
plt.subplot(311)
plt.xlabel('$t$')
plt.ylabel('$x_1(t)$')
plt.plot(t, x1, c='royalblue')

plt.subplot(312)
plt.xlabel('$t$')
plt.ylabel('$x_2(t)$')
plt.plot(t, x2, c='royalblue')

plt.subplot(313)
plt.xlabel('$t$')
plt.ylabel('$x_1(t-\\tau)$')
plt.plot(t1, x3, c='royalblue')
fig.set_size_inches(18.5, 10.5)
plt.show()

# subtract solution components
sol1 = dde.sample(tfinal-tfinal/10, tfinal, 0.1)
t = sol1['t']
x1 = sol1['x1']
x2 = sol1['x2']
sol3 = dde.sample(tfinal-tfinal/10-1.757290761249588, tfinal-1.757290761249588, 0.1)
t1 = sol3['t']
x3 = sol3['x1']

# save data to file for tikz
filename='/home/maikel/Dropbox/Documents/Uni/Papers/switching to nonhyperbolic cycles from codim 2 equilibria in DDEs/images/vdpo_torus1.dat'
np.savetxt(filename,np.transpose([t, x1, x2]),fmt='%.3f')
filename='/home/maikel/Dropbox/Documents/Uni/Papers/switching to nonhyperbolic cycles from codim 2 equilibria in DDEs/images/vdpo_torus2.dat'
np.savetxt(filename,np.transpose([t1, x3]),fmt='%.3f')

# subtract solution components
sol1 = dde.sample(tfinal-tfinal/10, tfinal, 0.001)
t = sol1['t']
x1 = sol1['x1']
x2 = sol1['x2']
sol3 = dde.sample(tfinal-tfinal/10-1.757290761249588, tfinal-1.757290761249588, 0.001)
t1 = sol3['t']
x3 = sol3['x1']


# Poincare section
def poincaresection(x1, x2, x3,x1_label, x2_label, val):
  zero_cross = np.where(np.diff(np.sign(x3-val)))[0]
#  zero_cross = zero_cross[1::100]
  fig = plt.figure()
  plt.figure(1)
  plt.xlabel(x1_label)
  plt.ylabel(x2_label)
  plt.plot(x1[zero_cross], x2[zero_cross],'.', c='royalblue')
  # save data to file for tikz
  filename='/home/maikel/Dropbox/Documents/Uni/Papers/switching to nonhyperbolic cycles from codim 2 equilibria in DDEs/images/vdpo_torus_poincare.dat'
  np.savetxt(filename,np.transpose([x1[zero_cross], x2[zero_cross]]),fmt='%.5f')
  plt.show()
  return;

x1_label='$x_1(t)$'
x2_label='$x_2(t)$'
poincaresection(x1,x3,x2,x1_label,x2_label,-0.009)

# plot the solution in phase-space
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.plot(x1, x2, x3, c='royalblue', alpha=0.2,zorder=-100)
m=len(x1)
ax.plot(x1[m-40000:m], x2[m-40000:m], x3[m-40000:m], c='indianred', alpha=1,zorder=-100)

X1,X2 = np.meshgrid(np.linspace(1.5*min(x1),1.5*max(x1),100),np.linspace(1.5*min(x3),1.5*max(x3),100))
Z = 0*(X1+X2)+-0.00

ax.plot_surface(X1, Z, X2, color='royalblue', alpha=1, linewidth=0,zorder=100)
ax.view_init(24,-16)

zero_cross = np.where(np.diff(np.sign(x2)))[0]
ax.plot(x1[zero_cross], 0*x1[zero_cross],x3[zero_cross],'.',c='indianred')

ax.legend()

ax.set_xlabel('$x_1(t)$')
ax.set_ylabel('$x_2(t)$')
ax.set_zlabel('$x_2(t-\\tau_0)$')

mpl.rcParams['legend.fontsize'] = 10

plt.show()

