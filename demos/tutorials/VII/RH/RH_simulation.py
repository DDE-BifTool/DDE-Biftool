import numpy as np
from pydelay import dde23
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import colors

# Number of time units
tfinal_cycle = 10000
tfinal_torus = 20000

# Define DDE
eqns = {
    'x' : 'y-a*pow(x,3)+b*pow(x(t-tau),2)-c*z+Iapp',
    'y' : 'c-d*pow(x,2)-y',
    'z' : 'r*(S*(x-chi)-z)'
}

# Set parameters for torus
params_torus = {
    'a':1,'b':3,'c':1,'d':5.0,'chi':-1.6,'r':1.4,
    'tau':0.940246941050084,
    'Iapp':-18.902420391705071,
    'S':-8.045234985422740
}

# Set parameters for torus
params_cycle = params_torus.copy()
params_cycle['Iapp']=-18.886177304147466
params_cycle['S']=-8.044197084548104

# Set number of timesteps from the end to plot
timesteps_torus=300
timesteps_cycle=10

# Select periodic orbit or torus
tfinal,params,timesteps=tfinal_cycle,params_cycle,timesteps_cycle
#tfinal,params,timesteps=tfinal_torus,params_torus,timesteps_torus

# Solve DDE
dde = dde23(eqns=eqns, params=params)
dde.set_sim_params(tfinal=tfinal)
dde.set_sim_params(tfinal=tfinal, dtmax=0.001)
histfunc = {
    'x': lambda t: 1.097167540709727, 
    'y': lambda t: -5.018883061935152,
    'z': lambda t: -21.577340325677817
}
dde.hist_from_funcs(histfunc, 51)
dde.run()

# Subtract solution components
sol0 = dde.sample(tfinal-timesteps,tfinal, 0.001)
t = sol0['t']
x = sol0['x']
y = sol0['y']
z = sol0['z']

# Plot time series
fig = plt.figure()
plt.figure(1)
plt.subplot(311)
plt.xlabel('$t$')
plt.ylabel('$x(t)$')
plt.plot(t, x, c='royalblue')

plt.subplot(312)
plt.xlabel('$t$')
plt.ylabel('$y(t)$')
plt.plot(t, y, c='royalblue')

plt.subplot(313)
plt.xlabel('$t$')
plt.ylabel('$z(t)$')
plt.plot(t, z, c='royalblue')
fig.set_size_inches(18.5, 10.5)
plt.show()

# Plot the solution in phase-space
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.plot(x, y, z, c='royalblue')

ax.set_xlabel('$x(t)$')
ax.set_ylabel('$y(t)$')
ax.set_zlabel('$z(t)$')

# Fix to get z_label inside the figure
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

ax.view_init(9,-44)
plt.show()


# Subtract solution components for Poincare-section
sol0 = dde.sample(tfinal-1000, tfinal, 0.001)
t = sol0['t']
x = sol0['x']
y = sol0['y']
z = sol0['z']

# Poincar!{\color{comment}\'e}! section
def poincaresection(x1, x2, x3,x1_label, x2_label, val):
  zero_cross = np.where(np.diff(np.sign(x3-val)))
  plt.figure(1)
  plt.xlabel(x1_label)
  plt.ylabel(x2_label)
  plt.plot(x1[zero_cross], x2[zero_cross],'.', c='royalblue')
  plt.show()
  return

x1_label='$x(t)$'
x2_label='$y(t)$'
poincaresection(x,y,z,x1_label,x2_label,-21.75)
