import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import polymer

dt = 0.0005
T = 100
N = 8
k_harm = 15.0
k_F = 15.0
R0 = 2.0
eps = 1.0
sigma = 1.0
omega = 10.0

rg = np.zeros(T)
cm = np.zeros(T)
E = np.zeros(T)

pol = polymer.Polymer(N,k_harm,k_F,R0,eps,sigma,omega)
pol.randomwalk_configuration(False)

def update(i):
    for j in range(0,100):
        pol.leapfrog_update(dt)
    cm[i] = np.sqrt(pol.cmx**2 + pol.cmy**2)    
    rg[i] = np.sqrt(pol.rgx + pol.rgy)
    E[i] = pol.U_tot() 
    ax.clear()
    ax.set_xlim(-N/2,N/2)
    ax.set_ylim(-N/2,N/2)
    ax.plot(pol.x, pol.y, 'ro-')
    print(i, ': ', E[i])

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
a = anim.FuncAnimation(fig, update, frames=T, repeat=False)
a.save("output/dynamics.gif",writer="imagemagick")

fig2 = plt.figure()
ax2 = fig2.add_subplot(1,1,1)
ax2.plot(np.arange(0,T)*dt, rg)
fig2.savefig('output/rg.png')

fig3 = plt.figure()
ax3 = fig3.add_subplot(1,1,1)
ax3.plot(np.arange(0,T)*dt, cm)
fig3.savefig('output/cm.png')

fig3 = plt.figure()
ax3 = fig3.add_subplot(1,1,1)
ax3.plot(np.arange(0,T)*dt, E)
fig3.savefig('output/energy.png')