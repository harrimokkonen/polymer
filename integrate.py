import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import polymer

dt = 0.005
T = 100

pol = polymer.Polymer(10,1,1,1,5)
fig = plt.figure()
ax = fig.add_subplot(1,1,1)

def update_polymer():
    fx, fy = pol.F_tot()
    vhx = np.zeros(pol.N)
    vhy = np.zeros(pol.N)
    for t in range(0,100):
        for i in range(0,pol.N):
            vhx[i] = pol.vx[i] + fx[i]*dt/2.0
            vhy[i] = pol.vy[i] + fy[i]*dt/2.0
            pol.x[i] += vhx[i] * dt
            pol.y[i] += vhy[i] * dt
        
        fx, fy = pol.F_tot()
        for i in range(0,pol.N):
            pol.vx[i] = vhx[i] + fx[i]*dt/2.0
            pol.vy[i] = vhy[i] + fy[i]*dt/2.0    

def update(i):
    update_polymer()
    ax.clear()
    ax.set_xlim(-10,10)
    ax.set_ylim(-10,10)
    ax.plot(pol.x, pol.y, 'ro')
    print(i, ': ', pol.H)


a = anim.FuncAnimation(fig, update, frames=T, repeat=False)
plt.show()

