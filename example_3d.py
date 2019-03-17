import numpy as np
from mayavi import mlab

import polymer 

dt = 0.0001
T = 100
N = 32
M = 15
k_harm = 30.0
k_F = 150.0
R0 = 2.0
eps = 1.0
sigma = 1.0
omega = 0.001

mlab.figure(1, bgcolor=(0, 0, 0), size=(800, 800))
mlab.clf()

pol = polymer.Polymer(N,k_harm,k_F,R0,eps,sigma,omega)
pol.randomwalk_configuration(True)

p = mlab.points3d(pol.x,pol.y,pol.z,
                  scale_factor=1.2,
                  resolution=20,
                  color=(1, 0, 0),
                  scale_mode='none')                  

def V(x,y,z):
    return -0.5*omega*(x**2 + y**2 + z**2)*0.1

X,Y,Z = np.mgrid[-M:M:100j,-M:M:100j,-M:M:100j]
source = mlab.pipeline.scalar_field(X,Y,Z,V)

min = V(X,Y,Z).min()
max = V(X,Y,Z).max()
vol = mlab.pipeline.volume(source,  vmin=min + 0.8 * (max - min),
                                    vmax=min + 0.9 * (max - min))


@mlab.animate(delay=50)
def anim():
    i = 0
    while True:
        for j in range(0,200):
            pol.leapfrog_update(dt)
        p.mlab_source.x = pol.x
        p.mlab_source.y = pol.y
        p.mlab_source.z = pol.z
        mlab.savefig(filename='output/frames/frame{0}.png'.format(i))
        i += 1
        yield

anim()

#mlab.view(132, 54, 45, [21, 20, 21.5])
mlab.show()