import os
import glob 
import subprocess
import numpy as np
from mayavi import mlab

import polymer 

wrkdir = '../../_output/frames/'

dt = 0.001
Nfr = 100
N = 64
M = 20
k_harm = 30.0
k_F = 15.0
R0 = 2.0
eps = 1.0
sigma = 1.0
omega = 0.0001
gamma = 0.9
T = 1.0

print('clearing workdir...')
files = glob.glob(wrkdir + '/*')
for f in files:
    os.remove(f)

mlab.figure(1, bgcolor=(0, 0, 0), size=(800, 800))
mlab.clf()

pol = polymer.Polymer(N,k_harm,k_F,R0,eps,sigma,omega,gamma,T,dt)
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
vol = mlab.pipeline.volume(source,  vmin=min + 0.85 * (max - min),
                                    vmax=min + 0.9 * (max - min))

@mlab.animate(delay=50)
def anim():
    for i in range(0,Nfr):
        for _ in range(0,50):
            pol.bbk_update()
        p.mlab_source.x = pol.x
        p.mlab_source.y = pol.y
        p.mlab_source.z = pol.z
        mlab.savefig(filename='../../_output/frames/frame{0}.png'.format(i))
        print('Wrote frame{0}.png'.format(i))
        yield

print('Creating animation frames...')
anim()
mlab.show()
print('Animation complete, ffmpeg rendering...')
proc = subprocess.Popen('ffmpeg -f image2 -r 10 -i frame%0d.png -vcodec mpeg4 -y out.mp4', cwd=wrkdir)
proc.wait()
print('Complete!')