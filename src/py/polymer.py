import numpy as np

class Polymer(object):

    def __init__(self, N, k_harm, k_F, R0, eps, sigma, omega, gamma, T, dt):
        self.N = N
        self.x = np.arange(0,N,dtype='float')
        self.y = np.ones(N)*N
        self.z = np.zeros(N)
        self.vx = np.zeros(N)
        self.vy = np.zeros(N)
        self.vz = np.zeros(N)
        self.vhx = np.zeros(N)
        self.vhy = np.zeros(N)
        self.vhz = np.zeros(N)
        self.fx = np.zeros(N)
        self.fy = np.zeros(N)
        self.fz = np.zeros(N)
        self.k_harm = k_harm
        self.k_F = k_F
        self.R0 = R0
        self.eps = eps
        self.sigma2 = sigma**2
        self.cutoff2 = 1.1224620483093729**2
        self.omega = omega
        self.cmx = self.cmy = self.cmz = 0.0
        self.rgx = self.rgy = self.rgz = 0.0
               
        self.fx, self.fy, self.fz = self.F_tot()
        self.update_cm()
        self.update_rg()  

        self.dt = dt
        self.T = T 
        self.gamma = gamma
        self.mu = np.sqrt(2.0*gamma*T/self.dt)

        self.rx = np.random.normal(0,1,self.N)*self.mu
        self.ry = np.random.normal(0,1,self.N)*self.mu
        self.rz = np.random.normal(0,1,self.N)*self.mu

    def U_harm(self):
        U = 0.0
        for i in range(0,self.N-1):
            r2 = (self.x[i] - self.x[i+1])**2 + \
                    (self.y[i] - self.y[i+1])**2 + \
                    (self.z[i] - self.z[i+1])**2
            U += 0.5*self.k_harm*r2
        return U

    def U_LJ(self):
        U = 0.0
        for i in range(0,self.N):
            for j in range(i+1,self.N):
                dx = self.x[i] - self.x[j]
                dy = self.y[i] - self.y[j]
                dz = self.z[i] - self.z[j]
                r2 = dx**2 + dy**2 + dz**2
                if r2 < self.cutoff2:
                    fr2 = self.sigma2 / r2
                    fr6 = fr2**3
                    U += 4.0*self.eps*fr6*(fr6 - 1.0) + self.eps
        return U

    def U_ext(self):
        U = 0.0
        for i in range(0,self.N):
            U += 0.5*self.omega*(self.x[i]**2 + self.y[i]**2 \
             + self.z[i]**2)
        return U

    def U_FENE(self):
        U = 0.0        
        for i in range(0,self.N-1):
            r2 = (self.x[i] - self.x[i+1])**2 +  \
                    (self.y[i] - self.y[i+1])**2 + \
                    (self.z[i] - self.z[i+1])**2
            U += 0.5*self.k_F*self.R0*np.log(1-r2/self.R0**2)
        return U

    def F_harm(self):
        F_x = np.zeros(self.N)
        F_y = np.zeros(self.N)
        F_z = np.zeros(self.N)

        if self.N == 1:
            return F_x, F_y, F_z
        
        F_x[0] = self.k_harm*(self.x[1] - self.x[0])
        F_y[0] = self.k_harm*(self.y[1] - self.y[0])
        F_z[0] = self.k_harm*(self.z[1] - self.z[0])

        F_x[self.N-1] = -self.k_harm*(self.x[self.N-1] - \
                self.x[self.N-2])
        F_y[self.N-1] = -self.k_harm*(self.y[self.N-1] - \
                self.y[self.N-2])
        F_z[self.N-1] = -self.k_harm*(self.z[self.N-1] - \
                self.z[self.N-2])

        for i in range(1,self.N-1):
            F_x[i] = self.k_harm*(self.x[i-1] - \
                     2.0*self.x[i] + self.x[i+1])
            F_y[i] = self.k_harm*(self.y[i-1] - \
                     2.0*self.y[i] + self.y[i+1])
            F_z[i] = self.k_harm*(self.z[i-1] - \
                     2.0*self.z[i] + self.z[i+1])
        
        return F_x, F_y, F_z

    def F_LJ(self):
        F_x = np.zeros(self.N)
        F_y = np.zeros(self.N)
        F_z = np.zeros(self.N)

        for i in range(0,self.N-1):
            for j in range(i+1,self.N):
                dx = self.x[i] - self.x[j]
                dy = self.y[i] - self.y[j]
                dz = self.z[i] - self.z[j]
                r2 = dx**2 + dy**2 + dz**2
                if r2 < self.cutoff2:
                    fr2 = self.sigma2 / r2
                    fr6 = fr2**3
                    fpr = 48.0*self.eps*(fr6*(fr6-0.5))/r2
                    F_x[i] += fpr*dx
                    F_y[i] += fpr*dy
                    F_z[i] += fpr*dz
                    F_x[j] -= fpr*dx
                    F_y[j] -= fpr*dy
                    F_z[i] -= fpr*dz
        return F_x, F_y, F_z
    
    def F_ext(self):
        F_x = np.zeros(self.N)
        F_y = np.zeros(self.N)
        F_z = np.zeros(self.N)

        for i in range(0,self.N):
            F_x[i] = -self.omega*self.x[i]
            F_y[i] = -self.omega*self.y[i]
            F_z[i] = -self.omega*self.z[i]

        return F_x, F_y, F_z

    def F_FENE(self):
        F_x = np.zeros(self.N)
        F_y = np.zeros(self.N)
        F_z = np.zeros(self.N)

        for i in range(0,self.N-1):
            dx = self.x[i+1] - self.x[i]
            dy = self.y[i+1] - self.y[i]
            dz = self.z[i+1] - self.z[i]
            r2 = dx**2 + dy**2 + dz**2
            fac = self.k_F/(1-r2/self.R0**2)
            F_x[i] += fac*dx
            F_y[i] += fac*dy
            F_z[i] += fac*dz
            F_x[i+1] -= fac*dx
            F_y[i+1] -= fac*dy
            F_z[i+1] -= fac*dz
        
        return F_x, F_y, F_z

    def F_tot(self):
        Fx_ext, Fy_ext, Fz_ext = self.F_ext()
        Fx_LJ, Fy_LJ, Fz_LJ = self.F_LJ()
        Fx_FENE, Fy_FENE, Fz_FENE = self.F_FENE()
        #Fx_harm, Fy_harm, Fz_harm = self.F_harm()
        return Fx_ext+Fx_LJ+Fx_FENE,Fy_ext+Fy_LJ+Fy_FENE, \
                Fz_ext+Fz_LJ+Fz_FENE
        #return Fx_harm+Fx_ext, Fy_harm+Fy_ext, Fz_harm+Fz_ext 

    def U_tot(self):
        U_kin = 0.0
        for i in range(0, self.N):
            U_kin += 0.5*(self.vx[i]**2 + self.vy[i]**2 + \
             self.vz[i]**2) 
        return self.U_LJ() +  self.U_FENE() + self.U_ext() + U_kin
        #return self.U_harm() + self.U_ext() + U_kin

    def update_cm(self):
        self.cmx = 0.0
        self.cmy = 0.0
        self.cmz = 0.0
        for i in range(0,self.N):
            self.cmx += self.x[i]
            self.cmy += self.y[i]
            self.cmz += self.z[i]
        self.cmx /= self.N
        self.cmy /= self.N
        self.cmz /= self.N

    def update_rg(self):
        self.update_cm()
        self.rgx = 0.0
        self.rgy = 0.0
        self.rgz = 0.0
        for i in range(0,self.N):
            self.rgx += (self.cmx - self.x[i])**2
            self.rgy += (self.cmy - self.y[i])**2
            self.rgz += (self.cmz - self.z[i])**2
        self.rgx /= self.N**2
        self.rgy /= self.N**2
        self.rgz /= self.N**2

    def leapfrog_update(self):
        for i in range(0,self.N):
            self.vhx[i] = self.vx[i] + self.fx[i]*self.dt/2.0
            self.vhy[i] = self.vy[i] + self.fy[i]*self.dt/2.0
            self.vhz[i] = self.vz[i] + self.fz[i]*self.dt/2.0
            self.x[i] += self.vhx[i] * self.dt
            self.y[i] += self.vhy[i] * self.dt
            self.z[i] += self.vhz[i] * self.dt
        
        self.fx, self.fy, self.fz = self.F_tot()
        for i in range(0,self.N):
            self.vx[i] = self.vhx[i] + self.fx[i]*self.dt/2.0
            self.vy[i] = self.vhy[i] + self.fy[i]*self.dt/2.0
            self.vz[i] = self.vhz[i] + self.fz[i]*self.dt/2.0

    def bbk_update(self):

        for i in range(0,self.N):
            self.vhx[i] = (1.0-0.5*self.gamma*self.dt)*self.vx[i] + \
                (self.fx[i] + self.rx[i]) * self.dt / 2.0
            self.vhy[i] = (1.0-0.5*self.gamma*self.dt)*self.vy[i] + \
                (self.fy[i] + self.ry[i]) * self.dt  / 2.0
            self.vhz[i] = (1.0-0.5*self.gamma*self.dt)*self.vz[i] + \
                (self.fz[i] + self.rz[i]) * self.dt  / 2.0

            self.x[i] += self.vhx[i] * self.dt
            self.y[i] += self.vhy[i] * self.dt
            self.z[i] += self.vhz[i] * self.dt

        self.fx, self.fy, self.fz = self.F_tot()        
        self.rx = np.random.normal(0,1,self.N)*self.mu
        self.ry = np.random.normal(0,1,self.N)*self.mu
        self.rz = np.random.normal(0,1,self.N)*self.mu

        for i in range(0,self.N):
            self.vx[i] = (self.vhx[i] + 0.5*(self.fx[i] + self.rx[i])\
                *self.dt)/(1.0+0.5*self.gamma*self.dt)
            self.vy[i] = (self.vhy[i] + 0.5*(self.fy[i] + self.ry[i])\
                *self.dt)/(1.0+0.5*self.gamma*self.dt)
            self.vz[i] = (self.vhz[i] + 0.5*(self.fz[i] + self.rz[i])\
                *self.dt)/(1.0+0.5*self.gamma*self.dt)


    def randomwalk_configuration(self, center):
        self.x[0] = self.y[0] = self.z[0] = 0.0
        for i in range(1,self.N):
            succeeded = False
            while not succeeded:
                theta = 1.0*np.pi*np.random.rand()
                phi = 2.0*np.pi*np.random.rand()
                self.x[i] = self.x[i-1] + np.sin(theta)*\
                                          np.cos(phi)
                self.y[i] = self.y[i-1] + np.sin(theta)*\
                                          np.sin(phi)
                self.z[i] = self.z[i-1] + np.cos(theta)

                succeeded = True 
                for j in range(0,i-1):
                    r2 = (self.x[i] - self.x[j])**2 +  \
                         (self.y[i] - self.y[j])**2 + \
                         (self.z[i] - self.z[j])**2
                    if r2 < self.cutoff2:
                        succeeded = False
        
        self.update_rg()
        self.update_cm()
        
        if center:
            for i in range(0,self.N):
                self.x[i] -= self.cmx
                self.y[i] -= self.cmy
                self.z[i] -= self.cmz
            