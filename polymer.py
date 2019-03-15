import numpy as np

class Polymer(object):

    def __init__(self, N, k_harm, eps, sigma, omega):
        self.N = N
        self.x = np.arange(0,N,dtype='float')
        self.y = np.ones(N)*N#np.arange(0,N,dtype='float')
        self.vx = np.zeros(N)
        self.vy = np.zeros(N)
        self.k_harm = k_harm
        self.eps = eps
        self.sigma2 = sigma**2
        self.cutoff = 1.1224620483093729
        self.omega = omega

        self.H = self.U_tot()

    def U_harm(self):
        U = 0.0
        for i in range(0,self.N-1):
            r2 = (self.x[i] - self.x[i+1])**2 + (self.y[i] - self.y[i+1])**2
            U += 0.5*self.k_harm*r2
        return U

    def U_LJ(self):
        U = 0.0
        for i in range(0,self.N):
            for j in range(i+1,self.N):
                dx = self.x[i] - self.x[j]
                dy = self.y[i] - self.y[j]
                r2 = dx**2 + dy**2
                if r2 < self.cutoff:
                    fr2 = self.sigma2 / r2
                    fr6 = fr2**3
                    U += 4*self.eps*fr6*(fr6 - 1.0)
        return U

    def U_ext(self):
        U = 0.0
        for i in range(0,self.N):
            U += 0.5*self.omega*(self.x[i]**2 + self.y[i]**2)
        return U

    def F_harm(self):
        F_x = np.zeros(self.N)
        F_y = np.zeros(self.N)
        
        F_x[0] = -self.k_harm*(self.x[0] - self.x[1])
        F_y[0] = -self.k_harm*(self.y[0] - self.y[1])

        F_x[self.N-1] = self.k_harm*(self.x[self.N-2] - self.x[self.N-1])
        F_y[self.N-1] = self.k_harm*(self.y[self.N-2] - self.y[self.N-1])

        for i in range(1,self.N-1):
            F_x[i] = self.k_harm*(self.x[i-1] -2*self.x[i] + self.x[i+1])
            F_y[i] = self.k_harm*(self.y[i-1] -2*self.y[i] + self.y[i+1])
        return F_x, F_y  

    def F_LJ(self):
        F_x = np.zeros(self.N)
        F_y = np.zeros(self.N)

        for i in range(0,self.N):
            for j in range(i+1,self.N):
                dx = self.x[i] - self.x[j]
                dy = self.y[i] - self.y[j]
                r2 = dx**2 + dy**2
                if r2 < self.cutoff:
                    fr2 = self.sigma2 / r2
                    fr6 = fr2**3
                    fpr = 48.0 * self.eps * fr6 * (fr6 - 0.5) / r2
                    F_x[i] += fpr*dx
                    F_y[i] += fpr*dy
                    F_x[j] -= fpr*dx
                    F_y[j] -= fpr*dy
        return F_x, F_y
    
    def F_ext(self):
        F_x = np.zeros(self.N)
        F_y = np.zeros(self.N)

        for i in range(0,self.N):
            F_x[i] = -self.omega*self.x[i]
            F_y[i] = -self.omega*self.y[i]

        return F_x, F_y

    def F_tot(self):
        Fx_ext, Fy_ext = self.F_ext()
        Fx_LJ, Fy_LJ = self.F_LJ()
        Fx_harm, Fy_harm = self.F_harm()
        #return Fx_ext+Fx_harm, Fy_ext+Fy_harm
        return Fx_ext+Fx_LJ+Fx_harm, Fy_ext+Fy_LJ+Fy_harm

    def U_tot(self):
        return self.U_LJ() +  self.U_harm() + self.U_ext()
