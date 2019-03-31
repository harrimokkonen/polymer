#include "polymer.hpp"
#include <iostream>
Polymer::Polymer(PolymerModel model, int N, float k_harm, float k_F, float R0,
                        float eps, float sigma, float omega,
                        float gamma, float T, float dt)
    : model(model), N(N), k_F(k_F), R0(R0), eps(eps),
      sigma(sigma), omega(omega), gamma(gamma),
      T(T), dt(dt),
      x(N, 0.0), y(N, 0.0), z(N, 0.0),
      vx(N, 0.0), vy(N, 0.0), vz(N, 0.0),
      vhx(N, 0.0), vhy(N, 0.0), vhz(N, 0.0),
      fx(N, 0.0), fy(N, 0.0), fz(N, 0.0),
      fx_harm(N, 0.0), fy_harm(N, 0.0), fz_harm(N, 0.0),
      fx_LJ(N, 0.0), fy_LJ(N, 0.0), fz_LJ(N, 0.0),
      fx_FENE(N, 0.0), fy_FENE(N, 0.0), fz_FENE(N, 0.0),
      fx_ext(N, 0.0), fy_ext(N, 0.0), fz_ext(N, 0.0),
      cmx(0.0), cmy(0.0), cmz(0.0),
      rgx(0.0), rgy(0.0), rgz(0.0),
      normal(0.0, 1.0), uniform(0.0, 1.0)
{
  mu = std::sqrt(2.0 * gamma * T / dt);
  PI = std::atan(2.0);
  sigma2 = sigma * sigma;
  cutoff2 = 1.25992104989 * sigma2;
  b1 = 1.0-0.5*gamma*dt;
  b2 = 1.0+0.5*gamma*dt;
  UpdateF();
};

void Polymer::ConfigureSAW(bool center)
{
  x[0] = 0.0;
  y[0] = 0.0;
  z[0] = 0.0;
  float theta, phi, dx, dy, dz, r2;

  for (int i = 1; i < N; ++i)
  {
    bool succeeded = false;
    while (!succeeded)
    {
      theta = PI * uniform(gen);
      phi = 2.0 * PI * uniform(gen);
      x[i] = x[i - 1] + std::sin(theta) * std::cos(phi);
      y[i] = y[i - 1] + std::sin(theta) * std::sin(phi);
      z[i] = z[i - 1] + std::cos(theta);

      succeeded = true;
      for (int j = 0; j < i - 1; ++j)
      {
        dx = x[i] - x[j];
        dy = y[i] - y[j];
        dz = z[i] - z[j];
        r2 = dx * dx + dy * dy + dz * dz;
        if (r2 < cutoff2)
          succeeded = false;
      }
    }
  }
  UpdateF();
  if (center)
  {
    Cm();
    for (int i = 0; i < N; ++i)
    {
      x[i] -= cmx;
      y[i] -= cmy;
      z[i] -= cmz;
    }
  }
};

void Polymer::Cm()
{
  cmx, cmy, cmz = 0.0;
  for (int i = 0; i < N; ++i)
  {
    cmx += x[i];
    cmy += y[i];
    cmz += z[i];
  }
  cmx /= N;
  cmy /= N;
  cmy /= N;
};

void Polymer::Rg()
{
  Cm();
  rgx, rgy, rgz = 0.0;
  for (int i = 0; i < N; ++i)
  {
    rgx += (cmx - x[i]) * (cmx - x[i]);
    rgy += (cmy - y[i]) * (cmy - y[i]);
    rgz += (cmz - z[i]) * (cmz - z[i]);
  }
  rgx /= N * N;
  rgy /= N * N;
  rgy /= N * N;
};

float Polymer::U_harm()
{
  float U = 0.0;
  float r2;
  for (int i = 0; i < N - 1; ++i)
  {
    r2 = (x[i] - x[i + 1]) * (x[i] - x[i + 1]) 
       + (y[i] - y[i + 1]) * (y[i] - y[i + 1]) 
       + (z[i] - z[i + 1]) * (z[i] - z[i + 1]);
    U += 0.5 * k_harm * r2;
  }
  return U;
};

float Polymer::U_LJ()
{
  float U = 0.0;
  for (int i = 0; i < N; ++i)
  {
    for (int j = i + 1; j < N; ++j)
    {
      r2 = (x[i] - x[i + 1]) * (x[i] - x[i + 1]) 
         + (y[i] - y[i + 1]) * (y[i] - y[i + 1]) 
         + (z[i] - z[i + 1]) * (z[i] - z[i + 1]);
      if (r2 < cutoff2)
      {
        fr6 = sigma2 * sigma2 * sigma2 / r2 / r2 / r2;
        U += 4.0 * eps * fr6 * (fr6 - 1.0) + eps;
      }
    }
  }
  return U;
};

float Polymer::U_FENE()
{
  float U = 0.0;
  for (int i = 0; i < N - 1; ++i)
  {
    r2 = (x[i] - x[i + 1]) * (x[i] - x[i + 1]) 
       + (y[i] - y[i + 1]) * (y[i] - y[i + 1]) 
       + (z[i] - z[i + 1]) * (z[i] - z[i + 1]);
    U += 0.5 * k_F * R0 * std::log(1 - r2 / R0 / R0);
  }
  return U;
};

float Polymer::U_ext()
{
  float U = 0.0;
  for (int i = 0; i < N; ++i)
  {
    U += 0.5 * omega * (x[i] * x[i] + y[i] * y[i] + z[i] * z[i]);
  }
  return U;
};

float Polymer::E_tot() {
  float E = 0.0;
  for (int i = 0; i < N; ++i) {
    E += 0.5*(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
  }
  switch (model)
  {
    case Harmonic:
      E += U_harm() + U_ext();
      break;
    case SelfAvoiding:
      E += U_FENE() + U_LJ() + U_ext();
      break;
    default:
      break;
  }
  return E;
}

void Polymer::UpdateF_harm()
{
  fx_harm[0] = k_harm * (x[1] - x[0]);
  fy_harm[0] = k_harm * (y[1] - y[0]);
  fz_harm[0] = k_harm * (z[1] - z[0]);

  fx_harm[N - 1] = -k_harm * (x[N - 1] - x[N - 2]);
  fy_harm[N - 1] = -k_harm * (y[N - 1] - y[N - 2]);
  fz_harm[N - 1] = -k_harm * (z[N - 1] - z[N - 2]);

  for (int i = 1; i < N - 1; i++)
  {
    fx_harm[i] = k_harm * (x[i - 1] - 2.0 * x[i] + x[i + 1]);
    fy_harm[i] = k_harm * (y[i - 1] - 2.0 * y[i] + y[i + 1]);
    fz_harm[i] = k_harm * (z[i - 1] - 2.0 * z[i] + z[i + 1]);
  }
};

void Polymer::UpdateF_LJ()
{
  std::fill(fx_LJ.begin(), fx_LJ.end(), (float)0.0);
  std::fill(fy_LJ.begin(), fy_LJ.end(), (float)0.0);
  std::fill(fz_LJ.begin(), fz_LJ.end(), (float)0.0);

  for (int i = 0; i < N - 1; ++i)
  {
    for (int j = i + 1; j < N; ++j)
    {
      dx = x[i] - x[j];
      dy = y[i] - y[j];
      dz = z[i] - z[j];
      r2 = dx*dx + dy*dy + dz*dz;

      if(r2 < cutoff2) {
        fr6 = sigma2*sigma2*sigma2/r2/r2/r2;
        fpr = 48.0*eps*(fr6*(fr6-0.5))/r2;
        fx_LJ[i] += fpr*dx;
        fy_LJ[i] += fpr*dy;
        fz_LJ[i] += fpr*dz;        
        fx_LJ[j] -= fpr*dx;
        fy_LJ[j] -= fpr*dy;
        fz_LJ[j] -= fpr*dz;
      }
    }
  }
};

void Polymer::UpdateF_FENE() {
  std::fill(fx_FENE.begin(), fx_FENE.end(), (float)0.0);
  std::fill(fy_FENE.begin(), fy_FENE.end(), (float)0.0);
  std::fill(fz_FENE.begin(), fz_FENE.end(), (float)0.0);

  for (int i = 0; i < N - 1; ++i) {
      dx = x[i] - x[i+1];
      dy = y[i] - y[i+1];
      dz = z[i] - z[i+1];
      r2 = dx*dx + dy*dy + dz*dz;
      fac = k_F/(1-r2/R0/R0);
      fx_FENE[i] -= fac*dx;
      fy_FENE[i] -= fac*dy;
      fz_FENE[i] -= fac*dz;
      fx_FENE[i+1] += fac*dx;
      fy_FENE[i+1] += fac*dy;
      fz_FENE[i+1] += fac*dz;
  }
};

void Polymer::UpdateF_ext() {
    for (int i = 0; i < N; ++i) {
      fx_ext[i] = -omega*x[i];
      fy_ext[i] = -omega*y[i];
      fz_ext[i] = -omega*z[i];
    }
};

void Polymer::UpdateF() {
  switch (model)
  {
    case Harmonic:
      UpdateF_ext();
      UpdateF_harm();
      for(int i = 0; i < N; ++i) {
        fx[i] = fx_harm[i] + fx_ext[i];
        fy[i] = fy_harm[i] + fy_ext[i];
        fz[i] = fz_harm[i] + fz_ext[i];
      }
      break;
    case SelfAvoiding:
      UpdateF_ext();
      UpdateF_FENE();
      UpdateF_LJ();
      for(int i = 0; i < N; ++i) {
        fx[i] = fx_LJ[i] + fx_FENE[i] + fx_ext[i];
        fy[i] = fy_LJ[i] + fy_FENE[i] + fy_ext[i];
        fz[i] = fz_LJ[i] + fz_FENE[i] + fz_ext[i];
      }
      break;
    default:
      break;
  }

  if(T > 0.0) {
    for(int i = 0; i < N; ++i) {
      fx[i] += mu*normal(gen);
      fy[i] += mu*normal(gen);
      fz[i] += mu*normal(gen);
    }
  }
  
};

void Polymer::UpdateLeapFrog() {
  for(int i = 0; i < N; ++i) {
    vhx[i] = vx[i] + fx[i]*dt/2.0;
    vhy[i] = vy[i] + fy[i]*dt/2.0;
    vhz[i] = vz[i] + fz[i]*dt/2.0;
    x[i] += vhx[i]*dt;
    y[i] += vhy[i]*dt;
    z[i] += vhz[i]*dt;
  }
  UpdateF();
  for(int i = 0; i < N; ++i) {
    vx[i] = vhx[i] + fx[i]*dt/2.0;
    vy[i] = vhy[i] + fy[i]*dt/2.0;
    vz[i] = vhz[i] + fz[i]*dt/2.0;
  }
};

void Polymer::UpdateBBK() {

  for(int i = 0; i < N; ++i) {
    vhx[i] = b1*vx[i] + fx[i]*dt/2.0;
    vhy[i] = b1*vy[i] + fy[i]*dt/2.0;
    vhz[i] = b1*vz[i] + fz[i]*dt/2.0;
    
    x[i] += vhx[i]*dt;
    y[i] += vhy[i]*dt;
    z[i] += vhz[i]*dt;
  }
  UpdateF();
  for(int i = 0; i < N; ++i) {
    vx[i] = (vhx[i] + 0.5*fx[i]*dt)/b2;
    vy[i] = (vhy[i] + 0.5*fy[i]*dt)/b2;
    vz[i] = (vhz[i] + 0.5*fz[i]*dt)/b2;
  }
};  