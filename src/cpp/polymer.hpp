#ifndef H_POLYMER
#define H_POLYMER

#include <vector>
#include <random>
#include <tuple>

enum PolymerModel {
  Harmonic,
  SelfAvoiding
};

class Polymer
{
private:
  std::default_random_engine gen;
  std::normal_distribution<float> normal;
  std::uniform_real_distribution<float> uniform;

  PolymerModel model;
  float k_harm, k_F, R0, eps, sigma, omega, gamma, T, dt;
  std::vector<float> vhx, vhy, vhz;
  std::vector<float> fx_harm, fy_harm, fz_harm;
  std::vector<float> fx_LJ, fy_LJ, fz_LJ;
  std::vector<float> fx_FENE, fy_FENE, fz_FENE;
  std::vector<float> fx_ext, fy_ext, fz_ext;

  float mu,PI,sigma2,cutoff2;

  //helper variables for numeric performance
  float r2,fac,fr6,fpr,dx,dy,dz,b1,b2;

public:
  int N;

  std::vector<float> x, y, z;
  std::vector<float> vx, vy, vz;
  std::vector<float> fx, fy, fz;

  float cmx, cmy, cmz;
  float rgx, rgy, rgz;

  Polymer(PolymerModel model, int N, float k_harm, float k_F, 
          float R0, float eps, float sigma, float omega,
          float gamma, float T, float dt);
  
  
  /*
   * Updates the center of mass
   */ 
  void Cm();
  
  /*
   * Updates the radius of gyration
   */ 
  void Rg();

  /*
   * Creates a self avoiding random walk configuarion for the polymer
   * parameters:
   * bool center: if the polymer is to be replaced to the center
   */
  void ConfigureSAW(bool center);

  /*
   * Harmonic energy of the polymer
   */
  float U_harm();

  /*
   * FENE energy of the polymer
   */
  float U_FENE();

  /*
   * Lennard-Jones energy of the polymer
   */
  float U_LJ();

  /*
   * External energy of the polymer
   */
  float U_ext();

  /*
   * Total potential energy of the polymer
   */
  float U_tot();

  /*
   * Total energy of polymer
   */
  float E_tot();

  /*
   * Update harmonic forces on polymer
   */
  void UpdateF_harm();

  /*
   * Update harmonic forces on polymer
   */
  void UpdateF_LJ();
  
  /*
   * Update FENE forces on polymer
   */
  void UpdateF_FENE();

  /*
   * Update external forces on polymer
   */
  void UpdateF_ext();

  /*
   * Update total forces of polymer
   */
  void UpdateF();
  
  /*
   * Integrate with Velocity Verlet type integration
   */
  void UpdateLeapFrog();
  
  /*
   * Integrate with BBK integration
   */
  void UpdateBBK();
};

#endif