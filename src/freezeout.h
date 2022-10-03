//using the freezeout-hypersurface, the relevant quantities to calculate 
//the local and global polarization vectors are calulated

#ifndef SRC_FREEZEOUT_H_
#define SRC_FREEZEOUT_H_
#include<cstring>
#include<iostream>
#include<cstdio>
#include<fstream>
#include<sstream>
#include<cmath>
#include <vector>

#define NY 200           // size of arry for storage of the y-spectrum (maximum)
#define NPT 100          // size of arry for storage of the pt-spectrum(maximum)
#define NPHI 100         // size of arry for storage of the phi-spectrum(maximum)
#define Fourindices 4
#define ParticleMax 320

typedef struct particle {           //imported from MUSIC
    int number;
    char name[50];
    double mass;
    double width;
    int degeneracy;
    int baryon;
    int strange;
    int charm;
    int bottom;
    int isospin;
    double charge;
    int decays;
    int stable;
    int ny;
    int npt;
    int nphi;
    double phimin;
    double phimax;
    double ymax;
    double deltaY;
    double dNdydptdphi[NY][NPT][NPHI+1];
    double slope;           // assymtotic slope of pt-spectrum
    double muAtFreezeOut;   // the chemical potential at freeze-out
                            // for the partial chemical equilibrium EoS
} Particle;


typedef struct de {      //imported from MUSIC
    int number;
    int  reso;          // Montecarlo number of decaying resonance
    int  numpart;       // number of daughter particles after decay
    double branch;      // branching ratio
    int    part[5];     // array of daughter particles Montecarlo number
} de;

typedef struct surfaceElement {             //imported from MUSIC and edited
    double x[4];                            // position in (tau, x, y, eta)
    double sinh_eta_s;                      // caching the sinh and cosh of eta_s
    double cosh_eta_s;
    double s[4];                            // 3-surface vector in (tau, x, y, eta)
    double u[4];                            // flow velocity in (tau, x, y, eta)
    double W[4][4];                         // W^{\mu\nu}
    double q[4];                            // baryon diffusion current
    double d_mu_betanu[4][4];               //gradient of inverse-temperature field 
    double pi_b;                            // bulk pressure
    double rho_B;                           // net baryon density

    double epsilon_f;
    double T_f;
    double mu_B;
    double mu_S;
    double mu_C;
    double eps_plus_p_over_T_FO;             // (energy_density+pressure)/temperature
} SurfaceElement;



class freezeout
{
private:

    //momentum grid parameters specifications
    double pt_min, pt_max;
    int N_steps_pt;
    double phi_min = 0, phi_max = 2*M_PI, dphi;
    int N_steps_phi;
    double y_max, dy;
    int N_steps_y;
    
    //invariant tensors
    int levicivita_4indices_sup[4][4][4][4];
    int gmunu_tilde[4][4];

    //hypersurface calculations
    int N_hypercells;
    std::vector<SurfaceElement> hypersurface_cells;
    Particle* Particle_list;
    Particle* my_Particle_list;
    std :: vector<int> PIDS;
    std::ostringstream surface_dat_file, pdg_dat_file, input_particles;
    de decay[2000];

    //calculating thermal spectrum of this particular particle at a time
    Particle particle_of_interest;

    // net baryon diffusion delta f 
    int deltaf_qmu_coeff_table_length_T;
    int deltaf_qmu_coeff_table_length_mu;
    double delta_qmu_coeff_table_T0, delta_qmu_coeff_table_mu0;
    double delta_qmu_coeff_table_dT, delta_qmu_coeff_table_dmu;
    double **deltaf_qmu_coeff_tb;

    int deltaf_coeff_table_14mom_length_T;
    int deltaf_coeff_table_14mom_length_mu;
    double delta_coeff_table_14mom_T0, delta_coeff_table_14mom_mu0;
    double delta_coeff_table_14mom_dT, delta_coeff_table_14mom_dmu;
    double **deltaf_coeff_tb_14mom_DPi, **deltaf_coeff_tb_14mom_BPi;
    double **deltaf_coeff_tb_14mom_BPitilde;
    double **deltaf_coeff_tb_14mom_BV, **deltaf_coeff_tb_14mom_DV;
    double **deltaf_coeff_tb_14mom_Bpi_shear;

    //quantities for polarization
    double diff_P_rest[NPT][NPHI][NY][Fourindices - 1]; 
    
public:
    freezeout();
    ~freezeout();
    void prepare_data_files(std :: string, std :: string);
    void count_number_of_hypercells();
    void fill_hypersurface_elements();
    void read_particles_info();
    //void test_function();
    void set_momentum_grid_params(double, double, int, int, double, int);
    double PseudoRap(double y, double pt, double m);
    double Rap(double eta, double pt, double m); 



    //functions to perform the phase-space integration 
    double single_point_integration(double pt, double phi, double y, int Flag, int index_mu = 0);
    void  phase_space_distribution_integration();
    void calc_polarization();
    void calc_pol_related_observables();
    double dydeta(double eta, double pt, double m);
   


    //functions to calculate final state observables (setting the benchmark)
    //void calculate_vn_pt_at_y_0(int n);
    //void calculate_pt_spectra_at_y_0();


    //imported from music
    void load_deltaf_qmu_coeff_table(std::string filename);
    void load_deltaf_qmu_coeff_table_14mom(std::string filename);
    double get_deltaf_qmu_coeff(double T, double muB);
    double get_deltaf_coeff_14moments(double T, double muB, double type);

};





#endif
