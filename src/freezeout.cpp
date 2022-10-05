#include "freezeout.h"
#include "util.h"
#include <omp.h>

#define current_max  320
#define hbarc 0.19733

using namespace std;

freezeout::freezeout(){
  load_deltaf_qmu_coeff_table_14mom( "tables/deltaf_coefficients_14moments.dat");
  load_deltaf_qmu_coeff_table("tables/Coefficients_RTA_diffusion.dat");
  int a[4];
  
  //assigning the invariant tensors
  //4-indices levi civita
  for(int ia1 = 0; ia1 < 4; ia1++) 
    for(int ia2 = 0; ia2 < 4; ia2++) 
      for(int ia3 = 0; ia3 < 4; ia3++) 
	for(int ia4 = 0; ia4 < 4; ia4++){
	  a[0] = ia1; a[1] = ia2; a[2] = ia3; a[3] = ia4;
	  int val = 1;
	  for (int i = 0; i < 3; i++){
	    for (int j = 3; j > i; j--){
	      val *= Sgn(a[j] - a[i]);
	    }
	  }
	  levicivita_4indices_sup[ia1][ia2][ia3][ia4] = val;
	}
  //gmunu-tilde (same value for both contravariant and covarinat cases)
  for (int ia1 = 0; ia1 < 4; ia1++) 
    for (int ia2 = 0; ia2 < 4; ia2++){
      gmunu_tilde[ia1][ia2] = 0;
    }
  gmunu_tilde[0][0] =  1; 
  gmunu_tilde[1][1] = -1;
  gmunu_tilde[2][2] = -1;
  gmunu_tilde[3][3] = -1;   
}

freezeout::~freezeout(){
  cout<<"This is the destructor called!"<<endl;
  hypersurface_cells.clear();
}


//imported from MUSIC
void freezeout::load_deltaf_qmu_coeff_table(string filename){
  ifstream table(filename.c_str());
  deltaf_qmu_coeff_table_length_T = 150;
  deltaf_qmu_coeff_table_length_mu = 100;
  delta_qmu_coeff_table_T0 = 0.05;
  delta_qmu_coeff_table_mu0 = 0.0;
  delta_qmu_coeff_table_dT = 0.001;
  delta_qmu_coeff_table_dmu = 0.007892;
  deltaf_qmu_coeff_tb = new double* [deltaf_qmu_coeff_table_length_T];
  for(int i = 0; i < deltaf_qmu_coeff_table_length_T; i++) {
    deltaf_qmu_coeff_tb[i] = new double [deltaf_qmu_coeff_table_length_mu];
  }
  
  double dummy;
  for(int j = 0; j < deltaf_qmu_coeff_table_length_mu; j++){
    for (int i = 0; i < deltaf_qmu_coeff_table_length_T; i++){
      table >> dummy >> dummy >> deltaf_qmu_coeff_tb[i][j];
    }
  }
  table.close();
}

void freezeout::load_deltaf_qmu_coeff_table_14mom(string filename){
  ifstream table(filename.c_str());
  deltaf_coeff_table_14mom_length_T = 190;
  deltaf_coeff_table_14mom_length_mu = 160;
  delta_coeff_table_14mom_T0 = 0.01;
  delta_coeff_table_14mom_mu0 = 0.0;
  delta_coeff_table_14mom_dT = 0.001;
  delta_coeff_table_14mom_dmu = 0.005;
  
  deltaf_coeff_tb_14mom_DPi = new double* [deltaf_coeff_table_14mom_length_T];
  deltaf_coeff_tb_14mom_BPi = new double* [deltaf_coeff_table_14mom_length_T];
  deltaf_coeff_tb_14mom_BPitilde = 
    new double* [deltaf_coeff_table_14mom_length_T];
  deltaf_coeff_tb_14mom_DV = new double* [deltaf_coeff_table_14mom_length_T];
  deltaf_coeff_tb_14mom_BV = new double* [deltaf_coeff_table_14mom_length_T];
  deltaf_coeff_tb_14mom_Bpi_shear =
    new double* [deltaf_coeff_table_14mom_length_T];
  for (int i = 0; i < deltaf_coeff_table_14mom_length_T; i++){
    deltaf_coeff_tb_14mom_DPi[i] =
      new double[deltaf_coeff_table_14mom_length_mu];
    deltaf_coeff_tb_14mom_BPi[i] =
      new double[deltaf_coeff_table_14mom_length_mu];
    deltaf_coeff_tb_14mom_BPitilde[i] =
      new double[deltaf_coeff_table_14mom_length_mu];
    deltaf_coeff_tb_14mom_DV[i] =
      new double[deltaf_coeff_table_14mom_length_mu];
    deltaf_coeff_tb_14mom_BV[i] =
      new double[deltaf_coeff_table_14mom_length_mu];
    deltaf_coeff_tb_14mom_Bpi_shear[i] =
      new double[deltaf_coeff_table_14mom_length_mu];
  }
  
  double dummy;
  for(int i = 0; i < deltaf_coeff_table_14mom_length_T; i++){
    for(int j = 0; j < deltaf_coeff_table_14mom_length_mu; j++){
      table >> dummy >> dummy >> deltaf_coeff_tb_14mom_DPi[i][j]
	    >> deltaf_coeff_tb_14mom_BPi[i][j] 
	    >> deltaf_coeff_tb_14mom_BPitilde[i][j]
	    >> deltaf_coeff_tb_14mom_DV[i][j] 
	    >> deltaf_coeff_tb_14mom_BV[i][j]
	    >> deltaf_coeff_tb_14mom_Bpi_shear[i][j];
    }
  }
  table.close();
  
  // convert units
  double hbarc3 = hbarc*hbarc*hbarc;
  double hbarc4 = hbarc3*hbarc;
  for(int i = 0; i < deltaf_coeff_table_14mom_length_T; i++){
    for(int j = 0; j < deltaf_coeff_table_14mom_length_mu; j++){
      deltaf_coeff_tb_14mom_DPi[i][j] *= hbarc4;          // fm^4/GeV
      deltaf_coeff_tb_14mom_BPi[i][j] *= hbarc4;          // fm^4/(GeV^2)
      deltaf_coeff_tb_14mom_BPitilde[i][j] *= hbarc4;     // fm^4/(GeV^2)
      deltaf_coeff_tb_14mom_DV[i][j] *= hbarc3;           // fm^3/GeV
      deltaf_coeff_tb_14mom_BV[i][j] *= hbarc3;           // fm^3/(GeV^2)
      deltaf_coeff_tb_14mom_Bpi_shear[i][j] *= hbarc4;    // fm^4/(GeV^2)
    }
  }
}


double freezeout::get_deltaf_qmu_coeff(double T, double muB){
  int idx_T = static_cast<int>((T - delta_qmu_coeff_table_T0)
			       /delta_qmu_coeff_table_dT);
  int idx_mu = static_cast<int>((muB - delta_qmu_coeff_table_mu0)
				/delta_qmu_coeff_table_dmu);
  double x_fraction = ((T - delta_qmu_coeff_table_T0)
		       /delta_qmu_coeff_table_dT - idx_T);
  double y_fraction = ((muB - delta_qmu_coeff_table_mu0)
		       /delta_qmu_coeff_table_dmu - idx_mu);
  
  // avoid overflow: return a large number so that delta f = 0
  if (idx_mu > deltaf_qmu_coeff_table_length_mu - 2){
    return(1e30);
  }
  if (idx_T > deltaf_qmu_coeff_table_length_T - 2){
    return(1e30);
  }
  
  // avoid underflow: return a large number so that delta f = 0
  if (idx_mu < 0){
    return(1e30);
  }
  if (idx_T < 0) {
    return(1e30);
  }
  
  double f1 = deltaf_qmu_coeff_tb[idx_T][idx_mu];
  double f2 = deltaf_qmu_coeff_tb[idx_T][idx_mu+1];
  double f3 = deltaf_qmu_coeff_tb[idx_T+1][idx_mu+1];
  double f4 = deltaf_qmu_coeff_tb[idx_T+1][idx_mu];
  
  double coeff = f1*(1. - x_fraction)*(1. - y_fraction) 
    + f2*(1. - x_fraction)*y_fraction
    + f3*x_fraction*y_fraction
    + f4*x_fraction*(1. - y_fraction);
  return(coeff);
}



double freezeout::get_deltaf_coeff_14moments(double T, double muB, double type){
  int idx_T = static_cast<int>((T - delta_coeff_table_14mom_T0)
			       /delta_coeff_table_14mom_dT);
  int idx_mu = static_cast<int>((muB - delta_coeff_table_14mom_mu0)
				/delta_coeff_table_14mom_dmu);
  double x_fraction = ((T - delta_coeff_table_14mom_T0)
		       /delta_coeff_table_14mom_dT - idx_T);
  double y_fraction = ((muB - delta_coeff_table_14mom_mu0)
		       /delta_coeff_table_14mom_dmu - idx_mu);
  
  double **deltaf_table = NULL;
  if(type == 0){
    deltaf_table = deltaf_coeff_tb_14mom_DPi;
  } else if(type == 1){
    deltaf_table = deltaf_coeff_tb_14mom_BPi;
  } else if(type == 2){
    deltaf_table = deltaf_coeff_tb_14mom_BPitilde;
  } else if(type == 3){
    deltaf_table = deltaf_coeff_tb_14mom_DV;
  } else if(type == 4){
    deltaf_table = deltaf_coeff_tb_14mom_BV;
  } else if(type == 5){
    deltaf_table = deltaf_coeff_tb_14mom_Bpi_shear;
  } else{
    cout<< ("error")<<endl;
       exit(-1);
  }
  
  double f1 = deltaf_table[idx_T][idx_mu];
  double f2 = deltaf_table[idx_T][idx_mu+1];
  double f3 = deltaf_table[idx_T+1][idx_mu+1];
  double f4 = deltaf_table[idx_T+1][idx_mu];
  
  double coeff = f1*(1. - x_fraction)*(1. - y_fraction) 
    + f2*(1. - x_fraction)*y_fraction
    + f3*x_fraction*y_fraction
    + f4*x_fraction*(1. - y_fraction);
  return(coeff);
}

//my addition starts
void freezeout::prepare_data_files(string _hypersurface_dat_file, string expected_particles){
  surface_dat_file.str(_hypersurface_dat_file);
  pdg_dat_file.str("pdg/pdg05.dat");
  input_particles.str(expected_particles);
  count_number_of_hypercells();
  fill_hypersurface_elements();
  read_particles_info();
}

void freezeout::set_momentum_grid_params(double _ptmin, double _ptmax, int n_pt, int n_phi, double y, int n_y){
  pt_min = _ptmin, pt_max = _ptmax, N_steps_pt = n_pt;
  if(N_steps_pt > NPT) {
    cout<< "Fatal Error!!! Alloted size increased in Pt Bin "<<endl;
    exit(-1);
  }
  N_steps_phi =  n_phi;
  if(N_steps_phi > NPHI) {
    cout<< "Fatal Error!!! Alloted size increased in Phi Bin "<<endl;
    exit(-1);
  }
  dphi = (phi_max - phi_min)/(static_cast<double>(N_steps_phi) - 1);
  y_max = y, N_steps_y = n_y;
  if(N_steps_y > NY) {
    cout<< "Fatal Error!!! Alloted size increased in y Bin"<<endl;
    exit(-1);
  }
  dy = (2*y_max)/(static_cast<double>(N_steps_y - 1));
}


void freezeout::count_number_of_hypercells(){
  int count = 0;   
  std::string line;
  std::ifstream readsurface(surface_dat_file.str());
  if(readsurface.is_open()){
    while (!readsurface.eof()){
      std:: getline(readsurface, line);
      count++;
    }
  }
  else{
    std :: cout << "[warning] Error in reading the hypersurface file...exiting!" << std :: endl;
    exit(1);
  }
  readsurface.close();
  N_hypercells = count;
  cout << "Total Number of freezeout hypercells: "<< N_hypercells <<endl;
}   


void freezeout::fill_hypersurface_elements(){
  std::ifstream readsurface(surface_dat_file.str());
  int i = 0;
  SurfaceElement temp_cell;
  cout << "Reading the hypersurface file: " << surface_dat_file.str() << endl;
  while (i < N_hypercells) {
    {   // position in (tau, x, y, eta)
      readsurface >> temp_cell.x[0] >> temp_cell.x[1]
		  >> temp_cell.x[2] >> temp_cell.x[3];
      
      // hypersurface vector in (tau, x, y, eta)
      readsurface >> temp_cell.s[0] >> temp_cell.s[1]
		  >> temp_cell.s[2] >> temp_cell.s[3];
      
      // flow velocity in (tau, x, y, eta)
      readsurface >> temp_cell.u[0] >> temp_cell.u[1]
		  >> temp_cell.u[2] >> temp_cell.u[3];
      
      readsurface >> temp_cell.epsilon_f >> temp_cell.T_f
		  >> temp_cell.mu_B >> temp_cell.mu_S >> temp_cell.mu_C
		  >> temp_cell.eps_plus_p_over_T_FO;
      
      // freeze-out Wmunupseudofreeze 1
      readsurface >> temp_cell.W[0][0] >> temp_cell.W[0][1]
		  >> temp_cell.W[0][2] >> temp_cell.W[0][3]
		  >> temp_cell.W[1][1] >> temp_cell.W[1][2]
		  >> temp_cell.W[1][3] >> temp_cell.W[2][2]
		  >> temp_cell.W[2][3] >> temp_cell.W[3][3];
      readsurface >> temp_cell.rho_B;
      //bulk-viscous effect is off, hence no Pi 
      //baryon-diffusion is off
      //readsurface >> temp_cell.q[0] >> temp_cell.q[1] >> temp_cell.q[2] >> temp_cell.q[3];
      
      // freeze-out Gradients
      for (int mu = 0; mu < 4; mu++)
	readsurface >> temp_cell.d_mu_betanu[mu][0]  >> temp_cell.d_mu_betanu[mu][1]
		    >> temp_cell.d_mu_betanu[mu][2]  >> temp_cell.d_mu_betanu[mu][3];         
      
    }
    temp_cell.sinh_eta_s = sinh(temp_cell.x[3]);
    temp_cell.cosh_eta_s = cosh(temp_cell.x[3]);
    
    hypersurface_cells.push_back(temp_cell);
    i++;
  }
  readsurface.close();
  cout << "Hypersurface file read properly...proceeding further!" << endl;
}

/*void freezeout :: test_function(){
  ofstream write("test.dat");
  int i = 0;
  while (i< N_hypercells)
  {
  if(i%10 == 0)
  write << hypersurface_cells[i].x[0] <<  "\t\t" << hypersurface_cells[i].x[2] <<
  "\t\t" << hypersurface_cells[i].d_mu_betanu[0][0] <<
  "\t\t" <<  hypersurface_cells[i].T_f * hbarc << endl;
  i++;
  }
  write.close();
  
  }*/


void freezeout::read_particles_info()
{   
  std :: cout << "Started reading the particle data group and user input... " << std:: endl;
  //reading the particle IDs from the user input file 
  ifstream read_inputPIDS(input_particles.str());
  if (!read_inputPIDS.is_open()){
    cout << "error in reading.Exiting...."<<endl;
    exit(1);
  }
  
  int temp_PID, i = 0;
  while (!read_inputPIDS.eof()){
    read_inputPIDS >>temp_PID;
    PIDS.push_back(temp_PID);
  }
  read_inputPIDS.close();
  
  //reading and storing particles from particle data group
  Particle_list = (Particle *)malloc(current_max*sizeof(Particle));
  FILE *p_file;
  p_file = fopen(pdg_dat_file.str().c_str(), "r");
  int j = 0;
  
  while (i <  current_max){ 
    int temp;
    temp = fscanf(p_file, "%d",  &Particle_list[i].number );
    temp = fscanf(p_file, "%s",   Particle_list[i].name   );
    temp = fscanf(p_file, "%lf", &Particle_list[i].mass   );
    temp = fscanf(p_file, "%lf", &Particle_list[i].width  );
    temp = fscanf(p_file, "%d",  &Particle_list[i].degeneracy);
    temp = fscanf(p_file, "%d",  &Particle_list[i].baryon );
    temp = fscanf(p_file, "%d",  &Particle_list[i].strange);
    temp = fscanf(p_file, "%d",  &Particle_list[i].charm  );
    temp = fscanf(p_file, "%d",  &Particle_list[i].bottom );
    temp = fscanf(p_file, "%d",  &Particle_list[i].isospin);
    temp = fscanf(p_file, "%lf", &Particle_list[i].charge );
    temp = fscanf(p_file, "%d",  &Particle_list[i].decays );   
    
    int h;
    for(int k = 0; k < Particle_list[i].decays; k++){
      h = fscanf(p_file, "%i%i%lf%i%i%i%i%i",
		 &decay[j].reso, &decay[j].numpart, &decay[j].branch, 
		 &decay[j].part[0], &decay[j].part[1], &decay[j].part[2],
		 &decay[j].part[3], &decay[j].part[4]);
      if (h != 8){
	printf("Error in scanf decay \n");
	exit(0);
      }
      j++;
    }
    i++;
  }
  fclose(p_file);
  
  
  //storing only the particles needed for the analysis
  i = 0 , j = 0;
  my_Particle_list = (Particle *)malloc((PIDS.size())*sizeof(Particle));
  for( i = 0; i < PIDS.size(); i++){
    for( j = 0; j < current_max; j++){
      if(Particle_list[j].number == PIDS[i]){
	my_Particle_list[i].number = Particle_list[j].number ;
	strncpy(my_Particle_list[i].name, Particle_list[j].name, 50);
	my_Particle_list[i].mass = Particle_list[j].mass;
	my_Particle_list[i].width = Particle_list[j].width;
	my_Particle_list[i].degeneracy = Particle_list[j].degeneracy;
	my_Particle_list[i].baryon = Particle_list[j].baryon;
	my_Particle_list[i].strange = Particle_list[j].strange;
	my_Particle_list[i].charm= Particle_list[j].charm;
	my_Particle_list[i].bottom = Particle_list[j].bottom;
	my_Particle_list[i].isospin = Particle_list[j].isospin;
	my_Particle_list[i].charge = Particle_list[j].charge;
	my_Particle_list[i].decays = Particle_list[j].decays;
      }
    }
  }
  cout << "Particle Information read successfully, no. of input: " << PIDS.size() << endl;
  cout << "Particle Name: " << "\t\t" << "Charge: " << endl;
  
  for (int i = 0; i < PIDS.size() ; i++)
    cout << my_Particle_list[i].name <<  "\t\t" << my_Particle_list[i].charge <<endl;
}



double freezeout::single_point_integration(double pt, double phi, double y, int Flag, double* MSP){

  // MSP : Mean Spin Vector
  for(int ii=0; ii<4; ii++){
      MSP[ii] = 0. ; 
  }

  double y_minus_eta_cut = 10.0;
  double total_sum[4] = {0.,0.,0.,0.};
  
  //check here
  double prefactor;
  
  //summing over the hypercells using openmp
#pragma omp parallel 
  {
    
    double temp_sum[4] = {0.,0.,0.,0.};                    //thread local variable
    for(int ii=0; ii<4; ii++){
       temp_sum[ii] = 0. ; 
    }

    
#pragma omp for  
    for (int icell = 0; icell < N_hypercells; icell++) {
      double tau        = hypersurface_cells[icell].x[0];
      double eta_s      = hypersurface_cells[icell].x[3];
      if (fabs(y - eta_s) < y_minus_eta_cut ){
	double cosh_eta_s = hypersurface_cells[icell].cosh_eta_s;
	double sinh_eta_s = hypersurface_cells[icell].sinh_eta_s;
	double T   = hypersurface_cells[icell].T_f  * hbarc;              // GeV
	double muB = hypersurface_cells[icell].mu_B * hbarc;              // GeV
	double muC = hypersurface_cells[icell].mu_C * hbarc;              // GeV
	double muS = hypersurface_cells[icell].mu_S * hbarc;              // GeV
	double mu  = particle_of_interest.baryon * muB
	  + particle_of_interest.charge * muC
	  + particle_of_interest.strange * muS;           // GeV
	double sigma_mu[4]={0.,0.,0.,0.};
	double u_flow[4]={0.,0.,0.,0.};
	for(int ii = 0; ii < 4; ii++){
	  sigma_mu[ii] = hypersurface_cells[icell].s[ii];
	  u_flow[ii] = hypersurface_cells[icell].u[ii];
	}
	       

	//check the beautiful interplay of the units
	double eps_plus_P_over_T = hypersurface_cells[icell].eps_plus_p_over_T_FO;  // fm^{-3}
	double prefactor_shear = 1./(2.*eps_plus_P_over_T*T*T*T)*hbarc;             // fm^{4}/GeV^{2}

	
	//energy-momemtum 4-vector in Milne-tilde coordinates 
	double mt = sqrt(pow(particle_of_interest.mass, 2) + pow(pt, 2));    
	double ptau = mt*(cosh(y)*cosh_eta_s
			  - sinh(y)*sinh_eta_s); 
	double peta = mt*(sinh(y)*cosh_eta_s
			  - cosh(y)*sinh_eta_s);
	double px  = pt*cos(phi);
	double py  = pt*sin(phi);
	double sum = 0. ;
	
	// compute p^mu*dSigma_mu [fm^3*GeV]
	double pdSigma = tau*(ptau*sigma_mu[0]
			      + px*sigma_mu[1]
			      + py*sigma_mu[2]
			      + peta/tau*sigma_mu[3]);
	//note the structure of g_munu                                
	double E = (ptau*u_flow[0] - px*u_flow[1] - py*u_flow[2] - peta*u_flow[3]); //GeV
	
	int sign; //(whether boson or fermion)
	if (particle_of_interest.baryon == 0) sign = -1.;
	else  sign = 1.;
	
	// this is the equilibrium f, f_0:
	double f = 1./(exp(1./T*(E - mu)) + sign);
        if( std::isinf(f) || std::isnan(f) ){
            std::cout << "E = " << E << ",  T =" << T << std::endl ; 
            std::cout << "f is inf/nan ... " << std::endl ;
             exit(1); 
        }        


	double sum0 = 0;
	double sum1 = 0;
	double sum2 = 0;
	double sum3 = 0;

	//the final integrand : put the flags here 
	switch (Flag){

	case InvYieldNoVisc:
	  {   
	    prefactor = 1.0;
	    sum0 = (prefactor * f) * pdSigma ;
             // Don't delete, although redundant
	    sum1 = (prefactor * f) * pdSigma ;
	    sum2 = (prefactor * f) * pdSigma ;
	    sum3 = (prefactor * f) * pdSigma ;

            if( std::isinf(sum) || std::isnan(sum) ){       
                std::cout << "ptau = " << ptau << ", px = " << px << ", py = " << py << ",  peta =" << peta << std::endl ; 
                std::cout << "stau = " << sigma_mu[0] << ", sx = " << sigma_mu[1] << ", sy = " << sigma_mu[2] << ",  seta =" << sigma_mu[3] << std::endl ;
                std::cout << "tau = " << tau << std::endl ;  
                std::cout << "f = " << f << ",  pdSigma =" << pdSigma << std::endl ; 
                std::cout << "Sum for InvYieldNoVisc is inf/nan ... " << std::endl ;
                exit(1); 
            }   
	    break;
	  }


	// ========================================= //
	// Flag to calculate Local Spin Polarization //
	// ========================================= //
	case NumSLab:
	  {       
	    prefactor = 1.0;
	    double d_mu_betanu[4][4], d_mu_beta_nu[4][4];
	    double omega_munu[4][4];
	    double p_mu[4];
	    p_mu[0] = ptau;
	    p_mu[1] = -px;
	    p_mu[2] = -py;
	    p_mu[3] = -peta;
	    
	    for (int mu = 0; mu < 4; mu++) 
	      for (int nu = 0; nu < 4; nu++){
		d_mu_betanu[mu][nu] = 0.0;
	      }
	    
	    //from the hypersurface, one obtains \partial_{\mu} \beta^{\nu}
	    for (int mu = 0; mu < 4; mu++) 
	      for (int nu = 0; nu < 4; nu++){
		d_mu_betanu[mu][nu] = hypersurface_cells[icell].d_mu_betanu[mu][nu];
	      }
	    

            // Conversion of \partial_{\mu} \beta^{\nu}  to \partial_{\mu} \beta_{\nu}
	    // d_\mu \beta_{\nu} = d_\mu \beta^{\rho} * g_{\rho\nu}
            // g_{\mu\nu} = g^{\mu\nu}
	    for (int mu = 0; mu < 4; mu++) 
	      for (int nu = 0; nu < 4; nu++) {
		double sum_murho = 0.0;
		for (int rho = 0; rho < 4; rho++){
		  sum_murho += d_mu_betanu[mu][rho]*gmunu_tilde[rho][nu]; 
		}
		d_mu_beta_nu[mu][nu] = sum_murho;
	      }
	    
	    
            // Calculation of \omega_{\mu \nu}
            // The form of \omega_{\mu}{\nu} is in Eq. (7) of arXiv:2011.03740
            // Now the calculated \omega_{\mu}{\nu} is in Milne co-ordinate           
	    for (int mu = 0; mu < 4; mu++) 
	      for (int nu = 0; nu < 4; nu++){
		omega_munu[mu][nu] = - 0.5 * (d_mu_beta_nu[mu][nu] - d_mu_beta_nu[nu][mu]);
	      }

            
           // Calculation of S^{\mu} = [1-f(x,p)] \epsilon^{\mu \nu \rho \sigma} p_{\sigma} \omega_{\nu \rho}
           // S^{\mu} defined here is same as in Eq. (6) of arXiv:2011.03740 but without the term -1/(2*m)*(S*(S+1))/3
           // The term -1/(2*m)*(S*(S+1))/3 will be incorporated in the calc_polarization() function.
           double S[4] = {0.,0.,0.,0.} ; 
           for(int imu=0; imu<4; imu++){
             S[imu] = 0. ;
             for(int inu=0; inu<4; inu++){
               for(int irho=0; irho<4; irho++){
                 for(int isigma=0; isigma<4; isigma++){
                    S[imu] += (1.0 -(prefactor * f))*
                      levicivita_4indices_sup[imu][inu][irho][isigma]*
                        p_mu[isigma] * omega_munu[inu][irho] ;
                 }
               }
             }
           }


           // Transformation of S^{\mu} from Milne to cartesian co-ordinate.
           double auxiliaryS[4] = {0.,0.,0.,0.} ; 
	    for(int ii = 0; ii < 4; ii++){
	       auxiliaryS[ii] = S[ii];   
            }
            S[1] = auxiliaryS[1] ;
            S[2] = auxiliaryS[2] ;
            S[0] = cosh_eta_s * auxiliaryS[0] + sinh_eta_s * auxiliaryS[3] ; 
            S[3] = sinh_eta_s * auxiliaryS[0] + cosh_eta_s * auxiliaryS[3] ;

            // numerator of Eq.(9) in arXiv:2011.03740
	    sum0 = pdSigma * (prefactor*f) * S[0];   
	    sum1 = pdSigma * (prefactor*f) * S[1];   
	    sum2 = pdSigma * (prefactor*f) * S[2];   
	    sum3 = pdSigma * (prefactor*f) * S[3];   
            	      
	    break;
	  }
          
	default:
	  std :: cout << "No matching flag was found! exiting ......" << std :: endl;
	  exit(1);
	  break;
	}

        
        // checks for unphysical instances. //
	if(sum > 10000){
	  std::cout << "sum>10000 in summation. sum = " << sum << ", f = " <<
                        f << ", pdSigma = " <<
                        pdSigma << ", T=" << T  << ", E = " << E << ", mu = " << mu;
	}
	if(f < 0.){
	  std::cout << " f_eq < 0.! f_eq = " << f << ", T = " << T << " GeV, mu = " <<
                        mu << " GeV, E = " << E << " GeV";
	}

	//integrating             
	temp_sum[0] += sum0;
	temp_sum[1] += sum1;
	temp_sum[2] += sum2;
	temp_sum[3] += sum3;

      } // if (y-etas) < y_minus_eta_cut
    } // loop over cell
#pragma omp critical
    {
      for(int ii=0; ii<4; ii++){
        total_sum[ii] += temp_sum[ii] ; 
      }
    }   // critical sum
    
  } // openpm thread sharing
  

  // required to group and deliver the mean spin polarisation vectors
  MSP[0] = total_sum[0];
  MSP[1] = total_sum[1];
  MSP[2] = total_sum[2];
  MSP[3] = total_sum[3];

  // required for denominator 
  return total_sum[0] ;
  
}


/*
void freezeout::phase_space_distribution_integration(){
  
  std:: ofstream writefile("yptphi_sprectra_all_particles.dat");
  double a_pt, b_phi, c_y, tempc_y;
  
  for (int a = 0; a < N_steps_pt; a++){
    a_pt = (pt_min+ (pt_max - pt_min)
	    *pow(static_cast<double>(a), 2.)
	    /pow(static_cast<double>(N_steps_pt - 1), 2.));
    for (int b = 0; b < N_steps_phi; b++){
      b_phi = b * dphi;
      for (int c = 0; c < N_steps_y; c++){
	tempc_y = - y_max + c*dy;
	writefile << tempc_y << " " << a_pt << " " << b_phi << " ";
	for (int i = 0; i < PIDS.size(); i++){
	  double result = 0;
	  particle_of_interest.baryon = my_Particle_list[i].baryon;
	  particle_of_interest.mass = my_Particle_list[i].mass;
	  particle_of_interest.number = my_Particle_list[i].number;
	  particle_of_interest.degeneracy = my_Particle_list[i].degeneracy;
          
	  c_y = Rap(tempc_y, a_pt, particle_of_interest.mass);
	  result = single_point_integration(a_pt, b_phi, c_y, MusicSameParam);
	  writefile << result  << " "; 
	  my_Particle_list[i].dNdydptdphi[c][a][b] = result;
	}
	writefile << endl;
      }
    }
  }
  writefile.close();
}

*/

void freezeout::calc_polarization(){
  std::cout << "Calculating polarization ...  " << std::endl;
  double eta_cut = 4.0;                                                               
  double a_pt, b_phi, c_y, tempc_y;
  //spin-of the particle is hard coded for time-being
  double spin = 0.5;            
  std::ofstream writefile("yptphi_Prest_all_particles.dat");
  
  for(int a = 0; a < N_steps_pt; a++){
    a_pt = (pt_min+ (pt_max - pt_min)
	    *pow(static_cast<double>(a), 2.)
	    / pow(static_cast<double>(N_steps_pt - 1), 2.));
    for (int b = 0; b < N_steps_phi; b++){
      b_phi = b * dphi;
      for (int c = 0; c < N_steps_y; c++){
	tempc_y = -y_max + c*dy;   //pseudorapidity
	// pseudorapidity cut  
	if (fabs(tempc_y) > eta_cut ) continue;
	writefile << tempc_y << " " << a_pt << " " << b_phi << " "; 
        //conversion to rapidity
	c_y = Rap(tempc_y, a_pt, particle_of_interest.mass); 
	
	for (int i = 0; i < PIDS.size(); i++){   //particle loop  
	  double val_num   = 0 ;
          double val_denom = 0 ;
	  double S_lab[4] = {0.,0.,0.,0.};
	  particle_of_interest.baryon     =  my_Particle_list[i].baryon;
	  particle_of_interest.mass       =  my_Particle_list[i].mass;
	  particle_of_interest.number     =  my_Particle_list[i].number;
	  particle_of_interest.degeneracy =  my_Particle_list[i].degeneracy;
	  double mass_h = particle_of_interest.mass;
	  
          double Wrapper[4] = {0.,0.,0.,0.} ; 

	  val_denom = single_point_integration(a_pt, b_phi, c_y, InvYieldNoVisc, Wrapper);
          val_num   = single_point_integration(a_pt, b_phi, c_y, NumSLab, Wrapper);


          // The term -1/(2*m)*(S*(S+1))/3 is the scalar co-efficient of S^{\mu} in Eq.(6) of arXiv:2011.03740.
          // Again another S is devided to get the average polarization in Lab frame.    
	  for (int mu = 0; mu < 4; mu++){
            S_lab[mu] = (-1./(2.*mass_h)) * (spin+1.) / 3. * Wrapper[mu] / val_denom ;
            if( std::isinf(S_lab[mu]) || std::isnan(S_lab[mu]) ){
              std::cout << "================ inside polarization calculation ===============" << std::endl ; 
              std::cout << "Numertor = " << val_num << ",  denominator =" << val_denom << std::endl ; 
              std::cout << "S_Lab is inf/nan ... " << std::endl ;
              exit(1); 
            }

          }
          

          // ================================================ //
	  // Lorentz transformation to particle's rest frame. //
          // ================================================ //

          // PolLab[3] contains 3 components of contravariant polarization tensor,  P^{x}, P^{y} and P^{z}.
	  double PolLab[3], MomentumLab[3]; 
	  PolLab[0]  = S_lab[1];
	  PolLab[1]  = S_lab[2];
	  //PolLab[2]  = S_lab[0]*sinh(c_y) + S_lab[3]*cosh(c_y);
	  PolLab[2]  = S_lab[3];
	  
          // MomentumLab[3] contains 3 components of contravariant p^{\mu} tensor, p^{x}, p^{y}, p^{z}.
	  MomentumLab[0] = a_pt*cos(b_phi);
	  MomentumLab[1] = a_pt*sin(b_phi);
	  MomentumLab[2] = pow (pow(a_pt, 2) + pow(particle_of_interest.mass, 2), 0.5) * sinh(c_y);
	  
	  double mod_p_sq = 0 ; // p^{2}
          double inner_Pp = 0 ; // inner product of particle momentum and polarization (p.P) 
          double E        = 0 ; // energy at specified momentum.
	  for (int i = 0; i < 3; i++){
	    mod_p_sq += pow(MomentumLab[i], 2);
	    inner_Pp += PolLab[i] * MomentumLab[i];
	  }
	  E = sqrt (pow(mass_h, 2) + mod_p_sq);

          // Storring space like components of Polarization(in particle rest frame) in diff_P_rest[pT][phi][eta][0]
          // diff_P_rest[pT][phi][eta][1] and diff_P_rest[pT][phi][eta][2].
          // Not storing diff_P_rest[pT][phi][y][t] = 0 
	  for (int i_mu = 1; i_mu < 4; i_mu++){
            // Eq.(40) of arXiv:2106.08125 or Eq.(11) of arXiv:2011.03740
	    double result =  PolLab[i_mu-1] - (inner_Pp)/(E * (E + mass_h))*MomentumLab[i_mu-1];
	    writefile << result  << " "; 
	    diff_P_rest[a][b][c][i_mu-1] = result;
	  }
	} // loop over PID
	writefile << endl;


      } // rapidity loop
    } // phi loop
    std::cout << "calculated upto pt: " << a_pt << std::endl; 
  } // pT loop
  writefile.close();
}


// returns pseudorapidity for a given rapidity, transverse momentum, and mass
double freezeout::PseudoRap(double y, double pt, double m){
  double eta = acosh(2*m*m/pt/pt*sinh(y)*sinh(y) + cosh(2*y))/2.;
  if (y < 0)
    eta *= -1.;
  return eta;
}


// returns rapidity for a given pseudorapidity, transverse momentum, and mass
double freezeout::Rap(double eta, double pt, double m) {
  double y = log((sqrt(m*m + pt*pt*cosh(eta)*cosh(eta)) + pt*sinh(eta))
		 /sqrt(m*m+pt*pt));
  return y;
} 

