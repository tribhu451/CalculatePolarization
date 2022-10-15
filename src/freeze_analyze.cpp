#include "freezeout.h"
#include <omp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


#define M_PI_ 3.1415927

void freezeout::calc_pol_related_observables(){
   double minrap = -1.0, maxrap = 1.0; 
    double minpt = 0.5, maxpt = 3.0; 
   std::stringstream output_filename;
   output_filename.str("");
   output_filename << "MinusPy_percent_vs_pT_" << particle_of_interest.number <<"_ycut_";
   output_filename << minrap << "_" << maxrap ;
   output_filename << ".dat";
   std:: ofstream writefile(output_filename.str().c_str());
   output_filename.str("");
   output_filename << "MinusPy_percent_" << particle_of_interest.number <<"_ycut_";
   output_filename << minrap << "_" << maxrap << "_ptcut_" ;
   output_filename << minpt << "_" << maxpt << ".dat";
   std:: ofstream writefile_all_integ(output_filename.str().c_str());

  //std:: ofstream writefile3("test_gsl_coarse.dat");
  std :: cout << "Calculating Polarization observables for the PID " << particle_of_interest.number << " ..." << std ::endl; 
  double pt_grid[NPT], pt_diff_P_rest[NPT], denominator_pt_diff[NPT];

 
  for (int i = 0; i < N_steps_pt; i++){
    pt_grid[i] = (pt_min + (pt_max - pt_min)*pow(static_cast<double>(i), 2.)
		  /pow(static_cast<double>(N_steps_pt - 1), 2.));
  }
  

  for(int i = 0; i < PIDS.size(); i++){  
    double mass_h = my_Particle_list[i].mass;
    for(int ipt = 0; ipt < N_steps_pt; ipt++){
      double sum = 0;
      double sum_denom = 0;
      double pt = pt_grid[ipt];
      // Integrate over phi using trapezoid rule
      for(int iphi = 0; iphi < N_steps_phi; iphi++){
	double ylist[NY] = {0};
	double etalist[NY] = {0};
	
	// Integrate over pseudorapidity using gsl
	double Py_rapidity_grid[NY] = {0};
	double denom_rapidity_grid[NY] = {0};
	for(int ieta = 0; ieta < N_steps_y; ieta++){
	  double rap_local = - y_max + ieta*dy;
	  double eta_local, y_local;
          
	  {
	    eta_local = rap_local;
	    etalist[ieta] = eta_local;
	    y_local = Rap(eta_local, pt, mass_h);    // get rapidity
	    ylist[ieta] = y_local;                  // get rapidity
	  } 
          
	  Py_rapidity_grid[ieta]    = pt * pt * cosh(eta_local)/sqrt(mass_h*mass_h + pt*pt*cosh(eta_local)*cosh(eta_local)) * diff_P_rest[ipt][iphi][ieta][2];
	  denom_rapidity_grid[ieta] = pt * pt * cosh(eta_local)/sqrt(mass_h*mass_h + pt*pt*cosh(eta_local)*cosh(eta_local)) * denom_dndydptdphi[ipt][iphi][ieta];
	  //writefile3 << ylist[ieta] << " " << Py_rapidity_grid[ieta] << std :: endl;
	  //Py_rapidity_grid[ieta] = diff_P_rest[ipt][iphi][ieta][1];
	}
	//writefile3.close();
	
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *spline    = gsl_spline_alloc (gsl_interp_linear, N_steps_y);
	gsl_spline_init(spline, etalist, Py_rapidity_grid, N_steps_y);
	
	gsl_interp_accel *accd = gsl_interp_accel_alloc ();
	gsl_spline *spline_denom = gsl_spline_alloc (gsl_interp_linear, N_steps_y);
	gsl_spline_init(spline_denom, etalist, denom_rapidity_grid, N_steps_y);
	
	/*
	  double test_min =  ylist[0];
	  double test_max =  ylist[N_steps_y-1];
	  std :: cout << test_min << "----->" << test_max << std :: endl;
	  
	  const int test_array_len =  10000;
	  double test_deta = (test_max - test_min)/(double)test_array_len;
	  double local_x, local_x_temp;
	  
	  std:: ofstream writefile2("test_gsl.dat");
	  for (int itest = 0; itest < test_array_len; itest++)
	  {
	  local_x_temp = test_min + itest*test_deta;
	  writefile2 << local_x_temp << " " << gsl_spline_eval(spline, local_x_temp, acc ) << std ::endl ;
	  }
	  writefile2.close();
	  exit(1);*/
        
	double Py_pseudorapidity_integrated;
	double denom_pseudorapidity_integrated;
	double normalize_pseudo_rap = fabs(maxrap - minrap);    
	// integrated from minrap to maxrap
	Py_pseudorapidity_integrated = gsl_spline_eval_integ(spline, minrap, maxrap, acc) ;
	denom_pseudorapidity_integrated = gsl_spline_eval_integ(spline_denom, minrap, maxrap, accd) ;
	//trapezoid
	sum += Py_pseudorapidity_integrated*(2. * M_PI_)/(N_steps_phi-1);
	sum_denom += denom_pseudorapidity_integrated*(2. * M_PI_)/(N_steps_phi-1);
	gsl_spline_free (spline);
	gsl_spline_free (spline_denom);
	gsl_interp_accel_free (acc);
	gsl_interp_accel_free (accd);
        
      }// phi loop
      pt_diff_P_rest[ipt] = sum; 
      denominator_pt_diff[ipt] = sum_denom ;    
      writefile << pt << "  " << - sum*(100) << "  " << sum_denom << "  " << - sum*(100) / sum_denom  << std :: endl;    
    } //pt-loop  


    
    gsl_interp_accel *acc2 = gsl_interp_accel_alloc ();
    gsl_spline *spline2 = gsl_spline_alloc (gsl_interp_linear, N_steps_pt );
    gsl_spline_init (spline2, pt_grid, pt_diff_P_rest, N_steps_pt);
    gsl_spline *spline2_denom = gsl_spline_alloc (gsl_interp_linear, N_steps_pt );
    gsl_spline_init (spline2_denom, pt_grid, denominator_pt_diff, N_steps_pt);
    
    double Py_total_integrated, denom_total_integrated ;
    double normalize_pt = fabs(maxpt - minpt);
    Py_total_integrated = gsl_spline_eval_integ(spline2, minpt, maxpt, acc2);
    denom_total_integrated = gsl_spline_eval_integ(spline2_denom, minpt, maxpt, acc2);
    writefile_all_integ << denom_total_integrated << "  " << -100*Py_total_integrated / denom_total_integrated << std::endl;
    
    gsl_spline_free (spline2);
    gsl_spline_free (spline2_denom);
    gsl_interp_accel_free (acc2);
      
  } //particle loop
  
}





void freezeout::calculate_pseudo_rapidity_differential_polarization(){

   double minpt = 0.5, maxpt = 3.0; 
   std::stringstream output_filename;
   output_filename.str("");
   output_filename << "MinusPy_percent_vs_eta_" << particle_of_interest.number <<"_ptcut_";
   output_filename << minpt << "_" << maxpt ;
   output_filename << ".dat";
   std:: ofstream writefile(output_filename.str().c_str());
  //std:: ofstream writefile3("test_gsl_coarse.dat");
  std :: cout << "Calculating Rapidity differential Polarization for the PID " << particle_of_interest.number << " ..." << std ::endl; 
  double y_grid[NY],  y_diff_P_rest[NY], denominator_y_diff[NY];

 
  for(int ieta = 0; ieta < N_steps_y; ieta++){
    y_grid[ieta] = - y_max + ieta*dy;
  }
  

  for(int i = 0; i < PIDS.size(); i++){  
    double mass_h = my_Particle_list[i].mass;
    for(int ieta = 0; ieta < N_steps_y; ieta++){
      double sum = 0;
      double sum_denom = 0;
      double eta = y_grid[ieta];
      // Integrate over phi using trapezoid rule
      for(int iphi = 0; iphi < N_steps_phi; iphi++){
	double ptlist[NPT] = {0};
	
	// Integrate over pt using gsl
	double Py_pt_grid[NPT] = {0};
	double denom_pt_grid[NPT] = {0};
	for(int ipt = 0; ipt < N_steps_pt; ipt++){
          double pt = (pt_min + (pt_max - pt_min)*pow(static_cast<double>(ipt), 2.)
		  /pow(static_cast<double>(N_steps_pt - 1), 2.)) ; 

           ptlist[ipt] = pt ; 

	  Py_pt_grid[ipt]    = pt * pt * cosh(eta)/sqrt(mass_h*mass_h + pt*pt*cosh(eta)*cosh(eta)) * diff_P_rest[ipt][iphi][ieta][2];
	  denom_pt_grid[ipt] = pt * pt * cosh(eta)/sqrt(mass_h*mass_h + pt*pt*cosh(eta)*cosh(eta)) * denom_dndydptdphi[ipt][iphi][ieta];
	  //writefile3 << ylist[ieta] << " " << Py_rapidity_grid[ieta] << std :: endl;
	  //Py_rapidity_grid[ieta] = diff_P_rest[ipt][iphi][ieta][1];
	}
	//writefile3.close();
	
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *splinept    = gsl_spline_alloc (gsl_interp_linear, N_steps_pt);
	gsl_spline_init(splinept, ptlist, Py_pt_grid, N_steps_pt);
	

	gsl_spline *spline_denompt = gsl_spline_alloc (gsl_interp_linear, N_steps_pt);
	gsl_spline_init(spline_denompt, ptlist, denom_pt_grid, N_steps_pt);
	
        
	double Py_pt_integrated;
	double denom_pt_integrated;


	// integrated from minpt to maxpt
	Py_pt_integrated = gsl_spline_eval_integ(splinept, minpt, maxpt, acc) ;
	denom_pt_integrated = gsl_spline_eval_integ(spline_denompt, minpt, maxpt, acc) ;
	//trapezoid
	sum += Py_pt_integrated*(2. * M_PI_)/(N_steps_phi-1);
	sum_denom += denom_pt_integrated*(2. * M_PI_)/(N_steps_phi-1);
	gsl_spline_free (splinept);
	gsl_spline_free (spline_denompt);
	gsl_interp_accel_free (acc);
        
      }// phi loop
      y_diff_P_rest[ieta] = sum; 
      denominator_y_diff[ieta] = sum_denom ;    
      writefile << eta << "  " << - sum*(100) << "  " << sum_denom << "  " << - sum*(100) / sum_denom  << std :: endl;    
    } //eta-loop  

      
  } //particle loop
  
}



double freezeout::dydeta(double eta, double pt, double m) {
  return pt*cosh(eta)/sqrt(m*m + pt*pt*cosh(eta)*cosh(eta));
}









