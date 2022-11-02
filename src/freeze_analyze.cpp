#include "freezeout.h"
#include <omp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


#define M_PI_ 3.1415927

void freezeout::calc_pol_related_observables(int yflag, double rap_low, double rap_up ){
   // if yflag > 0  rapidity
   // else          pseudo-rapidity


   double minrap = rap_low, maxrap = rap_up;
   double minpt = 0.5, maxpt = 3.0; 
   std::stringstream output_filename;
   output_filename.str("");
   output_filename << "results/MinusP_percent_vs_pT_" << particle_of_interest.number ;
   if(yflag > 0){
     output_filename <<"_ycut_" ; 
   }
   else{
     output_filename <<"_etacut_" ; 
   }
   output_filename << minrap << "_" << maxrap ;
   output_filename << ".dat";
   std:: ofstream writefile(output_filename.str().c_str());
   output_filename.str("");
   output_filename << "results/MinusPy_percent_" << particle_of_interest.number ;
   if(yflag > 0){
     output_filename <<"_ycut_" ; 
   }
   else{
     output_filename <<"_etacut_" ; 
   }
   output_filename << minrap << "_" << maxrap << "_ptcut_" ;
   output_filename << minpt << "_" << maxpt << ".dat";
   std:: ofstream writefile_all_integ(output_filename.str().c_str());

  //std:: ofstream writefile3("test_gsl_coarse.dat");
  std :: cout << "Calculating Polarization observables for the PID " << particle_of_interest.number << " ..." << std ::endl; 
  double pt_grid[NPT], pt_diff_Py_rest[NPT], denominator_pt_diff[NPT];

 
  for (int i = 0; i < N_steps_pt; i++){
    pt_grid[i] = (pt_min + (pt_max - pt_min)*pow(static_cast<double>(i), 2.)
		  /pow(static_cast<double>(N_steps_pt - 1), 2.));
  }
  

  for(int i = 0; i < PIDS.size(); i++){  
    double mass_h = my_Particle_list[i].mass;
    for(int ipt = 0; ipt < N_steps_pt; ipt++){
      double sum_pt = 0;
      double sum_px = 0;
      double sum_py = 0;
      double sum_pz = 0;
      double sum_denom = 0;
      double pt = pt_grid[ipt];
      // Integrate over phi using trapezoid rule
      for(int iphi = 0; iphi < N_steps_phi; iphi++){
	double ylist[NY] = {0};	
	// Integrate over pseudorapidity using gsl
	double Pt_rapidity_grid[NY] = {0};
	double Px_rapidity_grid[NY] = {0};
	double Py_rapidity_grid[NY] = {0};
	double Pz_rapidity_grid[NY] = {0};
	double denom_rapidity_grid[NY] = {0};
	for(int ieta = 0; ieta < N_steps_y; ieta++){
	  double eta_local = - y_max + ieta*dy;
	  double y_local = Rap(eta_local, pt, mass_h);    // get rapidity

          if( yflag > 0 ){
            if (rap_low < Rap(-y_max, pt, mass_h)){
                 std::cout << " Error in calculating pt differential Polarization with y_cut ..." << std::endl; 
                 std::cout << " min_rap = " <<  Rap(-y_max, pt, mass_h) << " min_rap_cut = " << rap_low << std::endl; 
               }  
            if (rap_up > Rap(y_max, pt, mass_h)){
                 std::cout << " Error in calculating pt differential Polarization with y_cut ..." << std::endl ;
                 std::cout << " max_rap = " <<  Rap(y_max, pt, mass_h) << " max_rap_cut = " << rap_up << std::endl;  
               }  
          }

          if(yflag > 0){ 
	    ylist[ieta] = y_local;
          }
          else{
	    ylist[ieta] = eta_local;
          }

          if(yflag > 0){ 
	    Pt_rapidity_grid[ieta]    = pt * diff_P_rest[ipt][iphi][ieta][0];
	    Px_rapidity_grid[ieta]    = pt * diff_P_rest[ipt][iphi][ieta][1];   
	    Py_rapidity_grid[ieta]    = pt * diff_P_rest[ipt][iphi][ieta][2];
	    Pz_rapidity_grid[ieta]    = pt * diff_P_rest[ipt][iphi][ieta][3];
	    denom_rapidity_grid[ieta] = pt * denom_dndydptdphi[ipt][iphi][ieta];                     
          }
          else{
	    Pt_rapidity_grid[ieta]    = pt * pt * cosh(eta_local)/sqrt(mass_h*mass_h + pt*pt*cosh(eta_local)*cosh(eta_local)) * diff_P_rest[ipt][iphi][ieta][0];
	    Px_rapidity_grid[ieta]    = pt * pt * cosh(eta_local)/sqrt(mass_h*mass_h + pt*pt*cosh(eta_local)*cosh(eta_local)) * diff_P_rest[ipt][iphi][ieta][1];
	    Py_rapidity_grid[ieta]    = pt * pt * cosh(eta_local)/sqrt(mass_h*mass_h + pt*pt*cosh(eta_local)*cosh(eta_local)) * diff_P_rest[ipt][iphi][ieta][2];
	    Pz_rapidity_grid[ieta]    = pt * pt * cosh(eta_local)/sqrt(mass_h*mass_h + pt*pt*cosh(eta_local)*cosh(eta_local)) * diff_P_rest[ipt][iphi][ieta][3];
	    denom_rapidity_grid[ieta] = pt * pt * cosh(eta_local)/sqrt(mass_h*mass_h + pt*pt*cosh(eta_local)*cosh(eta_local)) * denom_dndydptdphi[ipt][iphi][ieta];
          }
	}

	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *splinept    = gsl_spline_alloc (gsl_interp_linear, N_steps_y);
	gsl_spline_init(splinept, ylist, Pt_rapidity_grid, N_steps_y);
	gsl_spline *splinepx    = gsl_spline_alloc (gsl_interp_linear, N_steps_y);
	gsl_spline_init(splinepx, ylist, Px_rapidity_grid, N_steps_y);
	gsl_spline *splinepy    = gsl_spline_alloc (gsl_interp_linear, N_steps_y);
	gsl_spline_init(splinepy, ylist, Py_rapidity_grid, N_steps_y);
	gsl_spline *splinepz    = gsl_spline_alloc (gsl_interp_linear, N_steps_y);
	gsl_spline_init(splinepz, ylist, Pz_rapidity_grid, N_steps_y);
	
	gsl_interp_accel *accd = gsl_interp_accel_alloc ();
	gsl_spline *spline_denom = gsl_spline_alloc (gsl_interp_linear, N_steps_y);
	gsl_spline_init(spline_denom, ylist, denom_rapidity_grid, N_steps_y);
	
        
	double Pt_pseudorapidity_integrated ;
	double Px_pseudorapidity_integrated;
	double Py_pseudorapidity_integrated;
	double Pz_pseudorapidity_integrated;
	double denom_pseudorapidity_integrated;

	// integrated from minrap to maxrap
	Pt_pseudorapidity_integrated = gsl_spline_eval_integ(splinept, minrap, maxrap, acc) ;
	Px_pseudorapidity_integrated = gsl_spline_eval_integ(splinepx, minrap, maxrap, acc) ;
	Py_pseudorapidity_integrated = gsl_spline_eval_integ(splinepy, minrap, maxrap, acc) ;
	Pz_pseudorapidity_integrated = gsl_spline_eval_integ(splinepz, minrap, maxrap, acc) ;
	denom_pseudorapidity_integrated = gsl_spline_eval_integ(spline_denom, minrap, maxrap, accd) ;
	//trapezoid
	sum_pt += Pt_pseudorapidity_integrated*(2. * M_PI_)/(N_steps_phi-1);
	sum_px += Px_pseudorapidity_integrated*(2. * M_PI_)/(N_steps_phi-1);
	sum_py += Py_pseudorapidity_integrated*(2. * M_PI_)/(N_steps_phi-1);
	sum_pz += Pz_pseudorapidity_integrated*(2. * M_PI_)/(N_steps_phi-1);
	sum_denom += denom_pseudorapidity_integrated*(2. * M_PI_)/(N_steps_phi-1);
	gsl_spline_free (splinept);
	gsl_spline_free (splinepx);
	gsl_spline_free (splinepy);
	gsl_spline_free (splinepz);
	gsl_spline_free (spline_denom);
	gsl_interp_accel_free (acc);
	gsl_interp_accel_free (accd);
        
      }// phi loop
      pt_diff_Py_rest[ipt] = sum_py; 
      denominator_pt_diff[ipt] = sum_denom ;    
      writefile << pt << "  " << sum_denom << "  " << -sum_pt*(100)/sum_denom << "  " << -sum_px*(100)/sum_denom 
                                           << "  " << -sum_py*(100)/sum_denom << "  " << -sum_pz*(100)/sum_denom << std :: endl;    
    } //pt-loop  


    
    gsl_interp_accel *acc2 = gsl_interp_accel_alloc ();
    gsl_spline *spline2 = gsl_spline_alloc (gsl_interp_linear, N_steps_pt );
    gsl_spline_init (spline2, pt_grid, pt_diff_Py_rest, N_steps_pt);
    gsl_spline *spline2_denom = gsl_spline_alloc (gsl_interp_linear, N_steps_pt );
    gsl_spline_init (spline2_denom, pt_grid, denominator_pt_diff, N_steps_pt);
    
    double Py_total_integrated, denom_total_integrated ;
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
   output_filename << "results/MinusPy_percent_vs_eta_" << particle_of_interest.number <<"_ptcut_";
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









