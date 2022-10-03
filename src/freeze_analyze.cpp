#include "freezeout.h"
#include <omp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

void freezeout :: calc_pol_related_observables()
{   
    std:: ofstream writefile("Py_vs_pt.dat");

    //std:: ofstream writefile3("test_gsl_coarse.dat");
    
    std :: cout << "Calculating Polarization observables for the PID " << particle_of_interest.number << std ::endl; 
    double pt_grid[NPT], pt_diff_P_rest[NPT];

    for (int i = 0; i < N_steps_pt; i++) pt_grid[i] = (pt_min + (pt_max - pt_min)*pow(static_cast<double>(i), 2.)
                        /pow(static_cast<double>(N_steps_pt - 1), 2.));

    for (int i = 0; i < PIDS.size(); i++)  
    {   
        double mass_h = my_Particle_list[i].mass;
                       
        for (int ipt = 0; ipt < N_steps_pt; ipt++) {
            double sum = 0;
            double pt = pt_grid[ipt];
            // Integrate over phi using trapezoid rule
            for (int iphi = 0; iphi < N_steps_phi; iphi++) {
                double ylist[NY] = {0};
                double etalist[NY] = {0};

                // Integrate over pseudorapidity using gsl
                double Py_rapidity_grid[NY] = {0};
                for (int ieta = 0; ieta < N_steps_y; ieta++) {
                    double rap_local = - y_max + ieta*dy;
                    double eta_local, y_local;
                
                    {
                        eta_local = rap_local;
                        etalist[ieta] = eta_local;
                        y_local = Rap(eta_local, pt, mass_h);    // get rapidity
                        ylist[ieta] = y_local;                  // get rapidity
                    } 
                
                    Py_rapidity_grid[ieta] = pt*diff_P_rest[ipt][iphi][ieta][1];
		            //writefile3 << ylist[ieta] << " " << Py_rapidity_grid[ieta] << std :: endl;
                    //Py_rapidity_grid[ieta] = diff_P_rest[ipt][iphi][ieta][1];
                }
		        //writefile3.close();
		
		gsl_interp_accel *acc = gsl_interp_accel_alloc ();
                gsl_spline *spline    = gsl_spline_alloc (gsl_interp_linear, N_steps_y);
                
                gsl_spline_init (spline, etalist, Py_rapidity_grid, N_steps_y);
		        
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
                double minrap = -1.0, maxrap = 1.0; 
                double normalize_pseudo_rap = fabs(maxrap - minrap);
                {     
                // integrated from minrap to maxrap
                Py_pseudorapidity_integrated = gsl_spline_eval_integ(spline, minrap, maxrap, acc)/normalize_pseudo_rap;
                }
                //trapezoid
                sum += Py_pseudorapidity_integrated*(1.0)/(double)(N_steps_phi - 1);
                gsl_spline_free (spline);
                gsl_interp_accel_free (acc);
                
            }// phi loop
            pt_diff_P_rest[ipt] = sum;    
            writefile << pt << " " << - sum*(100) << std :: endl;    
        } //pt-loop  

        gsl_interp_accel *acc2 = gsl_interp_accel_alloc ();
        gsl_spline *spline2 = gsl_spline_alloc (gsl_interp_linear, N_steps_pt );
        gsl_spline_init (spline2, pt_grid, pt_diff_P_rest, N_steps_pt);
		        
        double Py_total_integrated;
        double minpt = 0.5, maxpt = 3.0; 
        double normalize_pt = fabs(maxpt - minpt);
        Py_total_integrated = gsl_spline_eval_integ(spline2, minpt, maxpt, acc2)/normalize_pt;
        std :: cout << "Final Integrated result is : " << -100*Py_total_integrated << std :: endl;
                    
        gsl_spline_free (spline2);
        gsl_interp_accel_free (acc2);

   

    } //particle loop

}
    

double freezeout :: dydeta(double eta, double pt, double m) {
    return pt*cosh(eta)/sqrt(m*m + pt*pt*cosh(eta)*cosh(eta));
}









