
// ===========================  Information for Koustav  ============================= // 
// * We only need to calculate the polarization of specific baryons only. 
//   So, there is no need to read the properties of all particles from PDG. 
// =================================================================================== //







#include"./freezeout.h"
#include<fstream>
#include<iostream>
#include<chrono>

using namespace std :: chrono;
using namespace std;
int main(int argc, char* argv[]){
    auto start = high_resolution_clock::now();
    
    double pt_min = 0.495, pt_max = 3.005, y_max = 1.05;
    int n_pt  = 16 ; // length of the total pt_array
    int n_y   = 33 ; // array length in pseudo-rapidity 
    int n_phi = 24 ;
    int number_of_particles = 1;
   
    freezeout *frzout;
    frzout = new freezeout();

    frzout->prepare_data_files(argv[1], argv[2]);
    // description about above function //
    // argument-1 : path to hypersurface file
    // argument-2 : path to the file which contains PID of the particles
    //              of interest. Polarization of those specific partilces
    //              will be calculated in the code.


    frzout->set_momentum_grid_params(pt_min, pt_max, n_pt, n_phi, y_max, n_y);
    //frzout->phase_space_distribution_integration();

    frzout->calc_polarization();
    frzout->calc_pol_related_observables();


    delete frzout;

    auto stop = high_resolution_clock :: now();
    auto duration = duration_cast<microseconds>(stop - start);
    std :: cout << "Execution time :" << duration.count()/pow(10, 6) << " seconds" << std:: endl;
    return EXIT_SUCCESS;
}

