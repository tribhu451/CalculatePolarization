// Copyright 2012 Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_EVOLVE_H_
#define SRC_EVOLVE_H_

#include <memory>
#include <vector>
#include "util.h"
#include "data.h"
#include "cell.h"
#include "grid.h"
#include "grid_info.h"
#include "eos.h"
#include "advance.h"
#include "hydro_source_base.h"
#include "u_derivative.h"
#include "pretty_ostream.h"
#include "HydroinfoMUSIC.h"

// this is a control class for the hydrodynamic evolution
class Evolve {
 private:
    const EOS &eos;        // declare EOS object
    const InitData &DATA;
    std::weak_ptr<HydroSourceBase> hydro_source_terms_ptr;

    Cell_info grid_info;
    Advance advance;
    U_derivative u_derivative;
    pretty_ostream music_message;


    // simulation information
    int rk_order;

    int facTau;

    // information about freeze-out surface
    // (only used when freezeout_method == 4)
    int n_freeze_surf;
    std::vector<double> epsFO_list;

    typedef std::unique_ptr<SCGrid, void(*)(SCGrid*)> GridPointer;

 public:
    Evolve(const EOS &eos, const InitData &DATA_in,
           std::shared_ptr<HydroSourceBase> hydro_source_ptr_in);
    int EvolveIt(SCGrid &arena_prev, SCGrid &arena_current,
                 SCGrid &arena_future, HydroinfoMUSIC &hydro_info_ptr);

    void AdvanceRK(double tau, GridPointer &arena_prev, GridPointer &arena_current, GridPointer &arena_future);

    int FreezeOut_equal_tau_Surface(double tau, SCGrid &arena_current);
    void FreezeOut_equal_tau_Surface_XY(double tau,
                                        int ieta, SCGrid &arena_current,
                                        int thread_id, double epsFO);
   //functions modified by Koustav                                   
    int FindFreezeOutSurface_Cornelius(double tau,
                                       SCGrid &arena_current,
                                       SCGrid &arena_freezeout,
                                       SCGrid &arena_prev,
                                       SCGrid &arena_just_before_freezeout
                                       );
    int FindFreezeOutSurface_Cornelius_XY(double tau, int ieta,
                                          SCGrid &arena_current,
                                          SCGrid &arena_freezeout, 
                                          SCGrid &arena_prev,
                                          SCGrid &arena_just_before_freezeout,
                                          int thread_id, double epsFO);
    int FindFreezeOutSurface_boostinvariant_Cornelius(
                double tau, SCGrid &arena_current, SCGrid &arena_freezeout);
   
   //functions modified by Koustav
   void store_previous_step_for_freezeout(SCGrid &arena_current,
                                                          SCGrid &arena_freezeout,
                                                            SCGrid &arena_prev,
                                                               SCGrid &arena_just_before_freezeout);
    void regulate_qmu(const double u[], const double q[], double q_regulated[]) const;
    void regulate_Wmunu(const double u[], const double Wmunu[4][4], double Wmunu_regulated[4][4]) const;

    void initialize_freezeout_surface_info();
    
};

#endif  // SRC_EVOLVE_H_

