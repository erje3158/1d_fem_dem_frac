//
//  disp_functions.h
//  Jensen_code
//
//  Created by Erik Jensen 9/20/2017.
//  Copyright �� 2017 Erik Jensen. All rights reserved.
//

#include "armadillo"

using namespace arma;

void finiteAppliedDisp(vec & time, vec & disp, vec & eps, int nsteps, double time_tot, double strainrate, double h);    

void shpbAppliedDisp(vec & time, vec & disp, vec & eps, int nsteps, double time_tot, double strainrate, double h);