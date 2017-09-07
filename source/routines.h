//
//  routines.h
//  Jensen_code
//
//  Created by Christopher Kung on 1/20/16.
//  Copyright © 2016 Christopher Kung. All rights reserved.
//

#include "userInput.h"

using namespace arma;
using namespace std;

#define PI 3.1415926535897932384626433832795

mat el_md_g1int(vec coordsx,
                vec params);

vec el_f_g4int(vec coordsx,
               vec params);

void copyfile (string inputFile,
               string outputFile);

void el_stress_ellip3d(const char * outputDir,
                       rowvec coordsx,
                       rowvec d,
                       rowvec d_last,
                       vec params,
                       int n_save,
					             int n_stop,
                       int n,
                       int el,
                       int ip,
                       cube & stress_el,
                       cube & isv_el,
                       double dt,
                       demInput demParams);

void main_ellip3d(float disp_top,
                  float disp_bot,
                  int num_runs,
                  int num_threads,
                  string dirName,
                  double dt,
		  demInput demParams);

void el_kd_g2int_ellip3d(const char * outputDir,
                         rowvec coordsx,
                         rowvec d,
                         vec params,
                         int n,
                         int el,
                         mat & stiff_el);

string getLastLine(std::ifstream& in);

vec el_f_g2int(rowvec coordsx,
               mat stress,
               vec params);

void timestepping(const char * outputDir,
                  int n0,
                  int n1,
                  double dt,
                  double multi,
                  vec dispfun_disp,
                  mat dsolve,
                  mat vsolve,
                  mat asolve,
                  mat tsolve,
                  field<cube> stress_solve,
                  field<cube> isv_solve,
                  mat F_Ssolve,
                  vec &d,
                  vec &v,
                  vec &a,
                  double * t,
                  cube & stress_el,
                  cube & isv_el,
                  vec & F_S,
                  int rank,
                  int numtasks);

void el_stress_isv(rowvec coordsx,
                   rowvec d,
                   vec params,
                   int el,
                   int ip,
                   cube & stress_el,
                   cube & isv_el);

void el_kd_g2int(rowvec coordsx,
                 rowvec d,
                 vec params,
                 int el,
                 mat & stiff_el);
    
