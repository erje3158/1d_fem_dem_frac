//
//  el_stress_ellip3d.cpp
//  Jensen_code
//
//  Created by Christopher Kung on 2/3/16.
//  Copyright �� 2016 Christopher Kung. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <string>
#include "armadillo"
#include "routines.h"

using namespace arma;

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
                       demInput demParams) {
    double x1, x2, el_length, numips, nstress, Area, nel;
    double const0, xi, lambda, mu, nisv, h0, el_mid;
    double sig11, mass, disp, trash;
    int ii;
    int num_threads;
    
    string dirName;
    
    string::size_type sz, addin;
    
    vec f_e;
    vec dudX, dudX1, E11, e11, le11;
    vec dudX_n, dudX1_n, le11_n;
    vec dtop_n, dtop, dbot;
    vec xi_vect(2);
    vec weight(2);
    vec gauss_loc(2);
    rowvec Bu;
    rowvec N(2);
    
    char cCurrentPath[FILENAME_MAX];
    
    getcwd(cCurrentPath, sizeof(cCurrentPath));
    
#pragma omp parallel
    num_threads = omp_get_num_threads();
#pragma end omp parallel

    cout << "Number of threads = " << num_threads << endl;

    sig11 = -1.0;
    mass  = -1.0;
    disp  = -1.0;
    sz = -1;
    
    x1 = coordsx(0);
    x2 = coordsx(1);
    el_length = x2-x1;
    
    lambda = params(0);
    mu = params(1);
    
    numips = round(params(4));
    nstress = round(params(5));
    nisv = round(params(6));
    
    Area = params(7);
    h0 = params(8);
    
    nel = params(9);
    
    
    // Initialize stress at Gauss points
    stress_el.zeros(nstress,numips,nel);
    isv_el.zeros(nisv,numips,nel);
    
    const0 = sqrt(1.0/3.0);
    xi_vect(0) = -const0;
    xi_vect(1) =  const0;
    el_mid = x1 + (el_length/2.0);
    gauss_loc(0) = el_mid + xi_vect(0) * (el_length/2.0);
    gauss_loc(1) = el_mid + xi_vect(1) * (el_length/2.0);
        
    cout << "el_stress_ellip3d: Iteration " << ip+1 << " of " << numips << endl;
        
    //cout << "OMP_NUM_THREADS " << omp_get_num_threads() << endl;
        
    xi = xi_vect(ip);
    Bu.ones(2);
    Bu(0) = -1/el_length;
    Bu(1) =  1/el_length;
        
    dudX = Bu*d(span(0,1)).t();
    dudX1 = dudX + 1;
    E11 = (dudX1 % dudX1 - 1.0)/2.0;
    e11 =  (1.0 - (1.0/(dudX1%dudX1)))/2.0;
    le11 = log(dudX1);
        
    dudX_n = Bu * d_last(span(0,1)).t();
    dudX1_n = 1 + dudX_n;
        
    le11_n = log(dudX1_n);
    dtop_n = h0 * (1-exp(le11_n));
    dtop  = (h0 * (1-exp(le11))) - dtop_n;
    dbot.zeros(1);
        
    dirName = string(outputDir) + "/el"+ to_string(el+1) + "_ip" + to_string(ip+1);
        
    cout << "dirName = " << dirName << endl;
        
    // Change to the output directory path so that ellips3d can run
    chdir(dirName.c_str());

    main_ellip3d(float(dtop(0)), float(dbot(0)), n_save, num_threads, dirName, dt, demParams);
        
    // Change back to the directory from where this code is being run so that file inputs can be referenced properly
    chdir(cCurrentPath);

    //Routines to look at comp_progress
    std:ifstream infile(string(outputDir) + "/el"+ to_string(el+1) + "_ip" + to_string(ip+1) + "/comp_progress");

    cout << "Here 1" << endl;

    //This routine retrieves the last line of comp_progress
    string lastLine = getLastLine(infile);

    cout << "Here 2" << endl;  

    // This loop goes through the last line of comp_progress and parses out the approrpriate parts
    for(ii = 0; ii < 131; ii++) {
        if(ii == 0)
            trash = stod(lastLine, &sz);
        else {
            trash = stod(lastLine.substr(sz), &addin);
            sz = sz + addin;
        }
            
        if(ii == 25) sig11 = trash;
        if(ii == 128) disp = trash;
    }

    cout << "Here 3" << endl;

    stress_el(0,ip,el) = sig11;
    stress_el(1,ip,el) = sig11;
    stress_el(2,ip,el) = sig11/dudX1(0);
    stress_el(3,ip,el) = e11(0);
    stress_el(4,ip,el) = E11(0);
    stress_el(5,ip,el) = le11(0);
        
    isv_el(0,ip,el) = mass;
    isv_el(1,ip,el) = disp;
        
    infile.close();
}
