//
//  el_kd_g2int_ellip3d.cpp
//  Jensen_code
//
//  Created by Christopher Kung on 2/4/16.
//  Copyright © 2016 Christopher Kung. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "armadillo"
#include "routines.h"

using namespace arma;

void el_kd_g2int_ellip3d(const char * outputDir, 
                         rowvec coordsx, 
			 rowvec d, 
			 vec params, 
			 int n, 
			 int el, 
			 mat & stiff_el) {
    double x1, x2, el_length, numips, Area;
    double const0, xi, lambda, mu, jj;
    double rho, grav, D11, trash;
    int ip, ii;
    
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
    
    x1 = coordsx(0);
    x2 = coordsx(1);
    el_length = x2-x1;
    
    lambda = params(0);
    mu = params(1);
    rho = params(2);
    grav = params(3);
    
    numips = round(params(4));
    
    Area = params(7);
    
    // Initialze Target Matrix
    stiff_el.zeros(2,2);
    
    // Set Gauss point coordinates in xi space
    const0 = 1/sqrt(3.0);
    xi_vect(0) = -const0;
    xi_vect(1) =  const0;
    weight(0)  = 1.0;
    weight(1)  = 1.0;
    
    //Loop through the Gauss Points
    for(ip = 0; ip < numips; ip++) {
        dirName = string(outputDir) + "/el" + to_string(el+1) + "_ip" + to_string(ip+1);
        
        // Code to read comp_prgress to get value for D11
        std:ifstream infile(dirName + "/comp_progress");
        string lastLine = getLastLine(infile);
        for(ii = 0; ii < 130; ii++) {
            if(ii == 0)
                trash = stod(lastLine, &sz);
            else {
                trash = stod(lastLine.substr(sz), &addin);
                sz = sz + addin;
            }
            if(ii == 129) D11 = trash;
        }
        infile.close();
        
        xi = xi_vect(ip);
        Bu.ones(2);
        Bu(0) = -1/el_length;
        Bu(1) =  1/el_length;
        
        dudX = Bu*d(span(0,1)).t();
        dudX1 = dudX + 1;
        
        jj = el_length / 2.0;
        
        stiff_el = stiff_el + Bu.t() * Bu*D11 * jj * weight(ip) * Area;
    }
}
