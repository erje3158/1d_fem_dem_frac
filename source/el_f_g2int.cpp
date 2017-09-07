//
//  el_f_g2int.cpp
//  Jensen_code
//
//  Created by Christopher Kung on 2/8/16.
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

vec el_f_g2int(rowvec coordsx, mat stress, vec params) {
    double x1, x2, el_length, rho, grav, numips, Area;
    double const0, xi, P11;
    int ii;
    
    vec f_e;
    vec xi_vect(2);
    vec weight(2);
    rowvec Bu(2);
    
    x1 = coordsx(0);
    x2 = coordsx(1);
    el_length = x2-x1;
    
    rho = params(2);
    grav = params(3);
    
    numips = round(params(4));
    
    Area = params(7);
    
    // Initialize Force Vector
    f_e.zeros(2,1);
    // Set Gauss point coordinates in xi space
    const0 = 1.0/sqrt(3.0);
    xi_vect(0) = -const0;
    xi_vect(1) = const0;
    weight(0) = 1.0;
    weight(1) = 1.0;
    
    for(ii = 0; ii < int(numips); ii++) {
        
        // Current Gauss point location
        xi = xi_vect(ii);
        
        // Calculate jacobian
        
        // shape function vector;
        Bu(0) = -1.0/el_length;
        Bu(1) =  1.0/el_length;
        
        P11 = stress(1,ii);
        
        f_e = f_e + Bu.t() * P11 * (el_length/2.0) * weight(ii) * Area;
    }
    
    return f_e;
}