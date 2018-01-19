//
//  el_kd_g2int.cpp
//  Jensen_code
//
//  Created by Christopher Kung on 2/11/16.
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

void el_kd_g2int(rowvec coordsx, 
		rowvec d, 
		vec params,
                int el,
		mat & stiff_el) {
    double x1, x2, el_length, numips, Area;
    double const0, lambda, mu, jj;
    double xi, rho, grav;
    double D11;
    int ii;
    
    string dirName;
    
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
    el_length = x2 - x1;
    
    lambda = params(0);
    mu = params(1);
    
    rho = params(2);
    grav = params(3);
    
    numips = round(params(4));
    Area = params(7);
    
    stiff_el.zeros(2,2);
    
    const0 = 1.0/sqrt(3.0);
    xi_vect(0) = -const0;
    xi_vect(1) =  const0;
    weight(0) = 1.0;
    weight(1) = 1.0;
    
    for(ii = 0; ii < numips; ii++) {
        xi = xi_vect(ii);
        
        Bu.ones(2);
        Bu(0) = -1.0/el_length;
        Bu(1) =  1.0/el_length;
        
        dudX = Bu*d(span(0,1)).t();
        dudX1 = 1 + dudX;
        
        D11 = mu+1.0/(dudX1(0) * dudX1(0) * (lambda + mu - lambda*log(dudX1(0))));
        
        jj = el_length/2.0;
        
        stiff_el = stiff_el + Bu.t() * Bu * D11 * jj * weight(ii) * Area;
    }
}
