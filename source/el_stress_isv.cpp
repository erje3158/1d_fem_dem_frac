//
//  el_stress_isv.cpp
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
#include <sstream>
#include <string>
#include "armadillo"
#include "routines.h"

using namespace arma;

void el_stress_isv(rowvec coordsx,
                   rowvec d,
                   vec params,
                   int el,
                   int ip,
                   cube & stress_el,
                   cube & isv_el) {
    
    double x1, x2, el_length, numips, Area;
    double lambda, mu;
    double const0, xi;
    double nisv, nstress, nel, dudX1, S11, P11, sig11;
    double E11, e11, le11;
    double trS, devS;
    
    vec f_e, dudX;
    vec xi_vect(2);
    vec weight(2);
    rowvec Bu(2);
    
    x1 = coordsx(0);
    x2 = coordsx(1);
    el_length = x2 - x1;
    
    lambda = params(0);
    mu = params(1);
    
    numips = round(params(4));
    nstress = round(params(5));
    nisv = round(params(6));
    
    Area = params(7);

    nel = params(9);
    
    stress_el.zeros(nstress,numips,nel);
    isv_el.zeros(nisv,numips,nel);
    
    const0 = sqrt(1.0/3.0);
    xi_vect(0) = -const0;
    xi_vect(1) = const0;
    
    weight(0) = 1.0;
    weight(1) = 1.0;
    
    xi = xi_vect(ip);
        
    Bu(0) = -1.0/el_length;
    Bu(1) = 1.0 /el_length;
        
    dudX = Bu * d(span(0,1)).t();
    dudX1 = 1.0 + dudX(0);
        
    S11 = mu + (lambda * log(dudX1) - mu) / (dudX1 * dudX1);
    P11 = dudX1 * S11;
    sig11 = P11;
        
    E11 = (dudX1 * dudX1 - 1.0)/2.0;
    e11 = (1.0 - (1.0/(dudX1*dudX1)))/2.0;
    le11 = log(dudX1);
        
    stress_el(0,ip,el) = sig11;
    stress_el(1,ip,el) = P11;
    stress_el(2,ip,el) = S11;
    stress_el(3,ip,el) = e11;
    stress_el(4,ip,el) = E11;
    stress_el(5,ip,el) = le11;
      
    trS = 0.0;
    devS = 0.0;
        
    isv_el(0,ip,el) = trS;
    isv_el(1,ip,el) = devS;
}
