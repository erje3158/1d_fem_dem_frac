//
//  timestepping.cpp
//  Jensen_code
//
//  Created by Christopher Kung on 2/9/16.
//  Copyright © 2016 Christopher Kung. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <sys/stat.h>

#include "armadillo"
#include "routines.h"

#include "mpi.h"

using namespace arma;

void timestepping(const char * outputDir,
                  int n_save,
                  int n_stop,
                  double dt,
                  double multi,
                  vec dispfun_disp,
                  mat d_stop,
                  mat v_stop,
                  mat a_stop,
                  mat t_stop,
                  field<cube> stress_stop,
                  field<cube> isv_stop,
                  mat F_S_stop,
                  vec &d,
                  vec &v,
                  vec &a,
                  double * t,
                  cube & stress_el,
                  cube & isv_el,
                  vec & F_S,
                  int rank,
                  int numtasks) {
    
    double rho, grav, numips, Area, dd;
    double lambda, mu, r, LDratio;
    double alphaM, h;
    double h_DEM, w_DEM, l_DEM, A_DEM;
    double time_tot, t_ramp, strainrate;
    double gd_n, gd, traction_max, traction;
    double F_F, beta, gamma;
    double tolr, tola;
    double d_1, d_2, d_3;
    double tract;
    double Rtol, normR;
    
    int ii, jj, kk, I, J, nstress, nisv, nel, neldof, ndof, nsteps, el, n,ip;
    int iter, iter_break;
    int nstart, nend;
    
    string dirName;
    
    cube mass_el, stiff_el;
    
    mat coords, LM, g;
    mat fg_el, fs_el, d_el, v_el, a_el, a_el_last, Delta_a_el;
    mat dsolve, vsolve, asolve;
    mat M, C, K, Fsolve, F_Ssolve;
    mat dR;
    mat g_n, d_el_last;
    mat K_el;
    
    field<cube> stress_solve, stress_solve_DEM, isv_solve;
    
    rowvec d_last, v_last, a_last;
    rowvec d_pred, v_pred;
    
    vec params, temp;
    vec tsolve, F, F_G, R, R0;
    vec dispfun_disp_new;
    vec del_a, Delta_a;
    vec F_S_el;
    
    double K_temp, K_total, F_S_temp, F_S_total;

    demInput demParams;
    
    lambda = 3.13e8;
    mu = 2.8e8;
    rho = 1.54;
    rho = rho*1000;
    
    grav = 0.0;
    
    dd = 0.5;
    dd = dd * 0.0254;
    r = dd/2;
    
    LDratio = 0.81;
    Area = PI*pow(r,2);
    
    alphaM = 0.0;
    
    h = dd * LDratio;
    
    coords.zeros(2,2);
    coords(0,1) = h/2.0;
    coords(1,0) = h/2.0;
    coords(1,1) = h;
    
    h_DEM = 0.003;
    w_DEM = 0.003;
    l_DEM = 0.003;
    A_DEM = l_DEM * w_DEM;
    
    numips = 2;
    
    nstress = 6;
    nisv = 2;
    ndof = 1;
    nel = 2;
    neldof = 2;
    
    params.zeros(10);
    params(0) = lambda;
    params(1) = mu;
    params(2) = rho;
    params(3) = grav;
    params(4) = numips;
    params(5) = nstress;
    params(6) = nisv;
    params(7) = Area;
    params(8) = h_DEM;
    params(9) = nel;
    
    LM.zeros(2,2);
    LM(1,0) = 1.0;
    LM(0,1) = 1.0;
    
    *t = 0.0;
    dt = dt/multi;
    
    time_tot = 0.001;
    nsteps = round(time_tot/dt);
    
    t_ramp = time_tot;
    
    strainrate = 387.0;
    
    g.zeros(2,2);
    g(1,1) = dispfun_disp(0);
    
    gd_n = 0.0;
    gd = 0.0;
    
    traction_max = 0.0;
    traction = 0.0;
    F_F = traction * Area;
    
    stress_el.zeros(nstress, numips,nel);
    mass_el.zeros(neldof, neldof,nel);
    stiff_el.zeros(neldof,neldof,nel);
    fg_el.zeros(neldof, nel);
    fs_el.zeros(neldof, nel);
    isv_el.zeros(nisv, numips, nel);
    d_el.zeros(nel, neldof);
    v_el.zeros(nel, neldof);
    a_el.zeros(nel, neldof);
    a_el_last.zeros(nel,neldof);
    Delta_a_el.zeros(nel, neldof);
    
    dsolve.zeros(nsteps, ndof);
    vsolve.zeros(nsteps, ndof);
    asolve.zeros(nsteps, ndof);
    tsolve.zeros(nsteps);
    M.zeros(ndof, ndof);
    C.zeros(ndof, ndof);
    K.zeros(ndof, ndof);
    F.zeros(ndof);
    F_G.zeros(ndof);
    F_S.zeros(ndof);
    Fsolve.zeros(nsteps, ndof);
    F_Ssolve.zeros(nsteps, ndof);
    stress_solve.set_size(nsteps);
    stress_solve_DEM.set_size(nsteps);
    isv_solve.set_size(nsteps);
    dR.zeros(ndof, ndof);
    R.zeros(ndof);
    
    for(ii = 0; ii < nsteps; ii++) {
        stress_solve(ii).zeros(nstress, numips, nel);
        stress_solve_DEM(ii).zeros(nstress, numips, nel);
        isv_solve(ii).zeros(nisv, numips, nel);
    }
    
    beta = 1.0/4.0;
    gamma = 1.0/2.0;
    
    tolr = 1.0e-4;
    tola = 1.0e-4;
    if (dt==1.0e-7) {
        iter_break = 5.0;
    } else if (dt==1.0e-8) {
        iter_break = 15.0;
    } else {
        cout << "Invalid dt = " << to_string(dt) << endl;
    }
    
    d.zeros(ndof);
    v.zeros(ndof);
    a.zeros(ndof);
    
    for(el = 0; el < nel; el++) {
        mass_el.slice(el) = el_md_g1int(coords.row(el).t(), params);
        fg_el.col(el) = el_f_g4int(coords.row(el).t(), params);
    }

    for(el = 0; el < nel; el++) {
        temp = conv_to<vec>::from(LM.col(kk));
        for(ii = 0; ii < neldof; ii++) {
            I = int(temp(ii));
            if(I > 0) {
                F_G(I-1) = F_G(I-1) + fg_el(ii, el);
                for(jj = 0; jj < neldof; jj++) {
                    J = int(temp(jj));
                    if(J > 0) {
                        M(I-1, J-1) = M(I-1, J-1) + mass_el(ii,jj,el);
                        C(I-1, J-1) = C(I-1, J-1) + alphaM * mass_el(ii,jj,el);
                    }
                }
            }
        }
    }

    F = (-F_F) - F_G;
    a = solve(M, -F);

    dsolve.row(0) = d.t();
    vsolve.row(0) = v.t();
    asolve.row(0) = a.t();
    
    if (rank == 0) {
		for(jj = 1; jj <= nel; jj++) {
			for(ii = 1; ii <= numips; ii++) {
				dirName = string(outputDir) + "/el" + to_string(jj) + "_ip" + to_string(ii);
				if (dt==1.0e-7) {
					copyfile(dirName + "/input_boundary_file_" + to_string(n_save), dirName + "/input_boundary_file");
					copyfile(dirName + "/input_particle_file_" + to_string(n_save), dirName + "/input_particle_file");
				} else if (dt==1.0e-8) {
					copyfile(dirName + "/input_boundary_file_" + to_string(n_save) + "_" + to_string(n_stop), dirName + "/input_boundary_file");
					copyfile(dirName + "/input_particle_file_" + to_string(n_save) + "_" + to_string(n_stop), dirName + "/input_particle_file");				
				} else {
					cout << "Location 1: dt isn't correct: dt = " << dt << endl;
					exit(1);
				}
			}
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);

    dsolve.row(n_stop*int(multi)) = d_stop.row(n_stop);
    vsolve.row(n_stop*int(multi)) = v_stop.row(n_stop);
    asolve.row(n_stop*int(multi)) = a_stop.row(n_stop);
    tsolve.row(n_stop*int(multi)) = t_stop.row(n_stop);
    stress_solve(n_stop*int(multi)) = stress_stop(n_stop);
    isv_solve(n_stop*int(multi)) = isv_stop(n_stop);
    F_Ssolve.row(n_stop*int(multi)) = F_S_stop.row(n_stop);
    
    nstart = n_stop * int(multi);
    nend   = n_stop * int(multi) + (int(multi) - 1);
    
    d_2 = dispfun_disp(n_stop+1);
    d_3 = dispfun_disp(n_stop+2);
    
    if(n_stop==1) {
        d_1 = 0.0;
    }
    else {
        d_1 = dispfun_disp(n_stop);
        d_1 = d_2 - ((d_2 - d_1)/multi);
    }
    
    iter = 0;
    dispfun_disp_new = zeros(nsteps,1);
    dispfun_disp_new(nstart) = d_1;
    for (ii = nstart; ii <= nend; ii++) {
        dispfun_disp_new(ii+1) = d_2 + ((d_3-d_2)/multi)*iter;
        iter = iter + 1;
    }
    
    gd = d_1;
    
    //this is super hardcoded
	if (rank == 0) {
		el = 0;
		ip = 0;
	} else if (rank == 1) {
		el = 0;
		ip = 1;
	} else if (rank == 2) {
		el = 1;
		ip = 0;
	} else {
		el = 1;
		ip = 1;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
    
    for(n = nstart; n <= nend; n++) {
		if (rank == 0) {
			for(el = 1; el <= nel; el++) {
				for(ip = 1; ip <= numips; ip++) {
					if (dt==1.0e-7) {
						copyfile(string(outputDir) + "/el" + to_string(el) + "_ip" + to_string(ip) + "/input_boundary_file", 
								string(outputDir) + "/el" + to_string(el) + "_ip" + to_string(ip) + "/input_boundary_file_" + to_string(n_save) + "_" + to_string(n));
						copyfile(string(outputDir) + "/el" + to_string(el) + "_ip" + to_string(ip) + "/input_particle_file",
								string(outputDir) + "/el" + to_string(el) + "_ip" + to_string(ip) + "/input_particle_file_" + to_string(n_save) + "_" + to_string(n));					
					} else if (dt==1.0e-8) {
						copyfile(string(outputDir) + "/el" + to_string(el) + "_ip" + to_string(ip) + "/input_boundary_file", 
								string(outputDir) + "/el" + to_string(el) + "_ip" + to_string(ip) + "/input_boundary_file_" + to_string(n_save) + "_" + to_string(n_stop) + "_" + to_string(n));
						copyfile(string(outputDir) + "/el" + to_string(el) + "_ip" + to_string(ip) + "/input_particle_file",
								string(outputDir) + "/el" + to_string(el) + "_ip" + to_string(ip) + "/input_particle_file_" + to_string(n_save) + "_" + to_string(n_stop) + "_" + to_string(n));						
					} else {
						cout << "Location 2: dt isn't correct: dt = " << dt << endl;
						exit(1);
					}
				}
			}
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
		
        *t = *t+dt;
        gd_n = gd;
        g_n.zeros(2,2);
        g_n(1,1) = gd_n;
        if(*t < t_ramp) {
            tract = traction_max * (*t/t_ramp);
        }
        else {
            tract = traction_max;
        }
        gd = dispfun_disp_new(n-1);
        g.zeros(2,2);
        g(1,1) = gd;
        
        F_F = tract*Area;
        
        d_last = dsolve.row(n-1);
        v_last = vsolve.row(n-1);
        a_last = asolve.row(n-1);
        
        Rtol = 1;
        normR = 1;
        kk = 0;
        
        d_pred = d_last + dt*v_last * pow(dt,2) * (1.0-2.0*beta) * a_last/2.0;
        v_pred = v_last + dt*(1.0 - gamma) * a_last;
        
        while ((Rtol > tolr) && (normR > tola)) {
            kk = kk + 1;
            
            if(kk == 1) {
                del_a = zeros(ndof);
            }
            else {
                del_a = solve(dR, -R);
            }
            
            a = a + del_a;
            
            Delta_a = a - a_last.t();
            
            R.zeros(ndof);
            dR.zeros(ndof, ndof);
            F_S.zeros(ndof);
            K.zeros(ndof, ndof);
            
            d = d_pred.t() + pow(dt,2) * beta * a;
            v = v_pred.t() + dt * gamma * a;
            d_el_last.zeros(nel, neldof);
            temp = conv_to<vec>::from(LM.col(el));
            for(ii = 0; ii < neldof; ii++) {
                I = temp(ii);
                if(I > 0) {
                    d_el(el,ii) = d(I-1);
                    d_el_last(el,ii) = d_last(I-1);
                    v_el(el,ii) = v(I-1);
                    a_el(el,ii) = a(I-1);
                    a_el_last(el,ii) = a_last(I-1);
                    Delta_a_el(el,ii) = Delta_a(I-1);
                }
                else {
                    d_el(el,ii) = g(ii,el);
                    d_el_last(el,ii) = g_n(ii,el);
                }
            }
                
            if (dt==1.0e-7) {
				el_stress_ellip3d(outputDir, coords.row(el), d_el.row(el), d_el_last.row(el), params, n_save, -1, n, el, ip, stress_el, isv_el, dt, demParams);
			} else if (dt==1.0e-8) {
				el_stress_ellip3d(outputDir, coords.row(el), d_el.row(el), d_el_last.row(el), params, n_save, n_stop, n, el, ip, stress_el, isv_el, dt, demParams);
			} else {
				cout << "Location 3: dt isn't correct: dt = " << dt << endl;
				exit(1);
			}
			
			MPI_Barrier(MPI_COMM_WORLD);
	    // All sorts of needs fixing
            el_kd_g2int_ellip3d(outputDir, coords.row(el), d_el.row(el), params, n_save, el, stiff_el.slice(el));
            fs_el.col(el) = el_f_g2int(coords.row(el), stress_el.slice(el), params);

            F_S_el = fs_el.col(el);
            K_el = stiff_el.slice(el);
                
            temp = conv_to<vec>::from(LM.col(el));
            for(ii = 0; ii < neldof; ii++) {
                I = temp(ii);
                if(I > 0) {
                    F_S(I-1) = F_S(I-1) + F_S_el(ii);
                    for(jj = 0; jj < neldof; jj++) {
                        J = temp(jj);
                        if(J > 0) {
                            K(I-1, J-1) = K(I-1, J-1) + K_el(ii,jj);
                        }
                    }
                }
            }

            for(ii = 0; ii < K.n_rows; ii++) {
            	for(jj = 0; jj < K.n_cols; jj++) {
            		 MPI_Allreduce(&K(ii,jj),&K_total,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            		 K(ii,jj) = K_total;
            	}
            }

            for(ii = 0; ii < F_S.n_elem; ii++) {
            	MPI_Allreduce(&F_S(ii), &F_S_total, 1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
            	F_S(ii) = F_S_total;
            }

            R = M*a + C*v + F_S - F_F - F_G;
            dR = M + dt*gamma*C + (dt*dt) * beta * K;
            if(kk == 1) {
                R0 = R;
            }
            Rtol = norm(R)/norm(R0);
            normR = norm(R);

            if(kk == iter_break) {
                if(dt==1.0e-7) {
                    timestepping(outputDir, n_save, n, dt, multi, dispfun_disp_new, dsolve, vsolve, asolve, tsolve, stress_solve, isv_solve, F_Ssolve,
                                 d, v, a, t, stress_el, isv_el, F_S, rank, numtasks);
                    
                    Rtol = tolr/10.0;
                    normR = tola/10.0;
                }
                else if(dt == 1.0e-8) {
                    cout << "Max Iterations reached" << endl;
                    break;
                }
                else {
                    cout << "Screwy dt" << endl;
                    break;
                }
            }

			for(ii = 0; ii < stress_el.n_rows; ii++) {
				for(jj = 0; jj < stress_el.n_cols; jj++) {
					for(kk = 0; kk < stress_el.n_slices; kk++) {
						MPI_Allreduce(&stress_el(ii,jj,kk),&K_total,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
						stress_el(ii,jj,kk) = K_total;
					}
				}
			}

			for(ii = 0; ii < isv_el.n_rows; ii++) {
				for(jj = 0; jj < isv_el.n_cols; jj++) {
					for(kk = 0; kk < isv_el.n_slices; kk++) {
						MPI_Allreduce(&isv_el(ii,jj,kk),&K_total,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
						isv_el(ii,jj,kk) = K_total;
					}
				}
			}

        }

        dsolve.row(n) = d;
        vsolve.row(n) = v;
        asolve.row(n) = a;
        tsolve(n) = *t;
        stress_solve(n) = stress_el;
        isv_solve(n) = isv_el;
        F_Ssolve.row(n) = F_S;
    }
}
