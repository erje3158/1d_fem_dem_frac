//
//  main.cpp
//  Jensen_code
//
//  Created by Christopher Kung on 1/20/16.
//  Copyright �� 2016 Christopher Kung. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstring>
#include <fstream>

#include "armadillo"
#include "routines.h"
#include "meshTools.h"

#include "mpi.h"

using namespace std;
using namespace arma;

int main(int argc, char * argv[]) {

    double lambda, mu, rho, grav, d, r, LDratio, Area, alphaM;
    double h;
    double h_DEM, w_DEM, l_DEM, A_DEM;
    double t, dt, t_ramp, time_tot;
    double strainrate;
    double gd_n, gd;
    double traction_max, traction, F_F;
    double beta, gamma, tolr, tola, iter_break;
    double tract, normR;
    double K_temp, K_total, F_S_temp, F_S_total;
    
    int numips, nstress, nisv, ndof, nel, neldof;
    int nsteps, print_int, n_print;
    int ii, jj, kk, ll, n, el, ip;
    int I, J;
    int Rtol;
    int numtasks, rank, rc;
    
    string dirName, dirPrefix;
    
    mat coords, g, g_n;
    mat fg_el, fs_el, d_el, d_el_last, v_el, a_el, a_el_last, Delta_a_el;
    mat M, C, K, K_el;
    
    cube stress_el, mass_el, stiff_el, isv_el;
    
    vec params, dispfun_time, dispfun_disp, dispfun_eps;
    vec F, F_G, F_S, F_S_el;
    vec dd, v, a;
    vec temp;
    vec del_a, Delta_a;
    
    vec d_last, v_last, a_last;
    vec d_pred, v_pred;
    
    umat LM;
    
    struct stat info;
    
    const char * BI_file_Path, * PI_file_Path, * outputDir, * qdel_Path, * fem_inputs, * dem_inputs;
    
    ofstream d_file, v_file, a_file, t_file, stress_file, isv_file, F_S_file;
    char cCurrentPath[FILENAME_MAX];
    
    // Initialize MPI
    rc = MPI_Init(&argc, &argv);
    if(rc != MPI_SUCCESS) {
    	cout << "Error starting program with MPI..." << endl;
    	MPI_Abort(MPI_COMM_WORLD, rc);
    }

    if(argc > 7 || argc < 7) {
        cout << "Input Format: " << argv[0] << " <Path to Boundary Input File> <Path to Particle Input File> <Path to qdelauny> <Directory to Write Outputs> <Path to FEM Inputs> <Path to DEM Inputs>" << endl;
        return -1;
    }

    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	cout << "rank = " << rank << endl;
	cout << "numtasks = " << numtasks << endl;

    BI_file_Path = argv[1];
    PI_file_Path = argv[2];
    qdel_Path    = argv[3];
    outputDir    = argv[4];
    fem_inputs   = argv[5];
    dem_inputs   = argv[6];

    if(rank==0) {
    	cout << "Armadillo version: " << arma_version::as_string() << endl;
    	if(ifstream(BI_file_Path)) {
    		cout << "Boundary Input File found..." << endl;
    	}
    	else {
    		cout << "Boundary Input File not found...exiting..." << endl;
    		return -1;
    	}
    	if(ifstream(PI_file_Path)) {
    		cout << "Particle Input File found..." << endl;
    	}
    	else {
    		cout << "Particle Input File not found...exiting..." << endl;
    		return -1;
    	}
    	if(ifstream(qdel_Path)) {
    		cout << "qdelaunay binary found..." << endl;
    	}
    	else {
    		cout << "qdelaunay binary not found...exiting..." << endl;
    		return -1;
    	}
    	if(ifstream(fem_inputs)) {
    		cout << "FEM input file found..." << endl;
    	}
    	else {
    		cout << "FEM input file not found...exiting..." << endl;
    		return -1;
    	}
        if(ifstream(dem_inputs)) {
            cout << "DEM input file found..." << endl;
        }
        else {
            cout << "DEM input file not found...exiting..." << endl;
            return -1;
        }
    	// Get Current Directory Path where the programming is running
    	if (!getcwd(cCurrentPath, sizeof(cCurrentPath)))
    	{
    		return errno;
    	}

    	cout << "Current Directory is " << cCurrentPath << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //Read data from input file (fem_inputs)
    femInput femParams;
    femParams.readData(fem_inputs);

    //Read data from input file (dem_inputs)
    demInput demParams;
    demParams.readData(dem_inputs);

    //Elatic parameters taking from dry mason sand calibration effort
    lambda     = femParams.lambda;  // Pa
    mu         = femParams.mu;      // Pa
    rho        = femParams.rho;     // kg/m^3
    
    //Gravitational Acceleration
    grav       = femParams.grav;    // m/s^2
       
    //Geometry
    d          = femParams.d;       // m
    r          = d/2.0;            // m
    LDratio    = femParams.LDratio; // Guess based on Luo et al
    Area       = PI * pow(r,2.0);  // m^2
    h          = d*LDratio;

    // DEM Geometry
    h_DEM      = femParams.h_DEM;         // m
    w_DEM      = femParams.w_DEM;         // m
    l_DEM      = femParams.l_DEM;         // m
    A_DEM      = l_DEM * w_DEM; // m^2
    
    // FEM Constants
    numips     = femParams.numips;
    nstress    = femParams.nstress;
    nisv       = femParams.nisv;
    ndof       = femParams.ndof;
    nel        = femParams.nel;
    neldof     = femParams.neldof;

    // Time Parameters
    t          = femParams.t;
    dt         = femParams.dt;
    print_int  = femParams.print_int;
    n_print    = femParams.n_print;
    time_tot   = femParams.time_tot;
    nsteps     = round(time_tot/dt);
    t_ramp     = time_tot;

    // Boundary Conditions
    strainrate = femParams.strainrate;

    //Damping
    alphaM     = femParams.alphaM;

    params.set_size(12);
    params(0 ) = lambda;
    params(1 ) = mu;
    params(2 ) = rho;
    params(3 ) = grav;
    params(4 ) = numips;
    params(5 ) = nstress;
    params(6 ) = nisv;
    params(7 ) = Area;
    params(8 ) = h_DEM;
    params(9 ) = nel;
    params(10) = neldof;
    params(11) = ndof;

    createCoords(coords,params,h);
    createLM(LM,params);
    whichELIP(rank, el, ip);
    
    dispfun_time.zeros(nsteps);
    dispfun_disp.zeros(nsteps);
    dispfun_eps.zeros(nsteps);
    
    for(ii = 1; ii < nsteps; ii++) {
        dispfun_time(ii) = (ii-1.0) * time_tot/double(nsteps);
        dispfun_disp(ii) = -h * (exp(strainrate * dispfun_time(ii))-1.0);
        dispfun_eps(ii)  = log(1.0+dispfun_disp(ii)/h);
    }
    
    createG(g, dispfun_disp, params, 0);
    gd_n = 0.0;
    gd   = 0.0;

    // Concentrated External Force Vector
    traction_max = 0.0;
    traction     = 0.0;
    F_F          = traction*Area;

    // Initialization
    stress_el.zeros(nstress,numips,nel);
    mass_el.zeros(neldof,neldof,nel);
    stiff_el.zeros(neldof,neldof,nel);
    fg_el.zeros(neldof,nel);
    fs_el.zeros(neldof,nel);
    isv_el.zeros(nisv,numips,nel);
    d_el.zeros(nel,neldof);
    v_el.zeros(nel,neldof);
    a_el.zeros(nel,neldof);
    a_el_last.zeros(nel,neldof);
    Delta_a_el.zeros(nel,neldof);
    M.zeros(ndof,ndof);
    C.zeros(ndof,ndof);
    K.zeros(ndof,ndof);
    F.zeros(ndof);
    F_G.zeros(ndof);
    F_S.zeros(ndof);

    // Newmark Method Parameters - O(2), explicit
    beta = 0.0;
    gamma = 1.0/2.0;
    
    // Initial Conditions
    // Both M and C are constants

    for(ii = 0; ii < nel; ii++) {
        mass_el.slice(ii) = el_md_g1int((coords.row(ii)).t(), params);
        fg_el.col(ii)     = el_f_g4int( (coords.row(ii)).t(), params);
    }
    
    for(kk = 0; kk < nel; kk++) {
        temp = conv_to<vec>::from(LM.col(kk));
        for(ii = 0; ii < neldof; ii++) {
            I = int(temp(ii));
            if(I > 0) {
                //Note that indices in C++/Armadillo start at 0, not 1 as in Matlab
                F_G(I-1) = F_G(I-1) + fg_el(ii,kk);
                for(jj = 0; jj < neldof; jj++) {
                    J = int(temp(jj));
                    if(J > 0) {
                        M(I-1, J-1) = M(I-1, J-1) + mass_el(ii,jj,kk);
                        C(I-1, J-1) = C(I-1, J-1) + alphaM * mass_el(ii,jj,kk);
                    }
                }
            }
        }
    }

    F = (-F_F) - F_G;
    //Need to add damping to this initialization
    a = solve(M,-F);
    
    //Need to change if ICs are anything other than zero
    dd.zeros(ndof);
    v.zeros(ndof);
    a.zeros(ndof);
    
    // Make the necessary directories - only need 1 per IP per EL
    // DEM "snapshot" per call

    if(rank == 0) {

    	for(jj = 1; jj <= nel; jj++) {
    		for(ii = 1; ii <= numips; ii++) {
    			dirName = string(outputDir) + "/el" + to_string(jj) + "_ip" + to_string(ii);
    			cout << "*** Attemping to Creating Directory: " << dirName << endl;
    			if( stat( dirName.c_str(), &info ) != 0 ) {
    				printf( "*** Making directory: %s\n", dirName.c_str() );
    				mkdir(dirName.c_str(),0700);
    			}
    			else if( info.st_mode & S_IFDIR )  // S_ISDIR() doesn't exist on my windows
    				printf( "*** %s already exists ***\n", dirName.c_str() );
    			else
    				printf( "*** %s is no directory\n", dirName.c_str() );

    			copyfile(BI_file_Path, dirName+"/input_boundary_file");
    			copyfile(PI_file_Path, dirName+"/input_particle_file");

    			// ellip3D is now coupled to this code. no need to precompile versions,
    			// but need to update the fracture version

    			copyfile(qdel_Path, dirName+"/qdelaunay");
    			string qdel_Path = dirName + "/qdelaunay";
    			chmod(qdel_Path.c_str(), S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH);
    		}
    	}

    }
	
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0)
    {
        femParams.echoData();
        demParams.echoData();
        printMesh(coords, LM, g);
    }

    femParams.~femInput();

    MPI_Barrier(MPI_COMM_WORLD);

    for(n = 1; n <= nsteps; n++) {
    	if(rank == 0) {
    		cout << "*** Current step = " << n << " ***" << endl;
    	}

    	MPI_Barrier(MPI_COMM_WORLD);
    	
        t = t + dt;

        createG(g_n, dispfun_disp, params, n-1);

        if(t < t_ramp) {
            tract = traction_max * (t/t_ramp);
        }
        else {
            tract = traction_max;
        }
        
        createG(g, dispfun_disp, params, n);
        
        F_F = tract*Area;
        
        d_last = dd;
        v_last = v;
        a_last = a;
        
        //predictors
        d_pred = d_last + dt * v_last + pow(dt,2) * (1.0-2.0*beta) * a_last/2.0;
        v_pred = v_last + dt * (1.0-gamma) * a_last;
               
        F_S.zeros(ndof);
        K.zeros(ndof,ndof);
        
        dd = d_pred;
        v  = v_pred;
            
        d_el_last.zeros(nel,neldof);

        temp = conv_to<vec>::from(LM.col(el));
        for(ii = 0; ii < neldof; ii++) {
            I = temp(ii);
            if(I > 0) {
                d_el(el,ii) = dd(I-1);
                d_el_last(el,ii) = d_last(I-1);
            }
            else {
                d_el(el,ii) = g(ii,el);
                d_el_last(el,ii) = g_n(ii,el);
            }
        }

        if (n%print_int==0) {
            n_print++;
        }
        
        el_stress_ellip3d(outputDir, coords.row(el), d_el.row(el), d_el_last.row(el), params, n_print, -1, -1, el, ip, stress_el, isv_el, dt, demParams);
        //el_stress_isv(coords.row(el), d_el.row(el), params, el, ip, stress_el, isv_el);
            
        MPI_Barrier(MPI_COMM_WORLD);

        //does this need to be sent to each node?
        el_kd_g2int_ellip3d(outputDir, coords.row(el), d_el.row(el), params, n, el, stiff_el.slice(el));
        //el_kd_g2int(coords.row(el), d_el.row(el), params, el, stiff_el.slice(el));

        fs_el.col(el) = el_f_g2int(coords.row(el),stress_el.slice(el),params);
            
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
                        K(I-1,J-1) = K(I-1,J-1) + K_el(ii,jj);
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

        a = solve(M+gamma*dt*C,-(F_S-F_F-F_G)-C*v);
        v = v + dt*gamma*a;
        
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

         if(rank==0) {
            cout << "dd = " << dd << endl;
            cout << "v = " << v << endl;
            cout << "a = " << a << endl;
            cout << "t = " << t << endl;
            cout << "stress_el = " << stress_el << endl;
            cout << "isv_el = " << isv_el << endl;
            cout << "F_S = " << F_S << endl;
         }

         MPI_Barrier(MPI_COMM_WORLD);
        
        if(rank==0 && n%print_int==0) {
			chdir(cCurrentPath);
    	    d_file.open("d.txt",ofstream::out | ofstream::app);
    	    v_file.open("v.txt",ofstream::out | ofstream::app);
    	    a_file.open("a.txt",ofstream::out | ofstream::app);
    	    t_file.open("t.txt",ofstream::out | ofstream::app);
    	    stress_file.open("stress.txt",ofstream::out | ofstream::app);
    	    isv_file.open("isv_el.txt",ofstream::out | ofstream::app);
    	    F_S_file.open("F_S.txt",ofstream::out | ofstream::app);

            d_file.precision(5);
            d_file.setf(ios::scientific);
            v_file.precision(5);
    	    v_file.setf(ios::scientific);
            a_file.precision(5);
    	    a_file.setf(ios::scientific);
            stress_file.precision(5);
    	    stress_file.setf(ios::scientific);
            isv_file.precision(5);
    	    isv_file.setf(ios::scientific);
            F_S_file.precision(5);
    	    F_S_file.setf(ios::scientific);
	    
            dd.raw_print(d_file);
    	    v.raw_print(v_file);
    	    a.raw_print(a_file);
            F_S.raw_print(F_S_file);
    	    t_file << t << endl;
			for (ii=0; ii<nel; ii++){
				for (jj=0; jj<numips; jj++){
					for (kk=0; kk<nstress; kk++){
						stress_file << stress_el(kk,jj,ii) << " ";
					}
					for (ll=0; ll<nisv; ll++){
						isv_file << isv_el(ll,jj,ii) << " ";
					}
				}
			}
			stress_file << endl;
			isv_file << endl;

    	    d_file.close();
    	    v_file.close();
    	    a_file.close();
    	    t_file.close();
     	    stress_file.close();
    	    isv_file.close();
    	    F_S_file.close();
        }

        MPI_Barrier(MPI_COMM_WORLD);

        stress_el.zeros();
        isv_el.zeros();
    }
    
    MPI_Finalize();
    return 0;
}
