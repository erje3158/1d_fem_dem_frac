#include <iostream>
#include <fstream>
#include <string>

using namespace std;

class femInput
{
  public:
  
  double lambda;         // Pa
  double mu;             // Pa
  double rho;            // kg/m^3

  double grav;           // m/s^2

  double d;              // m
  double LDratio;        // m/m

  double h_DEM;          // m
  double w_DEM;          // m
  double l_DEM;          // m

  int numips;            // # of Integration Points
  int nstress;           // # of Calculated Stresses
  int nisv;              // # of Internal State Variables
  int ndof;              // # of Degrees of Freedom
  int nel;               // # of Elements
  int neldof;            // # of Element Degrees of Freedom

  double t;              // s
  double dt;             // s
  int print_int;         // # of Global Iterations per FEM Output
  int n_print;           // # of Global Iterations per DEM Snapshot
  double time_tot;       // s

  double strainrate;     // s^-1

  double alphaM;         // Mass Proportional Parameter

  int whichDisp;         // If 0 - error
                         // If 1 - applied finite displacement
                         // If 2 - "correct" SHPB displacement

  int whichConst;        // If 0 - error
                         // If 1 - hyperelasticity
                         // If 2 - ellip3d DEM

  femInput()
  {
    lambda       = 0.0; 
    mu           = 0.0; 
    rho          = 0.0; 

    grav         = 0.0;

    d            = 0.0; 
    LDratio      = 0.0; 

    h_DEM        = 0.0; 
    w_DEM        = 0.0; 
    l_DEM        = 0.0; 

    numips       = 0;
    nstress      = 0;
    nisv         = 0;
    ndof         = 0;
    nel          = 0;
    neldof       = 0;  

    t            = 0.0;
    dt           = 0.0;
    print_int    = 0;
    n_print      = 0;
    time_tot     = 0.0; 
    
    strainrate   = 0.0;

    alphaM       = 0.0; 

    whichDisp    = 0;

    whichConst   = 0;
  }

  ~femInput();

  void readData(const char * inputFile);
  void echoData();
  void checkData();

};

class demInput
{
public:
  double maxOverlap;

  double youngsMod;
  double poisRatio;

  double timestep;

  double damping;
  double friction;

  double sigmaComp;
  double tensileCrit;
  double fracTough;
  int    isFrac;

  int    whichSeed; 

  demInput()
  {
    maxOverlap   = 0.0;

    youngsMod    = 0.0;
    poisRatio    = 0.0;

    timestep     = 0.0;

    damping      = 0.0;
    friction     = 0.0;

    sigmaComp    = 0.0;
    tensileCrit  = 0.0;
    fracTough    = 0.0;
    isFrac       = 0;

    whichSeed    = 0;
  }

  ~demInput();

  void readData(const char * inputFile);
  void echoData();
  void checkData();

};
