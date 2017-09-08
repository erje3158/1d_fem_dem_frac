//
//  userinput.cpp
//  Jensen_code
//
//  Created by Erik Jensen 8/11/2017.
//  Copyright �� 2017 Erik Jensen. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <string>
#include "userInput.h"

using namespace std;

femInput::~femInput(void)
{
  //cout << endl << endl << "femInput is being destroyed" << endl << endl;
}

//Read data from input file
void femInput::readData(const char * inputFile)
{
  string line; ifstream input(inputFile);
  while ( getline(input,line) )
    {
      if ( line == "$FEM Constitutive"  ) { input >> this->lambda     ;
                                            input >> this->mu         ;
                                            input >> this->rho        ;}
      if ( line == "$Gravity"           ) { input >> this->grav       ;}
      if ( line == "$FEM Geometry"      ) { input >> this->d          ;
                                            input >> this->LDratio    ;}
      if ( line == "$DEM Geometry"      ) { input >> this->h_DEM      ;
                                            input >> this->w_DEM      ;
                                            input >> this->l_DEM      ;}
      if ( line == "$FEM Constants"     ) { input >> this->numips     ;
                                            input >> this->nstress    ;
                                            input >> this->nisv       ;
                                            input >> this->ndof       ;
                                            input >> this->nel        ;
                                            input >> this->neldof     ;}
      if ( line == "$Time Parameters"   ) { input >> this->t          ;
                                            input >> this->dt         ;
                                            input >> this->print_int  ;
                                            input >> this->n_print    ;
                                            input >> this->time_tot   ;}
      if ( line == "$Strain Rate"       ) { input >> this->strainrate ;}
      if ( line == "$Mass Damping"      ) { input >> this->alphaM     ;}
    }
  input.close();
}

//Print user inputs for review
void femInput::echoData()
{
  cout << endl << endl;
  cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
  cout << "SUMMARY OF INPUTS:"                        << endl;
  cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
  cout << endl;
  cout << "FEM Constitutive:"                         << endl;
  cout << "   lambda         = " << this->lambda      << endl;
  cout << "   mu             = " << this->mu          << endl;
  cout << "   rho            = " << this->rho         << endl;
  cout << endl;  
  cout << "Gravity:"                                  << endl;
  cout << "   grav           = " << this->grav        << endl;
  cout << endl;
  cout << "FEM Geometry:"                             << endl;
  cout << "   d              = " << this->d           << endl;
  cout << "   LDratio        = " << this->LDratio     << endl;
  cout << endl;
  cout << "DEM Geometry:"                             << endl;
  cout << "   h_DEM          = " << this->h_DEM       << endl;
  cout << "   w_DEM          = " << this->w_DEM       << endl;
  cout << "   l_DEM          = " << this->l_DEM       << endl;
  cout << endl;
  cout << "FEM Constants:"                            << endl;
  cout << "   numips         = " << this->numips      << endl;
  cout << "   nstress        = " << this->nstress     << endl;
  cout << "   nisv           = " << this->nisv        << endl;
  cout << "   ndof           = " << this->ndof        << endl;
  cout << "   nel            = " << this->nel         << endl;
  cout << "   neldof         = " << this->neldof      << endl;
  cout << endl;
  cout << "Time Parameters:"                          << endl;
  cout << "   t              = " << this->t           << endl;
  cout << "   dt             = " << this->dt          << endl;
  cout << "   print_int      = " << this->print_int   << endl;
  cout << "   n_print        = " << this->n_print     << endl;
  cout << "   time_tot       = " << this->time_tot    << endl;
  cout << endl;
  cout << "Strain Rate:"                              << endl;
  cout << "   strainrate     = " << this->strainrate  << endl;
  cout << endl;
  cout << "Mass Damping:"                             << endl;
  cout << "   alphaM         = " << this->alphaM      << endl;
  cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
  cout << endl << endl << endl;
}

demInput::~demInput(void)
{
  //cout << endl << endl << "demInput is being destroyed" << endl << endl;
}

//Read data from input file
void demInput::readData(const char * inputFile)
{
  string line; ifstream input(inputFile);
  while ( getline(input,line) )
    {
      if ( line == "$Particle Overlap"  ) { input >> this->maxOverlap ;}
      if ( line == "$DEM Constitutive"  ) { input >> this->youngsMod  ;
                                            input >> this->poisRatio  ;}
      if ( line == "$Time Parameters"   ) { input >> this->timestep   ;}
      if ( line == "$Particle Contact"  ) { input >> this->damping    ;
                                            input >> this->friction   ;}
      if ( line == "$Partcile Fracture" ) { input >> this->sigmaCrit  ;
                                          { input >> this->sigmaComp  ;
                                            input >> this->tensileCrit;
                                            input >> this->sigmaF     ;}
    }
  input.close();
}

//Print user inputs for review
void demInput::echoData()
{
  cout << endl << endl;
  cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
  cout << "SUMMARY OF INPUTS:"                        << endl;
  cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
  cout << endl;
  cout << "Particle Overlap:"                         << endl;
  cout << "   maxOverlap     = " << this->maxOverlap  << endl;
  cout << endl;  
  cout << "DEM Constitutive:"                         << endl;
  cout << "   youngsMod      = " << this->youngsMod   << endl;
  cout << "   poisRatio      = " << this->poisRatio   << endl;
  cout << endl;
  cout << "Time Parameters:"                          << endl;
  cout << "   timestep       = " << this->timestep    << endl;
  cout << endl;
  cout << "Particle Contact:"                         << endl;
  cout << "   damping        = " << this->damping     << endl;
  cout << "   friction       = " << this->friction    << endl;
  cout << endl;
  cout << "Particle Fracture:"                        << endl;
  cout << "   sigmaCrit      = " << this->sigmaCrit   << endl;
  cout << "   sigmaComp      = " << this->sigmaComp   << endl;
  cout << "   tensileCrit    = " << this->tensileCrit << endl;
  cout << "   sigmaF         = " << this->sigmaF      << endl;
  cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
  cout << endl << endl << endl;
}


