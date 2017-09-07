//
//  meshTools.cpp
//  Jensen_code
//
//  Created by Erik Jensen 8/25/2017.
//  Copyright �� 2017 Erik Jensen. All rights reserved.
//

#include "meshTools.h"

using namespace std;
using namespace arma;

void createCoords(mat & coords, vec params, double h)
{

	coords.set_size(params(9),params(10));

	int dof = 0;

	for (int i = 0; i < params(10); i++)
	{
		for (int j = 0; j < params(9); j++)
		{
			coords(j,i) = dof * h / params(9);
			dof++;
		}
		dof = 1;
	}

}

void createLM(umat & LM, vec params)
{

	LM.set_size(params(10),params(9));

	int dof  = 0;

	for (int i = 0; i < params(10); i++)
	{
		for (int j = 0; j < params(9); j++)
		{
			LM(i,j) = dof;
			dof++;
		}
		dof = 1;
	}

	//Boundary Condition
	LM(params(10)-1,params(9)-1) = 0;
}

void createG(mat & g, vec disp, vec params, int n)
{

	g.zeros(params(10),params(9));
	g(params(10)-1,params(9)-1) = disp(n);

}

void printMesh(mat coords, umat LM, mat g)
{
	cout << endl << endl;
 	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
  	cout << "SUMMARY OF MESH:"                          << endl;
  	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
  	cout << endl;
  	cout << "coords = " 								<< endl;
  	coords.print();
    cout << endl;
  	cout << "LM = " 									<< endl;
  	LM.print();
  	cout << endl;
  	cout << "g = " 										<< endl;
  	g.print();  
  	cout << endl;	
  	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	cout << endl << endl << endl;
}

void whichELIP(int rank, int & el, int & ip)
{

	el = rank / 2;

	int mod = rank % 2;

	if (mod == 0) 
	{
		ip = 0;
	}
	else if (mod == 1)
	{
		ip = 1;
	}
	else {
		cout << endl << endl << "ERROR: Check meshTools::whichELIP()" << endl << endl;
	}
}

void printELIP(int rank, int el, int ip)
{
	
	cout << endl;
	cout << "rank = " << rank << " :: el = " << el << " :: ip = " << ip;
	cout << endl;

}

