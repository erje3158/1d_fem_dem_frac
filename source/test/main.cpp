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

int main (int argc, char * argv[]) 
{
  
  const char * inputFile = argv[1];
  
  userInput myinput;
  myinput.readData(inputFile);
  myinput.echoData();

  double testVar = myinput.alphaM;
  cout << "testVar = " << testVar << endl;

}
