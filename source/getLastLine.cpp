//
//  getLastLine.cpp
//  Jensen_code
//
//  Created by Christopher Kung on 2/8/16.
//  Copyright © 2016 Christopher Kung. All rights reserved.
//

#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

string getLastLine(std::ifstream& in)
{
    std::string line;
    while (in >> std::ws && std::getline(in, line)) // skip empty lines
        ;
    
    return line;
}