//
//  copyfile.cpp
//  Jensen_code
//
//  Created by Christopher Kung on 1/29/16.
//  Copyright © 2016 Christopher Kung. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>

using namespace std;

void copyfile (string inputFile, string outputFile) {
    int c;
    FILE *source_flag, *destination_flag;
    
    source_flag = fopen(inputFile.c_str(), "r");
    destination_flag = fopen(outputFile.c_str(),"w");
    
    while((c=getc(source_flag))!= EOF) {
        putc(c,destination_flag);
    }
    
    fclose(source_flag);
    fclose(destination_flag);
}