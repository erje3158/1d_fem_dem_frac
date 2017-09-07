#!/bin/bash

# Compilation script for ubuntu
# gcc -Wall userInput.cpp -o testrun -lstdc++

# Compilatoin script for max
g++ -c main.cpp
g++ -c userInput.cpp
g++ -c meshTools.cpp
g++ -o testrun main.o userInput.o meshTools.o
