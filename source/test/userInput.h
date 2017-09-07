#include <iostream>
#include <fstream>
#include <string>

using namespace std;

class userInput
{
 public:
  
  double alphaM;

  userInput()
    {
      alphaM = 0.0;
    }

  void readData(const char * inputFile);
  void echoData();

};
