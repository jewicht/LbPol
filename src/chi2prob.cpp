#include "TMath.h"
#include <iostream>

int main(int argc, char** argv)
{
  std::cout << TMath::Prob(atof(argv[1]),atoi(argv[2])) << std::endl;
}
