#ifndef PROGRESSBAR
#define PROGRESSBAR

#include <iostream>
#include <string>
#include <stdio.h>

class progressbar
{
 public:
  progressbar(const std::string&);
  ~progressbar(){};
  void print(double);
  void finish();
 private:
  std::string prefix;
  int prevperc;
  time_t seconds_0;
  static const int barlength=20;
};


#endif
