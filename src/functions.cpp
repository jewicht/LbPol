#include "functions.h"


// void mygetchar(std::string file, int line)
// {
//   if (line==0) {
//     std::cout << "Waiting for keypress..." << std::endl;
//   } else {
//     std::cout << "Waiting for keypress... (" << file << ":" << line << ")" << std::endl;
//   }

//   getchar();
// }

TString istr(int x,TString fill, bool zerohack)
{
  if (zerohack and (x==0)) return "--";
  char buffer1[50];
  //  int n = 
  sprintf(buffer1, "%" + fill + "d", x);
  return TString(buffer1);
}

TString str(double x)
{
  char buffer1[50];
  //  int n = 
  sprintf(buffer1, "%f", x);
  return TString(buffer1);
}


TString nicedate()
{
  time_t rawtime;
  struct tm * timeinfo;
  char buffer[80];

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );

  strftime (buffer,80,"%Y%m%d-%H%M",timeinfo);
  return TString(buffer);
}



