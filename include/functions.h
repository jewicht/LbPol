#include "TString.h"

#ifndef JWFUNCTION
#define JWFUNCTION

#define mygetchar std::cout << "Press enter (" << __FILE__ << ":" << __LINE__ << " [" << __PRETTY_FUNCTION__ << "])" << std::endl;getchar();

//void mygetchar(std::string file="", int line=0);

//TString istr(int,const TString& fill="", bool zerohack=false);
//TString str(double);
//TString nicedate();

//#include "functions.h"


// void mygetchar(std::string file, int line)
// {
//   if (line==0) {
//     std::cout << "Waiting for keypress..." << std::endl;
//   } else {
//     std::cout << "Waiting for keypress... (" << file << ":" << line << ")" << std::endl;
//   }

//   getchar();
// }


inline int chtoint(char c)
{
  return int(c)-48;
}


inline TString istr(int x, const TString& fill="", bool zerohack=false)
{
  if (zerohack and (x==0)) return "--";
  char buffer1[50];
  //  int n = 
  sprintf(buffer1, "%" + fill + "d", x);
  return TString(buffer1);
}

inline TString str(double x)
{
  char buffer1[50];
  //  int n = 
  sprintf(buffer1, "%f", x);
  return TString(buffer1);
}

inline bool atob(const std::string& mystr)
{
  if (mystr=="true") {
    return true; 
  } else {
    return false;
  }
}

inline std::string btoa(bool blu)
{
  if (blu) return "true";
  return "false";
}

inline TString nicedate()
{
  time_t rawtime;
  struct tm * timeinfo;
  char buffer[80];

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );

  strftime (buffer,80,"%Y%m%d-%H%M",timeinfo);
  return TString(buffer);
}



#endif
