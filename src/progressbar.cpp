#include "progressbar.h"

progressbar::progressbar(const std::string& _prefix)
{
  prefix=_prefix;
  prevperc=-1;
  seconds_0 = time(NULL);
  print(0.);
}

void progressbar::print(double perc)
{
  const int iperc=int(perc);
  if (iperc!=prevperc) {
    //    [>>>>>>>>>>>>>>>>>......] (74%, time left: 15 sec)
    std::string pbar;
    pbar.append("[");
    for (int i=0;i<barlength;i++) if (int(100*i/barlength)<perc) {pbar.append(">"); } else {pbar.append(".");}
    pbar.append("]");

    char ipercstr[50];
    sprintf(ipercstr,"%2d%%", iperc);

    time_t seconds_1=time(NULL);
    //did perc in (seconds_1 - seconds_0)
    const double speed= seconds_1!=seconds_0?perc/(seconds_1 - seconds_0):1.;
    const int timeleft = int((100.-perc)/speed);

    std::cerr << prefix << ": " << pbar << " (" << ipercstr  << ", time left: " << timeleft  << " sec)      \r";
    prevperc=iperc;
  }
}

void progressbar::finish()
{
  std::cerr << prefix << ": done                                             " << std::endl;
}
