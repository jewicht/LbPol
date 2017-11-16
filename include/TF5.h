
class TF5
{
 public:
  TF3(const char *name, Double_t (*fcn)(const Double_t *, const Double_t *), Double_t xmin=0, Double_t xmax=1, Double_t ymin=0,
       Double_t ymax=1, Double_t zmin=0, Double_t zmax=1, Int_t npar=0);

  double Eval(double,double,double,double,double);
  
};
