#include "angles.h"

TLorentzVector myboost(const TLorentzVector& a, const TVector3& aboost)
{
  TLorentzVector tmp(a);
  tmp.Boost(aboost);
  return tmp;
}

void calculatethetaNphiN(Float_t& costhetaN, Float_t& phiN, const TVector3& n, const TLorentzVector& mother, const TLorentzVector& daug1)//, const TLorentzVector& daug2)
{
  const TVector3 z1=mother.Vect().Unit();
  const TVector3 y1=(n.Cross(z1)).Unit();
  const TVector3 x1=(y1.Cross(z1)).Unit();
  
  const TVector3 boostToMother = -mother.BoostVector();
  const TLorentzVector lv_daug1_Motherframe = myboost(daug1, boostToMother);
  const TVector3 xyz_daug1_Motherframe = lv_daug1_Motherframe.Vect().Unit();
  costhetaN = z1.Dot(xyz_daug1_Motherframe);
  
  const TVector3 xyz_daug1_Motherframe_projete = (z1.Cross(xyz_daug1_Motherframe.Cross(z1))).Unit();
  const double cosphiN = x1.Dot(xyz_daug1_Motherframe_projete);
  const double sinphiN = y1.Dot(xyz_daug1_Motherframe_projete);
  phiN = sinphiN>0.0 ? acos(cosphiN) : -acos(cosphiN);
}


void calculateangles(const TLorentzVector& lv_Lambdab, 
		     const TLorentzVector& lv_Lambda0, const TLorentzVector& lv_proton, const TLorentzVector& lv_pion, 
		     const TLorentzVector& lv_Jpsi, const TLorentzVector& lv_muonp, const TLorentzVector& lv_muonm,
		     Float_t& costheta, Float_t& costheta1, Float_t& costheta2, Float_t& phi1, Float_t& phi2, const TVector3& proddir)
{
  costheta=0.;
  costheta1=0.;
  costheta2=0.;
  phi1=0.;
  phi2=0.;

  //  const TLorentzVector lv_Lambdab(lv_proton+lv_pion+lv_muonp+lv_muonm);
  //  const TLorentzVector lv_Lambda0(lv_proton+lv_pion);
  //  const TLorentzVector lv_Jpsi(lv_muonp+lv_muonm);

  const TVector3 xyz_Lambdab = lv_Lambdab.Vect().Unit();
  const TVector3 n((proddir.Cross(xyz_Lambdab)).Unit());
  
  //first boost everything in Lb frame

  const TVector3  boostToLambdab = -lv_Lambdab.BoostVector();
  const TLorentzVector lv_Lambdab_Lambdabframe = myboost(lv_Lambdab, boostToLambdab);
  const TLorentzVector lv_Lambda0_Lambdabframe = myboost(lv_Lambda0, boostToLambdab);
  const TLorentzVector lv_proton_Lambdabframe  = myboost(lv_proton,  boostToLambdab);
  const TLorentzVector lv_pion_Lambdabframe    = myboost(lv_pion,    boostToLambdab);
  const TLorentzVector lv_Jpsi_Lambdabframe    = myboost(lv_Jpsi,    boostToLambdab);
  const TLorentzVector lv_muonp_Lambdabframe   = myboost(lv_muonp,   boostToLambdab);
  const TLorentzVector lv_muonm_Lambdabframe   = myboost(lv_muonm,   boostToLambdab);
  
  {
    //theta
    const TVector3 xyz_Lambda0_Lambdabframe = TVector3(lv_Lambda0_Lambdabframe.Vect()).Unit();
    const TVector3 xyz_Jpsi_Lambdabframe = TVector3(lv_Jpsi_Lambdabframe.Vect()).Unit();
    costheta = n.Dot(xyz_Lambda0_Lambdabframe);
  }

  calculatethetaNphiN(costheta1,phi1,n,lv_Lambda0_Lambdabframe,lv_proton_Lambdabframe);
  calculatethetaNphiN(costheta2,phi2,n,lv_Jpsi_Lambdabframe,lv_muonp_Lambdabframe);

}

void calculateangles(const TLorentzVector& lv_Lambdab, 
		     const TLorentzVector& lv_Lambda0, const TLorentzVector& lv_proton, const TLorentzVector& lv_pion, 
		     const TLorentzVector& lv_Jpsi, const TLorentzVector& lv_muonp, const TLorentzVector& lv_muonm,
		     Float_t& costheta, Float_t& costheta1, Float_t& costheta2, Float_t& phi1, Float_t& phi2)
{
  const TVector3 proddir(0.,0.,1.);
  calculateangles(lv_Lambdab, 
		  lv_Lambda0, lv_proton, lv_pion, 
		  lv_Jpsi, lv_muonp, lv_muonm,
		  costheta, costheta1, costheta2, phi1, phi2, proddir);
}

TString lvcout(const TLorentzVector& lv)
{
  return "(" + str(lv.Px()) + ", " + str(lv.Py()) + ", " + str(lv.Pz()) + ", " +str(lv.E()) + ")";
}

TString lvcout(const TVector3& lv)
{
  return "(" + str(lv.Px()) + ", " + str(lv.Py()) + ", " + str(lv.Pz()) +  ")";
}

