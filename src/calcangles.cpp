//g++ -o calcangles calcangles.C -O2 `root-config --cflags --libs` 

#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "functions.h"
#include "progressbar.h"
#include "angles.h"

//for stripping17
typedef Float_t vartype;

//for stripping19
//typedef Double_t vartype;

const bool debug=false;


double pseudorapidity(const TLorentzVector& lv_lambdab)
{
  return -log(tan(0.5* lv_lambdab.Vect().Theta() ));
}

double rapidity(const TLorentzVector& lv_lambdab)
{
  const double E = lv_lambdab.E();
  const double pz = lv_lambdab.Pz();
  return 0.5*log( (E+pz)/(E-pz) );
}

void ksangle(const TLorentzVector& lv_proton, const TLorentzVector& lv_pion, Float_t& phi)
{
  const TVector3 x1(1.,0.,0.);
  const TVector3 y1(0.,1.,0.);
  const TVector3 proddir(0.,0.,1.);
  
  const TVector3 xyz_pion = lv_pion.Vect().Unit();
  const TVector3 xyz_proton = lv_proton.Vect().Unit();

  const TVector3 n(xyz_pion.Cross(xyz_proton).Unit());
  const TVector3 n_projete = proddir.Cross(n.Cross(proddir)).Unit();
  
  const double cosphi1 = x1.Dot(n_projete);
  const double sinphi1 = y1.Dot(n_projete);
  //  phi=asin(n.Dot(proddir));    

  phi = sinphi1>0.0 ? acos(cosphi1) : -acos(cosphi1);
}

int main(int argc, char** argv)
{

  if (argc!=2) {
    std::cerr << argv[0] << " file" << std::endl;
    return 1;
  }

  const double pi=atan(1.)*4.;

  const TString filename(argv[1]);

  TFile myfile(filename);


  const bool isMC=filename.Contains("MC");

  vartype
    lambdab_PE,
    lambdab_PX,
    lambdab_PY,
    lambdab_PZ,

    jpsi1s_PE,
    jpsi1s_PX,
    jpsi1s_PY,
    jpsi1s_PZ,

    lambda0_PE,
    lambda0_PX,
    lambda0_PY,
    lambda0_PZ,

    muplus_PE,
    muplus_PX,
    muplus_PY,
    muplus_PZ,
    
    muminus_PE,
    muminus_PX,
    muminus_PY,
    muminus_PZ,
    
    piminus_PE,
    piminus_PX,
    piminus_PY,
    piminus_PZ,
    
    pplus_PE,
    pplus_PX,
    pplus_PY,
    pplus_PZ;


 vartype
    lambdab_TRUEP_E,
    lambdab_TRUEP_X,
    lambdab_TRUEP_Y,
    lambdab_TRUEP_Z,

    jpsi1s_TRUEP_E,
    jpsi1s_TRUEP_X,
    jpsi1s_TRUEP_Y,
    jpsi1s_TRUEP_Z,

    lambda0_TRUEP_E,
    lambda0_TRUEP_X,
    lambda0_TRUEP_Y,
    lambda0_TRUEP_Z,

    muplus_TRUEP_E,
    muplus_TRUEP_X,
    muplus_TRUEP_Y,
    muplus_TRUEP_Z,
    
    muminus_TRUEP_E,
    muminus_TRUEP_X,
    muminus_TRUEP_Y,
    muminus_TRUEP_Z,
    
    piminus_TRUEP_E,
    piminus_TRUEP_X,
    piminus_TRUEP_Y,
    piminus_TRUEP_Z,
    
    pplus_TRUEP_E,
    pplus_TRUEP_X,
    pplus_TRUEP_Y,
    pplus_TRUEP_Z;



  Int_t
    lambdab_ID,
    jpsi1s_ID,
    lambda0_ID,
    muplus_ID,
    muminus_ID,
    piminus_ID,
    pplus_ID;


  Int_t
    lambdab_TRUEID,
    jpsi1s_TRUEID,
    lambda0_TRUEID,
    muplus_TRUEID,
    muminus_TRUEID,
    piminus_TRUEID,
    pplus_TRUEID;



  Float_t costheta, costheta1, costheta2, phi1, phi2, ks0_phi,lambdab0_eta, lambda0_eta;
  Float_t true_costheta, true_costheta1, true_costheta2, true_phi1, true_phi2;
  Float_t lambdab0_y, lambda0_y, piminus_y, pplus_y;

  //  Float_t F[8];

  std::vector<TString> directorylist;
  directorylist.push_back("Tuple_JpsiL0_betas");
  directorylist.push_back("Tuple_JpsiL0_betas_refitted");
  directorylist.push_back("Tuple_JpsiL0_B2XMuMu");
  directorylist.push_back("Tuple_JpsiL0_B2XMuMu_refitted");
  directorylist.push_back("Tuple_JpsiKs_detached_betas");
  directorylist.push_back("Tuple_JpsiKs_prescaled_betas");
  directorylist.push_back("Tuple_JpsiKs_detached_betas_refitted");
  directorylist.push_back("Tuple_JpsiKs_B2XMuMu");
  directorylist.push_back("Tuple_JpsiKs_B2XMuMu_refitted");




  TRandom3 rnd(0);


  TString output(filename);
  output.ReplaceAll(".root", ".angles.root");
  TFile *f = NULL;
  if (!debug) f = new TFile(output,"recreate");

  for (unsigned int idl=0; idl<directorylist.size(); idl++) {

    TDirectory* mydirectory = (TDirectory*)myfile.Get(directorylist[idl]);//.c_str());
    if (!mydirectory) continue;
    TTree* mytree = (TTree*)mydirectory->Get("DecayTree");

    mytree->SetBranchAddress("jpsi1s_ID", &jpsi1s_ID);
    mytree->SetBranchAddress("jpsi1s_PX", &jpsi1s_PX);
    mytree->SetBranchAddress("jpsi1s_PY", &jpsi1s_PY);
    mytree->SetBranchAddress("jpsi1s_PZ", &jpsi1s_PZ);
    mytree->SetBranchAddress("jpsi1s_PE", &jpsi1s_PE);

    mytree->SetBranchAddress("muplus_ID", &muplus_ID);
    mytree->SetBranchAddress("muplus_PX", &muplus_PX);
    mytree->SetBranchAddress("muplus_PY", &muplus_PY);
    mytree->SetBranchAddress("muplus_PZ", &muplus_PZ);
    mytree->SetBranchAddress("muplus_PE", &muplus_PE);
    
    mytree->SetBranchAddress("muminus_ID", &muminus_ID);
    mytree->SetBranchAddress("muminus_PX", &muminus_PX);
    mytree->SetBranchAddress("muminus_PY", &muminus_PY);
    mytree->SetBranchAddress("muminus_PZ", &muminus_PZ);
    mytree->SetBranchAddress("muminus_PE", &muminus_PE);



    const bool isJpsiKs=directorylist[idl].Contains("JpsiKs");
    if (isJpsiKs) {
      mytree->SetBranchAddress("piplus_ID", &pplus_ID);
      mytree->SetBranchAddress("piplus_PX", &pplus_PX);
      mytree->SetBranchAddress("piplus_PY", &pplus_PY);
      mytree->SetBranchAddress("piplus_PZ", &pplus_PZ);
      mytree->SetBranchAddress("piplus_PE", &pplus_PE);
      
      mytree->SetBranchAddress("ks0_ID", &lambda0_ID);
      mytree->SetBranchAddress("ks0_PX", &lambda0_PX);
      mytree->SetBranchAddress("ks0_PY", &lambda0_PY);
      mytree->SetBranchAddress("ks0_PZ", &lambda0_PZ);
      mytree->SetBranchAddress("ks0_PE", &lambda0_PE);
      
      mytree->SetBranchAddress("b0_ID", &lambdab_ID);
      mytree->SetBranchAddress("b0_PX", &lambdab_PX);
      mytree->SetBranchAddress("b0_PY", &lambdab_PY);
      mytree->SetBranchAddress("b0_PZ", &lambdab_PZ);
      mytree->SetBranchAddress("b0_PE", &lambdab_PE);
      
    } else {
      mytree->SetBranchAddress("pplus_ID", &pplus_ID);
      mytree->SetBranchAddress("pplus_PX", &pplus_PX);
      mytree->SetBranchAddress("pplus_PY", &pplus_PY);
      mytree->SetBranchAddress("pplus_PZ", &pplus_PZ);
      mytree->SetBranchAddress("pplus_PE", &pplus_PE);

      mytree->SetBranchAddress("lambda0_ID", &lambda0_ID);
      mytree->SetBranchAddress("lambda0_PX", &lambda0_PX);
      mytree->SetBranchAddress("lambda0_PY", &lambda0_PY);
      mytree->SetBranchAddress("lambda0_PZ", &lambda0_PZ);
      mytree->SetBranchAddress("lambda0_PE", &lambda0_PE);
      
      mytree->SetBranchAddress("lambda_b0_ID", &lambdab_ID);
      mytree->SetBranchAddress("lambda_b0_PX", &lambdab_PX);
      mytree->SetBranchAddress("lambda_b0_PY", &lambdab_PY);
      mytree->SetBranchAddress("lambda_b0_PZ", &lambdab_PZ);
      mytree->SetBranchAddress("lambda_b0_PE", &lambdab_PE);
    }

    mytree->SetBranchAddress("piminus_ID", &piminus_ID);
    mytree->SetBranchAddress("piminus_PX", &piminus_PX);
    mytree->SetBranchAddress("piminus_PY", &piminus_PY);
    mytree->SetBranchAddress("piminus_PZ", &piminus_PZ);
    mytree->SetBranchAddress("piminus_PE", &piminus_PE);
    


    mytree->SetBranchAddress("jpsi1s_ID", &jpsi1s_ID);
    mytree->SetBranchAddress("jpsi1s_PX", &jpsi1s_PX);
    mytree->SetBranchAddress("jpsi1s_PY", &jpsi1s_PY);
    mytree->SetBranchAddress("jpsi1s_PZ", &jpsi1s_PZ);
    mytree->SetBranchAddress("jpsi1s_PE", &jpsi1s_PE);

    mytree->SetBranchAddress("muplus_ID", &muplus_ID);
    mytree->SetBranchAddress("muplus_PX", &muplus_PX);
    mytree->SetBranchAddress("muplus_PY", &muplus_PY);
    mytree->SetBranchAddress("muplus_PZ", &muplus_PZ);
    mytree->SetBranchAddress("muplus_PE", &muplus_PE);
    
    mytree->SetBranchAddress("muminus_ID", &muminus_ID);
    mytree->SetBranchAddress("muminus_PX", &muminus_PX);
    mytree->SetBranchAddress("muminus_PY", &muminus_PY);
    mytree->SetBranchAddress("muminus_PZ", &muminus_PZ);
    mytree->SetBranchAddress("muminus_PE", &muminus_PE);



    if (isMC) {
      if (isJpsiKs) {
	mytree->SetBranchAddress("piplus_TRUEID", &pplus_TRUEID);
	mytree->SetBranchAddress("piplus_TRUEP_X", &pplus_TRUEP_X);
	mytree->SetBranchAddress("piplus_TRUEP_Y", &pplus_TRUEP_Y);
	mytree->SetBranchAddress("piplus_TRUEP_Z", &pplus_TRUEP_Z);
	mytree->SetBranchAddress("piplus_TRUEP_E", &pplus_TRUEP_E);
	
	mytree->SetBranchAddress("ks0_TRUEID", &lambda0_TRUEID);
	mytree->SetBranchAddress("ks0_TRUEP_X", &lambda0_TRUEP_X);
	mytree->SetBranchAddress("ks0_TRUEP_Y", &lambda0_TRUEP_Y);
	mytree->SetBranchAddress("ks0_TRUEP_Z", &lambda0_TRUEP_Z);
	mytree->SetBranchAddress("ks0_TRUEP_E", &lambda0_TRUEP_E);
	
	mytree->SetBranchAddress("b0_TRUEID", &lambdab_TRUEID);
	mytree->SetBranchAddress("b0_TRUEP_X", &lambdab_TRUEP_X);
	mytree->SetBranchAddress("b0_TRUEP_Y", &lambdab_TRUEP_Y);
	mytree->SetBranchAddress("b0_TRUEP_Z", &lambdab_TRUEP_Z);
	mytree->SetBranchAddress("b0_TRUEP_E", &lambdab_TRUEP_E);
	
      } else {
	mytree->SetBranchAddress("pplus_TRUEID", &pplus_TRUEID);
	mytree->SetBranchAddress("pplus_TRUEP_X", &pplus_TRUEP_X);
	mytree->SetBranchAddress("pplus_TRUEP_Y", &pplus_TRUEP_Y);
	mytree->SetBranchAddress("pplus_TRUEP_Z", &pplus_TRUEP_Z);
	mytree->SetBranchAddress("pplus_TRUEP_E", &pplus_TRUEP_E);
	
	mytree->SetBranchAddress("lambda0_TRUEID", &lambda0_TRUEID);
	mytree->SetBranchAddress("lambda0_TRUEP_X", &lambda0_TRUEP_X);
	mytree->SetBranchAddress("lambda0_TRUEP_Y", &lambda0_TRUEP_Y);
	mytree->SetBranchAddress("lambda0_TRUEP_Z", &lambda0_TRUEP_Z);
	mytree->SetBranchAddress("lambda0_TRUEP_E", &lambda0_TRUEP_E);
	
	mytree->SetBranchAddress("lambda_b0_TRUEID", &lambdab_TRUEID);
	mytree->SetBranchAddress("lambda_b0_TRUEP_X", &lambdab_TRUEP_X);
	mytree->SetBranchAddress("lambda_b0_TRUEP_Y", &lambdab_TRUEP_Y);
	mytree->SetBranchAddress("lambda_b0_TRUEP_Z", &lambdab_TRUEP_Z);
	mytree->SetBranchAddress("lambda_b0_TRUEP_E", &lambdab_TRUEP_E);
      }
      
      mytree->SetBranchAddress("piminus_TRUEID", &piminus_TRUEID);
      mytree->SetBranchAddress("piminus_TRUEP_X", &piminus_TRUEP_X);
      mytree->SetBranchAddress("piminus_TRUEP_Y", &piminus_TRUEP_Y);
      mytree->SetBranchAddress("piminus_TRUEP_Z", &piminus_TRUEP_Z);
      mytree->SetBranchAddress("piminus_TRUEP_E", &piminus_TRUEP_E);

      mytree->SetBranchAddress("jpsi1s_TRUEID", &jpsi1s_TRUEID);
      mytree->SetBranchAddress("jpsi1s_TRUEP_X", &jpsi1s_TRUEP_X);
      mytree->SetBranchAddress("jpsi1s_TRUEP_Y", &jpsi1s_TRUEP_Y);
      mytree->SetBranchAddress("jpsi1s_TRUEP_Z", &jpsi1s_TRUEP_Z);
      mytree->SetBranchAddress("jpsi1s_TRUEP_E", &jpsi1s_TRUEP_E);
      
      mytree->SetBranchAddress("muplus_TRUEID", &muplus_TRUEID);
      mytree->SetBranchAddress("muplus_TRUEP_X", &muplus_TRUEP_X);
      mytree->SetBranchAddress("muplus_TRUEP_Y", &muplus_TRUEP_Y);
      mytree->SetBranchAddress("muplus_TRUEP_Z", &muplus_TRUEP_Z);
      mytree->SetBranchAddress("muplus_TRUEP_E", &muplus_TRUEP_E);
      
      mytree->SetBranchAddress("muminus_TRUEID", &muminus_TRUEID);
      mytree->SetBranchAddress("muminus_TRUEP_X", &muminus_TRUEP_X);
      mytree->SetBranchAddress("muminus_TRUEP_Y", &muminus_TRUEP_Y);
      mytree->SetBranchAddress("muminus_TRUEP_Z", &muminus_TRUEP_Z);
      mytree->SetBranchAddress("muminus_TRUEP_E", &muminus_TRUEP_E);

    }
    
    TTree* T;
    if (!debug) {
      T = new TTree(directorylist[idl],"test friend trees");
      T->Branch("phi1",&phi1,"phi1/F");
      T->Branch("phi2",&phi2,"phi2/F");
      T->Branch("costheta", &costheta, "costheta/F");
      T->Branch("costheta1",&costheta1,"costheta1/F");
      T->Branch("costheta2",&costheta2,"costheta2/F");

      T->Branch("true_phi1",&true_phi1,"true_phi1/F");
      T->Branch("true_phi2",&true_phi2,"true_phi2/F");
      T->Branch("true_costheta", &true_costheta, "true_costheta/F");
      T->Branch("true_costheta1",&true_costheta1,"true_costheta1/F");
      T->Branch("true_costheta2",&true_costheta2,"true_costheta2/F");
  
      // T->Branch("F1",&F[1],"F1/F");
      // T->Branch("F2",&F[2],"F2/F");
      // T->Branch("F3",&F[3],"F3/F");
      // T->Branch("F4",&F[4],"F4/F");
      // T->Branch("F5",&F[5],"F5/F");
      // T->Branch("F6",&F[6],"F6/F");
      // T->Branch("F7",&F[7],"F7/F");

      T->Branch("piminus_y",&piminus_y,"piminus_y/F");
	
      if (isJpsiKs) {
	T->Branch("ks0_phi",&ks0_phi,"ks0_phi/F");
	T->Branch("b0_eta",&lambdab0_eta,"b0_eta/F");
	T->Branch("ks0_eta",&lambda0_eta,"ks0_eta/F");
	T->Branch("b0_y",&lambdab0_y,"b0_y/F");
	T->Branch("ks0_y",&lambda0_y,"ks0_y/F");
	T->Branch("piplus_y",&pplus_y,"piplus_y/F");
      } else {
	T->Branch("lambda0_phi",&ks0_phi,"lambda0_phi/F");
	T->Branch("lambda_b0_eta",&lambdab0_eta,"lambda_b0_eta/F");
	T->Branch("lambda0_eta",&lambda0_eta,"lambda0_eta/F");
	T->Branch("lambda_b0_y",&lambdab0_y,"lambda_b0_y/F");
	T->Branch("lambda0_y",&lambda0_y,"lambda0_y/F");
	T->Branch("pplus_y",&pplus_y,"pplus_y/F");
      }
    }

    progressbar pbar(std::string("Calculating angles in " + filename + ":" + directorylist[idl] ));
    const Long64_t nentries=debug?4:mytree->GetEntries();

    TVector3 proddir(0.,0.,1.);
      
    for (Long64_t i=0; i<nentries; i++) {
      //  for (Int_t i=0; i<10; i++) {
      mytree->GetEvent(i);
      
      //      TLorentzVector lv_muplus, lv_muminus; 
      
      const TLorentzVector lv_lambdab( lambdab_PX,  lambdab_PY,  lambdab_PZ,  lambdab_PE);
      const TLorentzVector lv_lambda0( lambda0_PX,  lambda0_PY,  lambda0_PZ,  lambda0_PE);
      const TLorentzVector lv_jpsi1s( jpsi1s_PX,  jpsi1s_PY,  jpsi1s_PZ,  jpsi1s_PE);
      
      const TLorentzVector lv_muplus( muplus_PX,  muplus_PY,  muplus_PZ,  muplus_PE);
      const TLorentzVector lv_muminus(muminus_PX, muminus_PY, muminus_PZ, muminus_PE);
      const TLorentzVector lv_pplus(pplus_PX,   pplus_PY,   pplus_PZ,   pplus_PE);
      const TLorentzVector lv_piminus(  piminus_PX, piminus_PY, piminus_PZ, piminus_PE);
      

      const TLorentzVector lv_true_lambdab( lambdab_TRUEP_X,  lambdab_TRUEP_Y,  lambdab_TRUEP_Z,  lambdab_TRUEP_E);
      const TLorentzVector lv_true_lambda0( lambda0_TRUEP_X,  lambda0_TRUEP_Y,  lambda0_TRUEP_Z,  lambda0_TRUEP_E);
      const TLorentzVector lv_true_jpsi1s( jpsi1s_TRUEP_X,  jpsi1s_TRUEP_Y,  jpsi1s_TRUEP_Z,  jpsi1s_TRUEP_E);
      
      const TLorentzVector lv_true_muplus( muplus_TRUEP_X,  muplus_TRUEP_Y,  muplus_TRUEP_Z,  muplus_TRUEP_E);
      const TLorentzVector lv_true_muminus(muminus_TRUEP_X, muminus_TRUEP_Y, muminus_TRUEP_Z, muminus_TRUEP_E);
      const TLorentzVector lv_true_pplus(pplus_TRUEP_X,   pplus_TRUEP_Y,   pplus_TRUEP_Z,   pplus_TRUEP_E);
      const TLorentzVector lv_true_piminus(  piminus_TRUEP_X, piminus_TRUEP_Y, piminus_TRUEP_Z, piminus_TRUEP_E);

      //       if (muminus_ID>0) {
      // 	lv_muplus  = _lv_muplus;
      // 	lv_muminus = _lv_muminus;
      //       } else {
      // 	lv_muplus  = _lv_muminus;
      // 	lv_muminus = _lv_muplus;
      //       }
      
      //      calculateangles(lv_pplus, lv_piminus, lv_muplus, lv_muminus,
      //		      costheta, costheta1, costheta2, phi1, phi2);
      
	//      std::cout << "1: " << costheta << " " << costheta1 << " " << costheta2 << " " << phi1 << " " << phi2 << std::endl;
      if (!debug) {
	calculateangles(lv_lambdab, 
			  lv_lambda0, lv_pplus, lv_piminus,
			  lv_jpsi1s, lv_muplus, lv_muminus,
			  costheta, costheta1, costheta2, phi1, phi2, proddir);

	if (isMC) calculateangles(lv_true_lambdab, 
			lv_true_lambda0, lv_true_pplus, lv_true_piminus,
			lv_true_jpsi1s, lv_true_muplus, lv_true_muminus,
			true_costheta, true_costheta1, true_costheta2, true_phi1, true_phi2, proddir);
	


	// if (!isJpsiKs) {
	//   calculateangles(lv_lambdab, 
	// 		  lv_lambda0, lv_pplus, lv_piminus,
	// 		  lv_jpsi1s, lv_muplus, lv_muminus,
	// 		  costheta, costheta1, costheta2, phi1, phi2, proddir);
	// } else {
	//   //	  if (lv_pplus.P()>lv_piminus.P()) {
	//   calculateangles(lv_lambdab, 
	// 		  lv_lambda0, lv_pplus, lv_piminus,
	// 		  lv_jpsi1s, lv_muplus, lv_muminus,
	// 		  costheta, costheta1, costheta2, phi1, phi2, proddir);
	//   //	  } else {
	//   //	    calculateangles(lv_lambdab, 
	//   //			    lv_lambda0, lv_piminus, lv_pplus,
	//   //			    lv_jpsi1s, lv_muplus, lv_muminus,
	//   //			    costheta, costheta1, costheta2, phi1, phi2);
	//   //	  }
	// }
	
	
	phi1/=pi;
	phi2/=pi;

	true_phi1/=pi;
	true_phi2/=pi;
	
	ksangle(lv_pplus,lv_piminus,ks0_phi);
	ks0_phi/=pi;
	
	lambdab0_eta=pseudorapidity(lv_lambdab);
	lambda0_eta=pseudorapidity(lv_lambda0);
	
	lambdab0_y=rapidity(lv_lambdab);
	lambda0_y=rapidity(lv_lambda0);
	piminus_y=rapidity(lv_piminus);
	pplus_y=rapidity(lv_pplus);
	
	// F[1]=costheta;
	// F[2]=costheta1;
	// F[3]=costheta*costheta1;
	// F[4]=0.5*(3.*costheta2*costheta2-1.);
	// F[5]=0.5*(3.*costheta2*costheta2-1.)*costheta;
	// F[6]=0.5*(3.*costheta2*costheta2-1.)*costheta1;
	// F[7]=0.5*(3.*costheta2*costheta2-1.)*costheta*costheta1;
	
      } else {
	
	TVector3 proddir(0.,0.,1.);
	std::cout << "lv_lambdab   = " << lambdab_ID << " " << lvcout(lv_lambdab) << std::endl;
	std::cout << "lv_lambdab   = " << lambdab_ID << " " << lvcout(lv_lambdab.Vect()) << std::endl;
	calculateangles(lv_lambdab, 
			lv_lambda0, lv_pplus, lv_piminus,
			lv_jpsi1s, lv_muplus, lv_muminus,
			costheta, costheta1, costheta2, phi1, phi2, proddir);
	std::cout << "proddir      = " << lvcout(proddir) << std::endl;
	std::cout << "costheta     = " << costheta << std::endl;
	std::cout << "costheta1    = " << costheta1 << std::endl;    
	std::cout << "costheta2    = " << costheta2 << std::endl;
	
	proddir = TVector3(0.,0.,-1.);
	  
	calculateangles(lv_lambdab, 
			lv_lambda0, lv_pplus, lv_piminus,
			lv_jpsi1s, lv_muplus, lv_muminus,
			costheta, costheta1, costheta2, phi1, phi2, proddir);
	std::cout << "proddir      = " << lvcout(proddir) << std::endl;
	std::cout << "costheta     = " << costheta << std::endl;
	std::cout << "costheta1    = " << costheta1 << std::endl;    
	std::cout << "costheta2    = " << costheta2 << std::endl << std::endl;
	
	for (int k=0;k<-1;k++) {
	  proddir = TVector3(rnd.Rndm(),rnd.Rndm(),rnd.Rndm()).Unit();
	  calculateangles(lv_lambdab, 
			  lv_lambda0, lv_pplus, lv_piminus,
			  lv_jpsi1s, lv_muplus, lv_muminus,
			  costheta, costheta1, costheta2, phi1, phi2, proddir);
	  phi1/=pi;
	  phi2/=pi;
	  
	  std::cout << "directory    = " << directorylist[idl] << std::endl;
	  std::cout << "lv_lambdab   = " << lambdab_ID << " " << lvcout(lv_lambdab) << std::endl;
	  std::cout << "lv_lambda0   = " << lambda0_ID << " " << lvcout(lv_lambda0) << std::endl;
	  std::cout << "lv_jpsi      = " << jpsi1s_ID << " " << lvcout(lv_jpsi1s) << std::endl;
	  std::cout << "lv_pplus     = " << pplus_ID << " " << lvcout(lv_pplus) << std::endl;
	  std::cout << "lv_piminus   = " << piminus_ID << " " << lvcout(lv_piminus) << std::endl;
	  std::cout << "lv_muonp     = " << muplus_ID << " " << lvcout(lv_muplus) << std::endl;
	  std::cout << "lv_muonm     = " << muminus_ID << " " << lvcout(lv_muminus) << std::endl;
	  std::cout << "proddir      = " << lvcout(proddir) << std::endl;
	  std::cout << "costheta     = " << costheta << std::endl;
	  std::cout << "costheta1    = " << costheta1 << std::endl;    
	  
	  calculateangles(lv_lambdab, 
			  lv_lambda0, lv_piminus, lv_pplus,
			  lv_jpsi1s, lv_muplus, lv_muminus,
			  costheta, costheta1, costheta2, phi1, phi2, proddir);
	  phi1/=pi;
	  phi2/=pi;
	  std::cout << "costheta1(v2)= " << costheta1 << std::endl;    
	  
	  std::cout << "costheta2    = " << costheta2 << std::endl;
	  std::cout << "phi1         = " << phi1 << std::endl;
	  std::cout << "phi2         = " << phi2 << std::endl;
	  std::cout << "B mass       = " << lv_lambdab.M() << std::endl;
	  std::cout << "B mass (4)   = " << (lv_muplus+lv_muminus+lv_pplus+lv_piminus).M() << std::endl;
	  std::cout << "B mass (2)   = " << (lv_jpsi1s+lv_lambda0).M() << std::endl;
	  std::cout << "Jpsi mass    = " << lv_jpsi1s.M() << std::endl;
	  std::cout << "lambda0 mass = " << lv_lambda0.M() << std::endl;
	  std::cout << "muplus mass  = " << lv_muplus.M() << std::endl;
	  std::cout << "muminus mass = " << lv_muminus.M() << std::endl;
	  std::cout << "pplus mass   = " << lv_pplus.M() << std::endl;
	  std::cout << "piminus mass = " << lv_piminus.M() << std::endl;
	  std::cout << std::endl;
	  
	}
	
      }
      
      //      calculateangles2(lv_pplus, lv_piminus, lv_muplus, lv_muminus,
      //		      costheta, costheta1, costheta2, phi1, phi2);
      //      std::cout << "2: " << costheta << " " << costheta1 << " " << costheta2 << " " << phi1 << " " << phi2 << std::endl;
      
      if (!debug) pbar.print(100.*i/(nentries-1));
      if (!debug) T->Fill();
    }
    pbar.finish();
    //    if (!debug) T->Write()AutoSave();
    if (!debug) T->AutoSave();
  }
  if (!debug) f->Close();
  myfile.Close();
}




// void calculateangles(const TLorentzVector& lv_proton, const TLorentzVector& lv_pion, const TLorentzVector& lv_muonp, const TLorentzVector& lv_muonm,
// 		     Float_t& costheta, Float_t& costheta1, Float_t& costheta2, Float_t& phi1, Float_t& phi2, const TVector3& proddir)
// {
//   costheta=0.;
//   costheta1=0.;
//   costheta2=0.;
//   phi1=0.;
//   phi2=0.;

//   const TLorentzVector lv_Lambdab(lv_proton+lv_pion+lv_muonp+lv_muonm);
//   const TLorentzVector lv_Lambda0(lv_proton+lv_pion);
//   const TLorentzVector lv_Jpsi(lv_muonp+lv_muonm);

//   const TVector3 xyz_Lambdab = lv_Lambdab.Vect().Unit();
//   const TVector3 n((proddir.Cross(xyz_Lambdab)).Unit());
  
//   //first boost everything in Lb frame

//   const TVector3  boostToLambdab = -lv_Lambdab.BoostVector();
//   const TLorentzVector lv_Lambdab_Lambdabframe = myboost(lv_Lambdab, boostToLambdab);
//   const TLorentzVector lv_Lambda0_Lambdabframe = myboost(lv_Lambda0, boostToLambdab);
//   const TLorentzVector lv_proton_Lambdabframe  = myboost(lv_proton, boostToLambdab);
//   const TLorentzVector lv_pion_Lambdabframe    = myboost(lv_pion, boostToLambdab);
//   const TLorentzVector lv_Jpsi_Lambdabframe    = myboost(lv_Jpsi, boostToLambdab);
//   const TLorentzVector lv_muonp_Lambdabframe   = myboost(lv_muonp, boostToLambdab);
//   const TLorentzVector lv_muonm_Lambdabframe   = myboost(lv_muonm, boostToLambdab);
  
//   {
//     //theta
//     const TVector3 xyz_Lambda0_Lambdabframe = TVector3(lv_Lambda0_Lambdabframe.Vect()).Unit();
//     const TVector3 xyz_Jpsi_Lambdabframe = TVector3(lv_Jpsi_Lambdabframe.Vect()).Unit();
//     costheta = n.Dot(xyz_Lambda0_Lambdabframe);
//   }


//   {
//     //theta 1, phi 1 (Lambda0 stuff)
//     const TVector3 z1=lv_Lambda0_Lambdabframe.Vect().Unit();
//     const TVector3 y1=(n.Cross(z1)).Unit();
//     const TVector3 x1=(y1.Cross(z1)).Unit();
    
//     const TVector3 boostToLambda0 = -lv_Lambda0_Lambdabframe.BoostVector();
//     const TLorentzVector lv_proton_Lambda0frame = myboost(lv_proton_Lambdabframe, boostToLambda0);
//     const TVector3 xyz_proton_Lambda0frame = lv_proton_Lambda0frame.Vect().Unit();
//     costheta1 = z1.Dot(xyz_proton_Lambda0frame);

//     const TVector3 xyz_proton_Lambda0frame_projete = (z1.Cross(xyz_proton_Lambda0frame.Cross(z1))).Unit();
//     const double cosphi1 = x1.Dot(xyz_proton_Lambda0frame_projete);
//     const double sinphi1 = y1.Dot(xyz_proton_Lambda0frame_projete);
//     phi1 = sinphi1>0.0 ? acos(cosphi1) : -acos(cosphi1);
//   }

//   {
//     //theta 2, phi 2 (Jpsi stuff)
//     const TVector3 z2=lv_Jpsi_Lambdabframe.Vect().Unit();
//     const TVector3 y2=(n.Cross(z2)).Unit();
//     const TVector3 x2=(y2.Cross(z2)).Unit();
    
//     const TVector3 boostToJpsi= -lv_Jpsi_Lambdabframe.BoostVector();
//     const TLorentzVector lv_muonp_Jpsiframe = myboost(lv_muonp_Lambdabframe, boostToJpsi);
//     const TVector3 xyz_muonp_Jpsiframe = lv_muonp_Jpsiframe.Vect().Unit();
//     costheta2 = z2.Dot(xyz_muonp_Jpsiframe);
    
//     const TVector3 xyz_muonp_Jpsiframe_projete = (z2.Cross(xyz_muonp_Jpsiframe.Cross(z2))).Unit();
//     const double cosphi2 = x2.Dot(xyz_muonp_Jpsiframe_projete);
//     const double sinphi2 = y2.Dot(xyz_muonp_Jpsiframe_projete);
//     phi2 = sinphi2>0.0 ? acos(cosphi2) : -acos(cosphi2);
//   }
// }



// void calculateangles(const TLorentzVector& lv_proton, const TLorentzVector& lv_pion, const TLorentzVector& lv_muonp, const TLorentzVector& lv_muonm,
// 		     Float_t& costheta, Float_t& costheta1, Float_t& costheta2, Float_t& phi1, Float_t& phi2)
// {
//   const TVector3 proddir(0.,0.,1.);
//   calculateangles(lv_proton, lv_pion, lv_muonp, lv_muonm,
// 		  costheta, costheta1, costheta2, phi1, phi2, proddir);
// }


