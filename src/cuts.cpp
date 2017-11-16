#include "cuts.h"

double getLambda0mass(bool doBd2JpsiKS0)
{
  if (doBd2JpsiKS0) {
    return 497.614;
  } else {
    return 1115.683;
  }
}

void varnames(bool doBd2JpsiKS0, TString& Lambdabname, TString& Lambda0name, TString& Jpsiname)
{
  Lambdabname="lambda_b0";
  Lambda0name="lambda0";
  Jpsiname="jpsi1s";
  
  //B0 to J/psi KS0
  if (doBd2JpsiKS0) {
    Lambdabname="b0";
    Lambda0name="ks0";
  }
}

// TCut pidcut(bool doBd2JpsiKS0)
// {
//   //  TCut piminusid_ll("piminus_PIDK>-40.&&piminus_PIDp>-40.");
//   //  TCut piplusid_ll("piplus_PIDK>-40.&&piplus_PIDp>-40.");
//   TCut pplusid_ll("pplus_PIDp>10.");

//   //  TCut piminusid_dd("piminus_PIDK>-30.&&piminus_PIDp>-30.");
//   //  TCut piplusid_dd("piplus_PIDK>-30.&&piplus_PIDp>-30.");
//   TCut pplusid_dd("pplus_PIDK>-20.&&pplus_PIDp>0.");

//   if (doBd2JpsiKS0) {
//     return TCut("(1==1)");//Lambda0_DD or (Lambda0_LL + piminusid_ll + piplusid_ll); 
//   } else {
//     return (Lambda0_DD + pplusid_dd) or (Lambda0_LL  + pplusid_ll);
//   }
// }

TCut triggercut(bool doBd2JpsiKS0)
{
  TString Lambdabname,Lambda0name,Jpsiname;
  varnames(doBd2JpsiKS0, Lambdabname, Lambda0name, Jpsiname);
  
  //   const TCut L0cut("L0MuonDecision==1||L0DiMuonDecision==1");
  //   const TCut HLT1cut("Hlt1DiMuonHighMassDecision==1");
  //   const TCut HLT2cut("Hlt2DiMuonJPsiDecision==1.");
  //   //  TCut TriggerCut("L0DiMuonDecision==1&&Hlt1DiMuonHighMassDecision==1&&Hlt2DiMuonJPsiDecision==1.&&jpsi1sL0Global_TOS==1&&jpsi1sHlt1Global_TOS==1&&jpsi1sHlt2Global_TOS==1" );
  
  std::vector<TString> L0decisions;
  L0decisions.push_back("L0MuonDecision");
  L0decisions.push_back("L0DiMuonDecision");

  std::vector<TString> HLT1decisions;
  HLT1decisions.push_back("Hlt1TrackMuonDecision");
  HLT1decisions.push_back("Hlt1TrackAllL0Decision");
  HLT1decisions.push_back("Hlt1DiMuonHighMassDecision");

//   HLT1decisions.push_back("Hlt1DiMuonLowMassDecision");
//   HLT1decisions.push_back("Hlt1SingleMuonNoIPDecision");


  std::vector<TString> HLT2decisions;
  HLT2decisions.push_back("Hlt2DiMuonDetachedJPsiDecision");

  //  HLT2decisions.push_back("Hlt2DiMuonJPsiDecision");

  //  HLT2decisions.push_back("Hlt2TopoMu2BodyBBDTDecision");
  //   //  HLT2decisions.push_back("Hlt2TopoMu3BodyBBDTDecision");//
  //   //  HLT2decisions.push_back("Hlt2TopoMu4BodyBBDTDecision");//
  //  HLT2decisions.push_back("Hlt2Topo2BodyBBDTDecision");
  //   //  HLT2decisions.push_back("Hlt2Topo3BodyBBDTDecision");//
  //   //  HLT2decisions.push_back("Hlt2Topo4BodyBBDTDecision");//
  //  HLT2decisions.push_back("Hlt2SingleMuonDecision");
  //   HLT2decisions.push_back("Hlt2DiMuonDetachedDecision");
  //   HLT2decisions.push_back("Hlt2DiMuonDetachedHeavyDecision");
  //   HLT2decisions.push_back("Hlt2DiMuonJPsiHighPTDecision");
  
  //const TString TOSon(Lambdabname);
  const TString& TOSon = Jpsiname;

  TCut L0TosSelection;
  for (unsigned int i=0;i<L0decisions.size();i++) L0TosSelection= L0TosSelection or TCut(TOSon + L0decisions[i] + "_TOS==1.");
  TCut HLT1TosSelection;
  for (unsigned int i=0;i<HLT1decisions.size();i++) HLT1TosSelection= HLT1TosSelection or TCut(TOSon + HLT1decisions[i] + "_TOS==1.");
  TCut HLT2TosSelection;
  for (unsigned int i=0;i<HLT2decisions.size();i++) HLT2TosSelection= HLT2TosSelection or TCut(TOSon + HLT2decisions[i] + "_TOS==1.");


  return L0TosSelection + HLT1TosSelection + HLT2TosSelection;
  //  return TCut(L0TosSelection) + TCut(HLT1TosSelection) + TCut(HLT2TosSelection);
  //  return L0cut + HLT1cut + HLT2cut;
  //  return L0cut + HLT2cut;
}


TCut masscut(bool doBd2JpsiKS0, const RooRealVar& mass)
{
  TString Lambdabname, Lambda0name, Jpsiname;
  varnames(doBd2JpsiKS0,Lambdabname, Lambda0name, Jpsiname);
  return TCut(Lambdabname + "_M>" + str(mass.getMin()) +  "&&" + Lambdabname + "_M<" + str(mass.getMax()));
}


void simplecut2(bool doBd2JpsiKS0, TCut& LL_selection, TCut& DD_selection)
{

  TString Lambdabname, Lambda0name, Jpsiname;
  varnames(doBd2JpsiKS0,Lambdabname, Lambda0name, Jpsiname);

  //Jpsi -> mu mu
  const TCut muplus_PT("muplus_PT>500.");
  const TCut muminus_PT("muminus_PT>500.");
  const TCut Jpsi_MM("abs(" + Jpsiname + "_MM-" + str(Jpsipdgmass) + ")<55.");

  //Lambda0 -> p pi
  const TCut pplus_PT("pplus_PT>500.");
  const TCut piminus_PT("piminus_PT>100.");
  const TCut pplus_P("pplus_P>2000.");
  const TCut piminus_P("piminus_P>2000.");

  //  const TCut Lambda0_IPCHI2_LL(Lambda0name + "_IPCHI2_OWNPV<8.");
  //  const TCut Lambda0_IPCHI2_DD(Lambda0name + "_IPCHI2_OWNPV<7.");

  TCut Lambda0_M("abs(" + Lambda0name + "_M-" + str(getLambda0mass(doBd2JpsiKS0)) + ")<6.");

  const TCut Lambda0_PT(Lambda0name + "_PT>1000.");

  //Lambda_b0 -> ...

  
  const TCut Lambdab_IPCHI2(Lambdabname + "_IPCHI2_OWNPV<20.");
  
  const TCut Lambdab_TAU(Lambdabname + "_TAU>0.000250");

  const TCut commonsel = Lambdab_TAU + Lambdab_IPCHI2 + Lambda0_M + Lambda0_PT + piminus_P + pplus_P + piminus_PT + pplus_PT + Jpsi_MM + muminus_PT + muplus_PT;
  LL_selection = commonsel;
  DD_selection = commonsel;

}

void simplecut(bool doBd2JpsiKS0, TCut& LL_selection, TCut& DD_selection)
{

  TString Lambdabname, Lambda0name, Jpsiname;
  varnames(doBd2JpsiKS0,Lambdabname, Lambda0name, Jpsiname);

  //  const TCut Lambdab_BDIRA_LL(Lambdabname + "_DIRA_OWNPV>0.999968");
  //  const TCut Lambdab_BDIRA_DD(Lambdabname + "_DIRA_OWNPV>0.999978");

  const TCut Lambdab_BDIRA_LL(Lambdabname + "_DIRA_OWNPV>0.9999");
  const TCut Lambdab_BDIRA_DD(Lambdabname + "_DIRA_OWNPV>0.9999");

  const TCut Lambdab_IPCHI2_LL(Lambdabname + "_IPCHI2_OWNPV<8.");
  const TCut Lambdab_IPCHI2_DD(Lambdabname + "_IPCHI2_OWNPV<7.");

  //  const TCut Lambdab_FDCHI2_LL(Lambdabname + "_FDCHI2_OWNPV>130.");
  //  const TCut Lambdab_FDCHI2_DD(Lambdabname + "_FDCHI2_OWNPV>180.");

  const TCut Lambdab_FDCHI2_LL(Lambdabname + "_FDCHI2_OWNPV>130.");
  const TCut Lambdab_FDCHI2_DD(Lambdabname + "_FDCHI2_OWNPV>180.");

  const TCut Lambdab_VCHI2DOF_LL(Lambdabname + "_ENDVERTEX_CHI2/" + Lambdabname + "_ENDVERTEX_NDOF<5.");
  const TCut Lambdab_VCHI2DOF_DD(Lambdabname + "_ENDVERTEX_CHI2/" + Lambdabname + "_ENDVERTEX_NDOF<4.");

  //  const TCut Lambdab_TAU_LL(Lambdabname + "_TAU>0.0002");
  //  const TCut Lambdab_TAU_DD(Lambdabname + "_TAU>0.0002");
  const TCut Lambdab_TAU_LL(Lambdabname + "_TAU>0.0002");
  const TCut Lambdab_TAU_DD(Lambdabname + "_TAU>0.0002");


  //  TCut Lambda0_M_LL("abs(" + Lambda0name + "_M-" + str(getLambda0mass(doBd2JpsiKS0)) + ")<20.");
  //  TCut Lambda0_M_DD("abs(" + Lambda0name + "_M-" + str(getLambda0mass(doBd2JpsiKS0)) + ")<15.");
  TCut Lambda0_M_LL("abs(" + Lambda0name + "_M-" + str(getLambda0mass(doBd2JpsiKS0)) + ")<9.");
  TCut Lambda0_M_DD("abs(" + Lambda0name + "_M-" + str(getLambda0mass(doBd2JpsiKS0)) + ")<6.");
  if (doBd2JpsiKS0) {
    Lambda0_M_LL = TString("abs(" + Lambda0name + "_M-" + str(getLambda0mass(doBd2JpsiKS0)) + ")<21.");
    Lambda0_M_DD = TString("abs(" + Lambda0name + "_M-" + str(getLambda0mass(doBd2JpsiKS0)) + ")<12.");
  }


  const TCut Lambda0_PT_LL(Lambda0name + "_PT>800.");
  const TCut Lambda0_PT_DD(Lambda0name + "_PT>1000.");

  const TCut Jpsi_VCHI2DOF_LL(Jpsiname + "_ENDVERTEX_CHI2/" + Jpsiname + "_ENDVERTEX_NDOF<9.");
  const TCut Jpsi_VCHI2DOF_DD(Jpsiname + "_ENDVERTEX_CHI2/" + Jpsiname + "_ENDVERTEX_NDOF<8.");

  const TCut Jpsi_MM_LL("abs(" + Jpsiname + "_MM-" + str(Jpsipdgmass) + ")<40.");
  const TCut Jpsi_MM_DD("abs(" + Jpsiname + "_MM-" + str(Jpsipdgmass) + ")<40.");

  LL_selection = Lambda0_LL + Lambdab_BDIRA_LL +  Lambdab_TAU_LL + Lambdab_IPCHI2_LL + Lambdab_FDCHI2_LL + Lambdab_VCHI2DOF_LL  + Lambda0_M_LL  +  Lambda0_PT_LL + Jpsi_MM_LL + Jpsi_VCHI2DOF_LL;
  DD_selection = Lambda0_DD + Lambdab_BDIRA_DD +  Lambdab_TAU_DD + Lambdab_IPCHI2_DD + Lambdab_FDCHI2_DD + Lambdab_VCHI2DOF_DD  + Lambda0_M_DD  +  Lambda0_PT_DD + Jpsi_MM_DD + Jpsi_VCHI2DOF_DD;
}

void tmvarectcut(bool doBd2JpsiKS0, TCut& LL_selection, TCut& DD_selection)
{

  TString Lambdabname, Lambda0name, Jpsiname;
  varnames(doBd2JpsiKS0,Lambdabname, Lambda0name, Jpsiname);
  
  TCut precut(Lambdabname + "_DIRA_OWNPV>0.99&&" + Lambdabname + "_TAU>0.&&" + Lambdabname + "_IPCHI2_OWNPV<30.");
  

// --- CutsGA                   : -------------------------------------------------------------------------------------
// --- CutsGA                   : Cut values for requested signal efficiency: 0.96
// --- CutsGA                   : Corresponding background efficiency       : 0.214286
// --- CutsGA                   : Transformation applied to input variables : None
// --- CutsGA                   : -------------------------------------------------------------------------------------
// --- CutsGA                   : Cut[ 0]:   -7.01373 <             log10(1.0000001-lambda_b0_DIRA_OWNPV) <=   -3.38773
// --- CutsGA                   : Cut[ 1]:  -0.154572 <                            lambda_b0_IPCHI2_OWNPV <=    26.0728
// --- CutsGA                   : Cut[ 2]:    1.61694 <                     log10(lambda_b0_FDCHI2_OWNPV) <=    5.80261
// --- CutsGA                   : Cut[ 3]: -0.0710959 < lambda_b0_ENDVERTEX_CHI2/lambda_b0_ENDVERTEX_NDOF <=    8.15207
// --- CutsGA                   : Cut[ 4]: 9.05971e-05 <                                     lambda_b0_TAU <=  0.0125028
// --- CutsGA                   : Cut[ 5]:  -0.134712 <                        abs(lambda0_M-1115.683000) <=    5.97692
// --- CutsGA                   : Cut[ 6]:  -0.140075 <                        abs(jpsi1s_MM-3096.916000) <=    79.7509
// --- CutsGA                   : Cut[ 7]: -0.0141781 <       jpsi1s_ENDVERTEX_CHI2/jpsi1s_ENDVERTEX_NDOF <=    16.1238
// --- CutsGA                   : Cut[ 8]:   -1.00679 <                  piminus_ProbNNk-piminus_ProbNNpi <=   0.626011
// --- CutsGA                   : Cut[ 9]:   -1.00565 <                  piminus_ProbNNp-piminus_ProbNNpi <=   0.335731
// --- CutsGA                   : Cut[10]:   -0.84198 <                       pplus_ProbNNk-pplus_ProbNNp <=   0.532617
// --- CutsGA                   : Cut[11]:  -0.769958 <                      pplus_ProbNNpi-pplus_ProbNNp <=   0.618048
// --- CutsGA                   : -------------------------------------------------------------------------------------


  const TCut Lambdab_BDIRA_LL("log10(1.0000001-" + Lambdabname + "_DIRA_OWNPV)<-3.38773");
  const TCut Lambdab_IPCHI2_LL(Lambdabname + "_IPCHI2_OWNPV<26.0");
  const TCut Lambdab_FDCHI2_LL("log10(" + Lambdabname + "_FDCHI2_OWNPV)>1.61694");
  const TCut Lambdab_VCHI2DOF_LL(Lambdabname + "_ENDVERTEX_CHI2/" + Lambdabname + "_ENDVERTEX_NDOF<8.2");
  const TCut Lambdab_TAU_LL(Lambdabname + "_TAU>0.0002");
  TCut Lambda0_M_LL("abs(" + Lambda0name + "_M-" + str(getLambda0mass(doBd2JpsiKS0)) + ")<6.0");
  //  const TCut Lambda0_PT_LL(Lambda0name + "_PT>854.");
  //  const TCut Jpsi_MM_LL("abs(" + Jpsiname + "_MM-" + str(Jpsipdgmass) + ")<73.");
  const TCut Jpsi_VCHI2DOF_LL(Jpsiname + "_ENDVERTEX_CHI2/" + Jpsiname + "_ENDVERTEX_NDOF<16.1");

  const TCut pplus_kp_LL("pplus_ProbNNk-pplus_ProbNNp<0.5");
  const TCut pplus_pip_LL("pplus_ProbNNpi-pplus_ProbNNp<0.2");

// --- CutsGA                   : -------------------------------------------------------------------------------------
// --- CutsGA                   : Cut values for requested signal efficiency: 0.91
// --- CutsGA                   : Corresponding background efficiency       : 0.185348
// --- CutsGA                   : Transformation applied to input variables : None
// --- CutsGA                   : -------------------------------------------------------------------------------------
// --- CutsGA                   : Cut[ 0]:   -7.04315 <             log10(1.0000001-lambda_b0_DIRA_OWNPV) <=    -3.2713
// --- CutsGA                   : Cut[ 1]:  -0.237002 <                            lambda_b0_IPCHI2_OWNPV <=    11.9528
// --- CutsGA                   : Cut[ 2]:    1.57368 <                     log10(lambda_b0_FDCHI2_OWNPV) <=    6.36213
// --- CutsGA                   : Cut[ 3]:  -0.038209 < lambda_b0_ENDVERTEX_CHI2/lambda_b0_ENDVERTEX_NDOF <=    6.90436
// --- CutsGA                   : Cut[ 4]: 0.000151011 <                                     lambda_b0_TAU <=  0.0118273
// --- CutsGA                   : Cut[ 5]:   -0.11156 <                        abs(lambda0_M-1115.683000) <=    11.9438
// --- CutsGA                   : Cut[ 6]:  -0.728259 <                        abs(jpsi1s_MM-3096.916000) <=    63.4592
// --- CutsGA                   : Cut[ 7]:  -0.069051 <       jpsi1s_ENDVERTEX_CHI2/jpsi1s_ENDVERTEX_NDOF <=     15.841
// --- CutsGA                   : Cut[ 8]:  -0.987364 <                  piminus_ProbNNk-piminus_ProbNNpi <=   0.182473
// --- CutsGA                   : Cut[ 9]:  -0.979794 <                  piminus_ProbNNp-piminus_ProbNNpi <=   0.666177
// --- CutsGA                   : Cut[10]:  -0.756661 <                       pplus_ProbNNk-pplus_ProbNNp <= -0.0554169
// --- CutsGA                   : Cut[11]:  -0.703686 <                      pplus_ProbNNpi-pplus_ProbNNp <=   0.729207
// --- CutsGA                   : -------------------------------------------------------------------------------------


  const TCut Lambdab_BDIRA_DD("log10(1.0000001-" + Lambdabname + "_DIRA_OWNPV)<3.2713");
  const TCut Lambdab_IPCHI2_DD(Lambdabname + "_IPCHI2_OWNPV<12.58");
  const TCut Lambdab_FDCHI2_DD("log10(" + Lambdabname + "_FDCHI2_OWNPV)>1.57368");
  const TCut Lambdab_VCHI2DOF_DD(Lambdabname + "_ENDVERTEX_CHI2/" + Lambdabname + "_ENDVERTEX_NDOF<6.9");
  const TCut Lambdab_TAU_DD(Lambdabname + "_TAU>0.0002");
  TCut Lambda0_M_DD("abs(" + Lambda0name + "_M-" + str(getLambda0mass(doBd2JpsiKS0)) + ")<11.94");
//  const TCut Lambda0_PT_DD(Lambda0name + "_PT>686.");
//  const TCut Jpsi_MM_DD("abs(" + Jpsiname + "_MM-" + str(Jpsipdgmass) + ")<53.");
  const TCut Jpsi_VCHI2DOF_DD(Jpsiname + "_ENDVERTEX_CHI2/" + Jpsiname + "_ENDVERTEX_NDOF<15.8");


  const TCut pplus_kp_DD("pplus_ProbNNk-pplus_ProbNNp<-0.06");
  const TCut pplus_pip_DD("pplus_ProbNNpi-pplus_ProbNNp<0.2");


  if (doBd2JpsiKS0) Lambda0_M_LL = TString("abs(" + Lambda0name + "_M-" + str(getLambda0mass(doBd2JpsiKS0)) + ")<21.");
  if (doBd2JpsiKS0) Lambda0_M_DD = TString("abs(" + Lambda0name + "_M-" + str(getLambda0mass(doBd2JpsiKS0)) + ")<12.");

  LL_selection = precut + Lambda0_LL + Lambdab_BDIRA_LL +  Lambdab_TAU_LL + Lambdab_IPCHI2_LL + Lambdab_FDCHI2_LL + Lambdab_VCHI2DOF_LL  + Lambda0_M_LL  + Jpsi_VCHI2DOF_LL;
  DD_selection = precut + Lambda0_DD + Lambdab_BDIRA_DD +  Lambdab_TAU_DD + Lambdab_IPCHI2_DD + Lambdab_FDCHI2_DD + Lambdab_VCHI2DOF_DD  + Lambda0_M_DD  + Jpsi_VCHI2DOF_DD;

  if (!doBd2JpsiKS0) {
    DD_selection+=pplus_kp_DD;
    DD_selection+=pplus_pip_DD;
    LL_selection+=pplus_kp_LL;
    LL_selection+=pplus_pip_LL;
  }

}




void bdtsel(bool doBd2JpsiKS0, TCut& LL_selection, TCut& DD_selection)
{

  //TMVA-bkgrej-ll-20120719-1115.root TMVA-bkgrej-dd-20120719-1115.root
  //dd
  //   // --- ==================================================================================================
  //   // --- Classifier   (  #signal, #backgr.)  Optimal-cut  S/sqrt(S+B)      NSig      NBkg   EffSig   EffBkg
  //   // --- --------------------------------------------------------------------------------------------------
  //   // ---        BDT:  (     4900,     3570)      -0.0533      64.9666  4655.219  479.2944     0.95   0.1343
  //   // ---     CutsGA:  (     4900,     3570)       0.9150      62.0119  4435.513   680.574   0.9052   0.1906
  //   // --- --------------------------------------------------------------------------------------------------
  // --- ==================================================================================================
  // --- Classifier   (  #signal, #backgr.)  Optimal-cut  S/sqrt(S+B)      NSig      NBkg   EffSig   EffBkg
  // --- --------------------------------------------------------------------------------------------------
  // ---        BDT:  (     5751,     5169)      -0.0851      66.3104   5238.13   1001.93   0.9108   0.1938
  // --- --------------------------------------------------------------------------------------------------
  

  //ll
  //   // --- ==================================================================================================
  //   // --- Classifier   (  #signal, #backgr.)  Optimal-cut  S/sqrt(S+B)      NSig      NBkg   EffSig   EffBkg
  //   // --- --------------------------------------------------------------------------------------------------
  //   // ---        BDT:  (     1640,      805)       0.0084      38.4983  1560.775  82.83217   0.9517   0.1029
  //   // ---     CutsGA:  (     1640,      805)       0.9650      37.4665  1564.308  178.9336   0.9538   0.2223
  //   // --- --------------------------------------------------------------------------------------------------
  // --- ==================================================================================================
  // --- Classifier   (  #signal, #backgr.)  Optimal-cut  S/sqrt(S+B)      NSig      NBkg   EffSig   EffBkg
  // --- --------------------------------------------------------------------------------------------------
  // ---        BDT:  (     1876,     1036)       0.0055      40.4286  1729.892  100.9894   0.9221  0.09748
  // --- --------------------------------------------------------------------------------------------------
  

  TString Lambdabname, Lambda0name, Jpsiname;
  varnames(doBd2JpsiKS0,Lambdabname, Lambda0name, Jpsiname);

  TCut BDTPRECUT=TCut(Lambdabname + "_DIRA_OWNPV>0.99&&" + Lambdabname + "_TAU>0.&&" + Lambdabname + "_IPCHI2_OWNPV<30.");

  //  if (doBd2JpsiKS0) BDTPRECUT+=TCut("(piplus_PT/piminus_PT>2.||piminus_PT/piplus_PT>2.)");
  //  if (doBd2JpsiKS0) BDTPRECUT+=TCut("(piplus_P/piminus_P>2.||piminus_P/piplus_P>2.)");

  DD_selection = Lambda0_DD + BDTPRECUT + TCut("bdt>-0.0851");// + TCut(Lambda0name + "_PT>2000.");
  LL_selection = Lambda0_LL + BDTPRECUT + TCut("bdt>0.0055");// + TCut(Lambda0name + "_PT>2000.");


  //  DD_selection+="abs(Lambda0_M_pmisid-497.614)>8.";
  //  LL_selection+="abs(Lambda0_M_pmisid-497.614)>14.";

  //  DD_selection = Lambda0_DD + BDTPRECUT + TCut("bdt>-0.5");// + TCut(Lambda0name + "_PT>2000.");
  //  LL_selection = Lambda0_LL + BDTPRECUT + TCut("bdt>-0.5");// + TCut(Lambda0name + "_PT>2000.");
}

TCut truemcsel(bool doBd2JpsiKS0)
{  
  TString Lambdabname, Lambda0name, Jpsiname;
  varnames(doBd2JpsiKS0,Lambdabname, Lambda0name, Jpsiname);
  
  TString Lambdab_PDGID("5122");
  if (doBd2JpsiKS0) Lambdab_PDGID = "511";
  return TCut("abs(" + Lambdabname + "_TRUEID)==" + Lambdab_PDGID);
}
