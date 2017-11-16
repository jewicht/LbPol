

Double_t myfunction0(Double_t *x, Double_t *par)
{
  Double_t costheta  = x[0];
  Double_t costheta1 = x[1];
  Double_t costheta2 = x[2];
  return pdfval(costheta,costheta1,costheta2,dcoeff3);
}



void readacceptance3()
{
  std::ifstream filestr;
  filestr.open("coeff3.txt");
  
  //build TF3
  double dcoeff3[ordermax3i][ordermax3j][ordermax3k];
  //  RooRealVar* coeff3[ordermax3i][ordermax3j][ordermax3k];
  //  createRooLegendre3(costheta, costheta1, costheta2, "accpdf", accpdf, coeff3);
  //  createRooLegendre3(costheta, costheta1, costheta2, "accpdf", accpdf, coeff3);
  //  createRooLbtoJpsiL0wAccPDF3(costheta, costheta1, costheta2, P_b, alpha_b, alpha_lambda, r_0, r_1, "wxacc", fullpdf, coeff3);

  //  accpdf->forceNumInt();
  for (int i=0; i<ordermax3i;i++) {
    for (int j=0; j<ordermax3j;j++) {
      for (int k=0; k<ordermax3k;k++) {
	double value, error;
	bool isconstant;
	int ii,jj,kk;
        filestr >>  ii >> jj >> kk >> value >> error >> isconstant;// >> std::endl;
	std::cout << ii << " " << jj << " " << kk << " " << value <<  " " << error << " " << isconstant << std::endl;
	std::cout << "coeff3_" << ii << jj << kk << " is set" << std::endl; 
	dcoeff3[ii][jj][kk] = value;
	//coeff3[ii][jj][kk]->setVal(value);
	//	coeff3[ii][jj][kk]->setError(error);
	//	  coeff3[i][j][k]->setConstant(isconstant);
	//	coeff3[ii][jj][kk]->setConstant();

      }
    }
  }
  filestr.close();


  func0 = new TF3("func0",myfunction0, -1., 1., -1., 1., -1., 1.);

  double func0_min = findminimum(*func0,1000.,-1., 1., -1., 1., -1., 1.);
  double func0_max = findmaximum(*func0,   0.,-1., 1., -1., 1., -1., 1.);
  std::cout << "func0_min = " << func0_min << std::endl;
  std::cout << "func0_max = " << func0_max << std::endl;  

}


