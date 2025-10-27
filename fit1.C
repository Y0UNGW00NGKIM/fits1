#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TRandom2.h"
#include "TROOT.h"


// fit1.C
// entries is the number of random samples filled into the histogram
void fit1(int entries=1000, bool save=false) {
   //Simple histogram fitting examples
  gROOT->Reset();  // useful to reset ROOT to a cleaner state

  TFile *tf=0;
  if (save) tf=new TFile("histo.root","recreate");

  TH1F *randomHist1 = new TH1F("randomHist1", "Random Histogram;x;frequency", 100, 0, 100);
  TRandom2 *generator=new TRandom2(0);  // parameter == seed, 0->use clock

  for (int i=0 ; i<entries ; i++){
    randomHist1->Fill(generator->Gaus(50,10)); // params: mean, sigma
  }
  // simple fits may be performed automatically
  // gStyle->SetOptFit(111);  // show reduced chi2 and params
  gStyle->SetOptFit(1111); // show reduced chi2, probability, and params
  randomHist1->Fit("gaus");  
  randomHist1->DrawCopy("e");  // "e" shows bin errors
  // Using DrawCopy vs Draw allows us to delete the original histogram
  // without removing it from the display.  If we save the histogran to a
  // file and close the file, it will be deleted from memory.


  // Above we used a built in function, gaus, in the fit
  // This function will be associated with the histogram
  // and may be retrieved to get parameter information
  // Refer to http://root.cern.ch/root/html/TF1.html
  // for a complete list of TF1 methods

  TF1 *fitfunc = randomHist1->GetFunction("gaus");
  cout << "\nFit Params and errors" << endl;
  cout << fitfunc->GetParameter(0) << " +- " << fitfunc->GetParError(0) << endl;
  cout << fitfunc->GetParameter(1) << " +- " << fitfunc->GetParError(1) << endl;
  cout << fitfunc->GetParameter(2) << " +- " << fitfunc->GetParError(2) << endl;
  cout << "Fit Probability: " << fitfunc->GetProb() << endl; // returns chi^2 p-value

  if (save) {
    tf->Write();
    tf->Close();
  }
  cout << "Use .q to exit root" << endl;
}
