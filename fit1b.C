#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TRandom2.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include <vector>
#include <string>

void fit1b(const char* outpdf="result2.pdf", int ntrials=1000, int entries=10) {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  TRandom2 r(0);
  std::vector<double> v_mu_chi2, v_mu_nll;
  v_mu_chi2.reserve(ntrials); v_mu_nll.reserve(ntrials);

  for (int t=0; t<ntrials; ++t) {
    TH1F h("h","",100,0,100);
    for (int i=0; i<entries; ++i) h.Fill(r.Gaus(50,10));

    TH1F* h1 = (TH1F*)h.Clone("h_chi2");
    h1->Fit("gaus","Q");
    TF1* f1 = h1->GetFunction("gaus");
    v_mu_chi2.push_back(f1->GetParameter(1));
    delete h1;

    TH1F* h2 = (TH1F*)h.Clone("h_nll");
    h2->Fit("gaus","QL");
    TF1* f2 = h2->GetFunction("gaus");
    v_mu_nll.push_back(f2->GetParameter(1));
    delete h2;
  }

  TH1F h_mu_chi2("h_mu_chi2","Mean from #chi^{2} fits;Mean;Counts",60,40,60);
  TH1F h_mu_nll ("h_mu_nll","Mean from NLL fits;Mean;Counts",60,40,60);
  for (size_t i=0;i<v_mu_chi2.size();++i){ h_mu_chi2.Fill(v_mu_chi2[i]); h_mu_nll.Fill(v_mu_nll[i]); }

  TCanvas c("c","Exercise 2: Means",900,900); 
  c.Divide(1,2);

  c.cd(1); h_mu_chi2.SetLineWidth(2); h_mu_chi2.Draw("hist");
  TLegend* L1 = new TLegend(0.66,0.78,0.90,0.90);
  L1->AddEntry(&h_mu_chi2, ("trials N="+std::to_string(ntrials)+", entries="+std::to_string(entries)).c_str(), "l"); 
  L1->Draw();

  c.cd(2); h_mu_nll.SetLineWidth(2); h_mu_nll.Draw("hist");
  TLegend* L2 = new TLegend(0.66,0.78,0.90,0.90);
  L2->AddEntry(&h_mu_nll, ("trials N="+std::to_string(ntrials)+", entries="+std::to_string(entries)).c_str(), "l"); 
  L2->Draw();

  c.SaveAs(outpdf);
}
