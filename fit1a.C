#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TRandom2.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include <vector>
#include <string>

void fit1a(const char* outpdf="result1.pdf", int ntrials=1000, int entries=1000) {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  TRandom2 r(0);
  TH1F h("randomHist1","Random Histogram;x;Entries",100,0,100);

  std::vector<double> v_chi2red, v_prob, v_mu, v_errmu;
  v_chi2red.reserve(ntrials); v_prob.reserve(ntrials); v_mu.reserve(ntrials); v_errmu.reserve(ntrials);

  for (int t=0; t<ntrials; ++t) {
    h.Reset();
    for (int i=0; i<entries; ++i) h.Fill(r.Gaus(50,10));
    h.Fit("gaus","Q");
    TF1* f = h.GetFunction("gaus");
    double chi2 = f->GetChisquare();
    double ndf  = f->GetNDF();
    v_chi2red.push_back(ndf>0 ? chi2/ndf : 0.0);
    v_prob.push_back(f->GetProb());
    v_mu.push_back(f->GetParameter(1));
    v_errmu.push_back(f->GetParError(1));
  }

  TH1F h_chi2("h_chi2","Reduced #chi^{2};Reduced #chi^{2};Counts",60,0,3);
  TH1F h_prob("h_prob","#chi^{2} probability;p-value;Counts",60,0,1);
  TH1F h_mu("h_mu","Fitted mean;Mean;Counts",60,40,60);
  TH1F h_emu("h_emu","Error on mean;Error on mean;Counts",60,0,2);
  for (size_t i=0;i<v_chi2red.size();++i){
    h_chi2.Fill(v_chi2red[i]); h_prob.Fill(v_prob[i]); h_mu.Fill(v_mu[i]); h_emu.Fill(v_errmu[i]);
  }

  TCanvas c("c","Exercise 1: Distributions",1000,800); 
  c.Divide(2,2);

  c.cd(1); h_chi2.SetLineWidth(2); h_chi2.Draw("hist");
  TLegend* L1 = new TLegend(0.68,0.78,0.90,0.90); L1->AddEntry(&h_chi2, ("trials N="+std::to_string(ntrials)).c_str(), "l"); L1->Draw();

  c.cd(2); h_prob.SetLineWidth(2); h_prob.Draw("hist");
  TLegend* L2 = new TLegend(0.68,0.78,0.90,0.90); L2->AddEntry(&h_prob, ("trials N="+std::to_string(ntrials)).c_str(), "l"); L2->Draw();

  c.cd(3); h_mu.SetLineWidth(2);   h_mu.Draw("hist");
  TLegend* L3 = new TLegend(0.68,0.78,0.90,0.90); L3->AddEntry(&h_mu,  ("trials N="+std::to_string(ntrials)).c_str(), "l"); L3->Draw();

  c.cd(4); h_emu.SetLineWidth(2);  h_emu.Draw("hist");
  TLegend* L4 = new TLegend(0.68,0.78,0.90,0.90); L4->AddEntry(&h_emu, ("trials N="+std::to_string(ntrials)).c_str(), "l"); L4->Draw();

  c.SaveAs(outpdf);
}
