#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TROOT.h"
#include <cmath>

static double nll_poisson_hist_scaled(const TH1* h, TF1* f, double scale) {
  const double eps = 1e-12;
  double s = 0.0;
  for (int i=1;i<=h->GetNbinsX();++i) {
    double x1 = h->GetXaxis()->GetBinLowEdge(i);
    double x2 = h->GetXaxis()->GetBinUpEdge(i);
    double mu = scale * f->Integral(x1,x2);
    double y  = h->GetBinContent(i);
    if (mu < eps) mu = eps;
    s += 2.0 * (mu - y * std::log(mu));
  }
  return s;
}

void fit1c(const char* outpdf="result3.pdf", const char* filename="histo25.root", const char* histname="randomHist1", int ntoys=5000, unsigned seed=12345) {  // left seed in their just because assignment asked for p-value and i want to keep the distribution the same for grading process as well
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  TFile fin(filename,"READ");
  TH1F* hdata = nullptr;
  fin.GetObject(histname,hdata);
  if (!hdata) {
      printf("hist not found\n"); 
      return; 
  }

  TH1F* h = (TH1F*)hdata->Clone("h_work");
  const double xmin = h->GetXaxis()->GetXmin();
  const double xmax = h->GetXaxis()->GetXmax();

  TF1 g("gfit","gaus", xmin, xmax);
  g.SetParameters(h->GetMaximum(), h->GetMean(), h->GetRMS());
  h->Fit(&g,"QLR");
  g.SetRange(xmin, xmax);

  const double data_sum  = h->Integral();
  const double model_sum = g.Integral(xmin, xmax);
  const double scale = (model_sum>0) ? (data_sum/model_sum) : 1.0;

  const double nll_data = nll_poisson_hist_scaled(h, &g, scale);

  TRandom3 rng(seed);
  TH1F toy("toy","",h->GetNbinsX(), xmin, xmax);

  double vmin = 0, vmax = 0;
  for (int t=0; t<ntoys; ++t) {
    toy.Reset();
    for (int i=1;i<=toy.GetNbinsX();++i) {
      double x1 = toy.GetXaxis()->GetBinLowEdge(i);
      double x2 = toy.GetXaxis()->GetBinUpEdge(i);
      double mu = scale * g.Integral(x1,x2);
      int y = rng.Poisson(mu);
      toy.SetBinContent(i, y);
    }
    double v = nll_poisson_hist_scaled(&toy, &g, scale);
    if (t==0) {
        vmin=v; vmax=v; 
    } else {
        if (v<vmin) vmin=v; 
        if (v>vmax) vmax=v; 
    }
  }

  TH1F hdist("hdist","NLL from toys;NLL;Counts",120, vmin, vmax);

  rng.SetSeed(seed);
  int ge = 0;
  for (int t=0; t<ntoys; ++t) {
    toy.Reset();
    for (int i=1;i<=toy.GetNbinsX();++i) {
      double x1 = toy.GetXaxis()->GetBinLowEdge(i);
      double x2 = toy.GetXaxis()->GetBinUpEdge(i);
      double mu = scale * g.Integral(x1,x2);
      int y = rng.Poisson(mu);
      toy.SetBinContent(i, y);
    }
    double v = nll_poisson_hist_scaled(&toy, &g, scale);
    hdist.Fill(v);
    if (v >= nll_data) ++ge;
  }

  double pvalue = (double)ge / (double)ntoys;
  printf("NLL(data) = %.6f,  p-value = %.6f\n", nll_data, pvalue);

  TCanvas c("c","Exercise 3: NLL toys",900,700);
  hdist.SetLineWidth(2); hdist.Draw("hist");
  TLine* L = new TLine(nll_data, 0, nll_data, hdist.GetMaximum()*1.05);
  L->SetLineColor(kRed); L->SetLineWidth(3); L->Draw();
  TLegend* leg = new TLegend(0.58,0.78,0.90,0.90);
  leg->AddEntry(&hdist, Form("toys N=%d", ntoys), "l");
  leg->AddEntry(L, Form("data NLL = %.3f", nll_data), "l");
  leg->Draw();
  c.SaveAs(outpdf);
}
