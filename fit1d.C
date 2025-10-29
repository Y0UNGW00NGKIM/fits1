// fit1d.C
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TPad.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/DistFunc.h"
#include <vector>
#include <algorithm>
#include <cmath>

static inline double gint(double x1, double x2, double mu, double s) {
    return ROOT::Math::normal_cdf(x2, s, mu) - ROOT::Math::normal_cdf(x1, s, mu);
}

static double nll_hist(const TH1 *h, double a, double mu, double s) {
    const double eps = 1e-12;
    double v = 0.0;
    for (int i = 1; i <= h->GetNbinsX(); i += 1) {
        double m = a * gint(h->GetXaxis()->GetBinLowEdge(i), h->GetXaxis()->GetBinUpEdge(i), mu, s);
        if (m < eps) m = eps;
        double y = h->GetBinContent(i);
        v += 2.0 * (m - y * std::log(m));
    }
    return v;
}

static double x2_hist(const TH1 *h, double a, double mu, double s) {
    const double eps = 1e-12;
    double v = 0.0;
    for (int i = 1; i <= h->GetNbinsX(); i += 1) {
        double m = a * gint(h->GetXaxis()->GetBinLowEdge(i), h->GetXaxis()->GetBinUpEdge(i), mu, s);
        if (m < eps) m = eps;
        double y = h->GetBinContent(i);
        double d = y - m;
        v += d * d / std::max(m, 1.0);
    }
    return v;
}

struct one_fit { bool ok; double val; double a; double s; int st1; int st2; };

static one_fit minimize_as(const TH1 *h, bool use_nll, double mu, double a0, double s0) {
    double xmin = h->GetXaxis()->GetXmin(), xmax = h->GetXaxis()->GetXmax();
    double bw = (xmax - xmin) / h->GetNbinsX();
    double a_lo = 0.0, a_hi = 1e12, s_lo = 0.5 * bw, s_hi = 0.5 * (xmax - xmin);
    auto obj = [&](const double *p) {
        double a = p[0], s = p[1];
        if (a <= a_lo || a >= a_hi || s <= s_lo || s >= s_hi) return 1e300;
        return use_nll ? nll_hist(h, a, mu, s) : x2_hist(h, a, mu, s);
    };
    ROOT::Math::Functor F(obj, 2);
    std::unique_ptr<ROOT::Math::Minimizer> M1(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex"));
    M1->SetMaxFunctionCalls(20000); M1->SetMaxIterations(20000); M1->SetFunction(F);
    M1->SetLimitedVariable(0, "A", std::max(1.0, a0), 0.2 * std::max(1.0, a0), a_lo, a_hi);
    M1->SetLimitedVariable(1, "s", std::max(s0, s_lo), 0.2 * std::max(1.0, s0), s_lo, s_hi);
    M1->Minimize(); const double *x1 = M1->X(); double v1 = M1->MinValue();
    std::unique_ptr<ROOT::Math::Minimizer> M2(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
    M2->SetMaxFunctionCalls(20000); M2->SetMaxIterations(20000); M2->SetFunction(F);
    M2->SetLimitedVariable(0, "A", x1[0], 0.05 * std::max(1.0, x1[0]), a_lo, a_hi);
    M2->SetLimitedVariable(1, "s", x1[1], 0.05 * std::max(1.0, x1[1]), s_lo, s_hi);
    bool ok = M2->Minimize(); const double *x = M2->X(); double v = ok ? M2->MinValue() : v1;
    return { ok, v, x[0], x[1], M1->Status(), M2->Status() };
}

struct scan_out { double mu_hat; double mu_err; };

static scan_out scan_curve(TH1 *h, bool use_nll, int npts, double mu_center, double halfspan, TGraph &gr) {
    double xmin = h->GetXaxis()->GetXmin(), xmax = h->GetXaxis()->GetXmax();
    double lo = std::max(xmin, mu_center - halfspan), hi = std::min(xmax, mu_center + halfspan);
    double a_seed = std::max(1.0, h->Integral());
    double s_seed = std::max(h->GetRMS(), 1e-3);
    std::vector<double> xs, ys; xs.reserve(npts + 1); ys.reserve(npts + 1);
    for (int i = 0; i <= npts; ++i) {
        double mu = lo + (hi - lo) * i / npts;
        one_fit r = minimize_as(h, use_nll, mu, a_seed, s_seed);
        if (!r.ok) continue;
        xs.push_back(mu); ys.push_back(r.val);
        a_seed = r.a; s_seed = r.s;
    }
    if (xs.size() < 5) { gr.SetPoint(0, 0.0, 0.0); return { mu_center, std::numeric_limits<double>::quiet_NaN() }; }
    auto it = std::min_element(ys.begin(), ys.end());
    double vmin = *it; double muhat = xs[it - ys.begin()];
    for (size_t k = 0; k < xs.size(); ++k) gr.SetPoint((int)k, xs[k], ys[k] - vmin);
    double L = muhat, R = muhat;
    for (size_t k = 0; k < xs.size(); ++k) if (xs[k] < muhat && ys[k] - vmin <= 1.0) L = xs[k];
    for (size_t k = 0; k < xs.size(); ++k) if (xs[k] > muhat && ys[k] - vmin <= 1.0) { R = xs[k]; break; }
    return { muhat, 0.5 * (R - L) };
}

static void draw_panel(TGraph &gr, const char *title, const char *ytitle, double mu_hat, double emu, double lx1, double ly1, double lx2, double ly2) {
    gr.SetTitle(title);
    gr.SetLineWidth(2);
    gr.Draw("AL");
    gr.GetXaxis()->SetTitle("mean");
    gr.GetYaxis()->SetTitle(ytitle);
    double xlo = gr.GetXaxis()->GetXmin(), xhi = gr.GetXaxis()->GetXmax();
    TLine l1(xlo, 1.0, xhi, 1.0), l4(xlo, 4.0, xhi, 4.0);
    l1.SetLineStyle(2); l4.SetLineStyle(2);
    l1.Draw(); l4.Draw();
    double yhi = gPad ? gPad->GetUymax() : gr.GetYaxis()->GetXmax();
    TLine lv(mu_hat, 0.0, mu_hat, 0.98 * yhi); lv.SetLineColor(kRed); lv.SetLineWidth(3); lv.Draw();
    TLegend leg(lx1, ly1, lx2, ly2);
    leg.AddEntry((TObject*)0, Form("#mu = %.3f #pm %.3f", mu_hat, emu), "");
    leg.AddEntry((TObject*)0, "#Delta=1,4 guides", "");
    leg.Draw();
}

void fit1d(const char *outpdf = "result4.pdf",
           const char *file_ll = "histo25.root", const char *hist_ll = "randomHist1",
           const char *file_x2 = "histo1k.root", const char *hist_x2 = "randomHist1",
           int npts = 201, double span_sigma = 6.0) {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    TFile f1(file_ll, "READ"); TH1 *h_ll = nullptr; f1.GetObject(hist_ll, h_ll);
    TFile f2(file_x2, "READ"); TH1 *h_x2 = nullptr; f2.GetObject(hist_x2, h_x2);
    if (!h_ll || !h_x2) { printf("missing histogram(s)\n"); return; }
    double xmin = h_x2->GetXaxis()->GetXmin(), xmax = h_x2->GetXaxis()->GetXmax();
    double bw = (xmax - xmin) / h_x2->GetNbinsX();
    double mu_c = h_x2->GetMean();
    double s_mu = h_x2->GetRMS() / std::sqrt(std::max(1.0, h_x2->GetEntries()));
    double half = std::max(bw, s_mu) * span_sigma + 5.0 * bw;
    TGraph g_ll, g_x2;
    scan_out s_ll = scan_curve(h_ll, true,  npts, mu_c, half, g_ll);
    scan_out s_x2 = scan_curve(h_x2, false, npts, mu_c, half, g_x2);
    TCanvas c("c", "result4", 1000, 900);
    c.Divide(1, 2);
    c.cd(1); draw_panel(g_ll, "-2 ln L scan", "#Delta(-2 ln L)", s_ll.mu_hat, s_ll.mu_err, 0.14, 0.74, 0.48, 0.89);
    c.cd(2); draw_panel(g_x2, "#chi^{2} scan", "#Delta#chi^{2}", s_x2.mu_hat, s_x2.mu_err, 0.60, 0.74, 0.94, 0.89);
    c.SaveAs(outpdf);
}
