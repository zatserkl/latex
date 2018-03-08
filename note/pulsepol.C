#include <TROOT.h>
#include <TEnv.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <Math/GSLIntegrator.h>

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TSpline.h>
#include <TH2.h>

#include <iostream>
#include <fstream>
#include <cmath>

using std::cout;     using std::endl;

// utils.C stuff
void pic(TString pathname="", TString suffix="");
TH1* htemp(const char* name=0, const char* title=0);
TF1* ftemp(const char* name=0, const char* title=0);
TGraph* gtemp(const char* name=0, const char* title=0);
void fitgaus(Double_t xmin=0, Double_t xmax=0, const char* opt="", const char* gopt="", TH1* h=0, TGraph* gr=0);
void fitpol(Double_t xmin=0, Double_t xmax=0, Int_t power=1, const char* opt="", const char* gopt="", TH1* h=0, TGraph* gr=0);
void rightgaus();
void leftgaus();
void left();
void right();
Int_t nbi(Int_t n=0);
const char* addtit(const char* title, TCanvas* can=0);

Double_t fplog(Double_t *xx, Double_t *par)
{
   Double_t anorm = par[0];
   Double_t x0 = par[1];
   Double_t tau = par[2];

   Double_t x = *xx - x0;

   //-- Double_t log = TMath::Log(anorm * (1. - TMath::Exp(-x/tau)));
   Double_t log = anorm * (1. - TMath::Exp(-x/tau));
   Double_t fun = log > 0? log: 0;
   return fun;
}

Double_t fpdlog(Double_t *xx, Double_t *par)
{
   Double_t anorm = par[0];
   Double_t x0 = par[1];
   Double_t tau = par[2];
   Double_t taud = par[3];

   Double_t x = *xx - x0;

   //-- Double_t log = TMath::Log(anorm * (1. - TMath::Exp(-x/tau))) - x/taud;
   Double_t log = anorm * (1. - TMath::Exp(-x/tau)) - x/taud;
   Double_t fun = log > 0? log: 0;
   return fun;
}

// TGraph* pulsepol(const char *ifname="Co60_STM_2x2_1.xml.root", Int_t evt=8, Int_t ch=3)
TGraph* pulsepol(const char *ifname="Co60_STM_2x2_1.xml.root", Int_t evt=15, Int_t ch=3)
{
   Int_t color = 2;

   TFile *ifile = TFile::Open(ifname);
   if (!ifile) {
      cout<< "File not found: " << ifname <<endl;
      return 0;
   }
   TTree *tree = (TTree*) ifile->Get("t");
   if (!tree) {
      cout<< "Tree t was not found" <<endl;
      return 0;
   }

   new TCanvas;
   tree->Draw("-y", Form("evt==%d &&ch==%d &&x>20&&x<120", evt,ch));
   // fitgaus(); right();
   Double_t hmean = tree->GetHistogram()->GetMean();
   Double_t hRMS = tree->GetHistogram()->GetRMS();
   Double_t hmin = tree->GetHistogram()->GetXaxis()->GetXmin();
   Double_t hmax = tree->GetHistogram()->GetXaxis()->GetXmax();
   Double_t xmin = hmean - 1.*hRMS;
   if (xmin < hmin) xmin = hmin;
   Double_t xmax = hmean + 1.*hRMS;
   if (xmax > hmax) xmax = hmax;
   fitgaus(xmin,xmax); right();
   Double_t bkg_mean = htemp()->GetFunction("gaus")->GetParameter("Mean");
   Double_t bkg_sigma = TMath::Abs(htemp()->GetFunction("gaus")->GetParameter("Sigma"));

   new TCanvas;
   tree->Draw("-y:x", Form("evt==%d &&ch==%d", evt,ch));

   Double_t x[1024];
   Double_t y[1024];
   Double_t ex[1024];
   Double_t ey[1024];
   Int_t np = 0;

   for (int i=0; i<tree->GetSelectedRows(); ++i) {
      x[np] = tree->GetV2()[i];
      y[np] = tree->GetV1()[i] - bkg_mean;
      ex[np] = 0.;
      ey[np] = bkg_sigma;
      ++np;
   }

   TGraphErrors *gr = new TGraphErrors(np, x, y, ex, ey);
   gr->SetNameTitle("gr", Form("%s;ns",ifname));
   gr->SetMarkerStyle(7);
   color = 2;
   gr->SetMarkerColor(color);
   gr->SetLineColor(color);

   new TCanvas;
   gr->Draw("ap");

   // spline

   TSpline3 *spline = new TSpline3(Form("Spline for %s",ifname),gr);
   spline->SetLineColor(kRed);
   new TCanvas;
   spline->Draw();

   // // spline smoothing

   // Double_t xs[1024];
   // Double_t ys[1024];
   // Double_t exs[1024];
   // Double_t eys[1024];
   // Int_t nps = np;
   // Int_t nrange = 10;       // +-2
   // for (int i=0; i<nrange; ++i) {
   //    xs[i] = x[i];
   //    ys[i] = y[i];
   //    exs[i] = ex[i];
   //    eys[i] = ey[i];
   //    xs[np-1-i] = x[np-1-i];
   //    ys[np-1-i] = y[np-1-i];
   //    exs[np-1-i] = ex[np-1-i];
   //    eys[np-1-i] = ey[np-1-i];
   // }
   // for (int i=nrange; i<np-nrange; ++i) {
   //    TSpline3 spline("spline", &x[i-nrange], &y[i-nrange], 2*nrange+1);
   //    ys[i] = spline.Eval(x[i]);
   //    xs[i] = x[i];
   //    exs[i] = ex[i];
   //    eys[i] = ey[i];
   // }

   // TGraphErrors *grs = new TGraphErrors(nps, xs, ys, exs, eys);
   // grs->SetNameTitle("grs", Form("Spline smoothing for %s;ns",ifname));
   // grs->SetMarkerStyle(7);
   // color = 6;
   // grs->SetMarkerColor(color);
   // grs->SetLineColor(color);

   // new TCanvas;
   // grs->Draw("ap");

   // derivative

   Double_t xd[1024];
   Double_t yd[1024];
   Double_t exd[1024];
   Double_t eyd[1024];
   Int_t npd = np;
   xd[0] = x[0];
   yd[0] = bkg_sigma;
   exd[0] = 0;
   eyd[0] = 0;
   for (int i=1; i<np; ++i) {
      xd[i] = x[i];
      exd[i] = 0;
      yd[i] = (y[i] - y[i-1]) / (x[i] - x[i-1]);
      // eyd[i] = 0.5*(ey[i] + ey[i-1]);
      eyd[i] = 0.;
   }

   TGraphErrors *grd = new TGraphErrors(npd, xd, yd, exd, eyd);
   grd->SetNameTitle("grd", Form("Derivative for %s;ns",ifname));
   grd->SetMarkerStyle(7);
   color = 2;
   grd->SetMarkerColor(color);
   grd->SetLineColor(color);

   new TCanvas;
   grd->Draw("ap");

   // integral

   Double_t xi[1024];
   Double_t yi[1024];
   Double_t exi[1024];
   Double_t eyi[1024];
   Int_t npi = np;

   Double_t integral = 0;

   xi[0] = x[0];
   yi[0] = y[0];
   exi[0] = ex[0];
   eyi[0] = ey[0];
   for (int i=1; i<np; ++i) {
      integral += 0.5 * (y[i] + y[i-1]) * (x[i] - x[i-1]) / 50. * 10/1.6;  // in Npe
      xi[i] = x[i];
      yi[i] = integral;
      exi[i] = ex[i];
      eyi[i] = ey[i];
   }

   cout<< "integral = " << integral <<endl;

   TGraphErrors *gri = new TGraphErrors(npi, xi, yi, exi, eyi);
   gri->SetNameTitle("gri", Form("Integral for %s;ns;Npe",ifname));
   gri->SetMarkerStyle(7);
   color = 8;
   gri->SetMarkerColor(color);
   gri->SetLineColor(color);

   new TCanvas;
   gri->Draw("ap");

   Double_t xilog[1024];
   Double_t yilog[1024];
   Double_t exilog[1024];
   Double_t eyilog[1024];
   Int_t npilog = npi;

   for (int i=0; i<npilog; ++i) {
      xilog[i] = xi[i];
      // yilog[i] = yi[i] > 1.? TMath::Log(yi[i]/integral): 0;
      //-- yilog[i] = yi[i] > 0.? TMath::Log(yi[i]): 0;
      yilog[i] = yi[i] > bkg_sigma? TMath::Log(yi[i]): 0;
      exilog[i] = 0;
      Double_t dy = eyi[i];
      eyilog[i] = TMath::Abs(TMath::Log(dy/integral));
      eyilog[i] = 0.;
   }

   TGraphErrors *grilog = new TGraphErrors(npilog, xilog, yilog, exilog, eyilog);
   grilog->SetNameTitle("grilog", Form("Normalized log for %s;ns;log(Npe)",ifname));
   grilog->SetMarkerStyle(7);
   color = 4;
   grilog->SetMarkerColor(color);
   grilog->SetLineColor(color);
   
   Double_t par[10];
   Double_t anorm = 10;
   Double_t x0 = 25;
   Double_t tau = 10;

   par[0] = anorm;
   par[1] = x0;
   par[2] = tau;

   // Double_t xfit_min = 150;
   Double_t xfit_min = 100;
   Double_t xfit_max = 240;

   TF1* fun_fplog = new TF1("fun_fplog", fplog, 20,300, 3);

   fun_fplog->SetParameters(par);

   new TCanvas;
   grilog->Draw("ap");
   grilog->Fit(fun_fplog, "R", "", xfit_min,xfit_max);

   TGraphErrors *grilog1 = (TGraphErrors*) grilog->Clone("grilog1");
   color = 46;
   grilog1->SetMarkerColor(color);
   grilog1->SetLineColor(color);

   TF1* fun_fpdlog = new TF1("fun_fpdlog", fpdlog, 20,300, 4);
   
   Double_t taud = 30;

   par[0] = anorm;
   par[1] = x0;
   par[2] = tau;
   par[3] = taud;

   fun_fpdlog->SetParameters(par);

   new TCanvas;
   grilog1->Draw("ap");
   grilog1->Fit(fun_fpdlog, "R", "", xfit_min,xfit_max);

   // return gri;

   return 0;
}
