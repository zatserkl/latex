#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include <iostream>
#include <cmath>

using std::cout;     using std::endl;

bool debug = 0;
bool gdebug = 1;
int debug_itimes = 0;

// utils.C stuff
TGraph* gtempClone(const char* name=0, const char* title=0);

Double_t I(Double_t x, Double_t tau, Double_t sigma)
{
   // static const pi2 = TMath::Sqrt(2.);

   // this is I' = exp(-((sigma^2/tau - x) * erfc((x/tau - 2*(sigma/tau)^2)
   // I'd like to rewrite: I'(x) = exp(-(x/t - sigma*sigma/2/tau/tau))

   Double_t arg;

   // Double_t exp_const = 0;
   // arg = x/tau - 2.*sigma*sigma/tau/tau;
   // if (TMath::Abs(arg) < 500.) exp_const = TMath::Exp(-arg);
   Double_t exp_const = 0;
   arg = x/tau - sigma*sigma/2./tau/tau;
   if (TMath::Abs(arg) < 500.) exp_const = TMath::Exp(-arg);

   Double_t erfc = 0;
   arg = (sigma*sigma/tau - x) / (sigma/TMath::Sqrt(2.));
   if (TMath::Abs(arg) < 500) erfc = TMath::Erfc(arg);
   else erfc = (arg < 0)? 2.: 0;

   Double_t I = exp_const * erfc;
   return I;
}

Double_t fpulse(Double_t *xx, Double_t *par)
{
   // fpulse = 0.5 * (I'(x,tau2,sigma) - I'(x,tau12,sigma))    NB tau2 in the first term

   Double_t fpulse = 0;

   Double_t a = par[0];
   Double_t x0 = par[1];
   Double_t tau1 = par[2];
   Double_t tau2 = par[3];
   Double_t sigma = par[4];

   Double_t x = *xx - x0;

   //if (debug) cout<< "a " << a << " x0 " << x0 << " tau1 " << tau1 << " tau2 " << tau2 << " sigma " << sigma << " x " << x <<endl;

   //if (x <= 0) return fpulse;
   
   Double_t I1 = I(x, tau2, sigma);
   Double_t tau12 = (tau1*tau2) / (tau1+tau2);
   Double_t I2 = I(x, tau12, sigma);

   fpulse = a*0.5*(I1 - I2);
   //-- fpulse = a*0.5*I1;
   return fpulse;
}

//------------------- convolution with scintillator decay function

Double_t IT(Double_t x, Double_t tau, Double_t T)
{
   if (debug_itimes > 0) {
      cout<< ".. debug_itimes " << debug_itimes << " x " << x << " tau " << tau << " T " << T <<endl;
      --debug_itimes;
   }

   Double_t IT = 0;
   if (x < 0) return IT;

   Double_t arg;

   Double_t exp_const = 0;
   arg = x/tau;
   if (TMath::Abs(arg) < 500.) exp_const = TMath::Exp(-arg);

   Double_t tau_exp_x = 0;
   Double_t inv_dtau = 1/tau - 1/T;
   arg = inv_dtau*x;
   if (TMath::Abs(arg) < 500.) {
      if (TMath::Abs(arg) < 1e-7) tau_exp_x = x;
      //-- if (TMath::Abs(arg) < 1e-12) tau_exp_x = x;
      else {
         tau_exp_x = 1./inv_dtau * (TMath::Exp(arg) - 1.);
      }
   }

   IT = exp_const * tau_exp_x;
   return IT;
}

Double_t fIT(Double_t *xx, Double_t *par)
{
   Double_t a = par[0];
   Double_t x0 = par[1];
   Double_t tau = par[2];
   Double_t T = par[3];

   Double_t x = *xx - x0;

   // normalizing: pnorm * snorm
   Double_t pnorm = 1./tau;
   Double_t snorm = 1./T;
   Double_t spnorm = pnorm*snorm;

   return (a*spnorm)*IT(x, tau, T);
}

Double_t fIT12(Double_t *xx, Double_t *par)
{
   Int_t npar = 0;
   Double_t a = par[npar++];
   Double_t x0 = par[npar++];
   Double_t tau1 = par[npar++];
   Double_t tau2 = par[npar++];
   Double_t T = par[npar++];

   Double_t x = *xx - x0;

   Double_t tau12 = tau1*tau2/(tau1+tau2);

   // normalizing: pnorm * snorm
   Double_t pnorm = (tau1+tau2)/tau2/tau2;
   Double_t snorm = 1./T;
   Double_t spnorm = pnorm*snorm;

   return (a*spnorm)*(IT(x, tau2, T) - IT(x, tau12, T));
}

Double_t fpulseT(Double_t *xx, Double_t *par)
{
   Double_t fpulseT = 0;

   Double_t a = par[0];
   Double_t x0 = par[1];
   Double_t tau1 = par[2];
   Double_t tau2 = par[3];
   Double_t T = par[4];

   Double_t x = *xx - x0;

   if (debug) cout<< "a " << a << " x0 " << x0 << " tau1 " << tau1 << " tau2 " << tau2 << " T " << T << " x " << x <<endl;

   //if (x <= 0) return fpulse;
   
   Double_t IT1 = IT(x, tau2, T);
   Double_t tau12 = (tau1*tau2) / (tau1+tau2);
   Double_t IT2 = IT(x, tau12, T);

   //fpulseT = (a/T)*(IT1 - IT2);
   //fpulseT = (a/T)*(IT1 - 0);
   //fpulseT = (a/T)*(0 - IT2);

   // normalizing: pnorm * snorm
   Double_t pnorm = (tau1+tau2)/tau2/tau2;
   Double_t snorm = 1./T;
   Double_t spnorm = pnorm*snorm;
   fpulseT = (a*spnorm)*(IT1 - IT2);
   return fpulseT;
}

void pulse_function(
      // Double_t tau1=1.
      Double_t tau1=5.
      , Double_t tau2=30.
      // , Double_t tau2=10.
      , Double_t sigma=1.
      , Double_t a=100.
      , Double_t x0=30.

      , Double_t xmin=0.
      , Double_t xmax=200.

      , bool todebug=0
      )
{
   debug = todebug;

   Double_t par[10];
   Int_t npar = 0;

   par[npar++] = a;
   par[npar++] = x0;
   par[npar++] = tau1;
   par[npar++] = tau2;
   par[npar++] = sigma;

   TF1 *fun_pulse = new TF1("fun_pulse", fpulse, xmin, xmax, npar);
   fun_pulse->SetParameters(par);
   cout<< "\n--> npar = " << npar <<endl;

   new TCanvas;
   fun_pulse->Draw();

   TGraph *gr_test = new TGraph(10000);
   gr_test->SetNameTitle("gr_test", "gr_test");
   gr_test->SetMarkerStyle(7);
   gr_test->SetMarkerColor(2);
   Int_t npoints = 0;

   // some test
   Double_t x1 = x0-5*tau1;
   Double_t x2 = x0+10*tau1;
   Double_t dx = 0.1;
   Double_t x = x1;
   while (x < x2) {
      Double_t fun_x = fun_pulse->Eval(x);
      Double_t tau12 = (tau1*tau2) / (tau1+tau2);
      Double_t I1 = I(x,tau2,sigma);
      Double_t I2 = I(x,tau12,sigma);
      cout<< "x " << x << "\t fun_x = " << fun_x << "\t I1 = " << I1 << "\t I2 = " << I2 <<endl;
      gr_test->SetPoint(npoints++,x, fun_x);
      x += dx;
   }
   gr_test->Set(npoints);
   new TCanvas;
   gr_test->Draw("ap");

   //--------------- convolution with scintillator decay function

   Double_t T = 40.;

   Double_t parT[10];
   Int_t nparT = 0;

   parT[nparT++] = a;
   parT[nparT++] = x0;
   parT[nparT++] = tau1;
   parT[nparT++] = tau2;
   parT[nparT++] = T;

   TF1 *fun_pulseT = new TF1("fun_pulseT", fpulseT, xmin, xmax, nparT);
   // TF1 *fun_pulseT = new TF1("fun_pulseT", fpulseT, -50, 50, nparT);
   fun_pulseT->SetLineColor(8);
   fun_pulseT->SetTitle("fun_pulseT;x;pulseT");
   fun_pulseT->SetParameters(parT);
   cout<< "\n--> nparT = " << nparT <<endl;

   new TCanvas;
   fun_pulseT->Draw();
}

void test_IT(
        Double_t tau1=5.
      , Double_t tau2=20.
      , Double_t T=40.
      , Double_t a=100.
      , Double_t x0=30.

      , Double_t xmin=0.
      , Double_t xmax=200.

      , bool todebug=true
      )
{
   debug = todebug;

   cout<< "tau1 " << tau1 << " tau2 " << tau2 << " T " << T <<endl;

   Double_t tau = 0;

   // IT1: use tau = tau2

   tau = tau2;
   cout<< "tau = " << tau <<endl;

   Double_t par1[10];
   Int_t npar1 = 0;

   par1[npar1++] = a;
   par1[npar1++] = x0;
   par1[npar1++] = tau;
   par1[npar1++] = T;
   cout<< "\t npar1 = " << npar1 <<endl;

   TF1* fun_IT1 = new TF1("fun_IT1", fIT, xmin, xmax, npar1);
   fun_IT1->SetNpx(1000);
   fun_IT1->SetTitle("fun_IT1");
   fun_IT1->SetLineColor(2);
   fun_IT1->SetParameters(par1);

   //debug_itimes = 1;

   new TCanvas;
   fun_IT1->DrawCopy();

   // IT2: use tau = tau12

   Double_t tau12 = tau1*tau2/(tau1+tau2);
   tau = tau12;
   cout<< "tau = " << tau <<endl;

   Double_t par2[10];
   Int_t npar2 = 0;

   par2[npar2++] = a;
   par2[npar2++] = x0;
   par2[npar2++] = tau;
   par2[npar2++] = T;
   cout<< "\t npar2 = " << npar2 <<endl;

   TF1* fun_IT2 = new TF1("fun_IT2", fIT, xmin, xmax, npar2);
   fun_IT2->SetNpx(1000);
   fun_IT2->SetTitle("fun_IT2");
   fun_IT2->SetLineColor(4);
   fun_IT2->SetParameters(par2);

   //debug_itimes = 1;

   new TCanvas;
   fun_IT2->DrawCopy();

   // IT12

   Double_t par12[10];
   Int_t npar12 = 0;

   par12[npar12++] = a;
   par12[npar12++] = x0;
   par12[npar12++] = tau1;
   par12[npar12++] = tau2;
   par12[npar12++] = T;
   cout<< "\t npar12 = " << npar12 <<endl;

   TF1* fun_IT12 = new TF1("fun_IT12", fIT12, xmin, xmax, npar12);
   fun_IT12->SetNpx(1000);
   fun_IT12->SetTitle(Form("fun_IT12   #tau_{1}=%0.1f, #tau_{2}=%0.1f, T=%0.1f",tau1,tau2,T));
   fun_IT12->SetLineColor(8);
   fun_IT12->SetParameters(par12);

   //debug_itimes = 1;

   new TCanvas;
   fun_IT12->DrawCopy();
}

void fit(const char* ifname="Co60_STM_LSO2x2_10pF_Ortamp_split.xml.root"
      , Int_t evt=10
      , Int_t ch=3
      )
{
   TFile *ifile = TFile::Open(ifname);
   if (!ifile) {
      cout<< "File not found: " << ifname <<endl;
      return;
   }

   TTree *tree = (TTree*) ifile->Get("t");
   if (!tree) {
      cout<< "TTree \"t\" was not found" <<endl;
      return;
   }

   new TCanvas;
   tree->Draw("-y:x",Form("ch==%d && evt==%d", ch, evt));
   
   TGraph* gr = gtempClone(Form("gr_ch%d_evt_%d",ch,evt), Form("gr_ch%d_evt_%d",ch,evt));
   gr->SetMarkerStyle(7);
   gr->SetMarkerColor(2);
   gr->Draw("ap");

   Double_t bkgmin = 0;
   Double_t bkgmax = 40;
   Double_t sigmin = 40;
   Double_t sigmax = 120;

   TH1F *hbkg = new TH1F("hbkg", "hbkg", 400, -100, 100);

   // find the level of flat background

   Double_t ybkg[1024];
   Int_t np_bkg = 0;
   for (int i=0; i<1024; ++i) {
      if (gr->GetX()[i] < bkgmin) continue;
      if (gr->GetX()[i] > bkgmax) break;
      ybkg[np_bkg] = gr->GetY()[i];
      ++np_bkg;
   }
   // sort bkg array to eliminate possible USB spikes: three highest points
   Int_t index[1024];
   TMath::Sort(np_bkg, ybkg, index, kFALSE);
   hbkg->Reset();
   for (int i=0; i<np_bkg-3; ++i) hbkg->Fill(ybkg[index[i]]);
   hbkg->Fit("gaus", "L0", "goff");
   Double_t bkg_mean = hbkg->GetFunction("gaus")->GetParameter("Mean");
   Double_t bkg_sigma = hbkg->GetFunction("gaus")->GetParameter("Sigma");
   if (gdebug) {
      hbkg->SetNameTitle(Form("hbkg_evt_%d_ch_%d",evt,ch),Form("hbkg_evt_%d_ch_%d",evt,ch));
      new TCanvas;
      hbkg->DrawCopy();
   }
   cout<< "bkg_mean = " << bkg_mean << " bkg_sigma = " << bkg_sigma <<endl;

   // make graph
   Double_t x[1024], y[1024], ex[1024], ey[1024];
   Int_t np = 0;

   for (int i=0; i<gr->GetN(); ++i) {
      if (gr->GetX()[i] < sigmin) continue;
      if (gr->GetX()[i] > sigmax) break;
      x[np] = gr->GetX()[i];
      y[np] = gr->GetY()[i] - bkg_mean;
      ex[np] = 0;
      ey[np] = bkg_sigma;
      ++np;
   }

   TGraphErrors *grsig = new TGraphErrors(np, x, y, ex, ey);
   grsig->SetNameTitle(Form("grsig_evt_%d_ch_%d",evt,ch), Form("grsig_evt_%d_ch_%d",evt,ch));
   grsig->SetMarkerStyle(24);
   grsig->SetMarkerColor(8);
   grsig->SetLineColor(8);
   if (gdebug) {
      new TCanvas;
      grsig->Draw("ap");
   }

   //new TCanvas;
   //grsig->Draw("ap");

   Double_t integral = 0;
   for (int i=0; i<np-1; ++i) integral += 0.5*(grsig->GetY()[i] + grsig->GetY()[i+1]) * (grsig->GetX()[i+1]-grsig->GetX()[i]);

   // normalize to pC by division to R = 50 Hom
   integral /= 50.;
   cout<< ".. integral = " << integral << " pC" <<endl;

   // IT12

   Double_t a = 1000;
   Double_t x0 = 50;
   Double_t tau1 = 1;
   Double_t tau2 = 5;
   Double_t T = 10;

   Double_t par12[10];
   Int_t npar12 = 0;

   par12[npar12++] = a;
   par12[npar12++] = x0;
   par12[npar12++] = tau1;
   par12[npar12++] = tau2;
   par12[npar12++] = T;
   cout<< "\t npar12 = " << npar12 <<endl;

   TF1* fun_IT12 = new TF1("fun_IT12", fIT12, sigmin, sigmax, npar12);
   fun_IT12->SetNpx(1000);
   fun_IT12->SetTitle(Form("fun_IT12   #tau_{1}=%0.1f, #tau_{2}=%0.1f, T=%0.1f",tau1,tau2,T));
   fun_IT12->SetLineColor(2);

   fun_IT12->SetParameters(par12);
   fun_IT12->SetParName(0, "a");
   fun_IT12->SetParName(1, "x0");
   fun_IT12->SetParName(2, "tau1");
   fun_IT12->SetParName(3, "tau2");
   fun_IT12->SetParName(4, "T");

   fun_IT12->FixParameter(2,1);

   fun_IT12->SetParLimits(3, 1., 20.);
   fun_IT12->SetParLimits(4, 1., 60.);

   //new TCanvas;
   //fun_IT12->DrawCopy();

   grsig->Fit(fun_IT12, "R", "", sigmin,sigmax);

   Double_t fit_x0 = grsig->GetFunction("fun_IT12")->GetParameter("x0");
   Double_t fit_tau2 = grsig->GetFunction("fun_IT12")->GetParameter("tau2");

   // refit

   Double_t sigmax_refit = grsig->GetFunction("fun_IT12")->GetMaximumX();
   TGraphErrors *grsig_refit = (TGraphErrors*) grsig->Clone("grsig_refit");
   grsig_refit->SetTitle(Form("%s refit",grsig_refit->GetTitle()));
   new TCanvas;
   grsig_refit->Draw("ap");
   fun_IT12->SetParameters(grsig->GetFunction("fun_IT12")->GetParameters());
   for (int i=0; i<fun_IT12->GetNpar(); ++i) cout<< i << "\t " << fun_IT12->GetParName(i) << "\t " << fun_IT12->GetParameter(i) <<endl;
   //grsig_refit->GetFunction("fun_IT12")->FixParameter(4, grsig->GetFunction("fun_IT12")->GetParameter(4));
   fun_IT12->FixParameter(4, grsig->GetFunction("fun_IT12")->GetParameter(4));
   // grsig_refit->Fit(fun_IT12, "R", "", sigmin,fit_x0+2.*fit_tau2);
   grsig_refit->Fit(fun_IT12, "R", "", sigmin,(fit_x0+sigmax_refit)/2.);

   /// TF1* fun_IT12_refit = (TF1*) grsig_refit->GetFunction("fun_IT12")->Clone();
   /// fun_IT12_refit->SetParameters(fun_IT12->GetParameters());
   /// // Double_t par_T = grsig_refit->GetFunction("fun_IT12")->GetParameter(4);
   /// Double_t par_T = fun_IT12_refit->GetParameter(4);
   /// cout<< "par_T = " << par_T <<endl;
   /// // grsig_refit->GetFunction("fun_IT12")->SetParLimits(4, 0.8*par_T, 1.2*par_T);
   /// fun_IT12_refit->FixParameter(4, par_T);
   /// grsig_refit->Fit(fun_IT12_refit, "R", "", sigmin,sigmax_refit);
}
