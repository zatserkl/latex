#include "drs.C"

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
#include <TH2.h>

#include <iostream>
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

bool debug = false;
bool gdebug = false;

class OscFit: public TObject
{
public:
   Int_t evt;
   Float_t adc[4];
   Float_t t[4];
   Float_t t0[4];
   Float_t tau[4];
   Float_t sigma[4];
   Float_t v[4];
   Float_t bkg[4];          // flat background before the signal
   Float_t chi2[4];
public:
   void clear() {
      evt = 0;
      for (int i=0; i<4; ++i) {
         adc[i] = 0;
         t[i] = 0;
         t0[i] = 0;
         tau[i] = 0;
         sigma[i] = 0;
         bkg[i] = 0;
         v[i] = 0;
         chi2[i] = 0;
      }
   }
   OscFit(): TObject() {clear();}
   ClassDef(OscFit,4);
};

#ifdef __MAKECINT__
#pragma link C++ class OscFit;
#endif

Double_t fpulse3_erfc(Double_t *xx, Double_t *par)
{
   Double_t fspulse = 0;

   Double_t a = par[0];
   Double_t x0 = par[1];
   Double_t tau = par[2];
   Double_t sigma = tau;

   Double_t x = *xx - x0;

   if (debug) cout<< "a " << a << " x0 " << x0 << " tau " << tau << " sigma " << sigma << " x " << x <<endl;

   Double_t arg;

   Double_t exp_const1 = 0;
   arg = sigma*sigma/(2.*tau*tau);
   if (TMath::Abs(arg) < 500.) exp_const1 = TMath::Exp(arg);
   //else cout<< "exp_const1: arg = " << arg <<endl;
   // exp_const1 = TMath::Exp(arg);

   Double_t exp_x1 = 0;
   arg = -x/tau;
   if (TMath::Abs(arg) < 500.) exp_x1 = TMath::Exp(arg);
   //else cout<< "exp_x1: arg = " << arg <<endl;
   // exp_x1 = TMath::Exp(arg);

   Double_t erfc1 = 0.;
   arg = (sigma*sigma/tau - x) / (sigma*TMath::Sqrt(2.));
   if (TMath::Abs(arg) < 500.) erfc1 = TMath::Erfc(arg);
   else erfc1 = (arg < 0.)? 2.: 0.;
   // erfc1 = TMath::Erfc(arg);

   Double_t I1 = (exp_const1*exp_x1)*erfc1;

   Double_t exp_const2 = 0;
   arg = 2.*sigma*sigma/(tau*tau);
   if (TMath::Abs(arg) < 500.) exp_const2 = TMath::Exp(arg);
   //else cout<< "exp_const2: arg = " << arg <<endl;
   // exp_const2 = TMath::Exp(arg);

   Double_t exp_x2 = 0;
   arg = -2*x/tau;
   if (TMath::Abs(arg) < 500.) exp_x2 = TMath::Exp(arg);
   //else cout<< "exp_x2: arg = " << arg <<endl;
   // exp_x2 = TMath::Exp(arg);

   Double_t erfc2 = 0;
   arg = (2.*sigma*sigma/tau - x) / (sigma*TMath::Sqrt(2.));
   if (TMath::Abs(arg) < 500.) erfc2 = TMath::Erfc(arg);
   //else erfc2 = (arg < 0.)? 2.: 0.;
   // erfc2 = TMath::Erfc(arg);

   Double_t I2 = (exp_const2*exp_x2)*erfc2;

   fspulse = a*(I1 - I2);
   return fspulse;
}

Double_t fpulse_erfc(Double_t *xx, Double_t *par)
{
   Double_t fspulse = 0;

   Double_t a = par[0];
   Double_t x0 = par[1];
   Double_t tau = par[2];
   Double_t sigma = par[3];

   Double_t x = *xx - x0;

   if (debug) cout<< "a " << a << " x0 " << x0 << " tau " << tau << " sigma " << sigma << " x " << x <<endl;

   Double_t arg;

   Double_t exp_const1 = 0;
   arg = sigma*sigma/(2.*tau*tau);
   if (TMath::Abs(arg) < 500.) exp_const1 = TMath::Exp(arg);
   //else cout<< "exp_const1: arg = " << arg <<endl;
   // exp_const1 = TMath::Exp(arg);

   Double_t exp_x1 = 0;
   arg = -x/tau;
   if (TMath::Abs(arg) < 500.) exp_x1 = TMath::Exp(arg);
   //else cout<< "exp_x1: arg = " << arg <<endl;
   // exp_x1 = TMath::Exp(arg);

   Double_t erfc1 = 0.;
   arg = (sigma*sigma/tau - x) / (sigma*TMath::Sqrt(2.));
   if (TMath::Abs(arg) < 500.) erfc1 = TMath::Erfc(arg);
   else erfc1 = (arg < 0.)? 2.: 0.;
   // erfc1 = TMath::Erfc(arg);

   Double_t I1 = (exp_const1*exp_x1)*erfc1;

   Double_t exp_const2 = 0;
   arg = 2.*sigma*sigma/(tau*tau);
   if (TMath::Abs(arg) < 500.) exp_const2 = TMath::Exp(arg);
   //else cout<< "exp_const2: arg = " << arg <<endl;
   // exp_const2 = TMath::Exp(arg);

   Double_t exp_x2 = 0;
   arg = -2*x/tau;
   if (TMath::Abs(arg) < 500.) exp_x2 = TMath::Exp(arg);
   //else cout<< "exp_x2: arg = " << arg <<endl;
   // exp_x2 = TMath::Exp(arg);

   Double_t erfc2 = 0;
   arg = (2.*sigma*sigma/tau - x) / (sigma*TMath::Sqrt(2.));
   if (TMath::Abs(arg) < 500.) erfc2 = TMath::Erfc(arg);
   //else erfc2 = (arg < 0.)? 2.: 0.;
   // erfc2 = TMath::Erfc(arg);

   Double_t I2 = (exp_const2*exp_x2)*erfc2;

   fspulse = a*(I1 - I2);
   return fspulse;
}

Double_t fsig(Double_t* xx, Double_t* par)
{
   Double_t x = *xx;
   Double_t a = par[0];
   Double_t mean = par[1];
   Double_t sigma = par[2];

   Double_t fsig = 0;
   Double_t arg = std::pow( (x-mean)/2/sigma, 2);
   if (std::abs(arg) > 50) return fsig;

   fsig = a*exp(-arg);
   return fsig;
}

Double_t fbkg(Double_t* xx, Double_t* par)
{
   *xx = *xx;
   return par[0];
}

Double_t fsigbkg(Double_t* xx, Double_t* par)
{
   Double_t fsigbkg = fbkg(xx,par) + fsig(xx, &par[1]);
   return fsigbkg;
}

Double_t fsigbkg_erfc(Double_t* xx, Double_t* par)
{
   Double_t fsigbkg = fbkg(xx,par) + fpulse_erfc(xx, &par[1]);
   return fsigbkg;
}

OscEvent* drs4fit(const char* ifname="sipm_25pe_split.xml.root"
      , Float_t bkgmin=75, Float_t bkgmax=95, Float_t sigmin=95, Float_t sigmax=104
      , Int_t entry_first=0, Int_t entry_last=-1
      , bool setdebug=false, bool setgdebug=false
      )
{
   if (setdebug) debug = true;
   if (setgdebug) gdebug = true;

   Double_t thres = 20.;
   Int_t nthres_min = 5;                  // data are required to contain at least nthres_min points over the thres value
   Double_t ysaturation = 499;

   TFile* ifile = TFile::Open(ifname);
   if (!ifile) {
      cout<< "File not found: " << ifname <<endl;
      return 0;
   }
   cout<< "Processing file " << ifname <<endl;

   TTree* tree = (TTree*) ifile->Get("t");
   if (!tree) {
      cout<< "tree \"t\" was not found in file " << ifname <<endl;
      return 0;
   }
   tree->SetMarkerStyle(7);
   tree->SetMarkerColor(2);

   // connect with buffers
   cout<< "// connect with buffers" <<endl;
   // OscEvent* oscEvent = 0;
   OscEvent* oscEvent = new OscEvent;
   tree->SetBranchAddress("oscEvent",&oscEvent);

   cout<< "tree->GetEntries() = " << tree->GetEntries() <<endl;

   // if (gdebug) {
   //    new TCanvas;
   //    tree->Draw("-y:x","ch==2&&evt==1");
   // }

   // output (fit results) tree
   TFile* ofile = TFile::Open(Form("%s-out.root",ifname),"recreate");
   TTree* otree = new TTree("ft", "Fit result tree");
   //restree::book(otree);
   OscFit* oscFit = new OscFit;
   otree->Branch("oscFit", "OscFit", &oscFit);

   // the number of channels in the data
   tree->GetEntry(entry_first);
   Int_t nchannels = oscEvent->oscChannels->GetEntries();

   Int_t isigmin = 0;
   Int_t isigmax = 0;
   Int_t ibkgmin = 0;
   Int_t ibkgmax = 0;

   //return;

   // prepare template histograms for background
   TH1F *hbkg[4];
   Int_t nbins = gEnv->GetValue("Hist.Binning.1D.x", 0);
   gEnv->SetValue("Hist.Binning.1D.x", 400);
   for (int ich=0; ich<nchannels; ++ich) {
      tree->Draw("-y", Form("ch==%d &&x>%f&&x<%f",((OscChannel*) oscEvent->oscChannels->At(ich))->ch,bkgmin,bkgmax), "goff");
      //new TCanvas;
      //tree->Draw("-y", Form("ch==%d &&x>%f&&x<%f",ch[ich],bkgmin,bkgmax), "");
      hbkg[ich] = (TH1F*) tree->GetHistogram()->Clone(Form("hbkg%d",ich+1));
   }

   //-- fit function
   TF1* fitpulse[4];

   Double_t par0[10];
   Double_t par4[10];

   par0[0] = 100.;                  // ampl of the pulse
   par0[1] = (sigmin+sigmax)/2.;    // approx. position
   par0[2] = 50.;                    // tau
   par0[3] = 50.;                    // sigma
   for (int i=0; i<10; ++i) par4[i] = par0[i];

   //-- Int_t npar_sig = 4;
   Int_t npar_sig = 3;
   for (int ich=0; ich<nchannels; ++ich) {
      fitpulse[ich] = new TF1(Form("fitpulse%d",ich),fpulse3_erfc,sigmin,sigmax,npar_sig);
      fitpulse[ich]->SetParameters(par4);
      fitpulse[ich]->SetParName(0,"A");
      fitpulse[ich]->SetParName(1,"mean");
      fitpulse[ich]->SetParName(2,"tau");
      // fitpulse[ich]->SetParName(3,"sigma");

      fitpulse[ich]->SetParLimits(0, 0., 10000.);
      fitpulse[ich]->SetParLimits(1, sigmin, sigmax);
      fitpulse[ich]->SetParLimits(2, 10., 100.);
      //-- fitpulse[ich]->SetParLimits(2, 0., 10.);
      // fitpulse[ich]->SetParLimits(3, 0., 10.);
   }

   if (entry_last < 0) entry_last = tree->GetEntries() - 1;

   for (int jentry=entry_first; jentry<=entry_last; ++jentry)
   {
      tree->GetEntry(jentry);
      cout<< "\n---> jentry = " << jentry << " evt = " << oscEvent->evt <<endl;

      oscFit->clear();
      oscFit->evt = oscEvent->evt;

      cout<< ".. loop over osc channels. oscEvent->oscChannels = " << oscEvent->oscChannels->GetEntries() <<endl;
      for (int ich=0; ich<nchannels; ++ich)
      {
         /// cout<< "---> ich = " << ich << " ch[ich] = " << ch[ich] <<endl;
         //cout<< "oscEvent->oscChannels->GetEntries() = " << oscEvent->oscChannels->GetEntries() <<endl;
         //-- const OscChannel* oscChannel = oscEvent->oscChannel(ich);
         const OscChannel* oscChannel = (OscChannel*) oscEvent->oscChannels->At(ich);

         cout<< "oscChannel->ch = " << oscChannel->ch <<endl;

         // find the level of flat background

         // look through the data

         //. // create a histogram and fit gaussian
         //. TH1F *h = hbkg[ich];
         //. h->Reset();
         Int_t nthres = 0;                // the number of channels which exceed threshold
         Float_t a = 0;
         for (int i=0; i<1024; ++i) {
            if (ibkgmin == 0 && i > 0 && oscChannel->x[i-1] < bkgmin && oscChannel->x[i] >= bkgmin) ibkgmin = i;
            if (ibkgmax == 0 && i > 0 && oscChannel->x[i-1] < bkgmax && oscChannel->x[i] >= bkgmax) ibkgmax = i;
            if (isigmin == 0 && i > 0 && oscChannel->x[i-1] < sigmin && oscChannel->x[i] >= sigmin) isigmin = i;
            if (isigmax == 0 && i > 0 && oscChannel->x[i-1] < sigmax && oscChannel->x[i] >= sigmax) isigmax = i;
            // maximum in signal window NB: (-1)*y
            if (oscChannel->x[i] > sigmin && oscChannel->x[i] < sigmax) if (-1.*oscChannel->y[i] > thres) ++nthres;
            //. if (oscChannel->x[i] >= bkgmin && oscChannel->x[i] <= bkgmax) h->Fill(-oscChannel->y[ic]);
            if (oscChannel->x[i] > sigmin && oscChannel->x[i] < sigmax) if (-1.*oscChannel->y[i] > a) a = -1.*oscChannel->y[i];
         }

         // apply const threshold
         if (nthres < nthres_min) {
            cout<< "skip event evt " << oscEvent->evt << " channel " << oscChannel->ch << " the number of points which exceed threshold " << thres << " is too small: " << nthres <<endl;
            continue;
         }

         // skip events with saturation
         if (a > ysaturation) {
            cout<< "skip saturated event evt " << oscEvent->evt << " channel " << oscChannel->ch <<endl;
            continue;
         }

         // find mean and RMS of the flat background
         //Double_t bkgmean1 = -TMath::Mean(ibkgmax - ibkgmin + 1, oscChannel->y);
         //Double_t bkgRMS1 = -TMath::RMS(ibkgmax - ibkgmin + 1, oscChannel->y);
         /// cout<< "ibkgmin = " << ibkgmin << " ibkgmax = " << ibkgmax << " isigmin = " << isigmin << " isigmax = " << isigmax <<endl;
         /// cout<< "bkgmean1 = " << bkgmean1 << " bkgRMS1 = " << bkgRMS1 <<endl;

         // create a histogram and fit gaussian
         TH1F *h = hbkg[ich];
         h->Reset();
         for (int ic=1; ic<=h->GetNbinsX(); ++ic) h->Fill(-oscChannel->y[ic]);
         // cout<< "--> fit background" <<endl;
         // //----------

         // tree->Draw("-y", Form("ch==%d &&evt==%d &&x>%f&&x<%f",ch[ich],oscEvent->evt,bkgmin,bkgmax), "goff");
         // //new TCanvas;
         // //tree->Draw("-y", Form("ch==%d &&x>%f&&x<%f",ch[ich],bkgmin,bkgmax), "");
         // h = (TH1F*) tree->GetHistogram()->Clone(Form("hbkg%d",ch[ich]));

         Double_t hmean = h->GetMean();
         Double_t hRMS = h->GetRMS();
         Double_t h_xmin = hmean - 1.*hRMS;
         if (h_xmin < h->GetXaxis()->GetXmax()) h_xmin = h->GetXaxis()->GetXmax();
         Double_t h_xmax = hmean + 1.*hRMS;
         if (h_xmax > h->GetXaxis()->GetXmax()) h_xmax = h->GetXaxis()->GetXmax();
         if (gdebug) {
            new TCanvas;
            h->Draw();
            h->Fit("gaus","RL","", h_xmin,h_xmax);
         }
         else h->Fit("gaus","LQ0","goff", h_xmin,h_xmax);     // no graphics, no printout
         // h->Fit("gaus","L0","goff");            // no graphics, printout
         Double_t meanbkg = h->GetFunction("gaus")->GetParameter("Mean");
         Double_t sigmabkg = h->GetFunction("gaus")->GetParameter("Sigma");
         cout<< "meanbkg = " << meanbkg << " sigmabkg = " << sigmabkg <<endl;
         // // if (TMath::Abs(sigmabkg) > 40) {
         // //    TH1F* hbkg_clone = (TH1F*) h->Clone(Form("hbkg_ch%d_evt_%d",ch[ich],oscEvent->evt));
         // //    //new TCanvas;
         // //    //hbkg_clone->Draw();
         // // }

         // create a graph

         Double_t x[1000], y[1000], ex[1000], ey[1000];
         Int_t np = 0;
         Double_t integral = 0;
         for (int i=isigmin; i<=isigmax; ++i)
         {
            integral += 0.5*((-1.)*oscChannel->y[i] + (-1.)*oscChannel->y[i+1] - 2.*meanbkg)*(oscChannel->x[i+1]-oscChannel->x[i]);

            x[np] = oscChannel->x[i];
            //y[np] = -1.*oscChannel->y[i];
            y[np] = -1.*oscChannel->y[i] - meanbkg;
            ex[np] = 0;
            //ey[np] = 1.;
            ey[np] = 1.*TMath::Abs(sigmabkg);
            ++np;
         }
         // normalize to pC by division to R = 50 Hom
         integral /= 50.;
         cout<< ".. integral = " << integral <<endl;

         TGraphErrors *gr = new TGraphErrors(np, x, y, ex, ey);
         gr->SetNameTitle(Form("grsig_evt_%d_ch_%d",oscEvent->evt,oscChannel->ch), Form("grsig_evt_%d_ch_%d",oscEvent->evt,oscChannel->ch));
         gr->SetMarkerStyle(20);
         gr->SetMarkerColor(46);
         gr->SetLineColor(46);
         // TGraph *gr = new TGraph(np, x, y);
         // gr->SetNameTitle(Form("grsig%d",ch[ich]), Form("grsig%d",ch[ich]));
         // gr->SetMarkerStyle(20);
         // gr->SetMarkerColor(46);
         // gr->SetLineColor(46);

         // fit in range
         /// cout<< "... fit for channel " << ch[ich] <<endl;
         /////// new TCanvas;
         /// cout<< "gPad->GetName() = " << gPad->GetName() <<endl;
         /// gr->Draw("awp");
         // TFitResultPtr rfit = gr->Fit(fitpulse[ich], "RSV", "", sigmin, sigmax);
         cout<< "--> fit data" <<endl;
         if (gdebug) {
            new TCanvas;
            gr->Draw("awp");
         }
         TFitResultPtr rfit = gr->Fit(fitpulse[ich], "RS", "", sigmin, sigmax);
         /// cout<< "rfit->Chi2() = " << rfit->Chi2() <<endl;

         // fit the leading edge
         cout<< "// fit the leading edge" <<endl;
         Double_t pos_max = fitpulse[ich]->GetX(fitpulse[ich]->GetMaximum());
         //const double *params_fit = rfit->GetParams();
         //-- Double_t tau = params_fit[2];
         // Double_t sigma = params_fit[3];

         TF1* fun_est = gr->GetFunction(fitpulse[ich]->GetName());
         Double_t tau = fun_est->GetParameter("tau");
         //if (tau < 0.5) tau = 0.5;
         //if (tau > 1.0) tau = 1.0;

         Double_t edge_min = pos_max - 5.*tau;
         Double_t edge_max = pos_max + 1.*tau;
         if (edge_min < sigmin) edge_min = sigmin;
         if (edge_max > sigmax) edge_max = sigmax;

         cout<< "---------> pos_max = " << pos_max << " edge_min = " << edge_min << " edge_max = " << edge_max <<endl;
         TFitResultPtr refit = gr->Fit(fitpulse[ich], "RS", "", edge_min, edge_max);

         // find time of half of maximum
         // fun->GetX(fun->GetMaximum());
         //const double *params_refit = refit->GetParams();
         TF1* fun_refit = gr->GetFunction(fitpulse[ich]->GetName());
         Double_t ymax = fun_refit->GetMaximum();
         Double_t xmax = fun_refit->GetX(ymax);
         Double_t xHalfMax = fun_refit->GetX(ymax/2., sigmin, xmax);
         cout<< "==> xHalfMax = " << xHalfMax <<endl;
         // assign fit result tree
         Int_t ichannel = oscChannel->ch - 1;                     // to count array element from 0
         oscFit->adc[ichannel] = integral;
         oscFit->t[ichannel] = xHalfMax;
         oscFit->t0[ichannel] = fun_refit->GetParameter("mean");
         oscFit->tau[ichannel] = fun_refit->GetParameter("tau");
         oscFit->v[ichannel] = ymax;
         oscFit->bkg[ichannel] = meanbkg;
         oscFit->sigma[ichannel] = sigmabkg;
         oscFit->chi2[ichannel] = (fun_refit->GetNDF() > 0)? fun_refit->GetChisquare() / fun_refit->GetNDF(): 0;
         
         if (fun_refit->GetNDF() < 10 || fun_refit->GetChisquare() / fun_refit->GetNDF() > 100.) fitpulse[ich]->SetParameters(par0);
      }
      otree->Fill();
   }

   ofile->Write();
   gEnv->SetValue("Hist.Binning.1D.x", nbins);

   return oscEvent;
}

TTree* pint(TTree *tree, const char *ofname
      , Double_t bkgmin=0, Double_t bkgmax=15, Double_t sigmin=10, Double_t sigmax=200
      , Int_t entry_first=0, Int_t entry_last=-1
      , bool setdebug=false, bool setgdebug=false
      )
{
   if (setdebug) debug = true;
   if (setgdebug) gdebug = true;

   //-- Double_t thres = 20.;
   Double_t thres = 10.;
   Int_t nthres_min = 5;                  // data are required to contain at least nthres_min points over the thres value
   Double_t ysaturation = 499;

   tree->SetMarkerStyle(7);
   tree->SetMarkerColor(2);

   // connect with buffers
   cout<< "// connect with buffers" <<endl;
   // OscEvent* oscEvent = 0;
   OscEvent* oscEvent = new OscEvent;
   tree->SetBranchAddress("oscEvent",&oscEvent);

   cout<< "tree->GetEntries() = " << tree->GetEntries() <<endl;

   // output (fit results) tree
   //-- TFile* ofile = TFile::Open(Form("%s-out.root",ifname),"recreate");
   TFile* ofile = TFile::Open(ofname,"recreate");
   TTree* otree = new TTree("ft", "Fit result tree");
   //restree::book(otree);
   OscFit* oscFit = new OscFit;
   otree->Branch("oscFit", "OscFit", &oscFit);

   // the number of channels in the data
   tree->GetEntry(entry_first);
   Int_t nchannels = oscEvent->oscChannels->GetEntries();
   // the maximum of the x value
   //if (sigmax < sigmin) sigmin = oscEvent->oscChannel(0)->x[1023];
   //-- if (sigmax < sigmin) sigmin = ((OscChannel*) oscEvent->oscChannels->At(0))->x[1023];
   
   // histogram for background fit
   //-- TH1F *hbkg = new TH1F("hbkg", "hbkg", 400, -20, 20);
   TH1F *hbkg = new TH1F("hbkg", "hbkg", 400, -100, 100);

   if (entry_last < 0) entry_last = tree->GetEntries() - 1;

   for (int jentry=entry_first; jentry<=entry_last; ++jentry)
   {
      tree->GetEntry(jentry);
      cout<< "\n---> jentry = " << jentry << " evt = " << oscEvent->evt <<endl;

      oscFit->clear();
      oscFit->evt = oscEvent->evt;

      cout<< ".. loop over osc channels. oscEvent->oscChannels = " << oscEvent->oscChannels->GetEntries() <<endl;
      for (int ich=0; ich<nchannels; ++ich)
      {
         //-- const OscChannel* oscChannel = oscEvent->oscChannel(ich);
         const OscChannel* oscChannel = (OscChannel*) oscEvent->oscChannels->At(ich);

         cout<< "oscChannel->ch = " << oscChannel->ch <<endl;

         // find the level of flat background

         Double_t x[1024], y[1024], ex[1024], ey[1024];
         Int_t np = 0;

         Double_t ybkg[1024];
         Int_t np_bkg = 0;
         for (int i=0; i<1024; ++i) {
            if (oscChannel->x[i] < bkgmin) continue;
            if (oscChannel->x[i] > bkgmax) break;
            ybkg[np_bkg] = -1.*oscChannel->y[i];
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
            new TCanvas;
            hbkg->DrawCopy();
         }

         Int_t nthres = 0;                // the number of channels which exceed threshold
         Float_t a = 0;
         np = 0;
         for (int i=0; i<1024; ++i)
         {
            if (oscChannel->x[i] < sigmin) continue;
            if (oscChannel->x[i] > sigmax) break;

            // maximum in signal window NB: (-1)*y
            if (oscChannel->x[i] > sigmin && oscChannel->x[i] < sigmax) if (-1.*oscChannel->y[i] > thres) ++nthres;
            //. if (oscChannel->x[i] >= bkgmin && oscChannel->x[i] <= bkgmax) h->Fill(-oscChannel->y[ic]);
            if (oscChannel->x[i] > sigmin && oscChannel->x[i] < sigmax) if (-1.*oscChannel->y[i] > a) a = -1.*oscChannel->y[i];

            // fill the graph
            x[np] = oscChannel->x[i];
            y[np] = -1.*oscChannel->y[i] - bkg_mean;
            ex[np] = 0;
            ey[np] = bkg_sigma;
            ++np;
         }

         TGraphErrors *gr = new TGraphErrors(np, x, y, ex, ey);
         gr->SetNameTitle(Form("grsig_evt_%d_ch_%d",oscEvent->evt,oscChannel->ch), Form("grsig_evt_%d_ch_%d",oscEvent->evt,oscChannel->ch));
         gr->SetMarkerStyle(20);
         gr->SetMarkerColor(46);
         gr->SetLineColor(46);
         if (gdebug) {
            new TCanvas;
            gr->Draw("ap");
         }

         Double_t integral = 0;
         for (int i=0; i<np-1; ++i) integral += 0.5*(gr->GetY()[i] + gr->GetY()[i+1]) * (gr->GetX()[i+1]-gr->GetX()[i]);

         // normalize to pC by division to R = 50 Hom
         integral /= 50.;
         cout<< ".. integral = " << integral << " pC" <<endl;

         // apply const threshold
         if (nthres < nthres_min) {
            cout<< "skip event evt " << oscEvent->evt << " channel " << oscChannel->ch << " the number of points which exceed threshold " << thres << " is too small: " << nthres <<endl;
            continue;
         }

         // skip events with saturation
         if (a > ysaturation) {
            cout<< "skip saturated event evt " << oscEvent->evt << " channel " << oscChannel->ch <<endl;
            continue;
         }

         // assign fit result tree
         Int_t ichannel = oscChannel->ch - 1;                     // to count array element from 0
         oscFit->adc[ichannel] = integral;
         // oscFit->t[ichannel] = xHalfMax;
         // oscFit->t0[ichannel] = fun_refit->GetParameter("mean");
         // oscFit->tau[ichannel] = fun_refit->GetParameter("tau");
         // oscFit->v[ichannel] = ymax;
         oscFit->bkg[ichannel] = bkg_mean;
         oscFit->sigma[ichannel] = bkg_sigma;
         // oscFit->chi2[ichannel] = (fun_refit->GetNDF() > 0)? fun_refit->GetChisquare() / fun_refit->GetNDF(): 0;
      }
      otree->Fill();
   }

   ofile->Write();
   return otree;
}

TTree* pint(const char *ifname="Pilas_ham4x4.xml.root"
      , Double_t bkgmin=0, Double_t bkgmax=15, Double_t sigmin=10, Double_t sigmax=200
      , Int_t entry_first=0, Int_t entry_last=-1
      , bool setdebug=false, bool setgdebug=false
      )
{
   TFile* ifile = TFile::Open(ifname);
   if (!ifile) {
      cout<< "File not found: " << ifname <<endl;
      return 0;
   }
   cout<< "Processing file " << ifname <<endl;

   TTree* tree = (TTree*) ifile->Get("t");
   if (!tree) {
      cout<< "tree \"t\" was not found in file " << ifname <<endl;
      return 0;
   }
   return pint(tree,Form("%s-ft.root",ifname), bkgmin,bkgmax, sigmin,sigmax, entry_first,entry_last, setdebug,setgdebug);
}
