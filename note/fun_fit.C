#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>

#include <iostream>
#include <cmath>

using std::cout;     using std::endl;

bool debug = 0;

Double_t fpulse3_erfc(Double_t *xx, Double_t *par)
{
   Double_t fpulse = 0;

   Double_t a = par[0];
   Double_t x0 = par[1];
   Double_t tau1 = par[2];
   Double_t tau2 = par[3];
   Double_t sigma = tau1;

   Double_t x = *xx - x0;

   if (debug) cout<< "a " << a << " x0 " << x0 << " tau1 " << tau1 << " tau2 " << tau2 << " sigma " << sigma << " x " << x <<endl;

   Double_t arg;

   return fpulse;
}
