#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>

#include <iostream>
#include <cmath>

using std::cout;     using std::endl;

bool debug = 0;

Double_t fun_conv(Double_t *xx, Double_t *par)
{
   Double_t t = *xx;          // it's t in the paper
   Double_t xpar = par[0];
   Double_t norm = par[1];
   Double_t shift = par[2];
   Double_t tau = par[3];
   Double_t sigma = par[4];

   Double_t x = xpar - shift;

   Double_t pi = std::acos(-1.);
   Double_t gnorm = 1./std::sqrt(2.*pi)/sigma;
   Double_t gauss = gnorm * std::exp(-0.5*std::pow((t-x)/sigma,2));
   Double_t fun_conv = norm*(1.-std::exp(-t/tau))*std::exp(-t/tau) * gauss;
   cout<< "gnorm = " << gnorm << " gauss = " << gauss << " fun_conv = " << fun_conv <<endl;
   cout<< "(1.-std::exp(-t/tau)) = " << (1.-std::exp(-t/tau)) <<endl;
   cout<< "std::exp(-t/tau) = " << std::exp(-t/tau) <<endl;
   return fun_conv;
}

// Double_t fun_deconv(Double_t *xx, Double_t *par)
// {
//    Double_t t = *xx;          // it's t in the paper
//    Double_t xpar = par[0];
//    Double_t norm = par[1];
//    Double_t shift = par[2];
//    Double_t tau = par[3];
//    Double_t sigma = par[4];
// 
//    Double_t pi = std::acos(-1.);
//    Double_t gnorm = 1./std::sqrt(2.*pi)/sigma;
// 
//    Double_t fun_conv = 0;
//    return fun_conv;
// }

Double_t fpulse0(Double_t *xx, Double_t *par)
{
   Double_t fpulse = 0;

   Double_t a = par[0];
   Double_t x0 = par[1];
   Double_t tau = par[2];
   Double_t sigma = par[3];

   Double_t x = *xx - x0;

   if (debug) cout<< "a " << a << " x0 " << x0 << " tau " << tau << " sigma " << sigma << " x " << x <<endl;

   Double_t arg;

   Double_t exp_const1 = 0;
   arg = sigma*sigma/(2.*tau*tau);
   // if (TMath::Abs(arg) < 50.) exp_const1 = TMath::Exp(arg);
   exp_const1 = TMath::Exp(arg);

   Double_t exp_x1 = 0;
   arg = -x/tau;
   // if (TMath::Abs(arg) < 50) exp_x1 = TMath::Exp(arg);
   exp_x1 = TMath::Exp(arg);

   Double_t erf1 = 0.;
   arg = (x - sigma*sigma/tau) / (sigma*TMath::Sqrt(2.));
   // if (TMath::Abs(arg) < 50) erf1 = TMath::Erf(arg);
   // else erf1 = (arg < 0.)? -1.: 1.;
   erf1 = TMath::Erf(arg);

   Double_t I1 = exp_const1*exp_x1*(1. + erf1);

   Double_t exp_const2 = 0;
   arg = 2.*sigma*sigma/(tau*tau);
   // if (TMath::Abs(arg) < 50.) exp_const2 = TMath::Exp(arg);
   exp_const2 = TMath::Exp(arg);

   Double_t exp_x2 = 0;
   arg = -2*x/tau;
   // if (TMath::Abs(arg) < 50) exp_x2 = TMath::Exp(arg);
   exp_x2 = TMath::Exp(arg);

   Double_t erf2 = 0;
   arg = (x - 2.*sigma*sigma/tau) / (sigma*TMath::Sqrt(2.));
   // if (TMath::Abs(arg) < 50) erf2 = TMath::Erf(arg);
   // else erf2 = (arg < 0.)? -1.: 1.;
   erf2 = TMath::Erf(arg);

   Double_t I2 = exp_const2*exp_x2*(1. + erf2);

   fpulse = a*(I1 - I2);
   return fpulse;
}

Double_t fpulse_tau_sigma(Double_t *xx, Double_t *par)
{
   Double_t fpulse = 0;

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

   fpulse = a*(I1 - I2);
   return fpulse;
}

Double_t fpulse(Double_t *xx, Double_t *par)
{
   Double_t fpulse = 0;

   Double_t a = par[0];
   Double_t x0 = par[1];
   Double_t tau = par[2];
   Double_t sigma = tau;

   Double_t x = *xx - x0;

   if (debug) cout<< "a " << a << " x0 " << x0 << " tau " << tau << " x " << x <<endl;

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

   fpulse = a*(I1 - I2);
   return fpulse;
}

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

/*
// oscillator behaviour with fun_test(0.5,5.0, 108,120)
*/
void fun_test(Double_t tau=0.5, Double_t sigma=1, Double_t xmin=90, Double_t xmax=110)
{
   TF1* fun_spulse = new TF1("fun_spulse", fpulse, xmin, xmax, 4);
   Double_t a = 1000;
   Double_t x0 = 100;
   // Double_t tau = 1;
   // Double_t sigma = 1;
   Double_t par[10];
   par[0] = a;
   par[1] = x0;
   par[2] = tau;
   par[3] = sigma;
   fun_spulse->SetParameters(par);

   //new TCanvas;
   fun_spulse->Draw();

   // 3 parameter function
   xmin = -200;
   xmax = 1000;
   TF1* fun_pulse3_erfc = new TF1("fun_pulse3_erfc", fpulse3_erfc, xmin, xmax, 3);
   Double_t par3[10];
   par3[0] = a;
   par3[1] = x0;
   par3[2] = tau;
   fun_pulse3_erfc->SetParameters(par3);
   fun_pulse3_erfc->SetParName(0, "A");
   fun_pulse3_erfc->SetParName(1, "mean");
   fun_pulse3_erfc->SetParName(2, "tau");
   // play with parameters
   fun_pulse3_erfc->SetParameter("A", 100);
   fun_pulse3_erfc->SetParameter("tau", 10);

   new TCanvas;
   fun_pulse3_erfc->Draw();
}
