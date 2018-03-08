#include <TROOT.h>
#include <TH1.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <TLatex.h>
#include <TLine.h>
#include <iostream>

using std::cout;    using std::endl;

void hgaus()
{
  Double_t area = 100;

  TH1F* hgaus = new TH1F("hgaus", ";MeV;events / 0.5 MeV", 20, -5, 5);
  hgaus->SetFillStyle(3001);
  hgaus->SetFillColor(38);

  hgaus->FillRandom("gaus", 1000);
  new TCanvas;
  hgaus->Draw();
  hgaus->Fit("gaus");
  gPad->SaveAs("hgaus.eps");
  gPad->SaveAs("hgaus.png");

  // picture to show relation between the FWHM and sigma

  TF1* fgaus = new TF1("fgaus", "gaus", -5,5);
  // fgaus->SetTitle("Gaussian y(x) = #frac{1}{#sqrt{2#pi}#sigma}exp(- #frac{x^{2}}{2#sigma^{2}})");
  fgaus->SetTitle("");

  Double_t mean = 0;
  Double_t sigma = 1;
  // Double_t A = area/(TMath::Sqrt(2.*TMath::Pi()) * sigma);     // normalize the area to 100
  Double_t A = 1./(TMath::Sqrt(2.*TMath::Pi()) * sigma);     // normalize the area to 1

  fgaus->SetParameter(0, A);
  fgaus->SetParameter(1, mean);
  fgaus->SetParameter(2, sigma);

  Double_t halfmax_y = A / 2.;
  Double_t halfmax_x1 = fgaus->GetX(halfmax_y, mean, mean - 3.*sigma);
  Double_t halfmax_x2 = fgaus->GetX(halfmax_y, mean, mean + 3.*sigma);

  // FWHM line
  TMarker* marker = new TMarker(halfmax_x2, halfmax_y, 20);
  TLine* line_fwhm = new TLine(halfmax_x1, halfmax_y, halfmax_x2, halfmax_y);
  line_fwhm->SetLineWidth(3);
  line_fwhm->SetLineColor(4);
  TLatex* text_fwhm = new TLatex(0, halfmax_y+0.02*A, "FWHM");
  text_fwhm->SetTextAlign(21);
  TLatex* text_marker = new TLatex(halfmax_x2+0.1*sigma, halfmax_y, "(x_{A/2}, A/2)");
  text_marker->SetTextAlign(11);

  // sigma line
  Double_t sigma_x = sigma;
  Double_t sigma_y = fgaus->Eval(sigma_x);	// A*exp(-1/2) = A/sqrt(e)
  TLine* line_sigma = new TLine(0, sigma_y, sigma_x, sigma_y);
  line_sigma->SetLineWidth(3);
  line_sigma->SetLineColor(4);
  TLatex* text_sigma = new TLatex((0+sigma_x)/2, sigma_y+0.02*A, "#sigma");
  text_sigma->SetTextAlign(21);

  new TCanvas;
  fgaus->Draw();
  line_fwhm->Draw();
  marker->Draw();
  text_fwhm->Draw();
  text_marker->Draw();
  line_sigma->Draw();
  text_sigma->Draw();
  gPad->SaveAs("fgaus.eps");
  gPad->SaveAs("fgaus.png");

  // integral over gaussian

  cout<< "\nIntegral function" <<endl;

  TF1* fint_gaus = new TF1("fint_gaus", "[0] * 0.5*(1.+TMath::Erf((x-[1])/(TMath::Sqrt(2)*[2])))", -5,5);
  // fint_gaus->SetParameter(0, area);     // normalize the area to 100
  fint_gaus->SetParameter(0, 1);     // normalize the area to 1
  fint_gaus->SetParameter(1, mean);
  fint_gaus->SetParameter(2, sigma);
  // fint_gaus->SetTitle("Integral over Gaussian: #frac{1}{2}(1+erf(#frac{x}{#sqrt{2}#sigma}))");
  fint_gaus->SetTitle("");
  fint_gaus->GetYaxis()->SetTitle("value of the integral");
  new TCanvas;
  fint_gaus->Draw();
  Double_t offset = 0.3;
  marker->DrawMarker(halfmax_x1, fint_gaus->Eval(halfmax_x1));
  text_marker->DrawText(offset+halfmax_x1, fint_gaus->Eval(halfmax_x1), Form("(-0.5*FWHM, %0.2f)", fint_gaus->Eval(halfmax_x1)));
  marker->DrawMarker(halfmax_x2, fint_gaus->Eval(halfmax_x2));
  text_marker->DrawText(offset+halfmax_x2, fint_gaus->Eval(halfmax_x2), Form("(0.5*FWHM, %0.2f)", fint_gaus->Eval(halfmax_x2)));
  gPad->SaveAs("fint_gaus.eps");
  gPad->SaveAs("fint_gaus.png");
}
