void erferfc()
{
   gStyle->SetPadGridX(0);
   gStyle->SetPadGridY(0);

   // to use latex scale factor 0.3
   gStyle->SetLabelSize(.075,"xyz");;

   TF1* erf = new TF1("erf", "TMath::Erf(x)", -3,3); erf->SetTitle("");
   TF1* erfc = new TF1("erfc", "TMath::Erfc(x)", -3,3); erfc->SetTitle("");

   new TCanvas;
   erf->Draw();
   pic("erf");

   new TCanvas;
   erfc->Draw();
   pic("erfc");
}
