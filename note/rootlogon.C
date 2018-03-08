{
gROOT->Macro(gSystem->ExpandPathName("$(HOME)/macros/rootlogon.C"));
cout<< "*-- Local rootlogon" << endl;

// cout<< "Local rootlogon.C: Defined const char s[] = \"same\"" <<endl;
// const char s[] = "same";

// cout<< "Load FWLite" <<endl;
// gSystem->Load("libFWCoreFWLite.so");
// AutoLibraryLoader::enable();
// gSystem->Load("libDataFormatsFWLite.so");
// gROOT->ProcessLine("namespace edm {typedef edm::Wrapper<vector<float> > Wrapper<vector<float,allocator<float> > >; }");
// gROOT->ProcessLine("namespace edm {typedef edm::Wrapper<vector<double> > Wrapper<vector<double,allocator<double> > >; }");

   //gStyle->SetLabelSize(.05,"xyz");
   //gStyle->SetTitleSize(.05,"xyz");

cout<< "Load macro pulse.tex.C+" <<endl;
gROOT->ProcessLine(".L pulse.tex.C+");
}
