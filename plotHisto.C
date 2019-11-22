{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //
  TFile f("Spectrum.root");
  
  
  TCanvas* c3 = new TCanvas("c1", "  ");
  c3->SetLogy(1);
  c3>cd();
  c3->Update();
  
  TH1D* hist1 = (TH1D*)f.Get("1");
  hist1->SetDirectory(0);
  hist1->Draw("HIST");    
}  
