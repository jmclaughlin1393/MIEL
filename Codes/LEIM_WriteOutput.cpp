Int_t WriteOutput(TString options,TH1F *h){
  TString fname = "../rootfiles/Spectrum_"+options+".root";
  TFile *fout = TFile::Open(fname,"RECREATE");
  fout->cd();
  h->Write();
  fout->Close();
  cout << " - Wrote spectrum to file " << fname << endl;
  return 0;
}

Int_t WriteOutput(TString options,TGraphAsymmErrors *h){
  TString fname = "../rootfiles/Spectrum_"+options+".root";
  TFile *fout = TFile::Open(fname,"RECREATE");
  fout->cd();
  h->Write();
  fout->Close();
  cout << " - Wrote spectrum to file " << fname << endl;
  return 0;
}
