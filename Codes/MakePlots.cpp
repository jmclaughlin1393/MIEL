void MakePlots_sameSiPM(TString sipm_type){
  gROOT->SetStyle("Modern");
  gROOT->ForceStyle();
  gStyle->SetPalette(87);
  gStyle->SetNumberContours(99);
  gStyle->SetOptStat(0);
  TString short_name = "spectrum_"+sipm_type+"S";
  TString short_filename = "../rootfiles/Spectrum_"+sipm_type+"S.root";
  TFile *short_file = new TFile(short_filename);
  TString medium_name = "spectrum_"+sipm_type+"M";
  TString medium_filename = "../rootfiles/Spectrum_"+sipm_type+"M.root";
  TFile *medium_file = new TFile(medium_filename);
  TString long_name = "spectrum_"+sipm_type+"L";
  TString long_filename = "../rootfiles/Spectrum_"+sipm_type+"L.root";
  TFile *long_file = new TFile(long_filename);

  //TH1F *h_short = (TH1F*)short_file->Get(short_name);
  //TH1F *h_medium = (TH1F*)medium_file->Get(medium_name);
  //TH1F *h_long = (TH1F*)long_file->Get(long_name);
  TGraphAsymmErrors *h_short = (TGraphAsymmErrors*)short_file->Get(short_name);
  TGraphAsymmErrors *h_medium = (TGraphAsymmErrors*)medium_file->Get(medium_name);
  TGraphAsymmErrors *h_long = (TGraphAsymmErrors*)long_file->Get(long_name);
  double sum_short,sum_short_err_hi,sum_short_err_lo,sum_medium,sum_medium_err_hi,sum_medium_err_lo,sum_long,sum_long_err_hi,sum_long_err_lo;
  double over_voltage[3];
  double over_voltage_err[3] = {0.001,0.001,0.001};
  if (sipm_type == "F"){
    over_voltage[0] = 12.785;
    over_voltage[1] = 12.288;
    over_voltage[2] = 12.100;
  }
  if (sipm_type == "H"){
    over_voltage[0] = 10.98;
    over_voltage[1] = 10.79;
    over_voltage[2] = 10.67;
    for (int i = 0; i < 3; i++) over_voltage_err[i] = 0.01;
  }
  double x,y,xerrh,xerrl,yerrh,yerrl;
  for (int i = 1; i < h_short->GetN(); i++){    
    h_short->GetPoint(i,x,y);
    if (x > 1020 || x < 450) continue;
    sum_short_err_hi += TMath::Power(h_short->GetErrorYhigh(i)*4,2.0);
    sum_short_err_lo += TMath::Power(h_short->GetErrorYlow(i)*4,2.0);
    sum_short += y*4;
    h_medium->GetPoint(i,x,y);
    sum_medium_err_hi += TMath::Power(h_medium->GetErrorYhigh(i)*4,2.0);
    sum_medium_err_lo += TMath::Power(h_medium->GetErrorYlow(i)*4,2.0);
    sum_medium += y*4;
    h_long->GetPoint(i,x,y);
    sum_long_err_hi += TMath::Power(h_long->GetErrorYhigh(i)*4,2.0);
    sum_long_err_lo += TMath::Power(h_long->GetErrorYlow(i)*4,2.0);
    sum_long += y*4;
  }
  sum_short_err_hi = TMath::Sqrt(sum_short_err_hi);
  sum_short_err_lo = TMath::Sqrt(sum_short_err_lo);
  sum_medium_err_hi = TMath::Sqrt(sum_medium_err_hi);
  sum_medium_err_lo = TMath::Sqrt(sum_medium_err_lo);
  sum_long_err_hi = TMath::Sqrt(sum_long_err_hi);
  sum_long_err_lo = TMath::Sqrt(sum_long_err_lo);
  
  cout << "Short Exposure Yield:  " << sum_short << " +" << sum_short_err_hi << "/-" << sum_short_err_lo << endl;
  cout << "Medium Exposure Yield: " << sum_medium << " +" << sum_medium_err_hi << "/-" << sum_medium_err_lo << endl;
  cout << "Long Exposure Yield:   " << sum_long << " +" << sum_long_err_hi << "/-" << sum_long_err_lo << endl;

  TCanvas *cPlotAll = new TCanvas("cPlotAll","cPlotAll",950,700);
  cPlotAll->SetGridx();
  cPlotAll->SetGridy();
  cPlotAll->SetTickx();
  cPlotAll->SetTicky();
  //TH2F *hPlotAll = new TH2F("hPlotAll","",1000,400,1050,1000,min,max*1.1);
  TH2F *hPlotAll = new TH2F("hPlotAll","",1000,450,1020,1000,-2e-9,40e-9);
  hPlotAll->GetXaxis()->SetLabelSize(0.03);
  hPlotAll->GetYaxis()->SetLabelSize(0.03);
  hPlotAll->GetXaxis()->SetTickSize(0.02);
  hPlotAll->GetYaxis()->SetTickSize(0.02);
  if (sipm_type == "H") hPlotAll->SetTitle("HPK VUV4 Dark Noise Emmission Spectrum vs Over Voltage;Wavelength [nm];#gamma/e^{-}/nm");
  if (sipm_type == "F") hPlotAll->SetTitle("FBK VUV-HD3 Dark Noise Emmission Spectrum vs Over Voltage;Wavelength [nm];#gamma/e^{-}/nm");
  hPlotAll->Draw();
  h_short->SetLineWidth(1);
  h_short->SetLineColor(kOrange+7);
  h_short->SetMarkerColor(kOrange+7);
  h_short->Draw("SAME P");
  h_medium->SetLineWidth(1);
  h_medium->SetLineColor(kAzure+6);
  h_medium->SetMarkerColor(kAzure+6);
  h_medium->Draw("SAME P");
  h_long->SetLineWidth(1);
  h_long->SetLineColor(kAzure+3);
  h_long->SetMarkerColor(kAzure+3);
  h_long->Draw("SAME P");
  TLegend *leg = new TLegend(.13,.65,.35,.83);
  if (sipm_type == "F"){
    leg->AddEntry(h_short,"V_{ov} = 12.8 V");
    leg->AddEntry(h_medium,"V_{ov} = 12.3 V");
    leg->AddEntry(h_long,"V_{ov} = 12.1 V");
  } else {
    leg->AddEntry(h_short,"V_{ov} = 11.0 V");
    leg->AddEntry(h_medium,"V_{ov} = 10.8 V");
    leg->AddEntry(h_long,"V_{ov} = 10.7 V");
  }
  leg->Draw("SAME");
  if (sipm_type == "F") cPlotAll->SaveAs("~/Documents/Pictures/Plots/leim/FBKVUV-HD3_Spectra_corrected.pdf");
  if (sipm_type == "H") cPlotAll->SaveAs("~/Documents/Pictures/Plots/leim/HamamatsuVUV4_Spectra_corrected.pdf");
}

void MakePlots_sameExposure(TString exposure){
  gROOT->SetStyle("Modern");
  gROOT->ForceStyle();
  gStyle->SetPalette(87);
  gStyle->SetNumberContours(99);
  gStyle->SetOptStat(0);
  TString HPK_name = "spectrum_H"+exposure;
  TString HPK_filename = "../rootfiles/Spectrum_H"+exposure+".root";
  TFile *HPK_file = new TFile(HPK_filename);
  TString FBK_name = "spectrum_F"+exposure;
  TString FBK_filename = "../rootfiles/Spectrum_F"+exposure+".root";
  TFile *FBK_file = new TFile(FBK_filename);

  TH1F *h_HPK = (TH1F*)HPK_file->Get(HPK_name);
  TH1F *h_FBK = (TH1F*)FBK_file->Get(FBK_name);
  double max = h_HPK->GetMaximum();
  double min = h_HPK->GetMinimum();
  if (h_FBK->GetMaximum() > max) max = h_FBK->GetMaximum();
  if (h_FBK->GetMinimum() < min) min = h_FBK->GetMinimum();

  TCanvas *cPlotAll = new TCanvas("cPlotAll","cPlotAll",950,700);
  cPlotAll->SetGridx();
  cPlotAll->SetGridy();
  cPlotAll->SetTickx();
  cPlotAll->SetTicky();
  //TH2F *hPlotAll = new TH2F("hPlotAll","",1000,400,1020,1000,min,max*1.1);
  TH2F *hPlotAll = new TH2F("hPlotAll","",1000,400,1020,1000,-0.012e-9,0.12e-9);
  hPlotAll->GetXaxis()->SetLabelSize(0.03);
  hPlotAll->GetYaxis()->SetLabelSize(0.03);
  hPlotAll->GetXaxis()->SetTickSize(0.02);
  hPlotAll->GetYaxis()->SetTickSize(0.02);

  if (exposure == "S") hPlotAll->SetTitle("Short Exposure Dark Noise Emmission Spectra from HPK VUV4 and FBK VUV-HD3 SiPMs;Wavelength [nm];Photons/e^{-}/nm");
  if (exposure == "M") hPlotAll->SetTitle("Mid-Length Exposure Dark Noise Emmission Spectra from HPK VUV4 and FBK VUV-HD3 SiPMs;Wavelength [nm];Photons/e^{-}/nm");
  if (exposure == "L") hPlotAll->SetTitle("Long Exposure Dark Noise Emmission Spectra from HPK VUV4 and FBK VUV-HD3 SiPMs;Wavelength [nm];Photons/e^{-}/nm");

  hPlotAll->Draw();
  h_HPK->SetLineWidth(1);
  h_HPK->SetLineColor(kOrange+7);
  h_HPK->SetMarkerColor(kOrange+7);
  h_HPK->Draw("SAME");
  h_FBK->SetLineWidth(1);
  h_FBK->SetLineColor(kAzure+3);
  h_FBK->SetMarkerColor(kAzure+3);
  h_FBK->Draw("SAME");
  TLegend *leg = new TLegend(.13,.67,.35,.81);
  leg->AddEntry(h_HPK,"HPK VUV4");
  leg->AddEntry(h_FBK,"FBK VUV-HD3");
  leg->Draw("SAME");
  if (exposure == "S") cPlotAll->SaveAs("~/Documents/Pictures/Plots/leim/SiPMSpectra_ShortExposure.pdf");
  if (exposure == "M") cPlotAll->SaveAs("~/Documents/Pictures/Plots/leim/SiPMSpectra_MidExposure.pdf");
  if (exposure == "L") cPlotAll->SaveAs("~/Documents/Pictures/Plots/leim/SiPMSpectra_LongExposure.pdf");
}
