#include <iostream>
#include <fstream>
using namespace std;

std::vector<float> split(TString line,char delimiter=','){
  TString value = "";
  std::vector<float> vec;
  for (int i = 0; i < line.Length(); i++){
    if (line[i] == ' ') continue;
    if (line[i] == delimiter){
      vec.push_back(value.Atof());
      value = "";
      continue;
    }
    value += line[i];
  }
  vec.push_back(value.Atof());
  return vec;
}

std::vector<std::vector<double> > ReadData(TString option){
  ifstream data;
  if (option=="H") data.open("../Data/newHPKVUV4_IV_curve.txt");
  else if (option=="F") data.open("../Data/newFBKVUV-HD3_IV_curve.txt");
  else {cout << "Invalid choice of SiPM data." << endl; exit(1);}
  if (!data.is_open()){cout << "Failed to open data file" << endl; exit(1);} 
  double V,dV,I,dI;
  std::vector<double> row;
  std::vector<std::vector<double> > data_vec;
  bool acceptInput = !data.eof();
  while (acceptInput) {
    data >> V >> dV >> I >> dI;
    row.push_back(V);
    row.push_back(dV);
    row.push_back(I);
    row.push_back(dI);
    data_vec.push_back(row);
    row.clear();
    acceptInput = !data.eof();
  }
  data.close();
  return data_vec;
}

void DrawIVCurve(TString option){
  gROOT->SetStyle("Modern");
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  TFile *fout = new TFile();
  if (option == "H") fout = TFile::Open("../Rootfiles/HPKVUV4_IVCurve.root","RECREATE");
  if (option == "F") fout = TFile::Open("../Rootfiles/FBKVUV-HD3_IVCurve.root","RECREATE");
  std::vector<std::vector<double> > data = ReadData(option);
  TGraphErrors *gData = new TGraphErrors(data.size());
  double V,dV,I,dI;
  for (int i = 0; i < data.size(); i++){
    V = data[i][0];
    dV = data[i][1];
    I = 1000*data[i][2];
    dI = 1000*data[i][3];
    gData->SetPoint(i,V,I);
    gData->SetPointError(i,dV,dI);
  }
  gData->SetName("IVCurve");
  double max,derivative,dx,Vbd;
  if (option == "H") max = 60;
  if (option == "F") max = 40;
  derivative = 0;
  dx = 0.01;
  for (double x = 30; x <= max; x += dx){
    double test_derivative = (TMath::Log(gData->Eval(x+dx)) - TMath::Log(gData->Eval(x)))/dx;
    cout << x << " " << test_derivative << endl;
    if (derivative < test_derivative){
      derivative = test_derivative;
      Vbd = x;
    }
  }
  cout << "Breakdown voltage = " << Vbd << endl;
  if (option == "H") gData->SetTitle("Hamamatsu VUV4 MPPC Current vs. Reverse Bias Voltage;Reverse Bias Voltage [V];Current [mA]");
  if (option == "F") gData->SetTitle("FBK VUV-HD3 SiPM Current vs. Reverse Bias Voltage;Reverse Bias Voltage [V];Current [mA]");
  gData->SetMarkerStyle(8);
  gData->SetMarkerSize(0.7);
  gData->SetMarkerColor(kBlue+3);
  gData->SetLineColor(kOrange+7);
  gData->SetMinimum(1.e-7);
  gData->SetMaximum(2);
  gData->GetXaxis()->SetLabelSize(0.03);
  gData->GetXaxis()->SetTickSize(0.02);
  gData->GetXaxis()->SetRange(10,96);
  gData->GetYaxis()->SetLabelSize(0.03);
  gData->GetYaxis()->SetTickSize(0.02);
  TCanvas *c = new TCanvas("c","c",900,600);
  c->SetLogy();
  c->SetGridx();
  c->SetGridy();
  c->cd();
  gData->Draw("APL");
  if (option == "H") c->SaveAs("../Images/HPKVUV4_IVCurve.pdf");
  if (option == "F") c->SaveAs("../Images/FBKVUV-HD3_IVCurve.pdf");
  fout->cd();
  gData->Write();
  fout->Close();
}
