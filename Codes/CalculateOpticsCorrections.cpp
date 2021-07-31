#include <iostream>
#include <fstream>
#include <vector>
#include "TGraph.h"
#include "TSystem.h"
#include "TString.h"
#include "TFile.h"
using namespace std;
double NA = 0.45;
double XD = 0.8;
std::vector<float> wavelength,n_value,k_value;
TGraph *G_Si = new TGraph();
TGraph *G_SiO2 = new TGraph();
TGraph *G_Atm = new TGraph();
TGraph *G_Abs = new TGraph();

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

void MakeTGraphs(){
  cout << "Makeing refractive index TGraphs..." << endl;
  TGraph *G_Si_n = new TGraph();
  TGraph *G_Si_k = new TGraph();
  TGraph *G_SiO2_n = new TGraph();
  TGraph *G_SiO2_k = new TGraph();
  TGraph *G_Atm_n = new TGraph();
  TGraph *G_abs = new TGraph();
  TFile *out = new TFile("../rootfiles/RefractiveIndex_vs_Wavelength.root","RECREATE");
  bool acceptInput;
  TString line;
  std::vector<float> values;
  ifstream Si_data,SiO2_data,Atm_data,absorption_lengths;
  
  Si_data.open("../data/csv_files/RefractiveIndexINFO_Si.csv");
  SiO2_data.open("../data/csv_files/RefractiveIndexINFO_SiO2.csv");
  Atm_data.open("../data/csv_files/RefractiveIndexINFO_Atm.csv");
  absorption_lengths.open("../data/absorption_lengths.txt");

  cout << "... Silicon" << endl;
  acceptInput = !Si_data.eof();
  while (acceptInput){
    Si_data >> line;
    values = split(line);
    if (values[0] < 0.25){
      acceptInput = !Si_data.eof();
      continue;
    }
    wavelength.push_back(1000.*values[0]);
    n_value.push_back(values[1]);
    k_value.push_back(values[2]);
    acceptInput = !Si_data.eof();
  }
  G_Si_n->Set(wavelength.size());
  G_Si_n->SetNameTitle("G_Si_n","Silicon Refractive Index (real part);Wavelength [nm];Re(n)");  
  G_Si_k->Set(wavelength.size());
  G_Si_k->SetNameTitle("G_Si_k","Silicon Refractive Index (imaginary part);Wavelength [nm];Im(n)");
  for (unsigned int i = 0; i < wavelength.size(); i++){
    G_Si_n->SetPoint(i,wavelength[i],n_value[i]);
    G_Si_k->SetPoint(i,wavelength[i],k_value[i]);
  }
  wavelength.clear();
  n_value.clear();
  k_value.clear();

  cout << "... Silicon Dioxide" << endl;
  acceptInput = !SiO2_data.eof();
  SiO2_data >> line;
  while (acceptInput){
    SiO2_data >> line;
    values = split(line);
    wavelength.push_back(1000.*values[0]);
    n_value.push_back(values[1]);
    k_value.push_back(values[2]);
    acceptInput = !SiO2_data.eof();
  }
  G_SiO2_n->Set(wavelength.size());
  G_SiO2_n->SetNameTitle("G_SiO2_n","Silicon Dioxide Refractive Index (real part);Wavelength [nm];Re(n)");
  G_SiO2_k->Set(wavelength.size());
  G_SiO2_k->SetNameTitle("G_SiO2_k","Silicon Dioxide Refractive Index (imaginary part);Wavelength [nm];Im(n)");
  for (unsigned int i = 0; i < wavelength.size(); i++){
    G_SiO2_n->SetPoint(i,wavelength[i],n_value[i]);
    G_SiO2_k->SetPoint(i,wavelength[i],k_value[i]);
  }
  wavelength.clear();
  n_value.clear();
  k_value.clear();

  cout << "... Atmosphere" << endl;
  acceptInput = !Atm_data.eof();
  Atm_data >> line;
  while (acceptInput){
    Atm_data >> line;
    values = split(line);
    wavelength.push_back(1000.*values[0]);
    n_value.push_back(values[1]);
    acceptInput = !Atm_data.eof();
  }
  G_Atm_n->Set(wavelength.size());
  G_Atm_n->SetNameTitle("G_Atm_n","Atmosphere Refractive Index (real part);Wavelength [nm];Re(n)");
  for (unsigned int i = 0; i < wavelength.size(); i++) G_Atm_n->SetPoint(i,wavelength[i],n_value[i]);
  wavelength.clear();
  n_value.clear();

  cout << "... Absorption length" << endl;
  acceptInput = !absorption_lengths.eof();
  double wav,abs;
  while (acceptInput){
    absorption_lengths >> wav >> abs;
    //values = split(line);
    wavelength.push_back(wav);
    n_value.push_back(abs);
    acceptInput = !absorption_lengths.eof();
  }
  G_abs->Set(wavelength.size());
  G_abs->SetNameTitle("G_abs","Absorption Length in Silicon;Wavelength [nm];Absorption Length [#mum]");
  for (unsigned int i = 0; i < wavelength.size(); i++) G_abs->SetPoint(i,wavelength[i],n_value[i]);
  wavelength.clear();
  n_value.clear();

  cout << "... writing to output file" << endl;
  out->cd();
  G_Si_n->Write();
  G_Si_k->Write();
  G_SiO2_n->Write();
  G_SiO2_k->Write();
  G_Atm_n->Write();
  G_abs->Write();
  out->Close();
  cout << "Done!" << endl;
}

double reflactance(double n1, double n2, double theta){
  double answer,sPol,pPol;
  sPol  = n1*TMath::Cos(theta) - n2*TMath::Sqrt(1 - TMath::Power(n1*TMath::Sin(theta)/n2,2.0));
  sPol /= n1*TMath::Cos(theta) + n2*TMath::Sqrt(1 - TMath::Power(n1*TMath::Sin(theta)/n2,2.0));
  sPol *= sPol;

  pPol  = n1*TMath::Sqrt(1 - TMath::Power(n1*TMath::Sin(theta)/n2,2.0)) - n2*TMath::Cos(theta);
  pPol /= n1*TMath::Sqrt(1 - TMath::Power(n1*TMath::Sin(theta)/n2,2.0)) + n2*TMath::Cos(theta);
  pPol *= pPol;

  answer += 0.5*(sPol + pPol);
  return answer;
}

double f_lambda(double lambda){
  double answer;
  double ref1,ref2,n_Si,n_SiO2,n_Atm,theta_max,theta_prime,dtheta,absorption_length,absorption_coeff;
  n_Si = G_Si->Eval(lambda);
  n_SiO2 = G_SiO2->Eval(lambda);
  n_Atm = G_Atm->Eval(lambda);
  absorption_length = G_Abs->Eval(lambda);
  theta_max = TMath::ASin(NA/n_Si);
  dtheta = theta_max/5000;
  answer = 0;
  for (double theta = 0; theta <= theta_max; theta += dtheta){
    absorption_coeff = TMath::Exp(-XD/(TMath::Cos(theta)*absorption_length));
    ref1 = reflactance(n_Si,n_SiO2,theta);
    theta_prime = TMath::ASin(n_Si*TMath::Sin(theta)/n_SiO2);
    ref2 = reflactance(n_SiO2,n_Atm,theta_prime);
    answer += dtheta*TMath::Sin(theta)*(1 - ref1)*(1 - ref2)*absorption_coeff;
  }
  answer /= 2.;
  return answer;
}

void CalculateOpticsCorrections(TString mode="NIR",TString sipm="VUV4"){
  if (mode=="VIS") NA = 0.4;
  if (sipm=="VUVHD3") XD = 0.145;
  bool makeGraphs = gSystem->AccessPathName("../rootfiles/RefractiveIndex_vs_Wavelength.root");
  if (makeGraphs) MakeTGraphs();
  TFile *f_refIndex = new TFile("../rootfiles/RefractiveIndex_vs_Wavelength.root");
  G_Si = (TGraph*)f_refIndex->Get("G_Si_n");
  G_SiO2 = (TGraph*)f_refIndex->Get("G_SiO2_n");
  G_Atm = (TGraph*)f_refIndex->Get("G_Atm_n");
  G_Abs = (TGraph*)f_refIndex->Get("G_abs");
  TGraph *theta_Si = new TGraph(501);
  TGraph *G_correction = new TGraph(501);
  double lambda,correction,theta;
  for (int i = 0; i <= 500; i++){
    lambda = 400. + 0.002*i*700;
    correction = f_lambda(lambda);
    G_correction->SetPoint(i,lambda,correction);
    theta = TMath::ASin(NA/G_Si->Eval(lambda));
    theta_Si->SetPoint(i,lambda,theta);
  }
  G_correction->SetNameTitle("G_correction","Rescaling Integral vs. Wavelength;Wavelength [nm];Integral, I(#lambda)");
  theta_Si->SetNameTitle("theta_Si","Acceptance Angle #theta_{Si} vs Wavelength;Wavelength [nm];Acceptance Angle #theta_{Si}");
  TString filename = "../rootfiles/CorrectionVsWavelength_";
  filename += mode;
  filename += "mode_";
  filename += sipm;
  filename += ".root";
  TFile *fout = new TFile(filename,"RECREATE");
  fout->cd();
  G_correction->Write();
  theta_Si->Write();
  fout->Close();
}
