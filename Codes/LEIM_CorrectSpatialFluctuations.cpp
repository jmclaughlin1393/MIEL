#include <iostream>
#include <fstream>

using namespace std;
double variance;
double SiPM_Integral(TH2F *h,TString options="HG"){
  double answer = 0;
  double f_right[2] = {326./7.,-993040./7.};
  double f_top[2] = {-7./598.,79960./23.};
  double f_left[2] = {649./5.,-13399.};
  double f_bottom[2] = {-10./589.,136520./589.};
  if (options=="FG"){
    f_right[0] = 605./13.; f_right[1] = -3211900./13.;
    f_top[0] = -18./1045.; f_top[1] = 1294484./209;
    f_left[0] = 606./11.; f_left[1] = -62200./11.;
    f_bottom[0] = 16./1041.; f_bottom[1] = 45670./347.;
  } else if (options=="HL"){
    f_top[0] = 0.0; f_top[1] = 2125.;
    f_bottom[0] = 0.0; f_bottom[1] = 1720.;
    f_left[0] = 1.e+09; f_left[1] = -1165.e+09;
    f_right[0] = 1.e+09; f_right[1] = -1770.e+09;
  } else if (options=="FL"){
    f_top[0] = 0.0; f_top[1] = 2125.;
    f_bottom[0] = 0.0; f_bottom[1] = 1720.;
    f_left[0] = 1.e+09; f_left[1] = -1165.e+09;
    f_right[0] = 1.e+09; f_right[1] = -1770.e+09;
  }
  double y_right,y_top,y_left,y_bottom;
  double startx = 0;
  double starty = 0;
  double endx = 3500;
  double endy = 3500;
  double counter = 0;
  int xBin,yBin;
  if (options.Contains("F")){endx=5500; endy=6500;}  
  for (double x = startx; x < endx; x += 5.0){
    y_right = f_right[0]*x + f_right[1];
    y_top = f_top[0]*x + f_top[1];
    y_left = f_left[0]*x + f_left[1];
    y_bottom = f_bottom[0]*x + f_bottom[1];
    for (double y = starty; y < endy; y++){
      if (y >= y_left || y >= y_top) continue;
      if (y <= y_right || y <= y_bottom) continue;
      xBin = h->GetXaxis()->FindBin(x);
      yBin = h->GetYaxis()->FindBin(y);
      answer += h->GetBinContent(xBin,yBin);
      counter += 1.0;
    }
  }
  variance = answer*0.7 + counter*3.26*3.26;
  return answer;
}

double SiPM_Integral_zoom(TH2F *h,TString options){
  double answer = 0;
  double counter = 0;
  double f_left,f_right;
  if (options=="G"){
    f_left = h->GetXaxis()->GetBinCenter(0);
    f_right = h->GetXaxis()->GetBinCenter(h->GetNbinsX()+1);
  }
  else if (options=="L"){
    f_left = 343.25;
    f_right = 348.2;
  }
  for (int xBin = 1; xBin < h->GetNbinsX(); xBin++){
    if (h->GetXaxis()->GetBinCenter(xBin) <= f_left || h->GetXaxis()->GetBinCenter(xBin) >= f_right) continue;
    for (int yBin = 1; yBin < h->GetNbinsY(); yBin++) {
      answer += h->GetBinContent(xBin,yBin);
      counter += 1.0;
    }
  }
  variance = answer*0.7 + counter*3.26*3.26;
  return answer;
}

double *CorrectSpatialFluctuations(TString options){
  gROOT->SetStyle("Modern");
  gROOT->ForceStyle();
  gStyle->SetPalette(87);
  gStyle->SetNumberContours(99);
  gStyle->SetOptStat(0);
  
  TFile *fHPK_full = new TFile("../Rootfiles/HPKVUV4_SiPM_fullImage_CRR.root");
  TFile *fHPK_zoom = new TFile("../Rootfiles/HPKVUV4_SiPM_zoomed.root");
  TFile *fFBK_full = new TFile("../Rootfiles/FBKVUV-HD3_SiPM_fullImage.root");
  TFile *fFBK_zoom = new TFile("../Rootfiles/FBKVUV-HD3_SiPM_zoomed.root");
  TH2F *hHPK_full = (TH2F*)fHPK_full->Get("HPKVUV4_full");
  TH2F *hHPK_zoom = (TH2F*)fHPK_zoom->Get("HPKVUV4_zoomed");
  TH2F *hFBK_full = (TH2F*)fFBK_full->Get("FBKVUV-HD3_full");
  TH2F *hFBK_zoom = (TH2F*)fFBK_zoom->Get("FBKVUV-HD3_zoomed"); 
  double *values = new double[2];
  
  if (options.Contains("H")){
    double HPK_global = SiPM_Integral(hHPK_full,"HG");
    double HPK_global_var = variance;
    double HPK_local = SiPM_Integral(hHPK_full,"HL");
    double HPK_local_var = 0;
    double HPK_zoom = SiPM_Integral_zoom(hHPK_zoom,"G");
    double HPK_zoom_var = variance;
    double HPK_slit = SiPM_Integral_zoom(hHPK_zoom,"L");
    double HPK_slit_var = 0;
    values[0] = (HPK_zoom/HPK_slit)*(HPK_global/HPK_local);
    values[1] =  TMath::Power(values[0]/HPK_zoom, 2.0)*HPK_zoom_var;
    values[1] += TMath::Power(values[0]/HPK_slit, 2.0)*HPK_slit_var;
    values[1] += TMath::Power(values[0]/HPK_global, 2.0)*HPK_global_var;
    values[1] += TMath::Power(values[0]/HPK_local, 2.0)*HPK_local_var;
    values[1] += TMath::Power(values[0]*(3./140.), 2.0);//assuming I was correct to within 3 pixels along x in isolating the local region
    values[1] += TMath::Power(values[0]*(3./70.), 2.0);//assuming I was correct to within 3 pixels along y in isolating the local region
    values[1] = TMath::Sqrt(values[1]);
    cout << " - Hamamatsu VUV4 Conversion..................." << endl;
    cout << "   ............................................" << endl;
    cout << "   HPK_global = " << HPK_global << endl;
    cout << "   HPK_local = " << HPK_local << endl;
    cout << "   HPK_zoom = " << HPK_zoom << endl;
    cout << "   HPK_slit = " << HPK_slit << endl;
    cout << "   slit->zoom multiplier = " << HPK_zoom/HPK_slit << endl;
    cout << "   zoom->full multiplier = " << HPK_global/HPK_local << endl;
    cout << "   ............................................" << endl;
  }else {    
    double FBK_global = SiPM_Integral(hHPK_full,"FG");
    double FBK_global_var = variance;
    double FBK_local = SiPM_Integral(hHPK_full,"FL");
    double FBK_local_var = 0;
    double FBK_zoom = SiPM_Integral_zoom(hFBK_zoom,"G");
    double FBK_zoom_var = variance;
    double FBK_slit = SiPM_Integral_zoom(hFBK_zoom,"L");
    double FBK_slit_var = 0;
    values[0] = (FBK_zoom/FBK_slit)*(FBK_global/FBK_local);
    values[1] =  TMath::Power(values[0]/FBK_zoom, 2.0)*FBK_zoom_var;
    values[1] += TMath::Power(values[0]/FBK_slit, 2.0)*FBK_slit_var;
    values[1] += TMath::Power(values[0]/FBK_global, 2.0)*FBK_global_var;
    values[1] += TMath::Power(values[0]/FBK_local, 2.0)*FBK_local_var;
    values[1] += TMath::Power(values[0]*(3./140.), 2.0);//assuming I was correct to within 3 pixels along x in isolating the local region
    values[1] += TMath::Power(values[0]*(3./70.), 2.0);//assuming I was correct to within 3 pixels along y in isolating the local region
    values[1] =  TMath::Sqrt(values[1]);
    cout << " - FBK VUV-HD3 Conversion......................" << endl;
    cout << "   ............................................" << endl;
    cout << "   FBK_global = " << FBK_global << endl;
    cout << "   FBK_local = " << FBK_local << endl;
    cout << "   FBK_zoom = " << FBK_zoom << endl;
    cout << "   FBK_slit = " << FBK_slit << endl;
    cout << "   slit->zoom multiplier = " << FBK_zoom/FBK_slit << endl;
    cout << "   zoom->full multiplier = " << FBK_global/FBK_local << endl;
    cout << "   ............................................" << endl;
  }
  return values;
}
