#include "TGraph.h"
#include "TFile.h"
#include <vector>
TFile *fcorrections_ir = new TFile("../rootfiles/InfraredRangeEfficiencyCurve.root");
TFile *fcorrections_vis = new TFile("../rootfiles/VisibleRangeEfficiencyCurve.root");
TFile *ferrors = new TFile("../rootfiles/EfficiencyCurve_ErrorBands.root");
TFile *fToyMC = new TFile("../rootfiles/MIELEfficiencyError_ToyMC.root");
TFile *fToyMC_systematic = new TFile("../rootfiles/MIELEfficiencyError_ToyMC_Systematic.root");
TFile *fToyMC_lt1percent = new TFile("../rootfiles/MIELEfficiencyError_ToyMC_lt1percent.root");
std::vector<double> CorrectEfficiency(std::vector<double> wavelength, std::vector<double> intensity, int row, bool range_vis){
  TGraph *G_ir = (TGraph*)fcorrections_ir->Get(Form("EfficiencyCurve_row%d",row));
  TGraph *G_vis = (TGraph*)fcorrections_vis->Get(Form("EfficiencyCurve_row%d",row));
  if (range_vis){
    for (unsigned int i = 0; i < intensity.size(); i++){
      intensity[i] /= G_vis->Eval(wavelength[i]);
    }
  } else {
    for	(unsigned int i = 0; i < intensity.size(); i++){
      intensity[i] /= G_ir->Eval(wavelength[i]);
    }
  }
  return intensity;
}


// std::vector<double> GetError(std::vector<double> wavelength, std::vector<double> intensity, std::vector<double> intensity_err, int row, bool range_vis, bool upperErrorBars){
//   TGraph *G_ir = (TGraph*)fcorrections_ir->Get(Form("EfficiencyCurve_row%d",row));
//   TGraph *G_vis = (TGraph*)fcorrections_vis->Get(Form("EfficiencyCurve_row%d",row));
//   TGraph *G_error;
//   if (upperErrorBars) G_error = (TGraph*)ferrors->Get("G_full_errLO");
//   else G_error = (TGraph*)ferrors->Get("G_full_errUP");//Switched because the upper error bar in efficiency drives the lower error bar in intensity
//   double correction,correction_err;
//   if (range_vis){
//     if (!upperErrorBars) return intensity_err;
//     for (unsigned int i = 0; i < intensity_err.size(); i++){
//       correction = G_vis->Eval(wavelength[i]);
//       correction_err = G_error->Eval(wavelength[i])*0.01/correction - 1.0;
//       intensity_err[i] = TMath::Sqrt(TMath::Power(intensity_err[i]/correction,2.0) + TMath::Power(intensity[i]*correction_err/(correction*correction),2.0));
//     }
//   } else {
//     for (unsigned int i = 0; i < intensity_err.size(); i++){
//       correction = G_ir->Eval(wavelength[i]);
//       correction_err = G_error->Eval(wavelength[i])*0.01/correction - 1.0;
//       intensity_err[i] = TMath::Sqrt(TMath::Power(intensity_err[i]/correction,2.0) + TMath::Power(intensity[i]*correction_err/(correction*correction),2.0));
//       if (G_ir->Eval(wavelength[i+1]) <= 0.03) break;
//     }
//   }
//   return intensity_err;
// }

std::vector<double> GetError(std::vector<double> wavelength, std::vector<double> corrected_intensity, std::vector<double> intensity_err, int row, bool range_vis, bool upperErrorBars){
  TH2F *hStdDev = (TH2F*)fToyMC->Get("hStdDev");
  TH2F *hStdDev_systematic = (TH2F*)fToyMC_systematic->Get("hStdDev");
  TH2F *hStdDev_lt1percent = (TH2F*)fToyMC_lt1percent->Get("hStdDev");
  TGraph *G_ir = (TGraph*)fcorrections_ir->Get(Form("EfficiencyCurve_row%d",row));
  TGraph *G_vis = (TGraph*)fcorrections_vis->Get(Form("EfficiencyCurve_row%d",row));
  double correction,truthCounts;
  if (range_vis){
    double wav,counts;
    for (unsigned int i = 0; i < intensity_err.size(); i++){
      wav = wavelength[i];
      counts = corrected_intensity[i];
      if (!upperErrorBars){
	if (wavelength[i] <= 400) wav = 401;
	if (wavelength[i] >= 550) wav = 549;
	if (corrected_intensity[i] < 1) counts = 1;
	if (corrected_intensity[i] > 99) counts = 99;
	intensity_err[i] = hStdDev_systematic->Interpolate(counts,wav);
      }
      else {
	correction = G_vis->Eval(wavelength[i]);
	counts = corrected_intensity[i];
	if (100*correction < 0.1) correction = 0.001 + 1.e-6;
	if (counts > 3000) counts = 3000 - 1;
	if (100*correction > 40) correction = 0.40 - 1.e-6;
	if (counts < 0) counts = 1;
	intensity_err[i] = hStdDev->Interpolate(counts,100*correction);
      }
    }
  }
  else {
    for (unsigned int i = 0; i < intensity_err.size(); i++){
      //if (G_ir->Eval(wavelength[i]) > 0.02) continue;
      correction = G_ir->Eval(wavelength[i]);
      truthCounts = corrected_intensity[i];
      if (100*correction <= 0.1) correction = 0.001 + 1.e-6;
      if (truthCounts >= 3000) truthCounts = 2999;
      if (100*correction >= 40) correction = 0.40 - 1.e-6;
      if (truthCounts <= 0) truthCounts = 1;
      if (100*correction >= 1.0) {
	intensity_err[i] = hStdDev->Interpolate(truthCounts,100*correction);
	//if (row == 199) cout << truthCounts << " " << 100*correction << " " << intensity_err[i] << endl;
      }
      else {intensity_err[i] = hStdDev_lt1percent->Interpolate(truthCounts,100*correction);}
    }
  }
  return intensity_err;
}
