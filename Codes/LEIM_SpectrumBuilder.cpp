#include "LEIM_ReadInData.cpp"
#include "LEIM_CalculateCharge.cpp"
#include "LEIM_CorrectSpatialFluctuations.cpp"
#include "LEIM_RemoveCosmicRays.cpp"
#include "LEIM_SubtractBaseline.cpp"
#include "LEIM_CorrectEfficiency.cpp"
#include "LEIM_SumRescaleRebin.cpp"
#include "LEIM_Stitch.cpp"
#include "LEIM_WriteOutput.cpp"

using namespace std;

std::vector< std::vector<double> > wavelength_vis,wavelength_ir,intensity_vis,intensity_ir,baseline_vis,baseline_ir;
std::vector< std::vector<double> > intensity_vis_err_up,intensity_ir_err_up,intensity_vis_err_lo,intensity_ir_err_lo;
std::vector<double> timestamps;
std::vector< std::vector<double> > *wavelength_vis_ptr = &wavelength_vis;
std::vector< std::vector<double> > *wavelength_ir_ptr = &wavelength_ir;
std::vector< std::vector<double> > *intensity_vis_ptr = &intensity_vis;
std::vector< std::vector<double> > *intensity_ir_ptr = &intensity_ir;
std::vector< std::vector<double> > *baseline_vis_ptr = &baseline_vis;
std::vector< std::vector<double> > *baseline_ir_ptr = &baseline_ir;
std::vector<double> *timestamps_ptr = &timestamps;
// This is the spectrum builder for LEIM
// Steps:
//// 1. Read data from .csv files
//// 2. Calculate charge for the measurement
//// 3. Remove cosmic rays from data
//// 4. Subtract baseline from spectrum, calculate associated uncertainty
//// 5. Correct for optical efficiencies, calculate associated uncertainty
//// 6. Sum counts over y-pixels, convert counts to electrons per avalanche,
////    rebin histogram, calculate associated uncertainty
//// 7. Stitch visible and infrared components together
//// 8. Write to output root file
//
// Options:
//// "F" = FBK VUVHD3
//// "H" = Hamamatsu VUV4
//// "S" = Short exposure time
//// "M" = Medium exposure time
//// "L" = Long exposure time
void SpectrumBuilder(TString options){
  ofstream logfile;
  logfile.open(options+"_log.txt");
  logfile << "---------------------------------------------------------" << endl;
  logfile << "Building spectrum for " << options << " --------------------------------" << endl;
  logfile << endl;
  if (options.Contains("F") && options.Contains("H")){
    logfile << "ERROR: More than one SiPM type specified." << endl;
    exit(1);
  }
  if (!options.Contains("F") && !options.Contains("H")){
    logfile << "ERROR: SiPM type invalid or not specified." << endl;
    exit(1);
  }
  if ((options.Contains("S") && options.Contains("M")) ||
      (options.Contains("S") && options.Contains("L")) ||
      (options.Contains("M") && options.Contains("L"))){
    logfile << "ERROR: More than one exposure time specified." << endl;
    exit(1);
  }
  if (!options.Contains("S") && !options.Contains("M") && !options.Contains("L")){
    logfile << "ERROR: Exposure time not valid or not specified." << endl;
    exit(1);
  }

  gROOT->SetStyle("Modern");
  gROOT->ForceStyle();
  gStyle->SetPalette(51);
  gStyle->SetNumberContours(99);
  gStyle->SetOptStat(0);
    
  logfile << "Importing data" << endl;
  double charge_vis,charge_ir,charge_vis_err,charge_ir_err,avg_current_vis,avg_current_ir,avg_current_vis_err,avg_current_ir_err,begin_time_vis,end_time_vis,begin_time_ir,end_time_ir;
  double *charge_calcs_vis = new double[4];
  double *charge_calcs_ir = new double[4];
  ReadInData(wavelength_vis_ptr, wavelength_ir_ptr, intensity_vis_ptr, intensity_ir_ptr, baseline_vis_ptr, baseline_ir_ptr, timestamps_ptr, options);
  begin_time_vis = timestamps[0];
  end_time_vis = timestamps[1];
  begin_time_ir = timestamps[2];
  end_time_ir = timestamps[3];
  intensity_vis_err_up.resize(intensity_vis.size());
  intensity_vis_err_lo.resize(intensity_vis.size());
  intensity_ir_err_up.resize(intensity_ir.size());
  intensity_ir_err_lo.resize(intensity_ir.size());
  TString hdraw_title = "";
  if (options.Contains("H")) hdraw_title += "Hamamatsu ";
  else if (options.Contains("F")) hdraw_title += "FBK ";
  hdraw_title += "Raw Data with previous steps and Integrated, Re-binned, and Scaled for ";
  if (options.Contains("S")) hdraw_title += "Short Exposure Time;Wavelength [nm];Photons/e^{-}/4 nm bin";
  else if (options.Contains("M")) hdraw_title += "Medium Exposure Time;Wavelength [nm];Photons/e^{-}/4 nm bin";
  else if (options.Contains("L")) hdraw_title += "Long Exposure Time;Wavelength [nm];Photons/e^{-}/4 nm bin";


  logfile << "Calculating charge" << endl;
  charge_calcs_vis = CalculateCharge(options,begin_time_vis,end_time_vis,true);
  charge_vis = charge_calcs_vis[0];
  charge_vis_err = charge_calcs_vis[1];
  avg_current_vis = charge_calcs_vis[2];
  avg_current_vis_err = charge_calcs_vis[3];
  charge_calcs_ir = CalculateCharge(options,begin_time_ir,end_time_ir,false);
  charge_ir = charge_calcs_ir[0];
  charge_ir_err = charge_calcs_ir[1];
  avg_current_ir = charge_calcs_ir[2];
  avg_current_ir_err = charge_calcs_ir[3];
  logfile << "   total charge (vis) = " << charge_vis << " +/- " << charge_vis_err << endl;
  logfile << "   mean current (vis) = " << avg_current_vis << " +/- " << avg_current_vis_err << endl;
  logfile << "   total charge (ir)  = " << charge_ir << " +/- " << charge_ir_err << endl;
  logfile << "   mean current (ir)  = " << avg_current_ir << " +/- " << avg_current_ir_err << endl;
  logfile << endl;
  double SiPM_area,SPAD_area,exposed_area,integral_scale,integral_scale_err;
  double *integral_scale_info = new double[2];
  logfile << "Correcting for spatial emission fluctuations" << endl;
  if (options.Contains("F")){
    SiPM_area = 5.77e-3 * 5.77e-3;
    SPAD_area = 35.e-6 * 35.e-6;
    exposed_area = (9.*35./35.3125)*(400.*35/35.3125)*1.e-12;// = 3.5366e-9
    integral_scale_info = CorrectSpatialFluctuations("F");
    integral_scale = integral_scale_info[0];
    integral_scale_err = integral_scale_info[1];    
    //exposed_area = 10.5*35.e-6*0.2*35.e-6;// = 2.5725e-9 
  }
  else if (options.Contains("H")){
    SiPM_area = 3.e-3 * 3.e-3;
    SPAD_area = 50.e-6 * 50.e-6;
    exposed_area = (9.*50./50.9474)*(400.*50./50.9474)*1.e-12;
    integral_scale_info = CorrectSpatialFluctuations("H");
    integral_scale = integral_scale_info[0];
    integral_scale_err = integral_scale_info[1];    
    //exposed_area = 7.4 * 50.e-6 * 0.14 * 50.e-6;//need to measure
  }
  logfile << "   Photon flux scale factor = " << integral_scale << " +/- " << integral_scale_err << endl;
  logfile << endl;

  logfile << "Removing cosmic rays" << endl;
  intensity_vis = RemoveCosmicRays(intensity_vis,13,25);
  intensity_ir = RemoveCosmicRays(intensity_ir,13,25);
  baseline_vis = RemoveCosmicRays(baseline_vis,13,25);
  baseline_ir = RemoveCosmicRays(baseline_ir,13,25);

  logfile << "Subtracting baseline" << endl;
  double baseline_noise_ir,baseline_noise_vis;
  for (int i = 0; i < 400; i++){
    baseline_noise_ir = baseline_noise_vis = TMath::Power(3.26,2.0); //measured this independently using wtfIsGoingOn.cpp and is consistent with data sheet
    intensity_vis[i] = SubtractBaseline(intensity_vis[i],baseline_vis[i]);
    intensity_ir[i] = SubtractBaseline(intensity_ir[i],baseline_ir[i]);
    intensity_vis_err_up[i] = GetError(intensity_vis[i],baseline_noise_vis);
    intensity_ir_err_up[i] = GetError(intensity_ir[i],baseline_noise_ir);
    intensity_vis_err_lo[i] = intensity_vis_err_up[i];
    intensity_ir_err_lo[i] = intensity_ir_err_up[i];
  }

  logfile << "Correcting for optical efficiencies" << endl;
  for (int i = 0; i < 400; i++){
    intensity_vis[i] = CorrectEfficiency(wavelength_vis[i], intensity_vis[i], i, true);
    intensity_ir[i] = CorrectEfficiency(wavelength_ir[i], intensity_ir[i], i, false);

    intensity_vis_err_up[i] = GetError(wavelength_vis[i], intensity_vis[i], intensity_vis_err_up[i], i, true, true);
    intensity_ir_err_up[i] = GetError(wavelength_ir[i], intensity_ir[i], intensity_ir_err_up[i], i, false, true);
    intensity_vis_err_lo[i] = GetError(wavelength_vis[i], intensity_vis[i], intensity_vis_err_lo[i], i, true, false);
    intensity_ir_err_lo[i] = GetError(wavelength_ir[i], intensity_ir[i], intensity_ir_err_lo[i], i, false, false);
  }

  logfile << "Summing, rescaling, and rebinning" << endl;
  //std::vector< std::vector<double> > visible_rebinned = SumRescaleRebin(wavelength_vis,intensity_vis,intensity_vis_err,400,600,4,charge_vis,exposed_area,SiPM_area);
  //std::vector< std::vector<double> > infrared_rebinned = SumRescaleRebin(wavelength_ir,intensity_ir,intensity_ir_err,500,1050,4,charge_ir,exposed_area,SiPM_area);
  std::vector< std::vector<double> > visible_rebinned = SumRescaleRebin(wavelength_vis,intensity_vis,intensity_vis_err_up,intensity_vis_err_lo,400,600,4,charge_vis,charge_vis_err,integral_scale,integral_scale_err,options);
  std::vector< std::vector<double> > infrared_rebinned = SumRescaleRebin(wavelength_ir,intensity_ir,intensity_ir_err_up,intensity_ir_err_lo,500,1050,4,charge_ir,charge_ir_err,integral_scale,integral_scale_err,options);
  logfile << endl;
  logfile << "Stitching visible and infrared" << endl;
  std::vector< std::vector<double> > stitched_spectrum = Stitch(visible_rebinned[0],visible_rebinned[1],visible_rebinned[2],visible_rebinned[3],
								infrared_rebinned[0],infrared_rebinned[1],infrared_rebinned[2],infrared_rebinned[3]);
  
  //TH1F *hstitched = new TH1F("","",stitched_spectrum[0].size(),400,1050);
  TGraphAsymmErrors *hstitched = new TGraphAsymmErrors(stitched_spectrum[0].size());
  TString hname = "spectrum_";
  hname += options;
  TString htitle = "";
  if (options.Contains("H")) htitle += "Hamamatsu VUV4 ";
  else if (options.Contains("F")) htitle += "FBK VUVHD3 ";
  htitle += "Dark Noise Emission Spectrum at ";
  if (options.Contains("S")) htitle += "Short Exposure Time;Wavelength [nm];Photons/e^{-}/nm bin";
  else if (options.Contains("M")) htitle += "Medium Exposure Time;Wavelength [nm];Photons/e^{-}/nm bin";
  else if (options.Contains("L")) htitle += "Long Exposure Time;Wavelength [nm];Photons/e^{-}/nm bin";
  //hstitched->Sumw2();
  hstitched->SetName(hname);
  hstitched->SetTitle(htitle);
  for (int i = 0; i < stitched_spectrum[0].size(); i++){
    //hstitched->Fill(stitched_spectrum[0][i],0.25*stitched_spectrum[1][i]);
    //hstitched->SetBinError(i,0.25*stitched_spectrum[2][i]);
    hstitched->SetPoint(i,stitched_spectrum[0][i],0.25*stitched_spectrum[1][i]);
    hstitched->SetPointError(i,2,2,0.25*stitched_spectrum[2][i],0.25*stitched_spectrum[3][i]);
  }

  logfile << "Writing output" << endl;
  WriteOutput(options,hstitched);
  logfile << "done." << endl;
  logfile << "---------------------------------------------------------" << endl;
  logfile.close();
}

void DoAll(){
  SpectrumBuilder("HS");
  SpectrumBuilder("HM");
  SpectrumBuilder("HL");
  SpectrumBuilder("FS");
  SpectrumBuilder("FM");
  SpectrumBuilder("FL");
}
