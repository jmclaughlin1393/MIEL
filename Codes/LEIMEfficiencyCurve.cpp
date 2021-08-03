#include <iostream>
#include <fstream>

using namespace std;
ifstream gratingFile150,gratingFile300,mirrorFile,objectiveFile20XIR,objectiveFile20X,objectiveFile10X,cameraFile,longpassFile,dichroicFile,modFile;
ifstream sourceFile,measurementFileVis,baselineFileVis,measurementFileIR,baselineFileIR,measurementFile675,baselineFile675;
std::vector<double> grating_W,mirror_W,objective_W,objective_W_sorted,longpass_W,longpass_W_sorted,dichroic_W,dichroic_W_sorted,camera_W,mod_W,source_W;
std::vector<double> grating_WC,mirror_WC,objective_WC,objective_WC_sorted,longpass_WC,longpass_WC_sorted,dichroic_WC,dichroic_WC_sorted,camera_WC,source_WC;
std::vector<double> grating_I,mirror_I,objective_I,longpass_I,dichroic_I,camera_I,mod_I,source_I;
std::vector< std::vector<double> > measurement_W,baseline_W,spectrum_W;
std::vector< std::vector<double> > measurement_I,baseline_I,spectrum_I;
std::vector< std::vector<double> > measurement_X,baseline_X,spectrum_X;
std::vector< std::vector<double> > measurement_Y,baseline_Y,spectrum_Y;
double wavelength_microscope[6] = {300,400,500,700,950,1200};
double efficiency_microscope[6] = {0,0.86,0.95,0.95,0.87,0.7};
TSpline5 *sp_microscope = new TSpline5("sp_microscope",wavelength_microscope,efficiency_microscope,6);
/* double wavelength_microscope[8] = {400,449.5,570,593,754.2,844.8,950,1200}; */
/* double efficiency_microscope[8] = {0.25,0.364,0.778844,0.775842,0.847311,0.764709,0.87,0.7}; */
/* TSpline3 *sp_microscope = new TSpline3("sp_microscope",wavelength_microscope,efficiency_microscope,8); */
TFile *correctionFile = new TFile();

TGraph *Gcon;
TGraph *Glin;
TGraph *Gquad;
TF1 *func;
bool information_loaded = false;
bool acceptInput;
double minWavelength,maxWavelength;

TGraph *G_measurement = new TGraph();
TGraph *G_source = new TGraph();
TGraph *G_source_weighted = new TGraph();
TGraph *G_objective = new TGraph();
TGraph *G_mirror = new TGraph();
TGraph *G_grating = new TGraph();
TGraph *G_camera = new TGraph();
TGraph *G_mod = new TGraph();
TGraph *G_longpass = new TGraph();
TGraph *G_dichroic = new TGraph();
TGraph *G_overallEfficiency = new TGraph(1000);
TGraph *G_observedEff = new TGraph(1000);
double prv_range;

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
double CorrectWavelength(int range, double lambda, int row){
  if (range != prv_range){
    if (range == 800) correctionFile = TFile::Open("../Rootfiles/WavelengthCorrectionParameters_Infrared.root");
    else if (range == 500 || range == 675) correctionFile = TFile::Open("../Rootfiles/WavelengthCorrectionParameters_Visible.root");
    Gcon = (TGraph*)correctionFile->Get("Gcon");
    Glin = (TGraph*)correctionFile->Get("Glin");
    Gquad = (TGraph*)correctionFile->Get("Gquad");
  }
  prv_range = range;
  double con,lin,quad;
  func = Gcon->GetFunction("fexpo");
  con = func->Eval(double(row));
  func = Glin->GetFunction("fexpo");
  lin = func->Eval(double(row));
  func = Gquad->GetFunction("fexpo");
  quad = func->Eval(double(row));
  if (range == 500) return lambda;
  return lambda + con + lin*lambda + quad*lambda*lambda;
}

void ClearVectors(){
  grating_W.clear();
  mirror_W.clear();
  objective_W.clear();
  longpass_W.clear();
  dichroic_W.clear();
  mod_W.clear();
  camera_W.clear();
  source_W.clear();
  grating_I.clear();
  mirror_I.clear();
  objective_I.clear();
  longpass_I.clear();
  dichroic_I.clear();
  mod_I.clear();
  camera_I.clear();
  source_I.clear();
  measurement_W.clear();
  baseline_W.clear();
  spectrum_W.clear();
  measurement_I.clear();
  baseline_I.clear();
  spectrum_I.clear();
  measurement_X.clear();
  baseline_X.clear();
  spectrum_X.clear();
  measurement_Y.clear();
  baseline_Y.clear();
  spectrum_Y.clear();
  information_loaded = false;
}

void GetInformation(int range){
  cout << "Getting information" << endl;
  gratingFile150.open("../Data/150gpmm_transmissionCurve.csv");
  gratingFile300.open("../Data/300gpmm_transmissionCurve.csv");
  objectiveFile20XIR.open("../Data/LCPLN20XIR_transmissionCurve.csv");
  objectiveFile20X.open("../Data/LMPLFLN20X_transmissionCurve_v2.csv");
  objectiveFile10X.open("../Data/LMPLFLN10X_transmissionCurve.csv");
  cameraFile.open("../Data/PyLoN_BRexcelon_QEcurve.csv");
  mirrorFile.open("../Data/Mirror_reflectionCurve.csv");
  longpassFile.open("../Data/LongPassFilter_TransmissionCurve.csv");
  dichroicFile.open("../Data/DichroicFilter_TransmissionCurve.csv");
  //microscopeFile.open("../Data/OlympusMicroscope_LeftSidePort_TransmissionCurve.csv");
  modFile.open("../Data/TransmissionCurve_Modifier_v2.csv");
  //sourceFile.open("../Data/LSVN0583.csv");
  sourceFile.open("../Data/LSVN0620.csv");
  //sourceFile.open("../Data/LSVN0583_CrystalPenner.csv");
  measurementFileIR.open("../Data/calibration_data/IntensityCalSource_IRspectrum.csv");
  //measurementFileIR.open("../Data/calibration_data/LSVN0620_NIRspectrum_200ms_50kHz.csv");
  baselineFileIR.open("../Data/calibration_data/IntensityCalSource_IRbackground.csv");
  //baselineFileIR.open("../Data/calibration_data/LSVN0620_NIRbackground_200ms_50kHz.csv");
  measurementFileVis.open("../Data/calibration_data/IntensityCalSource_VISspectrum_500.csv");
  //measurementFileVis.open("../Data/calibration_data/LSVN0620_VISspectrum_200ms_50kHz.csv");
  baselineFileVis.open("../Data/calibration_data/IntensityCalSource_VISbackground.csv");
  //baselineFileVis.open("../Data/calibration_data/LSVN0620_VISbackground_200ms_50kHz.csv");
  measurementFile675.open("../Data/calibration_data/IntensityCalibrationSource_150lpmmGrating_675nmCentre.csv");
  baselineFile675.open("../Data/calibration_data/IntensityCalSource_VISbackground.csv");
  
  if (range == 800) {minWavelength = 500; maxWavelength = 1050;}
  else if (range == 675) {minWavelength = 400; maxWavelength = 700;}
  else {minWavelength = 400; maxWavelength = 640;}
  std::vector<double> w,i,x,y,b;
  std::vector<float> values;
  TString line, line1, line2, baseline;
  double yval = 0;
  
  //Measured calibration source information--------------
  //-------------------
  if (range == 800) {acceptInput = !measurementFileIR.eof(); measurementFileIR >> line; baselineFileIR >> baseline;}
  else if (range == 675) {acceptInput = !measurementFile675.eof(); measurementFile675 >> line; baselineFile675 >> baseline;}
  else {acceptInput = !measurementFileVis.eof(); measurementFileVis >> line; baselineFileVis >> baseline;}
  cout << "   measurement data..." << endl;
  if (measurementFileIR.is_open()) cout << "file opened" << endl;
  while (acceptInput){    
    if (range == 800) {measurementFileIR >> line; baselineFileIR >> baseline;}
    else if (range == 675) {measurementFile675 >> line; baselineFile675 >> baseline;}
    else {measurementFileVis >> line; baselineFileVis >> baseline;}
    values = split(line);    
    if (values[3] == yval){
      w.push_back(values[0]);
      i.push_back(values[1]);
      x.push_back(values[2]);
      y.push_back(values[3]);
      yval = values[3];
      values = split(baseline);
      b.push_back(values[1]);
    }
    else {
      measurement_W.push_back(w);
      w.clear();
      w.push_back(values[0]);
      measurement_I.push_back(i);
      i.clear();
      i.push_back(values[1]);
      measurement_X.push_back(x);
      x.clear();
      x.push_back(values[2]);
      measurement_Y.push_back(y);
      y.clear();
      y.push_back(values[3]);
      yval = values[3];
      baseline_I.push_back(b);
      b.clear();
      values = split(baseline);
      b.push_back(values[1]);
    }
    if (range == 800) acceptInput = !measurementFileIR.eof();
    else if (range == 675) acceptInput = !measurementFile675.eof();
    else acceptInput = !measurementFileVis.eof();
  }
  cout << "   grating data..." << endl;
  //Grating efficiency information----------------------
  //----------------
  if (range == 800 || range == 675) acceptInput = !gratingFile150.eof();
  else acceptInput = !gratingFile300.eof();
  while (acceptInput){
    if (range == 800 || range == 675) gratingFile150 >> line1 >> line2;
    else gratingFile300 >> line1 >> line2;
    line = line1 + line2;
    values = split(line);
    if (values[1]/100. > 1.e-05 && values[1]/100. < 1 &&
	values[0] >= minWavelength && values[0] <= maxWavelength){
      grating_W.push_back(values[0]);
      grating_I.push_back(values[1]/100.);
    }
    if (range == 800 || range == 675) acceptInput = !gratingFile150.eof();
    else acceptInput = !gratingFile300.eof();
  }

  cout << "   objective data..." << endl;
  //Objective efficiency information-------------
  //-------------------
  if (range == 800) acceptInput = !objectiveFile20XIR.eof();
  else acceptInput = !objectiveFile20X.eof();
  while (acceptInput){
    if (range == 800) objectiveFile20XIR >> line1 >> line2;
    else objectiveFile20X >> line1 >> line2;
    line = line1 + line2;
    values = split(line);
    if (values[1]/100. > 1.e-05 && values[1]/100. < 1 &&
        values[0] >= minWavelength && values[0] <= maxWavelength){
      objective_W.push_back(values[0]);
      objective_W_sorted.push_back(values[0]);
      objective_I.push_back(values[1]/100.);
    }
    if (range == 800) acceptInput = !objectiveFile20XIR.eof();
    else acceptInput = !objectiveFile20X.eof();
  }
  std::sort(objective_W_sorted.begin(),objective_W_sorted.end());

  cout << "   camera data..." << endl;
  //Camera efficiency information------------------
  //----------------
  acceptInput =	!cameraFile.eof();
  double prv;
  while (acceptInput){
    cameraFile >> line1 >> line2;
    line = line1 + line2;
    values = split(line);
    if (values[0] <= prv){
      acceptInput = !cameraFile.eof();
      continue;
    }
    if (values[1]/100. > 1.e-05 && values[1]/100. < 1 &&
        values[0] >= minWavelength && values[0] <= maxWavelength){
      camera_W.push_back(values[0]);
      camera_I.push_back(values[1]/100.);
    }
    prv = values[0];
    acceptInput = !cameraFile.eof();
  }

  cout << "   mirror data..." << endl;
  //Mirror efficiency information--------------------
  //----------------
  acceptInput = !mirrorFile.eof();
  while (acceptInput){
    mirrorFile >> line1 >> line2;
    line = line1 + line2;
    values = split(line);
    if (values[1]/100. > 1.e-05 && values[1]/100. < 1 &&
        values[0] >= minWavelength && values[0] <= maxWavelength){
      mirror_W.push_back(values[0]);
      mirror_I.push_back(values[1]/100.);
    }
    acceptInput = !mirrorFile.eof();
  }

  cout << "   longpass data..." << endl;
  //Longpass filter efficiency information----------
  //----------------
  acceptInput = !longpassFile.eof();
  while (acceptInput){
    longpassFile >> line1 >> line2;
    line = line1 + line2;
    values = split(line);
    if (values[1] > 1.e-05 && values[1] < 1 &&
        values[0] >= minWavelength && values[0] <= maxWavelength){
      longpass_W.push_back(values[0]);
      longpass_W_sorted.push_back(values[0]);
      longpass_I.push_back(values[1]);
    }
    acceptInput = !longpassFile.eof();
  }
  std::sort(longpass_W_sorted.begin(),longpass_W_sorted.end());
  
  cout << "   dichroic data..." << endl;
  //Dichroic filter efficiency information----------
  //---------------- 
  acceptInput = !dichroicFile.eof();
  while (acceptInput){
    dichroicFile >> line1 >> line2;
    line = line1 + line2;
    values = split(line);
    if (values[1] > 1.e-05 && values[1] < 1 &&
        values[0] >= minWavelength && values[0] <= maxWavelength){
      dichroic_W.push_back(values[0]);
      dichroic_W_sorted.push_back(values[0]);
      dichroic_I.push_back(values[1]);
    }
    acceptInput = !dichroicFile.eof();
  }
  std::sort(dichroic_W_sorted.begin(),dichroic_W_sorted.end());
  
  cout << "   lsp data..." << endl;
  //Microscope LSP efficiency information----------
  //---------------- 
  acceptInput = !modFile.eof();
  while (acceptInput){
    modFile >> line;
    values = split(line);
    if (values[1] > 1.e-05 && values[1] < 1 &&
        values[0] >= minWavelength && values[0] <= maxWavelength*1.2){
      mod_W.push_back(values[0]);
      mod_I.push_back(values[1]);
    }
    acceptInput = !modFile.eof();
  }

  cout << "   LSVN0583 source data..." << endl;
  //Expected calibration source information--------------
  //-------------------
  acceptInput = !sourceFile.eof();
  while (acceptInput){
    sourceFile >> line;
    values = split(line);
    if (values[0] >= 400 && values[0] <= 1100){ //[1] if using crystal's data
      source_W.push_back(values[0]);//[1] if using crystal's data
      source_I.push_back(values[1]);//[2] if using crystal's data
    }
    acceptInput = !sourceFile.eof();
  }

  cout << "   closing files." << endl;
  gratingFile150.close();
  gratingFile300.close();
  objectiveFile20XIR.close();
  objectiveFile20X.close();
  objectiveFile10X.close();
  cameraFile.close();
  mirrorFile.close();
  modFile.close();
  sourceFile.close();
  measurementFileIR.close();
  baselineFileIR.close();
  measurementFileVis.close();
  baselineFileVis.close();
  measurementFile675.close();
  baselineFile675.close();
  information_loaded = true;
}

void LEIMEfficiencyCurve(int range,int row,bool writeData=false){
  gROOT->SetStyle("Modern");
  gROOT->ForceStyle();
  gStyle->SetPalette(51);
  gStyle->SetNumberContours(99);
  gStyle->SetOptStat(0);

  if (!information_loaded) GetInformation(range);
  //Correcting Wavelengths--------------------------------------
  //--------------------
  grating_WC.clear();mirror_WC.clear();objective_WC.clear();objective_WC_sorted.clear();longpass_WC.clear();longpass_WC_sorted.clear();dichroic_WC.clear();dichroic_WC_sorted.clear();camera_WC.clear();source_WC.clear();
  for (int i = 0; i < grating_W.size(); i++){
    grating_WC.push_back(0);
    grating_WC[i] = CorrectWavelength(range,grating_W[i],row);
  }
  for (int i = 0; i < mirror_W.size(); i++){
    mirror_WC.push_back(0);
    mirror_WC[i] = CorrectWavelength(range,mirror_W[i],row);
  }
  for (int i = 0; i < objective_W.size(); i++){
    objective_WC.push_back(0);
    objective_WC[i] = CorrectWavelength(range,objective_W[i],row);
  }
  for (int i = 0; i < objective_W_sorted.size(); i++){
    objective_WC_sorted.push_back(0);
    objective_WC_sorted[i] = CorrectWavelength(range,objective_W_sorted[i],row);
  }
  for (int i = 0; i < camera_W.size(); i++){
    camera_WC.push_back(0);
    camera_WC[i] = CorrectWavelength(range,camera_W[i],row);
  }
  for (int i = 0; i < longpass_W.size(); i++){
    longpass_WC.push_back(0);
    longpass_WC[i] = CorrectWavelength(range,longpass_W[i],row);
  }
  for (int i = 0; i < longpass_W_sorted.size(); i++){
    longpass_WC_sorted.push_back(0);
    longpass_WC_sorted[i] = CorrectWavelength(range,longpass_W_sorted[i],row);
  }
  for (int i = 0; i < dichroic_W.size(); i++){
    dichroic_WC.push_back(0);
    dichroic_WC[i] = CorrectWavelength(range,dichroic_W[i],row);
  }
  for (int i = 0; i < dichroic_W_sorted.size(); i++){
    dichroic_WC_sorted.push_back(0);
    dichroic_WC_sorted[i] = CorrectWavelength(range,dichroic_W_sorted[i],row);
  }
  for (int i = 0; i < source_W.size(); i++){
    source_WC.push_back(0);
    source_WC[i] = CorrectWavelength(range,source_W[i],row) + 0.9; //wavelength offset
  }
  
  //Constructing expected individual component efficiency curves
  //--------------------
  G_grating->Set(grating_WC.size());
  for (int i = 0; i < grating_WC.size(); i++) G_grating->SetPoint(i,grating_WC[i],grating_I[i]);
  G_mirror->Set(mirror_WC.size());
  for (int i = 0; i < mirror_WC.size(); i++) G_mirror->SetPoint(i,mirror_WC[i],mirror_I[i]);
  G_objective->Set(objective_WC.size());
  for (int i = 0; i < objective_WC_sorted.size(); i++){
    for (int j = 0; j < objective_WC.size(); j++){
      if (objective_WC_sorted[i] != objective_WC[j]) continue;
      G_objective->SetPoint(i,objective_WC_sorted[i],objective_I[j]);
    }
  }
  G_camera->Set(camera_WC.size());
  for (int i = 0; i < camera_WC.size(); i++) G_camera->SetPoint(i,camera_WC[i],camera_I[i]);
  G_longpass->Set(longpass_WC.size());
  for (int i = 0; i < longpass_WC.size(); i++) {
    for (int j = 0; j < longpass_WC.size(); j++){
      if (longpass_WC_sorted[i] != longpass_WC[j]) continue;
      G_longpass->SetPoint(i,longpass_WC_sorted[i],longpass_I[j]);
    }
  }
  G_dichroic->Set(dichroic_WC.size());
  for (int i = 0; i < dichroic_WC.size(); i++) {
    for (int j = 0; j < dichroic_WC.size(); j++){
      if (dichroic_WC_sorted[i] != dichroic_WC[j]) continue;
      G_dichroic->SetPoint(i,dichroic_WC_sorted[i],dichroic_I[j]);
    }
  }
  G_mod->Set(mod_W.size());
  for (int i = 0; i < mod_W.size(); i++) G_mod->SetPoint(i,mod_W[i],mod_I[i]);
  
  //Constructing expected overall efficiency curve--------------
  //--------------------
  double wavelength,overall_efficiency;
  for (int i = 0; i < 1000; i++){
    wavelength = minWavelength + (maxWavelength - minWavelength)*i/1000.;
    overall_efficiency  = G_objective->Eval(wavelength);
    if (range == 800) {overall_efficiency *= G_longpass->Eval(wavelength);}
    if (range == 800) {overall_efficiency *= G_dichroic->Eval(wavelength);}
    overall_efficiency *= TMath::Power(G_mirror->Eval(wavelength),3.0);
    overall_efficiency *= G_grating->Eval(wavelength);
    overall_efficiency *= sp_microscope->Eval(wavelength);
    overall_efficiency *= G_mod->Eval(wavelength);
    overall_efficiency *= G_camera->Eval(wavelength);
    G_overallEfficiency->SetPoint(i,wavelength,overall_efficiency);
  }

  //Constructing measured and expected spectra--------------
  //--------------------
  G_measurement->Set(measurement_W[row].size());
  for (int i = 0; i < measurement_W[row].size(); i++){
    if (range == 800) G_measurement->SetPoint(i,measurement_W[row][i],measurement_I[row][i]-baseline_I[row][i]);
    if (range == 675) G_measurement->SetPoint(i,measurement_W[row][i],measurement_I[row][i]);
    if (range == 500) G_measurement->SetPoint(i,measurement_W[row][i],measurement_I[row][i]-baseline_I[row][i]);
    //if (range == 500) G_measurement->SetPoint(i,measurement_W[row][i],measurement_I[row][i]-135.);
  }
  G_source->Set(source_WC.size());
  for (int i = 0; i < source_WC.size(); i++) G_source->SetPoint(i,source_WC[i],source_I[i]);
  
  //Constructing measured efficiency curve--------------
  //--------------------
  double wavelength_observedEff;
  double efficiency_observedEff;
  for (int i = 0; i < 1000; i++){
    wavelength_observedEff = minWavelength + (maxWavelength - minWavelength)*i/1000.;
    //efficiency_observedEff = (0.35/45000.)*G_measurement->Eval(wavelength_observedEff)/G_source->Eval(wavelength_observedEff);
    //if (wavelength_observedEff > 970) efficiency_observedEff = G_overallEfficiency->Eval(wavelength_observedEff);
    efficiency_observedEff = G_overallEfficiency->Eval(wavelength_observedEff);
    G_observedEff->SetPoint(i,wavelength_observedEff,efficiency_observedEff);
  }
  if (range == 500) G_observedEff->SetTitle("Observed Efficiency Curve for Visible Spectrum Measurements;Wavelength [nm]; Relative Efficiency [a.u.]");
  else G_observedEff->SetTitle("Observed Efficiency Curve for Infrared Spectrum Measurements;Wavelength [nm]; Relative Efficiency [a.u.]");

  
  if (writeData){
    TFile *out = new TFile();
    if (range == 500) out = TFile::Open("../Rootfiles/VisibleRangeEfficiencyCurve.root","UPDATE");
    if (range == 800) out = TFile::Open("../Rootfiles/InfraredRangeEfficiencyCurve.root","UPDATE");  
    out->cd();
    G_observedEff->SetName(Form("EfficiencyCurve_row%d",row));
    G_observedEff->Write();
    out->Close();
  } else {
    double WL,IN;
    G_source_weighted->Set(source_WC.size());
    for (int i = 0; i < source_WC.size(); i++){
      WL = source_WC[i];
      IN = source_I[i]*G_overallEfficiency->Eval(WL);
      //IN = source_I[i]*G_observedEff->Eval(WL);
      G_source_weighted->SetPoint(i,WL,IN);
    }

    TCanvas *cEff = new TCanvas("cEff","cEff");
    cEff->cd();
    TString title = "Comparing Efficiency Curves;Wavelength [nm];Efficiency [a.u.]";
    TH2F *h1 = new TH2F("h1",title,1000,minWavelength,maxWavelength,105,0,0.5);
    h1->Draw();
    G_overallEfficiency->SetLineColor(1);
    G_overallEfficiency->SetLineWidth(2);    
    G_overallEfficiency->Draw("SAME");

    TCanvas *cSpec = new TCanvas("cSpec","cSpec");
    cSpec->cd();
    double measMax = 0;
    double newMax = 0;
    double maxVal = 0;
    double XX,YY;
    for (int i = 0; i < G_measurement->GetN(); i++) {      
      G_measurement->GetPoint(i,XX,YY);
      //if (XX > maxWavelength || XX < 500) continue;
      if (measMax < YY) measMax = YY;
    }
    //int nPoints = G_measurement->GetN();
    //G_measurement->GetPoint(nPoints - 5,XX,YY);
    //measMax = YY;
    //cout << measMax << endl;
    for (int i = 0; i < G_measurement->GetN(); i++) {
      G_measurement->GetPoint(i,XX,YY);
      YY /= measMax;
      G_measurement->SetPoint(i,XX,YY);
    }
    for (int i = 0; i < G_source_weighted->GetN(); i++){
      G_source_weighted->GetPoint(i,XX,YY);
    //   if (XX > maxWavelength || XX < 500) continue;
      if (YY > newMax) newMax = YY;
    }
    //G_measurement->GetPoint(nPoints - 5,XX,YY);
    //newMax = G_source_weighted->Eval(XX);
    for (int i = 0; i < G_source_weighted->GetN(); i++){
      G_source_weighted->GetPoint(i,XX,YY);
      YY /= newMax;//*0.44;
      G_source_weighted->SetPoint(i,XX,YY);
      if (maxVal < YY) maxVal = YY;
    }
    double chi2 = 0;
    double NDF = 0;
    for (double x = minWavelength; x < maxWavelength; x += 5.0){
      chi2 += TMath::Power(G_measurement->Eval(x) - G_source_weighted->Eval(x),2.0)/0.001;
      NDF += 1.0;
    }
    cout << "chi2 = " << chi2 << endl;
    cout << "NDF  = " << NDF << endl;
    cout << "chi2/NDF = " << chi2/NDF << endl;

    TH2F *h2 = new TH2F("h2","Intensity Calibration Spectrum, Expected vs Observed;Wavelength [nm];Counts[a.u.]",1000,minWavelength,maxWavelength,110*maxVal,0,1.1*maxVal);
    h2->Draw();
    G_measurement->SetLineColor(1);
    G_measurement->Draw("SAME");
    G_source_weighted->SetLineColor(2);
    G_source_weighted->Draw("SAME");
    TLegend *leg2 = new TLegend(.15,.7,.4,.83);
    leg2->AddEntry(G_measurement,"Observed Spectrum");
    leg2->AddEntry(G_source_weighted,"Expected Spectrum");
    leg2->Draw("SAME");
    
    // TCanvas *cSource = new TCanvas("cSource","cSource");
    // cSource->cd();
    // cSource->SetGridx();
    // cSource->SetGridy();
    // G_source->SetTitle("LSVN0583 Independently Measured Spectrum;Wavelength [nm];Relative Amplitude [a.u.]");
    // G_source->SetLineColor(kAzure+6);
    // G_source->SetLineWidth(2);
    // G_source->Draw("ACP");
  }
}

void DoAll(){
  for (int i = 0; i < 400; i++){
    cout << "Visible range, row " << i << endl;
    LEIMEfficiencyCurve(500,i,true);
  }
  
  ClearVectors();
  for (int i = 0; i < 400; i++){
    cout << "Infrared range, row " << i << endl;
    LEIMEfficiencyCurve(800,i,true);
  }
  gROOT->ProcessLine(".q");
}
