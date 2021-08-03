TGraph *gCurrent = new TGraph();
double *CalculateCharge(TString options,double begin_time,double end_time,bool vis){
  ifstream current_vs_time;
  double *charge_values = new double[4];
  if (vis){
    if (options.Contains("F")){
      if (options.Contains("S")) current_vs_time.open("../Data/current_vs_time/CurrentVsTime_newFBKVUV-HD3-VIS_3h20m00s.txt");
      else if (options.Contains("M")) current_vs_time.open("../Data/current_vs_time/CurrentVsTime_newFBKVUV-HD3-VIS_4h45m43s.txt");
      else if (options.Contains("L")) current_vs_time.open("../Data/current_vs_time/CurrentVsTime_newFBKVUV-HD3-VIS_8h20m00s.txt");
    }
    else if (options.Contains("H")){
      if (options.Contains("S")) current_vs_time.open("../Data/current_vs_time/CurrentVsTime_newHamamatsuVUV4-VIS_3h20m00s.txt");
      else if (options.Contains("M")) current_vs_time.open("../Data/current_vs_time/CurrentVsTime_newHamamatsuVUV4-VIS_4h45m43s.txt");
      else if (options.Contains("L")) current_vs_time.open("../Data/current_vs_time/CurrentVsTime_newHamamatsuVUV4-VIS_8h20m00s.txt");
    }
  } else {
    if (options.Contains("F")){
      if (options.Contains("S")) current_vs_time.open("../Data/current_vs_time/CurrentVsTime_newFBKVUV-HD3-IR_3h20m00s.txt");
      else if (options.Contains("M")) current_vs_time.open("../Data/current_vs_time/CurrentVsTime_newFBKVUV-HD3-IR_4h45m43s.txt");
      else if (options.Contains("L")) current_vs_time.open("../Data/current_vs_time/CurrentVsTime_newFBKVUV-HD3-IR_8h20m00s.txt");
    }
    if (options.Contains("H")){
      if (options.Contains("S")) current_vs_time.open("../Data/current_vs_time/CurrentVsTime_newHamamatsuVUV4-IR_3h20m00s.txt");
      else if (options.Contains("M")) current_vs_time.open("../Data/current_vs_time/CurrentVsTime_newHamamatsuVUV4-IR_4h45m43s.txt");
      else if (options.Contains("L")) current_vs_time.open("../Data/current_vs_time/CurrentVsTime_newHamamatsuVUV4-IR_8h20m00s.txt");
    }
  }
  
  std::vector<double> current,nHours,nMins,nSecs;
  std::vector<float> values;
  float prv_hour;
  bool newday = false;
  TString line;
  bool acceptInput = !current_vs_time.eof();
  while (acceptInput) {
    current_vs_time >> line;
    values = split(line);
    if (values[0] <= 0 || values[0] >= 2.6e-3){
      acceptInput = !current_vs_time.eof();
      continue;
    }
    current.push_back(values[0]);
    if (prv_hour == 23 && values[1] == 0 && !newday) newday = true;
    if (newday) values[1] += 24.;
    nHours.push_back(values[1]);
    nMins.push_back(values[2]);
    nSecs.push_back(values[3]);
    prv_hour = values[1];
    acceptInput = !current_vs_time.eof();
  }
  gCurrent->Set(current.size());
  double curr,time;
  for (int i = 0; i < current.size(); i++){
    curr = current[i];
    time = nHours[i]*3600. + nMins[i]*60. + nSecs[i] - begin_time;
    gCurrent->SetPoint(i,time,curr);
  }

  double charge = 0;
  double charge_err = 0;
  double average = 0;
  double average_err = 0;
  double counter = 0;
  for (double t = 0; t < end_time-begin_time; t += 0.1){
    if (options.Contains("HM") && vis && t >= 6720) break;
    charge += gCurrent->Eval(t)*0.1;
    //charge_err += 100.e-18;
    charge_err += 25.e-16*0.1;//4.e-14*0.1
    average += gCurrent->Eval(t)*0.1;
  }
  if (options.Contains("HM") && vis){
    average /= 6720.;
    for (double t = 6720; t < end_time-begin_time; t += 0.1){
      charge += average*0.1;
      //charge_err += 100.e-18;
      charge_err += 25.e-16*0.1;//4.e-14*0.1
    }
  }
  else average /= end_time-begin_time;
  charge_err = TMath::Sqrt(charge_err);
  for (double t = 0; t < end_time-begin_time; t += 0.1){
    if (options.Contains("HM") && vis && t >= 6720) break;
    average_err += TMath::Power(gCurrent->Eval(t) - average,2.0);
    counter += 1.0;
  }
  average_err /= counter;
  average_err = TMath::Sqrt(average_err);

  charge_values[0] = charge;
  charge_values[1] = TMath::Sqrt(charge_err);
  charge_values[2] = average;
  charge_values[3] = average_err;
  
  TString filename = "../Rootfiles/";
  filename += options;
  filename += "_currentVtime";
  if (vis) filename += "_vis.root";
  else filename += "_ir.root";
  TFile *fout = new TFile(filename,"RECREATE");
  fout->cd();
  gCurrent->Write();
  fout->Close();
  return charge_values;
}
