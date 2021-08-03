#include <iostream>
#include <fstream>

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

int ReadInData(std::vector< std::vector<double> >* wavelength_vis, std::vector< std::vector<double> >* wavelength_ir,
	       std::vector< std::vector<double> >* intensity_vis, std::vector< std::vector<double> >* intensity_ir,
	       std::vector< std::vector<double> >* base_vis, std::vector< std::vector<double> >* base_ir,
	       std::vector<double>* timestamps, TString options){
  
  wavelength_vis->clear();
  wavelength_ir->clear();
  intensity_vis->clear();
  intensity_ir->clear();
  base_vis->clear();
  base_ir->clear();
  timestamps->clear();
  
  ifstream baseline_vis,spectrum_vis,baseline_ir,spectrum_ir;
  std::vector<float> values;
  std::vector<double> w,i,b;
  double yval;
  bool acceptInput;
  TString line,baseline;
  
  if (options.Contains("H") && options.Contains("S")){
    cout << " - Loading Hamamatsu short exposure" << endl;    
    baseline_vis.open("../Data/REDO_Hamamatsu-VUV4_VISbaseline_200min50kHz_200uA_11_Dec_2020-raw.csv");
    spectrum_vis.open("../Data/newHamamatsuVUV4_VISspectrum_12000s-2021_02_02-16_23_21.csv");
    baseline_ir.open("../Data/REDO_Hamamatsu-VUV4_IRbaseline_200min50kHz_200uA-2020_Dec_14-raw.csv");
    spectrum_ir.open("../Data/newHamamatsuVUV4_IRspectrum_12000s-2021_02_30-14_33_30.csv");
    timestamps->push_back(16*3600. + 23*60. + 21.);
    timestamps->push_back(timestamps->at(0) + 12000.);
    timestamps->push_back(14*3600. + 33*60. + 30.);
    timestamps->push_back(timestamps->at(2) + 12000.);
  } else if (options.Contains("H") && options.Contains("M")){
    cout << " - Loading Hamamatsu medium exposure" << endl; 
    //baseline_vis.open("../Data/REDO_Hamamatsu-VUV4_VISbaseline_20min50kHz_2mA_11_Dec_2020-raw.csv");
    baseline_vis.open("../Data/newHamamatsuVUV4_VISbackground_17143s-2021_02_04-11_34_45.csv");
    spectrum_vis.open("../Data/newHamamatsuVUV4_VISspectrum_17143s-2021_02_03-08_35_45.csv");
    //baseline_ir.open("../data/REDO_Hamamatsu-VUV4_IRbaseline_20min50kHz_2mA_14_Dec_2020-raw.csv");
    baseline_ir.open("../Data/newHamamatsuVUV4_IRbackground_17143s-2021_02_01-19_03_09.csv");
    spectrum_ir.open("../Data/newHamamatsuVUV4_IRspectrum_17143s-2021_01_30-17_59_54.csv");
    timestamps->push_back( 8*3600. + 35*60. + 45.);
    timestamps->push_back(timestamps->at(0) + 17143.);
    timestamps->push_back(17*3600. + 59*60. + 54.);
    timestamps->push_back(timestamps->at(2) + 17143.);
  } else if (options.Contains("H") && options.Contains("L")){
    cout << " - Loading Hamamatsu long exposure" << endl; 
    baseline_vis.open("../Data/REDO_Hamamatsu-VUV4_VISbaseline_9h50kHz_74uA_11_Dec_2020-raw.csv");
    spectrum_vis.open("../Data/newHamamatsuVUV4_VISspectrum_30000s-2021_02_02-19_49_43.csv");
    baseline_ir.open("../Data/REDO_Hamamatsu-VUV4_IRbaseline_9h50kHz_74uA-2020_Dec_14-raw.csv");
    spectrum_ir.open("../Data/newHamamatsuVUV4_IRspectrum_30000s-2021_01_31-09_55_44.csv");
    timestamps->push_back(19*3600. + 49*60. + 43.);
    timestamps->push_back(timestamps->at(0) + 30000.);
    timestamps->push_back( 9*3600. + 55*60. + 44.);
    timestamps->push_back(timestamps->at(2) + 30000.);
  } else if (options.Contains("F") && options.Contains("S")){
    cout << " - Loading FBK short exposure" << endl;
    baseline_vis.open("../Data/REDO_FBK-VUVHD3_VISbaseline_200min50kHz_200uA_09_Dec_2020-raw.csv");
    spectrum_vis.open("../Data/newFBKVUV-HD3_VISspectrum_12000s-2021_02_04-17_40_30.csv");
    baseline_ir.open("../Data/REDO_FBK-VUVHD3_IRbaseline_200min50kHz_200uA_07_Dec_2020-raw.csv");
    spectrum_ir.open("../Data/newFBKVUV-HD3_IRspectrum_12000s-2021_01_27-17_51_59.csv");
    timestamps->push_back(17*3600. + 40*60. + 30.);
    timestamps->push_back(timestamps->at(0) + 12000.);
    timestamps->push_back(17*3600. + 51*60. + 59.);
    timestamps->push_back(timestamps->at(2) + 12000.);
  } else if (options.Contains("F") && options.Contains("M")){
    cout << " - Loading FBK medium exposure" << endl;
    //baseline_vis.open("../Data/REDO_FBK-VUVHD3_VISbaseline_200min50kHz_200uA_09_Dec_2020-raw.csv");
    baseline_vis.open("../Data/newFBKVUV-HD3_VISbackground_17143s-2021_02_05-14_13_26.csv");
    spectrum_vis.open("../Data/newFBKVUV-HD3_VISspectrum_17143s-2021_02_05-07_37_57.csv");
    //baseline_ir.open("../Data/REDO_FBK-VUVHD3_IRbaseline_200min50kHz_200uA_07_Dec_2020-raw.csv");
    baseline_ir.open("../Data/newFBKVUV-HD3_IRbackground_17143s-2021_01_28-20_22_52.csv");
    spectrum_ir.open("../Data/newFBKVUV-HD3_IRspectrum_17143s-2021_01_28-12_01_00.csv");
    timestamps->push_back( 7*3600. + 37*60. + 57.);
    timestamps->push_back(timestamps->at(0) + 17143.);
    timestamps->push_back(12*3600. +  1*60. +  0.);
    timestamps->push_back(timestamps->at(2) + 17143.);
  } else if (options.Contains("F") && options.Contains("L")){
    cout << " - Loading FBK long exposure" << endl;
    baseline_vis.open("../Data/REDO_FBK-VUVHD3_VISbaseline_9h50kHz_74uA_09_Dec_2020-raw.csv");
    spectrum_vis.open("../Data/newFBKVUV-HD3_VISspectrum_30000s-2021_02_04-22_05_28.csv");
    baseline_ir.open("../Data/REDO_FBK-VUVHD3_IRbaseline_9h50kHz_74uA_07_Dec_2020-raw.csv");
    spectrum_ir.open("../Data/newFBKVUV-HD3_IRspectrum_30000s-2021_01_27-21_38_57.csv");
    timestamps->push_back(22*3600. +  5*60. + 28.);
    timestamps->push_back(timestamps->at(0) + 30000.);
    timestamps->push_back(21*3600. + 38*60. +  57.);
    timestamps->push_back(timestamps->at(2) + 30000.);    
  }
  
  if (!baseline_vis.is_open()){ cout << " ! failed to open baseline_vis." << endl; exit(1); }
  else if (!spectrum_vis.is_open()){ cout << "!  failed to open spectrum_vis." << endl; exit(1); }
  else if (!baseline_ir.is_open()){ cout << " ! failed to open baseline_ir." << endl; exit(1); }
  else if (!spectrum_ir.is_open()){ cout << " ! failed to open spectrum_ir." << endl; exit(1); }
  yval = 0;
  spectrum_ir >> line;
  values = split(line);
  baseline_ir >> baseline;
  acceptInput = !spectrum_ir.eof();
  while (acceptInput){
    spectrum_ir >> line;
    baseline_ir >> baseline;
    values = split(line);
    if (values[3] == yval){
      w.push_back(double(values[0]));
      i.push_back(double(values[1]));
      yval = values[3];
      values = split(baseline);
      b.push_back(double(values[1]));
    }
    else {
      wavelength_ir->push_back(w);
      w.clear();
      w.push_back(double(values[0]));
      intensity_ir->push_back(i);
      i.clear();
      i.push_back(double(values[1]));
      yval = values[3];
      base_ir->push_back(b);
      b.clear();
      values = split(baseline);
      b.push_back(double(values[1]));
    }
    acceptInput = !spectrum_ir.eof();
  }

  w.clear(); i.clear(); b.clear();
  spectrum_vis >> line;
  values = split(line);
  baseline_vis >> baseline;
  yval = 0;
  acceptInput = !spectrum_vis.eof();
  while (acceptInput){
    spectrum_vis >> line;
    baseline_vis >> baseline;
    values = split(line);
    if (options.Contains("FS") && values[0] >= 426 && values[0] < 574)  values[1] -= 0.4;
    if (options.Contains("FL") && values[0] >= 426 && values[0] < 574)  values[1] -= 0.2;
    if (options.Contains("HS") && values[0] >= 426 && values[0] < 574)  values[1] -= 0.5;
    if (options.Contains("HL") && values[0] >= 426 && values[0] < 574)  values[1] -= 0.2;
    if (values[3] == yval){
      w.push_back(double(values[0]));
      i.push_back(double(values[1]));
      yval = values[3];
      values = split(baseline);
      if (options.Contains("FS") && values[0] >= 426 && values[0] < 574)  values[1] -= 0.45;//-0.55
      if (options.Contains("FL") && values[0] >= 426 && values[0] < 574)  values[1] -= 0.6;//+0.2
      if (options.Contains("HM") && values[0] >= 426 && values[0] < 574)  values[1] -= 0.3;
      if (options.Contains("HL") && values[0] >= 426 && values[0] < 574)  values[1] -= 0.4;
      if (options.Contains("FL")) values[1] = (values[1] - 601.3)*25./27. + 601.3;
      /* if (options.Contains("H74") && values[0] >= 426 && values[0] < 574)  values[1] -= 0.;//+0.25 */
      b.push_back(double(values[1]));
    }
    else {
      wavelength_vis->push_back(w);
      w.clear();
      w.push_back(double(values[0]));
      intensity_vis->push_back(i);
      i.clear();
      i.push_back(double(values[1]));
      yval = values[3];
      base_vis->push_back(b);
      b.clear();
      values = split(baseline);
      b.push_back(double(values[1]));
    }
    acceptInput = !spectrum_vis.eof();
  }
  baseline_vis.close();
  spectrum_vis.close();
  baseline_ir.close();
  spectrum_ir.close();
  return 0;
}
