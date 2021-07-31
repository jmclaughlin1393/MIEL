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

int GetImageData(std::vector< std::vector<double> >* frame1_x, std::vector< std::vector<double> >* frame1_i,
		 std::vector< std::vector<double> >* frame2_x, std::vector< std::vector<double> >* frame2_i,
		 TString options){
  
  frame1_x->clear();
  frame1_i->clear();
  frame2_x->clear();
  frame2_i->clear();
  
  ifstream frame1,frame2;
  std::vector<float> values;
  std::vector<double> x,i;
  double yval;
  bool acceptInput;
  TString line;
  
  if (options.Contains("H")){
    cout << "Loading Hamamatsu image data" << endl;    
    frame1.open("../data/DarkNoisePaper/images/newHamamatsuVUV4_image_lowerHalf_4x_1200s.csv");
    frame2.open("../data/DarkNoisePaper/images/newHamamatsuVUV4_image_upperHalf_4x_1200s.csv");
  }
  else if (options.Contains("F") && options.Contains("TL")){
    cout << "Loading FBK top-left" << endl;
    frame1.open("../data/DarkNoisePaper/images/newFBKVUV-HD3_image_topLeft_4x_300s.csv");
    frame2.open("../data/DarkNoisePaper/images/newFBKVUV-HD3_image_upperLeft_4x_300s.csv");
  } else if (options.Contains("F") && options.Contains("UL")){
    cout << "Loading FBK upper-left" << endl;
    frame1.open("../data/DarkNoisePaper/images/newFBKVUV-HD3_image_upperLeft_4x_300s.csv");
    frame2.open("../data/DarkNoisePaper/images/newFBKVUV-HD3_image_middleLeft_4x_300s.csv");
  } else if (options.Contains("F") && options.Contains("ML")){
    cout << "Loading FBK middle-left" << endl;
    frame1.open("../data/DarkNoisePaper/images/newFBKVUV-HD3_image_middleLeft_4x_300s.csv");
    frame2.open("../data/DarkNoisePaper/images/newFBKVUV-HD3_image_lowerLeft_4x_300s.csv");
  } else if (options.Contains("F") && options.Contains("LL")){
    cout << "Loading FBK lower-left" << endl;
    frame1.open("../data/DarkNoisePaper/images/newFBKVUV-HD3_image_lowerLeft_4x_300s.csv");
    frame2.open("../data/DarkNoisePaper/images/newFBKVUV-HD3_image_bottomLeft_4x_300s.csv");
  } else if (options.Contains("F") && options.Contains("TC")){
    cout << "Loading FBK top-centre" << endl;
    frame1.open("../data/DarkNoisePaper/images/newFBKVUV-HD3_image_topCentre_4x_300s.csv");
    frame2.open("../data/DarkNoisePaper/images/newFBKVUV-HD3_image_upperCentre_4x_300s.csv");
  } else if (options.Contains("F") && options.Contains("UC")){
    cout << "Loading FBK upper-centre" << endl;
    frame1.open("../data/DarkNoisePaper/images/newFBKVUV-HD3_image_upperCentre_4x_300s.csv");
    frame2.open("../data/DarkNoisePaper/images/newFBKVUV-HD3_image_middleCentre_4x_300s.csv");
  } else if (options.Contains("F") && options.Contains("MC")){
    cout << "Loading FBK middle-centre" << endl;
    frame1.open("../data/DarkNoisePaper/images/newFBKVUV-HD3_image_middleCentre_4x_300s.csv");
    frame2.open("../data/DarkNoisePaper/images/newFBKVUV-HD3_image_lowerCentre_4x_300s.csv");
  } else if (options.Contains("F") && options.Contains("LC")){
    cout << "Loading FBK lower-centre" << endl;
    frame1.open("../data/DarkNoisePaper/images/newFBKVUV-HD3_image_lowerCentre_4x_300s.csv");
    frame2.open("../data/DarkNoisePaper/images/newFBKVUV-HD3_image_bottomCentre_4x_300s.csv");
  } else if (options.Contains("F") && options.Contains("TR")){
    cout << "Loading FBK top-right" << endl;
    frame1.open("../data/DarkNoisePaper/images/newFBKVUV-HD3_image_topRight_4x_300s.csv");
    frame2.open("../data/DarkNoisePaper/images/newFBKVUV-HD3_image_upperRight_4x_300s.csv");
  } else if (options.Contains("F") && options.Contains("UR")){
    cout << "Loading FBK upper-right" << endl;
    frame1.open("../data/DarkNoisePaper/images/newFBKVUV-HD3_image_upperRight_4x_300s.csv");
    frame2.open("../data/DarkNoisePaper/images/newFBKVUV-HD3_image_middleRight_4x_300s.csv");
  } else if (options.Contains("F") && options.Contains("MR")){
    cout << "Loading FBK middle-right" << endl;
    frame1.open("../data/DarkNoisePaper/images/newFBKVUV-HD3_image_middleRight_4x_300s.csv");
    frame2.open("../data/DarkNoisePaper/images/newFBKVUV-HD3_image_lowerRight_4x_300s.csv");
  } else if (options.Contains("F") && options.Contains("LR")){
    cout << "Loading FBK lower-right" << endl;
    frame1.open("../data/DarkNoisePaper/images/newFBKVUV-HD3_image_lowerRight_4x_300s.csv");
    frame2.open("../data/DarkNoisePaper/images/newFBKVUV-HD3_image_bottomRight_4x_300s.csv");
  }  
  if (!frame1.is_open()){ cout << "failed to open frame 1." << endl; exit(1); }
  else if (!frame2.is_open()){ cout << "failed to open frame 2." << endl; exit(1); }
  yval = 0;
  frame1 >> line;
  acceptInput = !frame1.eof();
  while (acceptInput){
    frame1 >> line;
    values = split(line);
    if (values[3] == yval){
      i.push_back(double(values[1]));
      x.push_back(double(values[2]));
      yval = values[3];
    }
    else {
      frame1_i->push_back(i);
      i.clear();
      i.push_back(double(values[1]));
      frame1_x->push_back(x);
      x.clear();
      x.push_back(double(values[2]));
      yval = values[3];
    }
    acceptInput = !frame1.eof();
  }

  i.clear(); x.clear();
  yval = 0;
  frame2 >> line;
  acceptInput = !frame2.eof();
  while (acceptInput){
    frame2 >> line;
    values = split(line);
    if (values[3] == yval){
      i.push_back(double(values[1]));
      x.push_back(double(values[2]));
      yval = values[3];
    }
    else {
      frame2_i->push_back(i);
      i.clear();
      i.push_back(double(values[1]));
      frame2_x->push_back(x);
      x.clear();
      x.push_back(double(values[2]));
      yval = values[3];
    }
    acceptInput = !frame2.eof();
  }
  frame1.close();
  frame2.close();
  return 0;
}

int GetImageData(std::vector< std::vector<double> >* frame_x, std::vector< std::vector<double> >* frame_i,
                 TString file){

  frame_x->clear();
  frame_i->clear();

  
  ifstream frame;
  std::vector<float> values;
  std::vector<double> x,i;
  double yval;
  bool acceptInput;
  TString line,dir,path;
  dir = "../data/DarkNoisePaper/images/";
  path = dir+file;
  frame.open(path);
  
  acceptInput = !frame.eof();
  yval = 0;
  frame >> line;
  while (acceptInput){
    frame >> line;
    values = split(line);
    if (values[3] == yval){
      i.push_back(double(values[1]));
      x.push_back(double(values[2]));
      yval = values[3];
    }
    else {
      frame_i->push_back(i);
      i.clear();
      i.push_back(double(values[1]));
      frame_x->push_back(x);
      x.clear();
      x.push_back(double(values[2]));
      yval = values[3];
    }
    acceptInput = !frame.eof();
  }

  frame.close();
  return 0;
}
