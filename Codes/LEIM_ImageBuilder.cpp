#include "LEIM_GetImageData.cpp" //input args: string options,std::vector< std::vector<double> > *pointers
#include "LEIM_TrimBlackZones.cpp"
#include "LEIM_FindBoundary.cpp" //input args: std::vector< std::vector<double> > *pointer1, std::vector< std::vector<double> > *pointer2, bool isXaxis
#include "LEIM_RemoveCosmicRays.cpp"
using namespace std;
std::vector< std::vector<double> > frame1_x,frame1_i,frame2_x,frame2_i,master_frame;
std::vector< std::vector<double> > *frame1_x_ptr = &frame1_x;
std::vector< std::vector<double> > *frame1_i_ptr = &frame1_i;
std::vector< std::vector<double> > *frame2_x_ptr = &frame2_x;
std::vector< std::vector<double> > *frame2_i_ptr = &frame2_i;

//options include "H" and "F"
void ImageBuilder(TString options){
  gROOT->SetStyle("Modern");
  gROOT->ForceStyle();
  gStyle->SetPalette(52);//52 or 87
  gStyle->SetNumberContours(99);
  gStyle->SetOptStat(0);
  TFile *fout = new TFile();
  if (options == "H"){
    fout = TFile::Open("../rootfiles/HPKVUV4_SiPM_fullImage.root","RECREATE");
    GetImageData(frame1_x_ptr,frame1_i_ptr,frame2_x_ptr,frame2_i_ptr,options);
    //frame1_i = RemoveCosmicRays(frame1_i,7,7,25);
    //frame2_i = RemoveCosmicRays(frame2_i,7,7,25);
    TrimBlackZones(frame1_i_ptr);
    TrimBlackZones(frame2_i_ptr);
    int boundary = FindBoundary(frame2_i,frame1_i,true);
    TH2F *h = new TH2F("HPKVUV4_full","Full Image of Hamamatsu VUV4 SiPM;X [#mum];Y [#mum]",frame1_i[0].size(),0,5*frame1_i[0].size(),frame2_i.size()+boundary,0,5*(frame2_i.size()+boundary));
    for (int y = 0; y < frame2_i.size(); y++){
      for (int x = 0; x < frame2_i[y].size(); x++){
	if (x < 5) h->SetBinContent(x,400+boundary-y,0);
     	else h->SetBinContent(x,400+boundary-y,frame2_i[y][x-5]);
      }
    }
    int counter = -1;
    for (int y = frame1_i.size()-boundary-1; y < frame1_i.size(); y++){
      for (int x = 0; x < frame1_i[y].size(); x++){
     	h->SetBinContent(x,boundary-counter,frame1_i[y][x]);
      }
      counter += 1;
    }
    h->SetMaximum(9500);
    h->SetMinimum(-1);
    TCanvas *chpk = new TCanvas("chpk","cphk",900,800);
    chpk->cd();
    h->GetXaxis()->SetTickSize(0.01);
    h->GetXaxis()->SetLabelSize(0.02);
    h->GetYaxis()->SetTickSize(0.01);
    h->GetYaxis()->SetLabelSize(0.02);
    h->GetZaxis()->SetLabelSize(0.02);
    h->Draw("COLZ");
    chpk->SaveAs("~/Documents/Pictures/Plots/leim/HPKVUV4_SiPM_fullImage.pdf");
    fout->cd();
    h->Write();
    fout->Close();
  }
  else if (options=="F") {//fbk image
    fout = TFile::Open("../rootfiles/FBKVUV-HD3_SiPM_fullImage.root","RECREATE");
    std::vector< std::vector<double> > left_frame,right_frame;
    int boundary;
    double norm;
    int offset = 8;
    for (int i = 0; i < 4; i++){
      if (i == 0) {boundary = 188; norm = 1;}
      if (i == 1) {boundary = 221; norm = 1/0.8;}
      if (i == 2) {boundary = 221; norm = 1/0.8;}
      if (i == 3) {boundary = 221; norm = 1/0.85;}

      if (i == 0) GetImageData(frame1_x_ptr,frame1_i_ptr,frame2_x_ptr,frame2_i_ptr,"FTL");
      if (i == 1) GetImageData(frame1_x_ptr,frame1_i_ptr,frame2_x_ptr,frame2_i_ptr,"FUL");
      if (i == 2) GetImageData(frame1_x_ptr,frame1_i_ptr,frame2_x_ptr,frame2_i_ptr,"FML");
      if (i == 3) GetImageData(frame1_x_ptr,frame1_i_ptr,frame2_x_ptr,frame2_i_ptr,"FLL");
      TrimBlackZones(frame1_i_ptr);
      TrimBlackZones(frame2_i_ptr);
      for (int j = 0; j < frame1_i.size(); j++) {
	for (int k = 0; k < frame1_i[j].size(); k++) frame1_i[j][k] *= norm;
	left_frame.push_back(frame1_i[j]);
	int n = left_frame.back().size()-768+325;
	for (int k = 0; k < n; k++) left_frame[left_frame.size()-1].pop_back();
	if (j == boundary) break;
      }
      
      for (int j = left_frame.size()-1; j >= 0; j--){
	std::vector<double>::iterator it = left_frame[j].begin();	
	left_frame[j].insert(it,3,left_frame[j][0]);
      }
      
      if (i == 3){
	for (int j = 0; j < frame2_i.size(); j++) {
	  left_frame.push_back(frame2_i[j]);
	  for (int k = 0; k < frame2_i[j].size(); k++) frame2_i[j][k] *= norm;
	  int n = left_frame.back().size()-768+325;
	  for (int k = 0; k < n; k++) left_frame[left_frame.size()-1].pop_back();
	}
      }
      cout << "Loading right side..." << endl;
      if (i == 0) {norm = 1/0.85; GetImageData(frame1_x_ptr,frame1_i_ptr,frame2_x_ptr,frame2_i_ptr,"FTR");}
      if (i == 1) {norm = 1/0.85; GetImageData(frame1_x_ptr,frame1_i_ptr,frame2_x_ptr,frame2_i_ptr,"FUR");}
      if (i == 2) {norm = 1/0.95; GetImageData(frame1_x_ptr,frame1_i_ptr,frame2_x_ptr,frame2_i_ptr,"FMR");}
      if (i == 3) {norm = 1/0.9; GetImageData(frame1_x_ptr,frame1_i_ptr,frame2_x_ptr,frame2_i_ptr,"FLR");}
      TrimBlackZones(frame1_i_ptr);
      TrimBlackZones(frame2_i_ptr);
      std::vector<double> row;
      if (i == 0){
	for (int k = 363-325; k < frame1_i[0].size(); k++) row.push_back(frame1_i[0][k]);
	for (int k = 0; k < offset; k++) right_frame.push_back(row);
      }
      for (int j = 0; j < frame1_i.size(); j++) {
	row.clear();
	for (int k = 0; k < frame1_i[j].size(); k++) frame1_i[j][k] *= norm;
	int start = j > boundary - offset ? 363-325 - 3: 363-325;
	for (int k = start; k < frame1_i[j].size(); k++) row.push_back(frame1_i[j][k]);
        right_frame.push_back(row);
        if (j == boundary) break;
      }
      if (i == 3){
	norm = 1;
	for (int j = 0; j < frame2_i.size()-offset; j++) {
	  row.clear();
	  for (int k = 0; k < frame2_i[j].size(); k++) frame2_i[j][k] *= norm;
	  int start = j > boundary - offset ? 363-325: 363-325;
	  for (int k = start; k < frame2_i[j].size(); k++) row.push_back(frame2_i[j][k]);
          right_frame.push_back(row);
        }
      }
    }    
    TH2F *h = new TH2F("FBKVUVHD3_full","Full Image of FBK VUV-HD3 SiPM;X [#mum];Y [#mum]",1100,0,5500,left_frame.size(),0,5*left_frame.size());
    double multiplier_row,multiplier_column;
    for (int y = 0; y < left_frame.size(); y++){
      for (int x = 0; x < left_frame[y].size() + right_frame[y].size(); x++){
	if (x < left_frame[y].size()) h->SetBinContent(x,left_frame.size()-y,left_frame[y][x]);
	else h->SetBinContent(x,right_frame.size()-y,right_frame[y][x-left_frame[y].size()]);
      }
    }
    h->GetXaxis()->SetTickSize(0.01);
    h->GetXaxis()->SetLabelSize(0.02);
    h->GetYaxis()->SetTickSize(0.01);
    h->GetYaxis()->SetLabelSize(0.02);
    h->GetZaxis()->SetLabelSize(0.02);

    TCanvas *cfbk = new TCanvas("cfbk","cfbk",900,800);
    cfbk->cd();
    h->SetMaximum(2500);
    h->SetMinimum(0);
    h->Draw("COLZ");
    cfbk->SaveAs("~/Documents/Pictures/Plots/leim/FBKVUV-HD3_SiPM_fullImage.pdf");
    fout->cd();
    h->Write();
    fout->Close();
  }
  else {
    GetImageData(frame1_x_ptr,frame1_i_ptr,options);
    //TrimBlackZones(frame1_i_ptr);
    frame1_i = RemoveCosmicRays(frame1_i,13,25);
    TString name,title,filename,pdf_filename;
    if (options.Contains("Hamamatsu")){
      name = "HPKVUV4_zoomed";
      if (options.Contains("Closed")) name += "_slitClosed";
      title = "HPK VUV4 zoomed";
      if (options.Contains("Close")) title += "with slit closed";
      title += ";X [#mum];Y [#mum]";
      filename = "../rootfiles/HPKVUV4_SiPM_zoomed_old";
      if (options.Contains("Close")) filename += "_slitClosed";
      pdf_filename = filename;
      filename += ".root";
      pdf_filename += ".pdf";
    }
    if (options.Contains("FBK")){
      name = "FBKVUV-HD3_spec_withCRR";
      if (options.Contains("Closed")) name += "_slitClosed";
      title = "FBK VUV-HD3 Spectrum with with CRR";
      if (options.Contains("Close")) title += "with slit closed";
      title += ";X [#mum];Y [#mum]";
      filename = "../rootfiles/FBKVUV-HD3_SiPM_spec_withCRR";
      if (options.Contains("Close")) filename += "_slitClosed";
      pdf_filename = filename;
      filename += ".root";
      pdf_filename += ".pdf";
    }
    double pixel_to_um = 1./1.008;
    TH2F *h = new TH2F(name,title,frame1_i[0].size(),0,pixel_to_um*frame1_i[0].size(),frame1_i.size(),0,pixel_to_um*frame1_i.size());
    for (int i = 0; i < frame1_i.size(); i++){
      for (int j = 0; j < frame1_i[i].size(); j++){
	if (frame1_i[i][j] < 0) frame1_i[i][j] = 0;
	h->SetBinContent(j,400-i,frame1_i[i][j]);
      }
    }
    TCanvas *c = new TCanvas("c","c",1.5*frame1_i[0].size(),1.5*frame1_i.size());
    c->cd();
    h->SetMaximum(900);
    h->SetMinimum(560);
    h->GetXaxis()->SetTickSize(0.01);
    h->GetXaxis()->SetLabelSize(0.02);
    h->GetYaxis()->SetTickSize(0.01);
    h->GetYaxis()->SetLabelSize(0.02);
    h->GetZaxis()->SetLabelSize(0.02);
    h->Draw("COLZ");
    pdf_filename = "~/Documents/Pictures/Plots/leim/PDFs/"+pdf_filename;
    c->SaveAs(pdf_filename);
    TFile *fout = new TFile(filename,"RECREATE");
    fout->cd();
    h->Write();
    fout->Close();
  }
}
