#include "TSpectrum.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"

TSpectrum *s = new TSpectrum();

TH2F *RemoveCosmicRays(TH2F *image,int iterations,double threshold){
  TH2F *image_copy = new TH2F();
  image_copy = (TH2F*)image->Clone("image_copy");
  image_copy->SetDirectory(0);
  TH1F *row;
  TH1* row_bkg;
  for (int i = 1; i <= image_copy->GetNbinsY(); i++){
    row = (TH1F*)image_copy->ProjectionX("",i,i,"");
    row_bkg = s->Background(row,iterations);
    for (int j = 1; j <= image->GetNbinsX(); j++) image_copy->SetBinContent(j,i,row_bkg->GetBinContent(j));
  }
  for (int i = 1; i <= image->GetNbinsY(); i++){
    for (int j = 1; j <= image->GetNbinsX(); j++){
      if (TMath::Abs(image->GetBinContent(j,i) - image_copy->GetBinContent(j,i)) < threshold){
        image_copy->SetBinContent(j,i,image->GetBinContent(j,i));
      }
    }
  }
  return image_copy;
}

std::vector< std::vector<double> > RemoveCosmicRays(std::vector< std::vector<double> > image, int iterations, double threshold){
  Int_t direction = TSpectrum::kBackDecreasingWindow;
  Int_t filterOrder = TSpectrum::kBackOrder2;
  Bool_t smoothing = kTRUE;
  Int_t smoothWindow = TSpectrum::kBackSmoothing3;
  Bool_t compton = kFALSE;
  Int_t ssize = image[0].size();
  Double_t *row = new Double_t[ssize];
  for (unsigned int i = 0; i < image.size(); i++){
    for (unsigned int j = 0; j < image[i].size(); j++) row[j] = image[i][j];
    s->Background(row,ssize,iterations,direction,filterOrder,smoothing,smoothWindow,compton);
    for (unsigned int j = 0; j < image[i].size(); j++){
      if (TMath::Abs(image[i][j] - row[j]) < threshold) continue;
      image[i][j] = row[j];
    }
  }
  return image;
}

std::vector< std::vector<double> > RemoveCosmicRays(std::vector< std::vector<double> > image, int local_x, int local_y, double threshold){
  if (local_x%2 == 0) local_x += 1;
  if (local_y%2 == 0) local_y += 1;
  int start_x = local_x/2 + 1;
  int start_y = local_y/2 + 1;
  int median_index;
  double median;
  std::vector<double> list_of_counts;
  for (int y = start_y; y < image.size() - local_y/2; y++){
    for (int x = start_x; x < image[y].size() - local_x/2; x++){
      list_of_counts.clear();
      for (int ly = y - local_y/2; ly <= y + local_y/2; ly++){
	for (int lx = x - local_x/2; lx <= x + local_x/2; lx++){
	  list_of_counts.push_back(image[ly][lx]);
	}
      }
      std::sort(list_of_counts.begin(),list_of_counts.end());
      median_index = list_of_counts.size()/2;
      median = list_of_counts[median_index - 1];
      if (image[start_y][start_x] > median + threshold) image[start_y][start_x] = median;
    }
  }
  return image;
}
