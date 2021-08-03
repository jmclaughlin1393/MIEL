std::vector< std::vector<double> > SumRescaleRebin(std::vector< std::vector<double> > wavelength,
						   std::vector< std::vector<double> > intensity,
						   std::vector< std::vector<double> > intensity_err_up,
						   std::vector< std::vector<double> > intensity_err_lo,
						   double min,
						   double max,
						   double binWidth,
						   double charge,
						   double charge_err,
						   double integral_scale,
						   double integral_scale_err,
						   TString option){
  std::vector<double> wavelength_rebinned,intensity_rebinned,intensity_rebinned_err_up,intensity_rebinned_err_lo;
  std::vector< std::vector<double> > rebinned_data;
  unsigned int index;
  bool vis_mode;
  double middle_wavelength = wavelength[0][0] + wavelength[0][wavelength[0].size() - 1];
  middle_wavelength /= 2.0;
  if (middle_wavelength < 600) vis_mode = true;
  else vis_mode = false;
  double underflow_out,underflow_in,overflow_out,overflow_in,sum,index_bin_width_upper,index_bin_width_lower;
  double underflow_out_err_up,underflow_in_err_up,overflow_out_err_up,overflow_in_err_up,sum_err_up;
  double underflow_out_err_lo,underflow_in_err_lo,overflow_out_err_lo,overflow_in_err_lo,sum_err_lo;
  for (double l = min; l <= max; l += binWidth){
    wavelength_rebinned.push_back(l + 0.5*binWidth);
    intensity_rebinned.push_back(0);
    intensity_rebinned_err_up.push_back(0);
    intensity_rebinned_err_lo.push_back(0);
  }

  //Summing and rebinning
  for (unsigned int i = 0; i < wavelength.size(); i++){
    
    index = 0;
    while (wavelength[i][index] < wavelength_rebinned[0]-0.5*binWidth) index += 1;
    
    for (unsigned int j = 0; j < wavelength_rebinned.size(); j++){
      sum = underflow_out = overflow_out = underflow_in = overflow_in = 0;
      sum_err_up = underflow_out_err_up = overflow_out_err_up = underflow_in_err_up = overflow_in_err_up = 0;
      sum_err_lo = underflow_out_err_lo = overflow_out_err_lo = underflow_in_err_lo = overflow_in_err_lo = 0;
      while (wavelength[i][index] >= wavelength_rebinned[j]-0.5*binWidth &&
	     wavelength[i][index] < wavelength_rebinned[j]+0.5*binWidth){
	
	if (index == 0) index_bin_width_lower = index_bin_width_upper = 0.5*(wavelength[i][index+1] - wavelength[i][index]);
	else if (index == wavelength[i].size()-1) index_bin_width_upper = index_bin_width_lower = 0.5*(wavelength[i][index] - wavelength[i][index-1]);
	else {index_bin_width_upper = 0.5*(wavelength[i][index+1] - wavelength[i][index]); index_bin_width_lower = 0.5*(wavelength[i][index] - wavelength[i][index-1]);}

	if (wavelength[i][index] - index_bin_width_lower < wavelength_rebinned[j]-0.5*binWidth){
	  underflow_out = intensity[i][index];
	  underflow_out_err_up = TMath::Power(intensity_err_up[i][index],2.0);
	  underflow_out_err_lo = TMath::Power(intensity_err_lo[i][index],2.0);
	  underflow_out *= (wavelength_rebinned[j] - 0.5*binWidth) - (wavelength[i][index] - index_bin_width_lower);
	  underflow_out_err_up *= (wavelength_rebinned[j] - 0.5*binWidth) - (wavelength[i][index] - index_bin_width_lower);
	  underflow_out_err_lo *= (wavelength_rebinned[j] - 0.5*binWidth) - (wavelength[i][index] - index_bin_width_lower);
	  underflow_out /= index_bin_width_upper + index_bin_width_lower;
	  underflow_out_err_up /= index_bin_width_upper + index_bin_width_lower;
	  underflow_out_err_lo /= index_bin_width_upper + index_bin_width_lower;
	  //underflow_out_err_up *= underflow_out_err_up;
	  //underflow_out_err_lo *= underflow_out_err_lo;
	}
	else if (wavelength[i][index] + index_bin_width_upper > wavelength_rebinned[j]+0.5*binWidth){
	  overflow_out = intensity[i][index];
	  overflow_out_err_up = TMath::Power(intensity_err_up[i][index],2.0);
	  overflow_out_err_lo = TMath::Power(intensity_err_lo[i][index],2.0);
	  overflow_out *= (wavelength[i][index] + index_bin_width_upper) - (wavelength_rebinned[j] + 0.5*binWidth);
	  overflow_out_err_up *= (wavelength[i][index] + index_bin_width_upper) - (wavelength_rebinned[j] + 0.5*binWidth);
	  overflow_out_err_lo *= (wavelength[i][index] + index_bin_width_upper) - (wavelength_rebinned[j] + 0.5*binWidth);
	  overflow_out /= index_bin_width_upper + index_bin_width_lower;
	  overflow_out_err_up /= index_bin_width_upper + index_bin_width_lower;
	  overflow_out_err_lo /= index_bin_width_upper + index_bin_width_lower;
	  //overflow_out_err_up *= overflow_out_err_up;
	  //overflow_out_err_lo *= overflow_out_err_lo;
	}
	else if (index > 1 && wavelength[i][index-1] < wavelength_rebinned[j]-0.5*binWidth && wavelength[i][index] - index_bin_width_lower > wavelength_rebinned[j]-0.5*binWidth){
	  underflow_in = intensity[i][index-1];
	  underflow_in_err_up = TMath::Power(intensity_err_up[i][index-1],2.0);
	  underflow_in_err_lo = TMath::Power(intensity_err_lo[i][index-1],2.0);
	  underflow_in *= (wavelength[i][index] - index_bin_width_lower) - (wavelength_rebinned[j]-0.5*binWidth);
	  underflow_in_err_up *= (wavelength[i][index] - index_bin_width_lower) - (wavelength_rebinned[j]-0.5*binWidth);
	  underflow_in_err_lo *= (wavelength[i][index] - index_bin_width_lower) - (wavelength_rebinned[j]-0.5*binWidth);
	  underflow_in /= 0.5*(wavelength[i][index] - wavelength[i][index-1]) + 0.5*(wavelength[i][index-1] - wavelength[i][index-2]);
	  underflow_in_err_up /= 0.5*(wavelength[i][index] - wavelength[i][index-1]) + 0.5*(wavelength[i][index-1] - wavelength[i][index-2]);
	  underflow_in_err_lo /= 0.5*(wavelength[i][index] - wavelength[i][index-1]) + 0.5*(wavelength[i][index-1] - wavelength[i][index-2]);
	  //underflow_in_err_up *= underflow_in_err_up;
	  //underflow_in_err_lo *= underflow_in_err_lo;
	}
	else if (index < wavelength[i].size()-2 &&
	    wavelength[i][index+1] > wavelength_rebinned[j]+0.5*binWidth &&
	    wavelength[i][index] + index_bin_width_upper < wavelength_rebinned[j]+0.5*binWidth){
	  overflow_in = intensity[i][index+1];
	  overflow_in_err_up = TMath::Power(intensity_err_up[i][index+1],2.0);
	  overflow_in_err_lo = TMath::Power(intensity_err_lo[i][index+1],2.0);
	  overflow_in *= (wavelength_rebinned[j]+0.5*binWidth) - (wavelength[i][index] + index_bin_width_upper);
	  overflow_in_err_up *= (wavelength_rebinned[j]+0.5*binWidth) - (wavelength[i][index] + index_bin_width_upper);
	  overflow_in_err_lo *= (wavelength_rebinned[j]+0.5*binWidth) - (wavelength[i][index] + index_bin_width_upper);
	  overflow_in /= 0.5*(wavelength[i][index+1] - wavelength[i][index]) + 0.5*(wavelength[i][index+2] - wavelength[i][index+1]);
	  overflow_in_err_up /= 0.5*(wavelength[i][index+1] - wavelength[i][index]) + 0.5*(wavelength[i][index+2] - wavelength[i][index+1]);
	  overflow_in_err_lo /= 0.5*(wavelength[i][index+1] - wavelength[i][index]) + 0.5*(wavelength[i][index+2] - wavelength[i][index+1]);
	  //overflow_in_err_up *= overflow_in_err_up;
	  //overflow_in_err_lo *= overflow_in_err_lo;
	}
	sum += intensity[i][index];
	sum_err_up += TMath::Power(intensity_err_up[i][index],2.0);
	sum_err_lo += TMath::Power(intensity_err_lo[i][index],2.0);
	index += 1;
      }
      intensity_rebinned[j] += sum + (underflow_in + overflow_in) - (underflow_out + overflow_out);
      intensity_rebinned_err_up[j] += sum_err_up + underflow_in_err_up + overflow_in_err_up + underflow_out_err_up + overflow_out_err_up;
      intensity_rebinned_err_lo[j] += sum_err_lo + underflow_in_err_lo + overflow_in_err_lo + underflow_out_err_lo + overflow_out_err_lo;
      //if (j == int(intensity_rebinned.size()/2)) cout << "UP " << TMath::Sqrt(intensity_rebinned_err_up[j]) << " DOWN " << TMath::Sqrt(intensity_rebinned_err_lo[j]) << endl;
    }
  }
  cout << endl;
  for (int i = 0; i < intensity_rebinned_err_up.size(); i++){
    intensity_rebinned_err_up[i] = TMath::Sqrt(intensity_rebinned_err_up[i]);
    intensity_rebinned_err_lo[i] = TMath::Sqrt(intensity_rebinned_err_lo[i]);
  }
  
  //Rescaling
  cout << " - Rescaling with the following factors:" << endl;
  cout << "      integral scale = " << integral_scale << " +/- " << integral_scale_err << endl; 
  cout << "      charge         = " << charge << " +/- " <<  charge_err << " (Coul.)" << endl;
  TFile *fcorr = new TFile();
  if (vis_mode && option.Contains("F")) fcorr = TFile::Open("../Rootfiles/CorrectionVsWavelength_VISmode_VUVHD3.root");
  else if (!vis_mode && option.Contains("F")) fcorr = TFile::Open("../Rootfiles/CorrectionVsWavelength_NIRmode_VUVHD3.root");
  else if (vis_mode && option.Contains("H")) fcorr = TFile::Open("../Rootfiles/CorrectionVsWavelength_VISmode_VUV4.root");
  else if (!vis_mode && option.Contains("H")) fcorr = TFile::Open("../Rootfiles/CorrectionVsWavelength_NIRmode_VUV4.root");
  TGraph *correction = (TGraph*)fcorr->Get("G_correction");
  double conversion = (integral_scale*1.602e-19)/charge;
  double conversion_err = TMath::Power((integral_scale_err*1.602e-19)/charge,2.0);
  conversion_err += TMath::Power((integral_scale*charge_err*1.602e-19)/(charge*charge),2.0);
  conversion_err = TMath::Sqrt(conversion_err);
  cout << "      conversion   = " << conversion << " +/- " << conversion_err << " ( (Photons/electron)/count )"  << endl;
  for (unsigned int i = 0; i < intensity_rebinned.size(); i++){
    intensity_rebinned_err_up[i] = TMath::Sqrt(TMath::Power(intensity_rebinned_err_up[i]*conversion,2.0) + TMath::Power(intensity_rebinned[i]*conversion_err,2.0)) / correction->Eval(wavelength_rebinned[i]);
    //intensity_rebinned_err_up[i] = TMath::Sqrt(TMath::Power(intensity_rebinned_err_up[i]*conversion,2.0)) * correction->Eval(wavelength_rebinned[i]);
    intensity_rebinned_err_lo[i] = TMath::Sqrt(TMath::Power(intensity_rebinned_err_lo[i]*conversion,2.0) + TMath::Power(intensity_rebinned[i]*conversion_err,2.0)) / correction->Eval(wavelength_rebinned[i]);
    //intensity_rebinned_err_lo[i] = TMath::Sqrt(TMath::Power(intensity_rebinned_err_lo[i]*conversion,2.0)) * correction->Eval(wavelength_rebinned[i]);
    intensity_rebinned[i] *= conversion;
    intensity_rebinned[i] /= correction->Eval(wavelength_rebinned[i]);
    //cout << intensity_rebinned_err_up[i]/intensity_rebinned[i] << endl;
  }
  rebinned_data.push_back(wavelength_rebinned);
  rebinned_data.push_back(intensity_rebinned);
  rebinned_data.push_back(intensity_rebinned_err_lo);
  rebinned_data.push_back(intensity_rebinned_err_up);
  return rebinned_data;
}
