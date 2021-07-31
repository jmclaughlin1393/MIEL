#include "TMath.h"
std::vector<double> SubtractBaseline(std::vector<double> signal, std::vector<double> baseline){
  for (unsigned int i = 0; i < signal.size(); i++){
    if (i >= baseline.size()) continue;    
    signal[i] -= baseline[i];
    signal[i] *= 0.7;
  }
  return signal;
}

std::vector<double> GetError(std::vector<double> signal, double baseline_noise){
  std::vector<double> error_vec;
  double error;
  for (unsigned int i = 0; i < signal.size(); i++){
    error = signal[i] >= 0 ? TMath::Sqrt(signal[i] + 2*baseline_noise) : TMath::Sqrt(2*baseline_noise);
    error_vec.push_back(error);
  }
  return error_vec;
}
