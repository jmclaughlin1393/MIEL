std::vector< std::vector<double> > Stitch(std::vector<double> vis_wavelength,
					  std::vector<double> vis_intensity,
					  std::vector<double> vis_intensity_err_lo,
					  std::vector<double> vis_intensity_err_up,
					  std::vector<double> ir_wavelength,
					  std::vector<double> ir_intensity,
					  std::vector<double> ir_intensity_err_lo,
					  std::vector<double> ir_intensity_err_up){
  TGraph *Gvis = new TGraph(vis_wavelength.size());
  TGraph *Gir = new TGraph(ir_wavelength.size());
  for (int i = 0; i < vis_wavelength.size(); i++) Gvis->SetPoint(i,vis_wavelength[i],vis_intensity[i]);
  for (int i = 0; i < ir_wavelength.size(); i++) Gir->SetPoint(i,ir_wavelength[i],ir_intensity[i]);
  double ratio = 0;
  for (double x = 550; x < 620; x += 0.5){
    ratio += Gir->Eval(x)/Gvis->Eval(x);
  }
  ratio /= 140.;
  std::vector< std::vector<double> > stitchedSpectrum;
  std::vector<double> wavelength,intensity,error_up,error_lo;
  for (int i = 0; i < vis_wavelength.size(); i++){
    if (vis_wavelength[i] > 550) continue;
    wavelength.push_back(vis_wavelength[i]);
    intensity.push_back(ratio*vis_intensity[i]);
    error_up.push_back(vis_intensity_err_up[i]);
    error_lo.push_back(vis_intensity_err_lo[i]);
  }
  for (int i = 0; i < ir_wavelength.size(); i++){
    if (ir_wavelength[i] <= 550) continue;
    wavelength.push_back(ir_wavelength[i]);
    intensity.push_back(ir_intensity[i]);
    error_up.push_back(ir_intensity_err_up[i]);
    error_lo.push_back(ir_intensity_err_lo[i]);
  }
  stitchedSpectrum.push_back(wavelength);
  stitchedSpectrum.push_back(intensity);
  stitchedSpectrum.push_back(error_up);
  stitchedSpectrum.push_back(error_lo);
  delete Gvis;
  delete Gir;
  return stitchedSpectrum;
}
