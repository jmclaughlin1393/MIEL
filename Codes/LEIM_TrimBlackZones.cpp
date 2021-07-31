int TrimBlackZones(std::vector<std::vector<double> >* frame){
  std::vector<std::vector<double> > frame_copy = *frame;
  std::vector<std::vector<double> > frame_trim;
  std::vector<double> row_trim;
  double bkg = 0;//585
  for (int i = 0; i < frame_copy.size(); i++){
    row_trim.clear();
    bkg = 0;
    for (int j = 0; j < 325; j++){
      bkg += frame_copy[i][j];
    }
    bkg/=325.;
    for (int j = 0; j < frame_copy[i].size(); j++){
      if (j > 325 && j < 1025) row_trim.push_back(frame_copy[i][j] - bkg);
    }
    frame_trim.push_back(row_trim);
  }
  *frame = frame_trim;
  return 0;
}
