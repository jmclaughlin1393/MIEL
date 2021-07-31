int FindBoundary(std::vector< std::vector<double> > frame1, std::vector< std::vector<double> > frame2, bool row_symmetry){
  int boundary_index,n_terms,n_xpixels;
  double correlation,min_correlation;
  n_xpixels = frame1[0].size();
  boundary_index = 0;
  min_correlation = 1.e+50;
  if (row_symmetry){
    //Assume frame1 is the top and frame2 is the bottom
    for (int y1 = 399; y1 >= 0; y1--){
      correlation = 0;
      n_terms = 0;
      for (int y2 = 0; y2 <= 399 - y1; y2++){
	for (int x = 0; x < n_xpixels; x++){
	  correlation += TMath::Power(frame1[y1+y2][x] - frame2[y2][x],2.0);
	  n_terms += 1;
	}
	if (correlation/double(n_terms) < min_correlation){
	  min_correlation = correlation/double(n_terms);
	  boundary_index = y1;
	}
      }
    }
  }
  else{
    //Assume frame1 is the left and frame2 is the right
    for	(int x1 = n_xpixels-1; x1 >= 0; x1--){
      correlation = 0;
      n_terms =	0;
      for (int x2 = 0; x2 <= n_xpixels - x1; x2++){
        for (int y = 0; y < 400; y++){
          correlation += TMath::Power(frame1[y][x1+x2] - frame2[y][x2],2.0);
          n_terms += 1;
	}
        if (correlation/double(n_terms) < min_correlation){
          min_correlation = correlation/double(n_terms);
          boundary_index = x1;
	}
      } 
    }
  }
  return boundary_index;
}
