% nanopict(file, ROI, check_roi_only, codage, dr_in, pixelsize, Nimg_max, Nsub_img)
check_roi_only = 0; codage = 'alpha'; magnification = 100*1.5; pixelsize = 16/magnification; dr_in = 0.1;
file = 'cel1.tif';
ROI = [130 405 281-130   507-405];
nanopict(file, ROI, check_roi_only, codage, dr_in, pixelsize, 10, 10)
nanopict(file, ROI, check_roi_only, codage, dr_in, pixelsize, 100, 100)
nanopict(file, ROI, check_roi_only, codage, dr_in, pixelsize, 1000, 200)