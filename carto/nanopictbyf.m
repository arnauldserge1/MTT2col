function nanopictbyf(filename)

files = dir(filename);
Nfiles = length(files);
dir_nano = 'nanopicts';

for ifile = 1:Nfiles
    filei = files(ifile).name ;
    outfile = [dir_nano '/' filei '.png'];
    if ~isempty(dir(outfile)), disp ([outfile ' already done']), continue, end
    
    %%%%%%%%%%%%%%%%%%%%%%
    nanopict(filei)%, ROI, check_roi_only, codage, dr_in, pixelsize, Nimg_max, Nsub_img);
    %%%%%%%%%%%%%%%%%%%%%%
    
    if isempty(dir(dir_nano)), mkdir(dir_nano), end
    disp(['saving ' outfile])
    saveas(gcf, outfile, 'png')
    if ifile<Nfiles, close gcf, end
end