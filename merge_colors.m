function merge_colors(filename, max_only)

% function merge_colors(filename)
% merge images from 2 color acquisition, expecting (green and red, usually)
% images to be initially side by side.
% as defined by Metamorph for 2 cameras
% Hence 'green' = left & 'red' = right
% save RGB movie (with blue channel empty..) and max. projection
%
% filename = '*.tif' by default
%
% see also MTT_2_colors, reg_2_colors

if nargin < 1, filename = '*.tif'; end
if nargin < 2, max_only = 1; end

files = dir2(filename);
Nf = length(files);
if Nf == 0, disp('niente??'), return, end

if contains(cd, 'split') || Nf >= 20,  max_only  = 1; end % Nf < 20 added to limit for timelapse (~100 files => ~100 Go!!, duplicates folder size) 23/5/2017
 
merge_dir = 'merged files';
max_dir = 'max';
if ~isdir(merge_dir) && ~max_only, mkdir(merge_dir), end 
if ~isdir(max_dir), mkdir(max_dir), end

do_reg = 0;
if ~isempty(dir(['dic' filesep 'reg' filesep 'mean_tform.mat']))
    do_reg = 1;
    mean_tform = importdata(['dic' filesep 'reg' filesep 'mean_tform.mat']);
end

for nf = 1:Nf
    file = files(nf).name;
    
    file_merg = [merge_dir filesep file(1:end-4) '_merg'];
    file_max = [max_dir filesep file(1:end-4) '_max'];
    if ~isempty(dir([file_max '.tif'])), disp([file_max ' already done']), continue, end
    
    img1 = imread(file, 1) ;
    Nim = length(imfinfo(file));
    
    %% split green/red
    h = size(img1, 1); % 512, a priori
    middle = size(img1, 2)/2; % 512, a priori
    
    if ~max_only, im_merg = uint16(zeros([h, middle, 3])); end
    im_max = uint16(zeros([h, middle, 3]));
    
    if do_reg == 1, R = imref2d([h, middle]); end
    
    fprintf('merging %s    ', file)
    
    for nt = 1:Nim
        im = imread(file, nt);
        imG = im(:, 1:middle);
        imR = im(:, middle+1:end);
        
        if do_reg, imR = imwarp(imR, mean_tform, 'OutputView', R); end
        
        %% merging corrected red & green
        if ~max_only
            im_merg(:, :, 1) = imR;
            im_merg(:, :, 2) = imG;
        end
        
        im_max(:, :, 1) = max(im_max(:, :, 1), imR);
        im_max(:, :, 2) = max(im_max(:, :, 2), imG);
        
        %% save merged images
        if ~max_only 
            if nt == 1, writemd = 'overwrite'; else, writemd = 'append'; end
            imwrite(im_merg, [file_merg '_tmp.tif'], 'tiff', 'Compression', 'none', 'WriteMode', writemd)
        end
        
        fprintf('\b\b\b\b%4i', nt)
    end
    fprintf('\r')
    
    sat = 0.002 ; % saturation 0.2% min-max du contraste
    im_max_sat = imadjust(im_max, stretchlim(im_max, sat), [0 1]); % imshow(im_max_sat), pause(.1)
    imwrite(im_max_sat, [file_max '.tif'], 'tiff', 'Compression', 'none')
    
    if ~contains(cd, 'split') && Nf < 20 && ~max_only, movefile([file_merg '_tmp.tif'], [file_merg '.tif']), end % if fileattrib().UserWrite == 1??
end
%%%