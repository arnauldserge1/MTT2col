function [im_thr, threshold_out] = threshold_img(data_in, do_plot, med_filt_size, size_dil_erod, do_fill, threshold_in, relative_threshold, old_fashion, flatten_background, do_save, area_min, area_max)

%% function [im_thr, threshold_out] = threshold_img(data_in, do_plot, med_filt_size, size_dil_erod, do_fill, threshold_in, relative_threshold, old_fashion, flatten_background, do_save, area_min, area_max)
% data_in: image or image filename
% threshold_in = -1 => automatic computation,  with graythresh ; -2: with opthr
% def: im_thr = threshold_img(data_in, 0, 3, 8, 1, -1, 0, 0, 0, 0, 0, inf)

% Caution, default values used for detect_cell_contact


if nargin < 1, data_in = uigetfile('*.tif;*.stk;*.lsm*.jpg', 'select image file'); end
if nargin < 2, do_plot = 0; end
if nargin < 3, med_filt_size = 3; end % odd size is better
if nargin < 4, size_dil_erod = 8; end
if nargin < 5, do_fill = 1; end
if nargin < 6, threshold_in = -1; end % => -1: automatic computation, with graythresh (Otsu method), -2: opthr, -3: maxentropie, -4: Moment
if nargin < 7, relative_threshold = 0; end
if nargin < 8, old_fashion = 0; end % old way: dilate8/erode8; new: dil4/er8/dil4; new2: erode/dil (old_fashion = 2!!)
if nargin < 9, flatten_background = 0; end % 15?
if nargin < 10, do_save = 0; end
if nargin < 11, area_min = 0; end
if nargin < 12, area_max = inf; end


med_filt_size = med_filt_size * [1 1];
warning off images:imshow:magnificationMustBeFitForDockedFigure

if ischar(data_in), filename = data_in; im = imread(filename); else, im = data_in; end
if (do_plot == 1)
    figure('WindowStyle','docked')
    subplot(221) % if (do_plot == 1), subplot(221), elseif (do_plot == 2), subplot(121), end
    imagesc(im), axis image off, colormap(gray)
    if do_save, imwrite(im, 'threshold_img_illustration/im.tif'), end
    if ischar(data_in), title(filename); else, title('raw'), end
    drawnow, pause(.1)
end % figure('windowstyle','docked'), 

%% flatten_background
if flatten_background > 0
    background_im = imopen(im, strel('disk', flatten_background));
    im = im - background_im + mean(background_im(:));
end

%% filt
if med_filt_size > 0
    im_filt = medfilt2(im, med_filt_size);
    if (do_plot == 1), subplot(222), imagesc(im_filt), axis image off, title(sprintf('filtered over %ix%i pxl', med_filt_size)), pause(.1), end
    if do_save, imwrite(im_filt, 'threshold_img_illustration/im_filt.tif'), end
else
    im_filt = im;
end

%% threshold
if (threshold_in < 1)% || relative_threshold
    threshold_auto = graythresh(im_filt) * max(im_filt(:)); % Otsu method : Otsu, N., "A Threshold Selection Method from Gray-Level Histograms," IEEE Transactions on Systems, Man, and Cybernetics, Vol. 9, No. 1, 1979, pp. 62-66.
    if (threshold_in == -2), threshold_auto = opthr(double(im)); end % added 6/7/2016, code from Felix Toran Marti
    if (threshold_in == -3), threshold_auto = maxentropie(im); end % author F.Gargouri
    if (threshold_in == -4), threshold_auto = momentsthresh(uint8(im)); end % Author: Daniel Martín REF: Tsai W. (1985) Moment-preserving thresholding: a new approach, Computer Vision, Graphics, and Image Processing, vol. 29, pp. 377-393
    if (threshold_in == -5), threshold_auto = opthr_pos(double(im)); end % added 6/7/2016, code from Felix Toran Marti

    if (relative_threshold > 0), threshold_out = threshold_auto * relative_threshold; 
    else, threshold_out = threshold_auto;
    end
    fprintf('Automatic threshold at %g\r', threshold_out)
else
    threshold_out = threshold_in;
end 
im_thr = (im_filt > threshold_out); % invert???
if (do_plot == 1), subplot(223), imshow(im_thr), title(sprintf('thresholded @ %g', threshold_out)), pause(.1), end %, title(['threshold =' num2str(threshold)])
if do_save, imwrite(im_thr, 'threshold_img_illustration/im_thr.tif'), end

%% smooth (dilate/erode), fill
n = 8;
if size_dil_erod > 0 % 13/2/18
    if old_fashion == 1
        disk = strel('disk', size_dil_erod, n);
        im_thr = imdilate(im_thr, disk);
        im_thr = imerode(im_thr, disk);
    elseif old_fashion == 2
        disk = strel('disk', size_dil_erod, n);
        im_thr = imerode(im_thr, disk);
        im_thr = imdilate(im_thr, disk);
    else
        disk1 = strel('disk', round(size_dil_erod/2), n);
        disk2 = strel('disk', size_dil_erod, n);
        im_thr = imerode(im_thr, disk1);
        im_thr = imdilate(im_thr, disk2);
        im_thr = imerode(im_thr, disk1);
    end
end

if do_fill, im_thr = imfill(im_thr, 'holes'); end

% stats = regionprops(im_thr,'Area');
% area =stats.area;

if (area_min > 0), im_thr = bwareafilt(logical(im_thr),[area_min, area_max]); end

% if (do_plot == 1), subplot(224), elseif (do_plot == 2), subplot(122), end
% if (do_plot > 0), imshow(im_thr), drawnow, pause(.1), end
% if (do_plot == 1), title (sprintf('smoothed by %i pxl', size_dil_erod))
% elseif (do_plot == 2), title(sprintf('thresholded @ %g', threshold_out)), end
if (do_plot == 1) && ((size_dil_erod > 0) || do_fill || (area_min > 0))
    subplot(224)
    imshow(im_thr)
    title(sprintf('smoothed by %i pxl', size_dil_erod))
    drawnow, pause(.1)
end
if (do_plot == 2)
%     clf %     figure('WindowStyle','docked')
    imshowpair(im, im_thr)
    title(sprintf('thresholded @ %g', threshold_out))
    drawnow, pause(.1)
end
if (do_plot == 3)
    clf %     figure('WindowStyle','docked')
    title(sprintf('thresholded @ %g', threshold_out))
    pause(.1)
end

if do_save, imwrite(im_thr, 'threshold_img_illustration/im_thr2.tif'), end

% im_smooth = imerode(im_thr, disk);
% im_smooth = imdilate(im_smooth, disk);
% if do_fill, im_smooth = imfill(im_smooth, 'holes'); end
% if do_plot, subplot(224), imshow(im_smooth), drawnow, title ('smoothed Dil/Erod'), end

%%%