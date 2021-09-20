function img_max = image_max(file, tab_num, show_result, keep_class)

% function img_max = image_max(file, tab_num, show_result, keep_class)
% compute max for all pixels of a video file 
% => along z or t, according to stack type, movie in 3D or over time
%  for frames tab_num=n1:n2 (default all: 1:end)
% see also meanimage.m
%
% AS, 1/7/10

if strcmp(file(end-3:end), '.stk'), use_stk = 1; else, use_stk = 0; end

if nargin < 4
    keep_class = 0;
end

if nargin < 3
    show_result = 1;
end

if use_stk
    if nargin < 2
        [img1, Nimg] = tiffread(file,1); % AS 16/4/2020   [img1, Nimg] = tiffread(file,tab_num(1)); ???
        tab_num = 1:Nimg;
    else
        img1 = tiffread(file,tab_num(1));
    end
else
    if nargin < 2
        Nimg = length(imfinfo(file));
        tab_num = 1:Nimg;
    end
    img1 = imread(file,tab_num(1));
    if ~keep_class, img1 = double(img1); end
end



if show_result, fprintf('getting max for %s    ', file), end

img_max = img1;
for i = tab_num(2:end)
    if use_stk
        imgi = tiffread(file,i);
    else
        imgi = imread(file,i);
    end
    if ~keep_class, imgi = double(imgi); end
    img_max = max(img1, imgi);
    if show_result
        fprintf('\b\b\b\b%4i', i)
    end
end

if show_result
    figure('WindowStyle','docked')
    imagesc(img_max);
    axis image;
    colormap(1-gray);
    fprintf('\r')
end
