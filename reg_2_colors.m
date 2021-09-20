function [mean_tform, sd_tform] = reg_2_colors(file_in, optimizer, metric, do_plot, check_dic_folder)

% function [mean_tform, sd_tform] = reg_2_colors(file_in, optimizer, metric, do_plot, checK_dic_folder)
% default: [optimizer, metric] = imregconfig('multimodal'); reg_2_colors('*.tif', optimizer, metric, 0, 1)
% intended to run from DIC folder, to compute alignment matrix Tr from left
% & right images, then averaging to give mean_tform (& sd)

if nargin < 1, file_in = '*.tif'; end
if nargin < 2
    [optimizer, metric] = imregconfig('multimodal');
    optimizer.InitialRadius = optimizer.InitialRadius / 10;
%     optimizer.MaximumIterations = optimizer.MaximumIterations * 10;
end
if nargin < 4, do_plot = 1; end
if nargin < 5, check_dic_folder = 1; end

if ~contains(lower(cd), 'dic') && check_dic_folder, disp('! caution, not in dic folder, registration cancelled !'), return, end

if ~isempty(dir(['reg' filesep 'mean_tform.mat']))
    mean_tform = importdata(['reg' filesep 'mean_tform.mat']);
    sd_tform = [];
    return
end


warning off images:imshow:magnificationMustBeFitForDockedFigure

files = dir(file_in);
Nf = length(files);
if Nf == 0, disp('niente??'), return, end


Tr = zeros(3, 3, Nf);
%            [cos -sin  0]   [scale -sin  0]
% Tr = scale*[sin  cos  0] ~ [sin  scale  0] if scale ~ 1, sin ~ 0, cos ~ 1
%            [dx   dy   1]   [dx    dy    1]

if isempty(dir('reg')), mkdir('reg'), end

if do_plot, figure('windowstyle', 'docked'), figure(gcf), end

for nf = 1:Nf
    file = files(nf).name;
    if strcmp(file(1), '.'), continue, end % 24/5/2017
    im = imread(file,1);%    im = tiffread(file,1); 28/2/18
    
    %% split green/red
    middle = size(im, 2)/2; % 512, a priori
    imG = im(:, 1:middle);
    imR = im(:, middle+1:end); % % % imR = imR*mean(imG(:))/mean(imR(:));
    
    %% tranformation for registration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tform = imregtform(imR, imG, 'similarity', optimizer, metric); % 'translation' rigid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Tr(:, :, nf) = tform.T;
    
    scale = Tr(1, 1, nf);
    sin_theta = Tr(2, 1, nf);
    delta_x =  Tr(3, 1, nf);
    delta_y =  Tr(3, 2, nf);
    
    fprintf('file %s: theta = %g deg, delta_x = %g pxl, delta_y = %g pxl, scale = %g\r', file, asind(sin_theta), delta_x, delta_y, scale)
    
    if do_plot
        R = imref2d(size(imR)) ;
        imRreg = imwarp(imR, tform, 'OutputView', R); % [a, b] = size(imR); [A, B] = size(imRreg); a0 = round(Amax/2-a/2); b0 = round(Bmax/2-b/2);...
    
%         subplot(211), imagesc([imG imR imRreg])
%         axis image off
%         subplot(234), imshowpair(imG(5:end-5, 5:end-5), imG(5:end-5, 5:end-5))
%         subplot(235), imshowpair(imG(5:end-5, 5:end-5), imR(5:end-5, 5:end-5))
%         subplot(236), imshowpair(imG(5:end-5, 5:end-5), imRreg(5:end-5, 5:end-5))
%         drawnow, pause(.1)
        imshowpair(imG, imRreg)
        title(sprintf('%s: theta = %.2g deg, dx = %.2g pxl, dy = %.2g pxl, scale = %.3g', file, asind(sin_theta), delta_x, delta_y, scale), 'interpreter', 'none')
        pause(.1)
        saveas(gcf, ['reg' filesep file '_check_reg.png'])
    end
end

scale = mean(Tr(1, 1, :)); % cos_theta = Tr(1, 1, :); %Tr(2, 2, :)];% cos_theta = mean(cos_theta(:));
sin_theta = mean(Tr(2, 1, :)); % -Tr(1, 2, :)];
delta_x = mean(Tr(3, 1, :));
delta_y = mean(Tr(3, 2, :));

sd_scale = std(Tr(1, 1, :)); % cos_theta = Tr(1, 1, :); %Tr(2, 2, :)];% cos_theta = mean(cos_theta(:));
sd_sin_theta = std(Tr(2, 1, :)); % -Tr(1, 2, :)];
sd_delta_x = std(Tr(3, 1, :));
sd_delta_y = std(Tr(3, 2, :));
TF = true(size(Tr));
Tr2 = Tr;

if verLessThan('matlab','9.2') % matlab 9.2 = 2017a, intro of isoutlier function
    remove_outliers = 0;
else
    for i = 1:3, for j = 1:3, TF(i,j,:) = isoutlier(Tr(i,j,:)); end, end
    Tr2(TF) = NaN;
    remove_outliers = 1;
end
 
data_pos = 0; do_log = 0; use_sem = 1; clr = -1; x_offset = 0; 

if do_plot
    subplot(221), scatter_hist(squeeze(Tr(1, 1, :)), 'scale', '', data_pos, do_log, use_sem, clr, x_offset, remove_outliers)
    subplot(222), scatter_hist(squeeze(Tr(2, 1, :)*180/2/pi), 'sin theta (deg)', '', data_pos, do_log, use_sem, clr, x_offset, remove_outliers)
    subplot(223), scatter_hist(squeeze(Tr(3, 1, :)), 'delta x', '', data_pos, do_log, use_sem, clr, x_offset, remove_outliers)
    subplot(224), scatter_hist(squeeze(Tr(3, 2, :)), 'delta y', '', data_pos, do_log, use_sem, clr, x_offset, remove_outliers)
    saveas(gcf, ['reg' filesep 'reg_results.png'])
end

mean_tform = tform;
mean_tform.T = nanmean(Tr2, 3);
sd_tform = nanstd(Tr2, [], 3);

if nargin < 1 % as called from cartobyf during MTT
    fprintf('all files: theta = %.2g+/-%.2g deg, delta_x = %.2g+/-%.2g pxl, delta_y = %.2g+/-%.2g pxl, scale = %.3g+/-%.2g\r', asind(sin_theta), asind(sd_sin_theta), delta_x, sd_delta_x, delta_y, sd_delta_y, scale, sd_scale)
    save(['reg' filesep 'mean_tform.mat'], 'mean_tform');
    if do_plot, close(gcf), end
end

%%%