function crop_png(filename, im_type, ROI, sub_crop)
% function crop_png(filename, im_type, ROI, sub_crop)
% crop (recadre / ROI) images des videos MTT, 
% 'Cont' niveaux 2D, '3D' surf 3D, 'xyt' traj_xyt (...)

if nargin==0,
    filename = 'cell1c.stk';
    cd('M:\Arnauld Serge\2005-09-29 SPT specif\carto_v3++++\cell1c_3Dmovie')
    sub_crop = [70 40 60 90];
%     ROI = [];%%%2D iroi=350; jroi=390;   3D iroi=370; jroi=400;
end

if nargin<2, im_type = 'Cont'; end
if nargin<3, ROI = []; end % ROI = [iroi jroi];
if nargin<4, sub_crop = zeros(1,4); end

% if isdir(['carto_v3_output22\' filename(1:end-4) '_movie'])
%     cd (['carto_v3_output22\' filename(1:end-4) '_movie'])
%     chg_dir = 1;
% else
%     chg_dir = 0;
% %     disp('on est où?'), return
% end

if strcmp(im_type, 'Cont') %% 2D contour
    filename_out = [filename(1:end-4) '_Cont'];
    io = 68+sub_crop(1)*3; jo = 255+sub_crop(2)*3;
%     Lcrop = 720-sub_crop(3)*3; Hcrop = 720-sub_crop(4)*3;
    Lcrop = 733-sub_crop(3)*3; Hcrop = 733-sub_crop(4)*3;
%     if ~isempty(ROI), iroi=350; jroi=390; Hroi = 240; Lroi = 240; end
    if ~isempty(ROI), iroi=460; jroi=510; Hroi = 184; Lroi = 92; end
    sampling = 1; % 2; 
%     DD = 2; % dim, pour nom
elseif strcmp(im_type, '3D')  %% 3D (320x240 pxl, full ou ROI)
    filename_out = [filename(1:end-4) '_3D'];
    io = 120; jo = 140;
    Hcrop = 720; Lcrop = 960; 
    if ~isempty(ROI), iroi = 370; jroi = 400; Hroi = 240; Lroi = 320; end % ROI = [371 370+240 401 400+320];
    sampling = 1; % 3; 
%     DD = 3;
elseif strcmp(im_type, 'xyt') % trajs
    filename_out = [filename(1:end-4) '_xyt'];
    io = 68*2+sub_crop(1)*6; jo = 255*2+sub_crop(2)*6;
    Lcrop = 733*2-sub_crop(3)*6; Hcrop = 733*2-sub_crop(4)*6;
    if ~isempty(ROI), iroi=460*2; jroi=510*2; Hroi = 184*2; Lroi = 92*2; end
    sampling = 1; % 2; 
elseif strcmp(im_type, 'testpeak') % testpeaks
    filename_out = [filename(1:end-4) '_testpeak'];
    io = 68+sub_crop(1)*3; jo = 255+sub_crop(2)*3;
    Lcrop = 733-sub_crop(3)*3; Hcrop = 733-sub_crop(4)*3;
%     if ~isempty(ROI), iroi=460*2; jroi=510*2; Hroi = 184*2; Lroi = 92*2; end
    sampling = 1; % 2; 
end% if use_contour

% if isempty(ROI), if ~isdir([num2str(DD) 'D crop']), mkdir([num2str(DD) 'D crop']), end
% else if ~isdir([num2str(DD) 'D ROI']), mkdir([num2str(DD) 'D ROI']), end
% end

N = length(dir([filename_out '*.png'])); % N = 993%100%
% format_out = 'png'; % 'jpg'
% mask = imread('\\amenophis\UsersData$\serge\Bureau\video carto\video mask zoom scalebar.tif');
% mask = repmat(mask,[1 1 3]); % N&B => RVB

fprintf('it is started my kiki           ')

%% go through images
for n=1:N
    if strcmp(im_type, 'xyt') || strcmp(im_type, 'testpeak')
        nb = num2str(n,'%04.0f');
    else
        nb = num2str(n);
    end
    im = imread([filename_out nb '.png']); 
    imcrop = im(io+1:sampling:io+Hcrop, jo+1:sampling:jo+Lcrop, :); % imagesc(imcrop);axis image
%     imwrite(imcrop, [num2str(DD) 'D crop\' filename_out '_crop' num2str(n,'%04i') '.' format_out], format_out) %,'png','tif','Compression','none')...
    imwrite(imcrop, [filename_out '_combi.tif'],...
        'tif', 'WriteMode', 'append', 'compression', 'none')

    if ~isempty(ROI)
        imroi = imcrop(iroi+1:iroi+Hroi, jroi+1:jroi+Lroi, :);
%         imwrite(imroi, [num2str(DD) 'D ROI\' filename_out '_ROI' num2str(n,'%04i') '.' format_out], format_out)
    end
%     im_combi = [imresize(imcrop, [512,512]) imresize(imroi, [512,256])]; % 'nearest'
% %     im_combi = im_combi + mask;
% %     text(4,16,sprintf('%2.1d s',(n-1)*0.036),'font','arial','color',')
%     imwrite(im_combi, [filename_out '_combi.tif'],...
%         'tif', 'WriteMode', 'append', 'compression', 'none')
     fprintf([repmat('\b',1,9) '%4.0i/%4.0i'],n,N)
end

% if chg_dir, cd .. , cd .., end
fprintf('\r')

%% Dans ImageJ
% filtre/gaussblur 0.5, si forte compression requise
% image/type/8 bits color (16 col ok en 2D, 256 en 3D)
% plugin/stack/time stamp 0.036 s, 1 decimal (noir ou blanc taille 28)
% qtime sor3 low 14 fps
% ou gif anim (ou avi + VDub....)