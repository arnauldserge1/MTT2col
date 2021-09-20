function detect_cell_contact(side, do_plot)

% function detect_cell_contact(side, do_plot)
%
% detect contour of contact from MTT trajs of given side
% and save associated trajs for the other side
% (default, right, assuming JAM-C in red, hence selecting JAM-B trajs, in left)
%
% see also MTT_2_colors, split_params_left_right


global N_PARAM PARAM_I PARAM_J PARAM_ALPHA
if isempty(N_PARAM), MTTparams_def; end

warning off images:imshow:magnificationMustBeFitForDockedFigure

if nargin < 1, side = 'right'; end
if nargin < 2, do_plot = 1; end

% % % Nrows = 512; Ncolumns = 512;
all_files = '*.tif';
dirname = 'output23'; % pdef{4};
files = dir(all_files);
Nfiles = length(files);
if isempty(files), disp('No data... Check dir & filename !'), return, end

dir_mask = 'cell_contact';
% dir_out = 'out_of_contact';
if ~isdir([dir_mask filesep dirname]), mkdir([dir_mask filesep dirname]), end

figure('windowstyle', 'docked')

for nf = 1:Nfiles
    file = files(nf).name;
    filename_mask = [file(1:end-4) '_' dir_mask '.tif'];
    filename_out = [file(1:end-4) '_out_of_contact.tif'];
    if ~isempty(dir([dir_mask filesep filename_mask])), disp([dir_mask filesep filename_mask ' done']), continue, end
    % img1 = imread(file);% [Nrows, Ncolumns] = size(img1);
    
    %% get data from given side
    filename_full = [dirname filesep file '_tab_param.mat'];
    if isempty(dir(filename_full)), continue, end
    tab_param = importdata(filename_full);
    img1 = imread(file, 1);
    [Nrows, Ncolumns] = size(img1);
    middle = Ncolumns/2; % size(img1, 2)/2;
    
    [tab_param_left, tab_param_right] = split_params_left_right(tab_param, middle); % tab_param_one_side = split_params_left_right(tab_param, side);
    
    if strcmp(side, 'left'), tab_param_JAMC = tab_param_left; tab_param_JAMB = tab_param_right; JB_color = 'r'; JC_color = 'g';
    else, tab_param_JAMB = tab_param_left; tab_param_JAMC = tab_param_right; JB_color = 'g'; JC_color = 'r';
    end
    
    tab_iC = tab_param_JAMC(PARAM_I-1:N_PARAM:end, :);
    tab_jC = tab_param_JAMC(PARAM_J-1:N_PARAM:end, :);
    tab_iC_int = floor(tab_iC(tab_iC > 0))+1;
    tab_jC_int = floor(tab_jC(tab_jC > 0))+1;
    indC = sub2ind([Nrows, Ncolumns], tab_iC_int, tab_jC_int); % get pixel positions for red staining (right camera)
    
    %% get mask around data positions
    imC = zeros([Nrows, Ncolumns]);
    imC(indC) = 1;
    
    maskC = threshold_img(imC, 0, 3, 8, 1, -1, 0, 1); % def: im_thr = threshold_img(data_in, 0, 3, 8, 1, -1, 0, 0, 0) % keep only largest???
    mask_out = 1 - maskC;
    maskC(1,1) = 0; % for points at 0,0 (before begin of traces)(cannot induce artefact, since MTT do not detect in image borders)
    
    [~, L, N] = bwboundaries(maskC, 'noholes'); % imshow(maskC)
    areae = zeros(1, N);
    for k = 1:N
        stats = regionprops(L, 'Area', 'Centroid');
        areae(k) = stats(k).Area;
    end
    
    big_one = find(areae==max(areae), 1); % plusieurs trucs??
    if (N > 0), contact_center = stats(big_one).Centroid;
    else, contact_center = []; end
    
    if do_plot
        imshow(mask_out), axis ij % ....
        title(filename_mask, 'interpreter', 'none')
        hold on
    end
    
    %% get all JAM-B trajs
    tab_iB = tab_param_JAMB(PARAM_I-1:N_PARAM:end, :);
    tab_jB = tab_param_JAMB(PARAM_J-1:N_PARAM:end, :);% tab_iB = [1 20 300 310;1.1 21 300.1 309]; tab_jB = [1 10 80 70; 1.1 11 80.1 71];
    tab_iB_int = floor(tab_iB)+1;
    tab_jB_int = floor(tab_jB)+1;
    if do_plot, plot(tab_jB, tab_iB, [JB_color '.']), end % xy = ji!!
    
    %% get JAM-B in mask & out
    indB =sub2ind([Nrows, Ncolumns], tab_iB_int, tab_jB_int); % select_pk_in_mask
    indB_mask = maskC(indB); % JAM-B positions within JAM-C mask (= cell contact) % [iB_mask, jB_mask] =ind2sub([h, w], indB_mask);
    indB_out = mask_out(indB); % JAM-B positions out of JAM-C mask (= stromal, out of cell contact)
    
    tab_iB_mask = tab_iB .* indB_mask; %   tab_iB_mask = tab_iB(indB_mask);
    tab_jB_mask = tab_jB .* indB_mask; %   tab_jB_mask = tab_jB(indB_mask);
            
    tab_iB_out = tab_iB .* indB_out;
    tab_jB_out = tab_jB .* indB_out;
    
    %% plot everybody
    if do_plot
        plot(tab_jC, tab_iC, [JC_color '.'])
        plot(tab_jB_mask, tab_iB_mask, 'yo')
        plot(tab_jB_out, tab_iB_out, [JB_color '.'])
        
        %% center of contact: yellow cross
        if ~isempty(contact_center), plot(contact_center(1), contact_center(2), 'y+', 'markersize', 18), end
        drawnow
        
        saveas(gcf, [dir_mask filesep filename_mask])
        
        %% plot out only
        clf
        imshow(mask_out), axis ij, title(filename_out, 'interpreter', 'none'), hold on
        plot(tab_jB_out, tab_iB_out, [JB_color '.'])
        drawnow
        
        saveas(gcf, [dir_mask filesep filename_out])
    end
    
    %% sort & save in & out
    tab_param_JAMB_mask = tab_param_JAMB;
    tab_param_JAMB_out = tab_param_JAMB;
    for np = 2:N_PARAM % handling blink, keep ij after bleach ???
        tab_par_np = tab_param_JAMB(np:N_PARAM:end, :);
        tab_param_JAMB_mask(np:N_PARAM:end, :) = tab_par_np .* indB_mask;
        tab_param_JAMB_out(np:N_PARAM:end, :) = tab_par_np .* indB_out;
    end
    
    alpha_jB_mask = tab_param_JAMB_mask(PARAM_ALPHA-1:N_PARAM:end, :); % intensity of fitted SM
    TrcLen_mask = sum(alpha_jB_mask > 0);
    tab_param = tab_param_JAMB_mask(:, TrcLen_mask > 0); % NOT tab_param_mask, since ct4 and others want a variable named tab_param !!!
    % min trc length?????
    N_trc_jB_mask = size(tab_param, 2); % sum(TrcLen > 0);
    fprintf('%s: saving %i trajs for JAM-B in mask\r', file, N_trc_jB_mask)
    save([dir_mask filesep dirname filesep filename_mask  '_tab_param.mat'], 'tab_param')
    
    alpha_jB_out = tab_param_JAMB_out(PARAM_ALPHA-1:N_PARAM:end, :);
    TrcLen_out = sum(alpha_jB_out > 0);
    tab_param = tab_param_JAMB_out(:, TrcLen_out > 0);
    N_trc_jB_out = size(tab_param, 2);
    fprintf('%s: saving %i trajs for JAM-B out of mask\r', file, N_trc_jB_out)
    save([dir_mask filesep dirname filesep filename_out  '_tab_param.mat'], 'tab_param')

end % file

if do_plot, close gcf, end
[pxl_size, time_lag] = get_calib3;

cd(dir_mask)
ct4_by_file(['*' dir_mask '.tif'], dirname, 'coloc', [], pxl_size, time_lag);
ct4_by_file('*_out_of_contact.tif', dirname, 'coloc', [], pxl_size, time_lag);
ct5_by_file(['*' dir_mask '.tif'], dirname, [], pxl_size, time_lag, side);
ct5_by_file('*_out_of_contact.tif', dirname, [], pxl_size, time_lag, side);
% % % % % removed 19/2/2021
% % % % % amax = 40; disp('Sorting traces for D & v')
% % % % % fit_directed_by_file(['*' dir_mask '.tif'], amax, pxl_size, time_lag); % side = ''!!
% % % % % fit_directed_by_file('*_out_of_contact.tif', amax, pxl_size, time_lag);
% % % % % cd ..

%%%

%cd('carto_coloc')%'merged files') % if isdir('max'), cd('max'), dirmax = 1; else  dirmax = 0; end
% if strcmp(filename(1), '.'), continue, end
% im = imread(file); im = imresize(im, 512/size(im, 2)); im = im((size(im, 1) - 512 + 1):end, :, :); % remove coloc title!  imR = (im(:,:,1) > 254);
%     if do_plot, subplot(221), imshow(im_thr), end
%     im_thr_R = imdilate(imR, disk);    im_thr_R = imerode(im_thr_R, disk);    im_thr_R = imfill(im_thr_R, 'holes');    if do_plot, subplot(222), imshow(im_thr_R), end
%     im_thr_R = threshold_img(im(:,:,1)>254);    %     im_thr_G = threshold_img(im(:,:,2)>254);    %     im_thr_Y = threshold_img(im(:,:,2));
% if dirmax, cd .., end