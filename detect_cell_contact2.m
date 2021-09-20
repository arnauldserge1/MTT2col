function detect_cell_contact2

% function detect_cell_contact2
%
% detect contour of contact from dic_mask, with black holes on KG1 cells, drawn with imageJ
% and save associated trajs in/out (default, selecting JAM-B trajs, left)
% ! caution: assume JAM-C is right (and no used, as compared to detect_cell_contact), and select JAM-B, at left, from dic_mask
%
% see also MTT_2_colors, split_params_left_right

global N_PARAM PARAM_I PARAM_J PARAM_ALPHA redo
if isempty(N_PARAM), MTTparams_def; end

redo = 1;
dic2 = 'dic - Copy'; % 'dic_mask'

dirname = 'output23'; % pdef{4};
files = dir2('*.tif');
Nfiles = length(files);
if isempty(files), disp('No data... Check dir & filename !'), return, end

dir_mask1 = 'cell_contact';
dir_mask2 = 'cell_contact_ImageJ';
do_ana = 0;
if isdir(dir_mask1) &&  ~isdir(dir_mask1) % keep copy
    disp(['copying ' dir_mask1 '...'])
    copyfile(dir_mask1, dir_mask2) %if ~isdir([dir_mask2 filesep dirname]), mkdir([dir_mask2 filesep dirname]), end
    do_ana = 1;
end
str = {'in', 'out'};

figure('windowstyle', 'docked')

for nf = 1:Nfiles
    file = files(nf).name;
    
    filename_in = [file(1:end-4) '_cell_contact.tif'];
    filename_out = [file(1:end-4) '_out_of_contact.tif'];   
    file_out = [dir_mask1 filesep dirname filesep eval(['filename_' str{2}])  '_tab_param.mat'];
    if ~isempty(dir(file_out)), disp([file_out ' done']), continue, end
%     if ~isempty(dir(filename_in)), disp([filename_in ' done']), continue, end
    if isempty(dir([dic2 filesep file])), disp(['no dic image for ' dic2 filesep file ', skept']), continue, end
    
    %% get all JAM-B trajs
    filename_full = [dirname filesep file '_tab_param.mat'];
    if isempty(dir(filename_full)), continue, end
    tab_param = importdata(filename_full);
    img1 = imread(file, 1);
    [Nrows, Ncolumns] = size(img1);
    middle = Ncolumns/2; % size(img1, 2)/2;
    
    tab_param_JAMB = split_params_left_right(tab_param, middle, 'left');  % if strcmp(side, 'left'), tab_param_JAMC = tab_param_left; tab_param_JAMB = tab_param_right; JB_color = 'r'; JC_color = 'g';
    tab_iB = tab_param_JAMB(PARAM_I-1:N_PARAM:end, :);
    tab_jB = tab_param_JAMB(PARAM_J-1:N_PARAM:end, :);% tab_iB = [1 20 300 310;1.1 21 300.1 309]; tab_jB = [1 10 80 70; 1.1 11 80.1 71];
    tab_iB_int = floor(tab_iB)+1;
    tab_jB_int = floor(tab_jB)+1;
    
    %% get mask made with ImageJ
    im_mask = imread([dic2 filesep file]);
    im_mask = im_mask(:, 1:middle);
    mask_in = (im_mask == 0); % black area (0) = KG1 cell(s)
    mask_out = ~mask_in;
    
    %% get JAM-B in mask & out
    indB =sub2ind([Nrows, Ncolumns], tab_iB_int, tab_jB_int); % select_pk_in_mask
    indB_in = mask_in(indB); % JAM-B positions within mask (= cell contact) % [iB_mask, jB_mask] =ind2sub([h, w], indB_mask);
    indB_out = mask_out(indB); % JAM-B positions out of mask (= out of cell contact)
    
    %% sort & save in & out
    tab_param_JAMB_in = tab_param_JAMB;
    tab_param_JAMB_out = tab_param_JAMB;
    for np = 2:N_PARAM % handling blink, keep ij after bleach ???
        tab_par_np = tab_param_JAMB(np:N_PARAM:end, :);
        tab_param_JAMB_in(np:N_PARAM:end, :) = tab_par_np .* indB_in;
        tab_param_JAMB_out(np:N_PARAM:end, :) = tab_par_np .* indB_out;
    end
    
    %% save images (so ct4 can use them for '*in....tif' '*out....tif')
    clf, imagesc(im_mask), colormap('gray'), hold on
    plot(tab_param_JAMB_in(3:8:end,:), tab_param_JAMB_in(2:8:end,:), 'og'), pause(.1)
    saveas(gcf, [dir_mask2 filesep filename_in])
    
    clf, imagesc(mask_out), hold on
    plot(tab_param_JAMB_out(3:8:end,:), tab_param_JAMB_out(2:8:end,:), 'og'), pause(.1)
    saveas(gcf, [dir_mask2 filesep filename_out])

    data = {tab_param_JAMB_in, tab_param_JAMB_out};
    for side = 1:2
        alpha = data{side}(PARAM_ALPHA-1:N_PARAM:end, :); % intensity of fitted SM  % a=eval(['tab_param_JAMB_' str{side} '(PARAM_ALPHA-1:N_PARAM:end, :)']);
        TrcLen = sum(alpha > 0);    %% min trc length?????
        tab_param = data{side}(:, TrcLen > 0); % NOT tab_param_mask, since ct4 and others want a variable named tab_param !!!
        fprintf('%s: saving %i trajs for JAM-B %s mask\r', file, size(tab_param, 2), str{side})
        save(file_out, 'tab_param')
    end
end % file

close(gcf)

%% run analysis for in/out contact
if do_ana
    [pxl_size, time_lag] = get_calib3;
    amax = 40;
    
    cd(dir_mask2)
    files = {'*cell_contact.tif', '*_out_of_contact.tif'};
    disp('ana: ct4, ct5 & sorting traces for D & v')
    for side = 1:2
        ct4_by_file(files{side}, dirname, 'coloc', [], pxl_size, time_lag);
        ct5_by_file(files{side}, dirname, [], pxl_size, time_lag, 'left');
        fit_directed_by_file(files{side}, amax, pxl_size, time_lag); % side = ''!!
    end
    cd ..
end

redo = 0; % 7/5/2019

%%%

% cd('/Volumes/spinning/Oksana/data for article1/non blocking 19H36')
% cd('/Volumes/spinning/Oksana/data for article1/blocking 208206')
% manips = dir('2017*--ARTICLE--');
% for nm = 1:length(manips)
%     cd(manips(nm).name)
%     detect_cell_contact2
%     cd ..
% end