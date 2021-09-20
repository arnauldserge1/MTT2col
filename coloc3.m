function dist = coloc3(file_in, do_plot, do_all_files, string_reg)

% function dist = coloc3(file_in, do_plot, do_all_files, string_reg)
% compute distance of nearest green spot for each red spot, in each frame
% d = min(sqrt(di2 + dj2)) matrix of N_image by N_red_trace
% or cell if multiple files
% default: dist = coloc3('*.tif', 1, 1, 'reg')
%
% see also traj_xy_coloc, MTT23i

global N_PARAM PARAM_I PARAM_J PARAM_ALPHA
if isempty(N_PARAM), MTTparams_def; end

if nargin < 1, file_in = '*.tif'; end
if nargin < 2, do_plot = 1; end
if nargin < 3, do_all_files = 1; end
if nargin < 4, string_reg = 'reg'; end %_by_file

files = dir2(file_in);
Nf = length(files);
if Nf == 0, disp('niente??'), return, end
if Nf <= 1, do_all_files = 0; end
dist = cell(1, Nf);

if isempty(dir('Red-Green distances')), mkdir('Red-Green distances'), end

if do_plot, figure('windowstyle', 'docked'), figure(gcf), end

for nf = 1:Nf
    file = files(nf).name;
    
    name_out = ['Red-Green distances' filesep file(1:end-4) '_RGdist' string_reg];
    if ~isempty(dir([name_out '.mat'])), dist{nf} = importdata([name_out_all '.mat'], 'dd'); continue, end
    
    %% load i j for R & G
    tab_param = importdata(['output23' filesep file '_tab_param.mat']); % tab_param = fread_all_params(file);
    if isempty(tab_param), disp('no data..'), continue, end
    
    img1 = imread(file, 1);
    middle = size(img1, 2)/2;
    [tab_paramG, tab_paramR] = split_params_left_right(tab_param, middle);
    
    iiG = tab_paramG(PARAM_I-1:N_PARAM:end, :); % i=y
    jjG = tab_paramG(PARAM_J-1:N_PARAM:end, :); % j=x
    alphaG = tab_paramG(PARAM_ALPHA-1:N_PARAM:end, :);
    
    iiG(alphaG == 0) = nan;
    jjG(alphaG == 0) = nan;
    
    iiR = tab_paramR(PARAM_I-1:N_PARAM:end, :); % i=y
    jjR = tab_paramR(PARAM_J-1:N_PARAM:end, :); % j=x
    alphaR = tab_paramR(PARAM_ALPHA-1:N_PARAM:end, :);
    
    iiR(alphaR == 0) = nan;
    jjR(alphaR == 0) = nan;
    
    [Nimg, nR] = size(jjR);
    nG = size(jjG, 2);
    if nR == 0, disp('no red..'), continue, end
    if nG == 0, disp('no green..'), continue, end
    
    % % %     if ~isempty(dir(['DIC' filesep 'mean_tform.mat'])) && do_reg
    if do_plot, clf, subplot(1,3,[1 2]), plot(iiG, jjG, 'g.', iiR, jjR, 'm.'), axis tight off, hold on, end
    % % %         tform = importdata(['DIC' filesep 'mean_tform.mat']);
    
    %% reg (not by file)
    searchType = 'file';
    if ~isempty(dir('DIC')) && exist('imregtform', searchType) && ~isempty(string_reg)%builtin')
        cd('DIC'), tform = reg_2_colors(file); cd ..
        [jjRreg, iiRreg] = transformPointsForward(tform, jjR, iiR); % x,y == j,i and transformPointsForward expect x,y, hence j,i !!!
    else
        jjRreg = jjR;
        iiRreg = iiR;
    end
    if do_plot, plot(iiRreg, jjRreg, 'r.'), hold off, pause(.1), end
    % % %     end
    
    %% dist
    dd = zeros(Nimg, nR);
    
    fprintf('Computing red-green distances in each image:    ')
    for nt = 1:Nimg
        iRmat = repmat(iiRreg(nt, :), nG, 1);
        iGmat = repmat(iiG(nt, :)', 1, nR);
        jRmat = repmat(jjRreg(nt, :), nG, 1);
        jGmat = repmat(jjG(nt, :)', 1, nR);
        distances = sqrt((iRmat-iGmat).^2 + (jRmat-jGmat).^2); % matrix nG x nR, dist to all green = sqrt(di2 + dj2) 
        dd(nt, :) = min(distances); % min == nearest, for each red
        
        fprintf('\b\b\b\b%4i', nt)
    end
    fprintf('\r')
    
    %% show & save histos (lin & log)
    if do_plot
        N = round(2*sqrt(sum(~isnan(dd(:))))); % nan occurs when red blinks, bleaches or not yet detected
        subplot(233), hist(dd(:), N), xlabel('nearest distance red-green (pxl)'), title(file, 'interpreter', 'none')
        subplot(236), hist(log10(dd(:)), N), xlabel('nearest distance red-green (log10, pxl)')
        saveas(gcf, [name_out '.png'])
    end
    save([name_out '.mat'], 'dd')
    dist{nf} = dd;
end

%% all
if do_all_files
    name_out_all = ['Red-Green distances' filesep file_in '_RGdist_um' string_reg];
    name_out_all(strfind(name_out_all, '*')) = '@'; % all files: expecting a name such as '*.tif'
    if ~isempty(dir([name_out_all '.mat'])), close(gcf), return, end

    files = dir(file_in);
    Nf = length(files);
    dist = cell(1, Nf);
    Nt = zeros(1, Nf);
    ok = ones(1, Nf);
    pxl_size = get_calib3;
    
    for nf = 1:Nf
        file = files(nf).name;
        name_in = ['Red-Green distances' filesep file(1:end-4) '_RGdist' string_reg '.mat'];
        if ~isempty(dir(name_in))
            load(name_in, 'dd') % dd = importdata(name_in); 
            dist{nf} = dd*pxl_size; % in um
            Nt(nf) = size(dd,1);
        else
           ok(nf) = 0;
        end
    end
    
    for nf = 1:Nf
        [l, c] = size(dist{nf});
        dist{nf} = [dist{nf}; NaN(max(Nt)-l, c)]; % complete with NaN if lower number of frame for current file
    end
    
    dist = dist(ok>0);
    Nf = sum(ok);
    dist_all = cell2mat(dist);
    
    %% histo all: compil or m & sd
    if do_plot
        N = round(2*sqrt(sum(~isnan(dist_all(:)))));
        clf
        subplot(121), hist(dist_all(:), N), xlabel('nearest distance red-green (um)'), title(file_in, 'interpreter', 'none')
        subplot(122), hist(log10(dist_all(:)), N), xlabel('nearest distance red-green (log10, um)')
        saveas(gcf, [name_out_all '.png'])
        
        subplot(121), mean_hist(dist), xlabel('nearest distance red-green (um)'), title(file_in, 'interpreter', 'none')
        log_dist = dist; for nf = 1:Nf, log_dist{nf} = log10(dist{nf}); end
        subplot(122), mean_hist(log_dist), xlabel('nearest distance red-green (log10, um)')
        saveas(gcf, [name_out_all 'mean.png'])
    end
    save([name_out_all '.mat'], 'dist')
end

%%%