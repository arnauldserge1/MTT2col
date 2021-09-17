function carto_movie_2colors2(file)
% trace les trajs (xy) image apres image (t)


global PARAM_I PARAM_J PARAM_ALPHA N_PARAM
global coloc_dist_max coloc_time_min

if isempty(coloc_dist_max), coloc_dist_max = 2; end % pxl
if isempty(coloc_time_min), coloc_time_min = 2, end % frm
if isempty(N_PARAM), MTTparams_def; end


if nargin < 1
    [file, newdir] = uigetfile('*.tif','please select file for movie');
    cd(newdir)
    cd('/Users/serge/Documents/labo/articles/JAM nanodyn/videos/raw')
%     cd(\spinning\Oksana\non blocking 19H36\2017-06-21 MS5 JB KG1 JC+ non blocking\ stream3
    %      cd('\\SAUV02\spinning\jessica\data SPT\manip A647 KG1 JC+\2015-05-07 MS5 anti-rb A488 + KG1 anti-mouse A647 colo primaire puis secondaire'), file = 'Stream5.tif';
    %     cd('/Volumes/spinning/jessica/data SPT/manip A647 KG1 JC+/2015-05-07 MS5 anti-rb A488 + KG1 anti-mouse A647 colo primaire puis secondaire'), file = 'Stream5.tif';
%     if ispc, cd('D:\DATA_D\jessica\data SPT\manip A647 KG1 JC+\2015-05-07 MS5 anti-rb A488 + KG1 anti-mouse A647 reglage 491-642'),
%     else, cd('/Volumes/spinning/jessica/data SPT/manip A647 KG1 JC+/2015-05-07 MS5 anti-rb A488 + KG1 anti-mouse A647 reglage 491-642');     end
%     file = 'Stream8.tif';
    % % %     cd('D:\DATA_D\jessica\2015-02-26 MS5-JB-A488 KG1-JC-Qd655 SPT8'), file = 'Stream7.tif'; DRIFT!!!
end % files = dir('*.tif'); file = files(1).name;

dir_out = '/Users/serge/Documents/labo/articles/JAM nanodyn/videos/raw/';
file_out = [dir_out file(1:end-4) '_xyt.avi'];
if exist(file_out, 'file'), disp([file_out ' already existing! Erasing?? press any key or wait 3 seconds to continue or ctl+C to stop']); pause(3), end

[pxl_size, time_lag] = get_calib3;
fprintf('assuming coloc for distance below %g pxl = %g nm\r', coloc_dist_max, coloc_dist_max*pxl_size*1000) % cf. d_min = 0.6 lambda/NA = 0.6*580/1.45 = 240 nm
fprintf('assuming coloc for duration above %g frm = %g s\r', coloc_dist_max, coloc_dist_max*time_lag)

params_def = MTTparams_def; dirname = params_def{4};
string_reg = 'reg';

if isempty(file), disp('No data... Check dir & filename !'), return, end

%% data
filename_full = [dirname filesep file '_tab_param.mat'];
if isempty(dir(filename_full)), disp('no data??'), return, end
tab_param = importdata(filename_full);
Tmax = min(500, size(tab_param, 1)/N_PARAM);

%% split Red/Green
img1 = tiffread(file, 1);
middle = size(img1, 2)/2; % default: 512
[tab_paramG, tab_paramR] = split_params_left_right(tab_param, middle); % Green = Left = JAM-B / Red = Right = JAM-C

%% Green peaks
tab_iG = tab_paramG(PARAM_I-1:N_PARAM:end, :);
tab_jG = tab_paramG(PARAM_J-1:N_PARAM:end, :);
tab_alphaG = tab_paramG(PARAM_ALPHA-1:N_PARAM:end, :);
tab_iG(tab_alphaG == 0) = nan;
tab_jG(tab_alphaG == 0) = nan;
r2 = calcul_r2(tab_paramG); % r2: Nfrm-1 * Ntraj Green (JAM-B)
logDG = log10(r2*pxl_size^2/time_lag); % hist(logDG(:), 100) D by step

%% Red peaks
tab_iR = tab_paramR(PARAM_I-1:N_PARAM:end, :);
tab_jR = tab_paramR(PARAM_J-1:N_PARAM:end, :);
tab_alphaR = tab_paramR(PARAM_ALPHA-1:N_PARAM:end, :);
tab_iR(tab_alphaR == 0) = nan;
tab_jR(tab_alphaR == 0) = nan;
r2 = calcul_r2(tab_paramR);
logDR = log10(r2*pxl_size^2/time_lag); % compute log of diff for color coding 


%% compute or load dist for coloc
dist_file = ['Red-Green distances' filesep file(1:end-4) '_RGdist' string_reg '.mat' ];
if isempty(dir(dist_file)), RGdist = coloc3(file, 0, 0, string_reg); RGdist = RGdist{1};
else, RGdist = importdata(dist_file, 'dd'); end
coloc = (RGdist < coloc_dist_max); % Nfrm * Ntrc, logical

%% --- go through frames ---
writerObj = VideoWriter(file_out);
writerObj.FrameRate = 30; % 10 for real time or 30, 50... for accelerated 3x, 5x, thus shorter
writerObj.Quality = 100;
% writerObj.Height = 512; % or more for super res?
% writerObj.VideoCompressionMethod = 'H.264'; % as recommended by NatMet
open(writerObj);
disp('frame            ')
figure%('WindowStyle', 'docked')
set(gcf, 'position', [10 10 240 240])%set(gcf, 'position', [10 10 512 512])
for t = 1:Tmax-1
    DIC_name = dicname(file, 0);
    DIC_image(file, DIC_name, 0, 0, 0, '_right', 0, 0);
    title('')
    set(gca, 'position', [0 0 1 1])
    
    %%  blue shadow (plotted before, to appear under the traces)
    for t2 = 1:t
        tt = [t2, t2+1];
        
        ok = find(tab_alphaG(t2, :) > 0);
        for nt = 1:length(ok) % Ntrc ok
            plot(tab_jG(tt, ok(nt)), tab_iG(tt, ok(nt)), 'color', [0 0 0.5], 'linewidth', 2)
        end
        
        ok = find(tab_alphaR(t2, :) > 0);
        for nt = 1:length(ok)
            plot(tab_jR(tt, ok(nt)), tab_iR(tt, ok(nt)), 'color', [0 0 0.5], 'linewidth', 2)
        end
    end
    
    for t2 = 1:t
        tt = [t2, t2+1];
        
        %% Green
        ok = find(tab_alphaG(t2, :) > 0);
        for nt = 1:length(ok) % Ntrc ok
            logD = logDG(t2, ok(nt)); % coding for D @ current step (or fast/diff/lin/conf??)
            logD = sort([-4, logD, 0]); logD = logD(2); % Dmin = 10-4, Dmax = 1, log distrib
            if (logD > -2), clr = [0.5 1 0.5] - [1 0 1]*(logD/4); else, clr = [0 1.5 0] + [0 1 0]*(logD/4); end
            plot(tab_jG(tt, ok(nt)), tab_iG(tt, ok(nt)), 'color', clr, 'linewidth', 1.5)
        end
        
        %% Magenta
        ok = find(tab_alphaR(t2, :) > 0);
        for nt = 1:length(ok)
            logD = logDR(t2, ok(nt));
            logD = sort([-4, logD, 0]); logD = logD(2);
% %             if (logD > -2), clr = [1 0.5 0.5] - [0 1 1]*(logD/4); else, clr = [1.5 0 0] + [1 0 0]*(logD/4); end
            if (logD > -2), clr = [1 0.5 1] - [0 1 0]*(logD/4); else, clr = [1.5 0 1.5] + [1 0 1]*(logD/4); end
            plot(tab_jR(tt, ok(nt)), tab_iR(tt, ok(nt)), 'color', clr, 'linewidth', 1.5)
        end
        
        %% White
        ok = find(coloc(t2, :)); % coloc points for this frame
        h = zeros(size(ok));
        for nc = 1:length(ok) % Ncoloc @ t
            h(nc) = plot(tab_jR(t2, ok(nc)), tab_iR(t2, ok(nc)), 'wo');%, 'markersize', 16)
        end
    end % t2
    
    text(50, 12, [num2str(t*time_lag, '%.1f') ' s'], 'color', 'w', 'HorizontalAlignment', 'right')
    
    %% save pict
    frame = getframe;
    writeVideo(writerObj, frame);
    fprintf([repmat('\b', 1, 11) '%5i/%5i'], t, Tmax-1)
    clf
end % t

close(writerObj);
%%%