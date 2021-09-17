function [N_coloc, coloc_duration] = traj_xy_coloc(file, dirname, string_reg, codage)

% function [N_coloc, coloc_duration] = traj_xy_coloc(file, dirname, string_reg)
%
% plot all traces on first or transm. image
% color: magenta, green or white (coloc)
%
% see also coloc3, cartobyf


global PARAM_I PARAM_J PARAM_ALPHA N_PARAM
global coloc_dist_max coloc_time_min

if isempty(coloc_dist_max), coloc_dist_max = 2; end % pxl
if isempty(coloc_time_min), coloc_time_min = 2, end % frm
if isempty(N_PARAM), MTTparams_def; end

[pxl_size, time_lag] = get_calib3;
fprintf('assuming coloc for distance below %g pxl = %g nm\r', coloc_dist_max, coloc_dist_max*pxl_size*1000) % cf. d_min = 0.6 lambda/NA = 0.6*580/1.45 = 240 nm
fprintf('assuming coloc for duration above %g frm = %g s\r', coloc_dist_max, coloc_dist_max*time_lag)


if nargin < 1, files = dir('*.tif'); if isempty(files), files = dir('*.stk'); end, file = files(1).name; end
if nargin < 2, params_def = MTTparams_def; dirname = params_def{4}; end
if nargin < 3, string_reg = 'reg'; end
if nargin < 4, codage = 'diff'; end % or time

if isempty(file), disp('No data... Check dir & filename !'), return, end

%% data
filename_full = [dirname filesep file '_tab_param.mat'];
if isempty(dir(filename_full)), disp('no data??'), N_coloc = []; coloc_duration = []; return, end
tab_param = importdata(filename_full);
Tmax = size(tab_i, 1);

%% split Red/Green
img1 = imread(file, 1);
middle = size(img1, 2)/2; % default: 512
[tab_paramG, tab_paramR] = split_params_left_right(tab_param, middle); % Green = Left = JAM-B / Red = Right = JAM-C

tab_iG = tab_paramG(PARAM_I-1:N_PARAM:end,:);
tab_jG = tab_paramG(PARAM_J-1:N_PARAM:end,:);
tab_alphaG = tab_paramG(PARAM_ALPHA-1:N_PARAM:end,:);
r2 = calcul_r2(tab_paramG); % r2: (Nfrm-1) * Ntraj_Green
logDG = log10(r2*pxl_size^2/time_lag); % hist(logDG(:), 100) D by step

tab_iR = tab_paramR(PARAM_I-1:N_PARAM:end,:);
tab_jR = tab_paramR(PARAM_J-1:N_PARAM:end,:);
tab_alphaR = tab_paramR(PARAM_ALPHA-1:N_PARAM:end,:);
r2 = calcul_r2(tab_paramR); % r2: (Nfrm-1) * Ntraj_Red
logDR = log10(r2*pxl_size^2/time_lag); % hist(logDG(:), 100) D by step

%% apply reg with mean_tform to i,j Red => i,j Red reg (Rr) 11/5/2017
if ~isempty(dir(['dic' filesep 'reg' filesep 'mean_tform.mat']))
    mean_tform = importdata(['dic' filesep 'reg' filesep 'mean_tform.mat']);
    [tab_jRreg, tab_iRreg] = transformPointsForward(mean_tform, tab_jR, tab_iR); % x,y == j,i and transformPointsForward expect x,y, hence j,i !!!
else
    disp(['! Caution, couldn''t find reg' filesep 'mean_tform.mat, red traces not registered !'])
    tab_jRreg = tab_jR;
    tab_iRreg = tab_iR;
end

% if isdir('max'), show_max_image(file)
DIC_image(file, dicname(file), 0, 0, 0, '_right', 1, 0); % figure... imagesc(pict), axis image ij off, colormap('gray'), hold on

%% compute or load dist
dist_file = ['Red-Green distances' filesep file(1:end-4) '_RGdist' string_reg '.mat' ];
if isempty(dir(dist_file))
    RGdist = coloc3(file, 0, 0, string_reg); % reg
    RGdist = RGdist{1}; % close(gcf)
else
    RGdist = importdata(dist_file, 'dd');
end

%% --- go through traces ---

%% Green
NtrcG = size(tab_iG, 2);
disp('green traj:            ')
for ntrc = 1:NtrcG
    ok = find(tab_alphaG(:, ntrc) > 0);
    for nf = 1:length(ok)-1
        t = ok(nf);
        tt = [ok(nf) ok(nf+1)];
        if strcmp(codage, 'diff')
            logD = logDG(t, ok(nt));
            logD = sort([-4, logD, 0]);
            logD = logD(2);
            if (logD > -2), clr = [1 0.5 1] - [0 1 0]*(logD/4);
            else, clr = [1.5 0 1.5] + [1 0 1]*(logD/4);
            end
        else
            if t < Tmax/2, clr = [0.5 1 0.5] - [1 0 1]*(t/Tmax);
            else, clr = [0 1.5 0] - [0 1 0]*(t/Tmax);
            end
        end
        % clr = [1/2-tt/Tmax; 3/2-tt/Tmax; 1/2-tt/Tmax]; clr = min(max(clr, 0), 1);
        % patch('XData', [nan; tab_jG; nan], 'YData', [nan; tab_iG(t, 1); nan], 'CData', [nan; clr'; nan], 'FaceColor', 'interp', 'EdgeColor', 'interp')
        plot(tab_jG(tt, ntrc), tab_iG(tt, ntrc), 'color', clr)
        % plot(tab_jG(ok, ntrc), tab_iG(ok, ntrc), 'color', clr)
    end
    if mod(ntrc, 10) == 0, fprintf([repmat('\b', 1, 11) '%5i/%5i'], ntrc, NtrcG), drawnow expose, end
end
fprintf([repmat('\b', 1, 11) '%5i/%5i\r'], NtrcG, NtrcG)

%% Red (Magenta..)
NtrcR = size(tab_iRreg, 2);
disp('red traj:            ')
for ntrc = 1:NtrcR
    ok = find(tab_alphaR(:, ntrc) > 0);
    for nf = 1:length(ok)-1
        t = ok(nf);
        tt = [ok(nf) ok(nf+1)];
        if strcmp(codage, 'diff')
            logD = logDR(t, ok(nt));
            logD = sort([-4, logD, 0]);
            logD = logD(2);
            if (logD > -2), clr = [1 0.5 1] - [0 1 0]*(logD/4);
            else, clr = [1.5 0 1.5] + [1 0 1]*(logD/4);
            end
        else
            if t < Tmax/2, clr = [1 0.5 1] - [0 1 0]*(t/Tmax);
            else, clr = [1.5 0 1.5] - [1 0 1]*(t/Tmax);
            end
        end
        plot(tab_jRreg(tt, ntrc), tab_iRreg(tt, ntrc), 'color', clr)
    end
    if mod(ntrc, 10) == 0, fprintf([repmat('\b', 1, 11) '%5i/%5i'], ntrc, NtrcR), drawnow expose, end
end
fprintf([repmat('\b', 1, 11) '%5i/%5i\r'], NtrcR, NtrcR)

%% white
disp('white coloc event:            ')
coloc = (RGdist < coloc_dist_max); % Nfrm * Ntrc, logical
ntrcW = find(max(coloc) > 0); % index of coloc traces (at least 1 point < coloc_dist)
NtrcW = length(ntrcW);
coloc_enlarged = [zeros(1, size(coloc, 2)); coloc; zeros(1, size(coloc, 2))]; % add zeros for diff
start_coloc = (diff(coloc_enlarged) == 1); % Nfrm+1 * Ntrc, logical
end_coloc = (diff(coloc_enlarged) == -1);
coloc_duration = cell(NtrcW, 1);
N_coloc = zeros(NtrcW, 1);

for ntrc = 1:NtrcW
    coloc_duration{ntrc} = find(end_coloc(:, ntrcW(ntrc))) - find(start_coloc(:, ntrcW(ntrc)));
    N_coloc(ntrc) = length(coloc_duration{ntrc});
    ok = find(coloc(:, ntrcW(ntrc))); % coloc points for this trace
    for nW = 1:N_coloc(ntrc)
        % plot(tab_jR(ok, ntrcW(ntrc)), tab_iR(ok, ntrcW(ntrc)), 'y.', 'markersize', 28)
        if coloc_duration{ntrc}(nW) >= coloc_time_min
            % clr = [1 1 0] + [-0.5 -0.5 0.5]/coloc_duration{ntrc}(nW); % grey to white??
            % plot(tab_jR(ok(nW), ntrcW(ntrc)), tab_iR(ok(nW), ntrcW(ntrc)), '.', 'color', clr, 'markersize', 28)
            plot(tab_jRreg(ok(nW), ntrcW(ntrc)), tab_iRreg(ok(nW), ntrcW(ntrc)), 'o', 'color', [1 1 1]*0.99, 'markersize', coloc_duration{ntrc}(nW)*3)
        end
    end
    fprintf([repmat('\b', 1, 11) '%5i/%5i'], ntrc, NtrcW), drawnow expose
end
fprintf([repmat('\b', 1, 11) '%5i/%5i\r'], NtrcW, NtrcW)

coloc_duration = cell2mat(coloc_duration);
N_coloc = sum(N_coloc);

%%%