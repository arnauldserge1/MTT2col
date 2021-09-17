function [N_coloc, coloc_duration] = traj_xy_coloc(file, dirname, string_reg)

% function [N_coloc, coloc_duration] = traj_xy_coloc(file, dirname, string_reg)
%
% plot all traces on first or transm. image
% color red, green or yellow
%
% see also coloc3, cartobyf


global PARAM_I PARAM_J PARAM_ALPHA N_PARAM
global coloc_dist_max coloc_time_min

if isempty(coloc_dist_max), coloc_dist_max = 2; end % pxl
if isempty(coloc_time_min), coloc_time_min = 2; end % frm
if isempty(N_PARAM), MTTparams_def; end

[pxl_size, time_lag] = get_calib3;
fprintf('assuming coloc for distance below %g pxl = %g nm\r', coloc_dist_max, coloc_dist_max*pxl_size*1000) % cf. d_min = 0.6 lambda/NA = 0.6*580/1.45 = 240 nm
fprintf('assuming coloc for duration above %g frm = %g s\r', coloc_dist_max, coloc_dist_max*time_lag)


if nargin < 1, files = dir('*.tif'); if isempty(files), files = dir('*.stk'); end, file = files(1).name; end
if nargin < 2, params_def = MTTparams_def; dirname = params_def{4}; end
if nargin < 3, string_reg = 'reg'; end

if isempty(file), disp('No data... Check dir & filename !'), return, end

%% data
filename_full = [dirname filesep file '_tab_param.mat'];
if isempty(dir(filename_full)), disp('no data??'), N_coloc = []; coloc_duration = []; return, end
tab_param = importdata(filename_full);
tab_i = tab_param(PARAM_I-1:N_PARAM:end,:);
tab_j = tab_param(PARAM_J-1:N_PARAM:end,:);
tab_alpha = tab_param(PARAM_ALPHA-1:N_PARAM:end,:);
tab_i(tab_alpha == 0) = nan;
tab_j(tab_alpha == 0) = nan;
Tmax = size(tab_i, 1);

%% Red/Green
img1 = imread(file, 1); % img1 = tiffread(file, 1); 28/4/17
middle = size(img1, 2)/2; % 512, a priori

mj = nanmean(tab_j); % not Mickael Jackson!! ;-)
tab_iG = tab_i(:, mj <= middle);
tab_jG = tab_j(:, mj <= middle);
tab_alphaG = tab_alpha(:, mj <= middle);
tab_iR = tab_i(:, mj > middle);
tab_jR = tab_j(:, mj > middle) - middle;
tab_alphaR = tab_alpha(:, mj > middle);

% if isdir('max'), show_max_image(file)
DIC_image(file, dicname(file), 0, 0, 0, '_right'); % figure... imagesc(pict), axis image ij off, colormap('gray'), hold on

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
        if t < Tmax/2, clr = [0.5 1 0.5] - [1 0 1]*(t/Tmax);
        else, clr = [0 1.5 0] - [0 1 0]*(t/Tmax);
        end
        % clr = [1/2-tt/Tmax; 3/2-tt/Tmax; 1/2-tt/Tmax]; clr = min(max(clr, 0), 1);
        % patch('XData', [nan; tab_jG; nan], 'YData', [nan; tab_iG(t, 1); nan], 'CData', [nan; clr'; nan], 'FaceColor', 'interp', 'EdgeColor', 'interp')
        plot(tab_jG(tt, ntrc), tab_iG(tt, ntrc), 'color', clr)
        % plot(tab_jG(ok, ntrc), tab_iG(ok, ntrc), 'color', clr)
    end
    if mod(ntrc, 10) == 0, fprintf([repmat('\b', 1, 11) '%5i/%5i'], ntrc, NtrcG), drawnow expose, end
end
fprintf([repmat('\b', 1, 11) '%5i/%5i\r'], NtrcG, NtrcG)

%% Red
NtrcR = size(tab_iR, 2);
disp('red traj:            ')
for ntrc = 1:NtrcR
    ok = find(tab_alphaR(:, ntrc) > 0);
    for nf = 1:length(ok)-1
        t = ok(nf);
        tt = [ok(nf) ok(nf+1)];
        if t < Tmax/2, clr = [1 0.5 0.5] - [0 1 1]*(t/Tmax);
        else, clr = [1.5 0 0] - [1 0 0]*(t/Tmax);
        end
        plot(tab_jR(tt, ntrc), tab_iR(tt, ntrc), 'color', clr)
    end
    if mod(ntrc, 10) == 0, fprintf([repmat('\b', 1, 11) '%5i/%5i'], ntrc, NtrcR), drawnow expose, end
end
fprintf([repmat('\b', 1, 11) '%5i/%5i\r'], NtrcR, NtrcR)

%% Yellow
disp('yellow coloc event:            ')
coloc = (RGdist < coloc_dist_max); % Nfrm * Ntrc, logical
ntrcY = find(max(coloc) > 0); % index of coloc traces (at least 1 point < coloc_dist)
NtrcY = length(ntrcY);
coloc_enlarged = [zeros(1, size(coloc, 2)); coloc; zeros(1, size(coloc, 2))]; % add zeros for diff
start_coloc = (diff(coloc_enlarged) == 1); % Nfrm+1 * Ntrc, logical
end_coloc = (diff(coloc_enlarged) == -1);
coloc_duration = cell(NtrcY, 1);
N_coloc = zeros(NtrcY, 1);

for ntrc = 1:NtrcY
    coloc_duration{ntrc} = find(end_coloc(:, ntrcY(ntrc))) - find(start_coloc(:, ntrcY(ntrc)));
    N_coloc(ntrc) = length(coloc_duration{ntrc});
    ok = find(coloc(:, ntrcY(ntrc))); % coloc points for this trace
    for nY = 1:N_coloc(ntrc)
        % plot(tab_jR(ok, ntrcY(ntrc)), tab_iR(ok, ntrcY(ntrc)), 'y.', 'markersize', 28)
        if coloc_duration{ntrc}(nY) >= coloc_time_min
            % clr = [1 1 0] + [-0.5 -0.5 0.5]/coloc_duration{ntrc}(nY); % grey to yellow
            % plot(tab_jR(ok(nY), ntrcY(ntrc)), tab_iR(ok(nY), ntrcY(ntrc)), '.', 'color', clr, 'markersize', 28)
            plot(tab_jR(ok(nY), ntrcY(ntrc)), tab_iR(ok(nY), ntrcY(ntrc)), 'yo', 'markersize', coloc_duration{ntrc}(nY)*3)
        end
    end
    fprintf([repmat('\b', 1, 11) '%5i/%5i'], ntrc, NtrcY), drawnow expose
end
fprintf([repmat('\b', 1, 11) '%5i/%5i\r'], NtrcY, NtrcY)

coloc_duration = cell2mat(coloc_duration);
N_coloc = sum(N_coloc);

%%%