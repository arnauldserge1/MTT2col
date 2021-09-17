function [N_coloc, coloc_duration] = traj_xy_Dv_2colors(filename, D, v)%, codage, dirname, D, v)

% function ok = traj_xyc_2colors(filename, D, v)
% plot all traces on first or transm. image, color coded 
% for side (green/red) and motion (conf, Brownian, linear, mix)
% and coloc as white circles
%
% see also traj_xyc, coloc3, cartobyf


global PARAM_I PARAM_J PARAM_ALPHA N_PARAM
global coloc_dist_max coloc_time_min

if isempty(coloc_dist_max), coloc_dist_max = 2; end % pxl
if isempty(coloc_time_min), coloc_time_min = 2; end % frm
if isempty(N_PARAM), MTTparams_def; end
params_def = MTTparams_def; dirname = params_def{4};

if (nargin < 1)
    files = dir2('*.tif');
    if isempty(files), files = dir('*.stk'); end
    filename = files(1).name;
end
if (nargin < 2), [~, ~, D, v] = fit_directed_motion4(filename, 0); end

if isempty(filename), disp('No data... Check dir & filename !'), return, end

[~, ~, ~, ~, ~, v_threshold, D_threshold] = fit_directed_motion4('');
string_reg = 'reg';

%% data
filename_full = [dirname filesep filename '_tab_param.mat'] ;
tab_param = importdata(filename_full);% % ok = ~isempty(tab_param);

DIC_image(filename, dicname(filename), 0, 0, 0, 'left', 1, 0);

im1 = imread(filename, 1);
[tab_param_side{1}, tab_param_side{2}] = split_params_left_right(tab_param, size(im1,2)/2);
side = {'Left' 'Right'};

%%          conf,   Brownian, linear,  mix
color{1} = [0 0.2 0; 0 0.6 0; 0 1 0.2; 0 1 0.8]; % green variants
color{2} = [0.2 0 0; 0.6 0 0; 1 0 0.2; 1 0 0.8]; % red variants

%% left (green) & right (red)
for ns = 1:2
    tab_i = tab_param_side{ns}(PARAM_I-1:N_PARAM:end ,:);
    tab_j = tab_param_side{ns}(PARAM_J-1:N_PARAM:end, :);
    tab_alpha = tab_param_side{ns}(PARAM_ALPHA-1:N_PARAM:end, :);
    tab_i(tab_alpha==0) = nan;
    tab_j(tab_alpha==0) = nan;
    
    %% apply reg with mean_tform to i,j Red => i,j Red reg (Rr) 11/5/2017
    if (ns == 2)
        if ~isempty(dir(['dic' filesep 'reg' filesep 'mean_tform.mat']))
            mean_tform = importdata(['dic' filesep 'reg' filesep 'mean_tform.mat']);
            [tab_j, tab_i] = transformPointsForward(mean_tform, tab_j, tab_i); % x,y == j,i and transformPointsForward expect x,y, hence j,i !!!
        else
            disp(['! Caution, couldn''t find reg' filesep 'mean_tform.mat, red traces not registered !'])
        end
    end

    %% --- go through traces ---
    Ntrc = size(tab_param_side{ns}, 2);
    disp([' ' side{ns} ', traj:            '])

    for nt = 1:Ntrc
        ind = (tab_i(:, nt) > 0);
        %% *** color code, according to mean v & mean D ***
        if (v(nt) <= v_threshold) && (D(nt) <= D_threshold), c = 1; % conf
        elseif (v(nt) <= v_threshold) && (D(nt) > D_threshold), c = 2; % Brownian
        elseif (v(nt) > v_threshold) && (D(nt) <= D_threshold), c = 3; % lin
        elseif (v(nt) > v_threshold) && (D(nt) > D_threshold), c = 4; % mix
        end
        
        plot(tab_j(ind, nt), tab_i(ind, nt), 'color', color{ns}(c, :)) % j i c % with blink??
        
        if (mod(nt, 10) == 0)
            fprintf([repmat('\b', 1, 11) '%5i/%5i'], nt, Ntrc)
            drawnow expose
        end
    end % for itrc = 1:ntrc
    fprintf([repmat('\b', 1, 11) '%5i/%5i'], Ntrc, Ntrc)
end % side
axis image

%% compute or load dist for coloc
dist_file = ['Red-Green distances' filesep filename(1:end-4) '_RGdist' string_reg '.mat'];
if isempty(dir(dist_file))
    RGdist = coloc3(filename, 0, 0, string_reg); % reg
    RGdist = RGdist{1}; % close(gcf)
else
    RGdist = importdata(dist_file, 'dd');
end

%% coloc - white circles
disp(' Traj with coloc:            ')
coloc = (RGdist < coloc_dist_max); % Nfrm * Ntrc, logical
ntrcW = find(max(coloc) > 0); % index of coloc traces (at least 1 point < coloc_dist)
NtrcW = length(ntrcW);
coloc_enlarged = [zeros(1, size(coloc, 2)); coloc; zeros(1, size(coloc, 2))]; % add zeros for diff
start_coloc = (diff(coloc_enlarged) == 1); % Nfrm+1 * Ntrc, logical
end_coloc = (diff(coloc_enlarged) == -1);
coloc_duration = cell(NtrcW, 1);
N_coloc = zeros(NtrcW, 1); % nota: coloc defined for red positions, presumably registered (tab_iRreg = tab_i, see traj_xy_coloc)

for ntrc = 1:NtrcW
    coloc_duration{ntrc} = find(end_coloc(:, ntrcW(ntrc))) - find(start_coloc(:, ntrcW(ntrc)));
    N_coloc(ntrc) = length(coloc_duration{ntrc});
    ok = find(coloc(:, ntrcW(ntrc))); % coloc points for this trace
    for nW = 1:N_coloc(ntrc)
        if (coloc_duration{ntrc}(nW) >= coloc_time_min)
            plot(tab_j(ok(nW), ntrcW(ntrc)), tab_i(ok(nW), ntrcW(ntrc)), 'o', 'color', [1 1 1]*0.99, 'markersize', coloc_duration{ntrc}(nW)*3) % nota: scaling *3 is aritrary
        end
    end
    fprintf([repmat('\b', 1, 11) '%5i/%5i'], ntrc, NtrcW), drawnow expose
end

coloc_duration = cell2mat(coloc_duration);
N_coloc = sum(N_coloc);
fprintf(' - Found %i coloc events\r', N_coloc)
%%%