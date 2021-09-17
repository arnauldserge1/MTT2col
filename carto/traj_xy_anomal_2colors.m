function [N_coloc, coloc_duration] = traj_xy_anomal_2colors(filename, D, gamma)

% function ok = traj_xy_anomal_2colors(filename, D, gamma)
% plot all traces on first or transm. image, color coded 
% for side (green/red) and motion: conf, Brownian, linear, according to anomalous fit MSD = 4Dt^gamma + 2sig^2
% and coloc as white circles
%
% see also traj_xyc, coloc3, cartobyf


global PARAM_I PARAM_J PARAM_ALPHA N_PARAM
global coloc_dist_max coloc_time_min

if isempty(coloc_dist_max), coloc_dist_max = 2; end % pxl
if isempty(coloc_time_min), coloc_time_min = 2; end % frm
if isempty(N_PARAM), MTTparams_def; end
params_def = MTTparams_def; dirname = params_def{4};

%% data

filename_full = [dirname filesep filename '_tab_param.mat'] ;
tab_param = importdata(filename_full);% % ok = ~isempty(tab_param);

if (nargin < 1)
    files = dir2('*.tif');
    if isempty(files), files = dir('*.stk'); end
    filename = files(1).name;
end

if (nargin < 2) || isempty(D)
    msddata = msd(detect_reconnex_to_trc(tab_param)); % left, then right
    [D, gamma] = fit_anomal2(msddata,'',0);
    D = D'; gamma = gamma';
end

if isempty(filename), disp('No data... Check dir & filename !'), return, end

string_reg = 'reg';

DIC_image(filename, dicname(filename), 0, 0, 0, 'left', 1, 0);

im1 = imread(filename, 1);
[tab_param_side{1}, tab_param_side{2}] = split_params_left_right(tab_param, size(im1,2)/2);
side = {'Left' 'Right'};

%%          conf,   Brownian, linear
color{1} = [0 0.3 0; 0 1 0; 0.7 1 0.7]; % green variants
color{2} = [0.3 0 0.3; 1 0 1;1 0.7 1]; % magenta variants

%% left (green) & right (red)
Ntrc = zeros(1,2);
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
    
    %% load data from Brownian simul
    if ispc, load('\\10.36.1.25\spinning\Oksana\data for article1\non blocking 19H36\anomal_data\carto_data.mat', 'gamma_mean_bro', 'two_sigma_gamma_bro')
    else, load('/volumes/spinning/Oksana/data for article1/non blocking 19H36/anomal_data/carto_data.mat', 'gamma_mean_bro', 'two_sigma_gamma_bro')
    end
    
    %% param for selecting lin, bro, conf from gamma & correction for D
%     alpha = gamma - gamma_mean_bro - log10(D/Dmean_bro) * slope_bro ;

    %% --- go through traces ---
    Ntrc(ns) = size(tab_param_side{ns}, 2);
    disp([' ' side{ns} ', traj:            '])
    for nt = 1:Ntrc(ns)
        if (ns == 2), nt2 = nt + Ntrc(1); else, nt2 = nt; end % red trc numbers shifted!
        ind = (tab_i(:, nt) > 0);
        
        %% *** color code, according to alpha & alpha1_2 ***  
%         
% %         if alpha(nt) < -alpha1_2_bro, c = 1; % conf
%         if alpha(nt2) < -alpha1_2_bro_nomansland, c = 1; % && alpha_C > - alpha1_2_bro, c = 2; % low conf
%         elseif alpha(nt2) < alpha1_2_bro_nomansland && alpha(nt2) > -alpha1_2_bro_nomansland, c = 2; % Brownian
%         elseif alpha(nt2) > alpha1_2_bro_nomansland, c = 3; % c = 4; % low directed   
% %         elseif alpha(nt) > alpha1_2_bro, c = 5; % directed
%         end
%         plot(tab_j(ind, nt), tab_i(ind, nt), 'color', color{ns}(c, :)) % j i c % with blink??
%         if (mod(nt, 10) == 0) % fprintf([repmat('\b', 1, 11) '%5i/%5i'], nt, Ntrc(ns)) % drawnow expose % end
      
        %% *** color code, according to gamma ***
        
        if gamma(nt2) < gamma_mean_bro - two_sigma_gamma_bro, c = 1;
        elseif gamma(nt2) < gamma_mean_bro + two_sigma_gamma_bro && gamma(nt2) > gamma_mean_bro - two_sigma_gamma_bro, c = 2 ;
        elseif gamma(nt2) > gamma_mean_bro + two_sigma_gamma_bro, c = 3; 
        end 
        
        plot(tab_j(ind, nt), tab_i(ind, nt), 'color', color{ns}(c, :)) % j i c % with blink??
        
        if (mod(nt, 10) == 0)
            fprintf([repmat('\b', 1, 11) '%5i/%5i'], nt, Ntrc(ns))
            drawnow expose
        end

    end % for itrc = 1:ntrc
    fprintf([repmat('\b', 1, 11) '%5i/%5i'], Ntrc(ns), Ntrc(ns))
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