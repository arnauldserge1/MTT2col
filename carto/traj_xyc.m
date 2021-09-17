function ok = traj_xyc(filename, codage, dirname, D, v)

% function ok = traj_xyc(filename, codage, dirname)
%
% plot all traces on first or transm. image
% co:lor code by trace, not by step (cf traj_xy2)

global PARAM_I PARAM_J PARAM_ALPHA N_PARAM
if isempty(N_PARAM), MTTparams_def; end

if nargin<1, files = dir('*.tif'); if isempty(files), files = dir('*.stk'); end, filename = files(1).name; end
if nargin<2, codage = 'SCI'; end
if nargin<3, params_def = MTTparams_def; dirname = params_def{4}; end
if nargin<5, if strcmp(codage, 'Dv'), [~, ~, D, v] = fit_directed_motion4(filename); end, end

if isempty(filename), disp('No data... Check dir & filename !'), return, end

if strcmp(codage, 'Dv')
    [~, ~, ~, ~, ~, v_threshold, D_threshold] = fit_directed_motion4('');
end
SCI_threshold = 0.6;
v_min = 0.01;

%% data
filename_full = [dirname filesep filename '_tab_param.mat'] ;
tab_param = importdata(filename_full);
tab_i = tab_param(PARAM_I-1:N_PARAM:end ,:);
tab_j = tab_param(PARAM_J-1:N_PARAM:end, :);
tab_alpha = tab_param(PARAM_ALPHA-1:N_PARAM:end, :);
tab_i(tab_alpha==0) = 0;
tab_j(tab_alpha==0) = 0;

ok = ~isempty(tab_param);

DIC_image(filename, dicname(filename), 0, 0, 0, '', 1, 0);

%% --- go through traces ---
Ntrc = size(tab_param, 2);
disp('traj :            ')
color = [0 0 1; 0 1 0; 1 0 0; 1 .5 0]; % B G R O


for nt = 1:Ntrc
    
    ind = (tab_i(:, nt) > 0);
    %% *** codage couleur ***
    if strcmp(codage, 'Dv')
        if (v(nt) > v_threshold) && (D(nt) <= D_threshold), c = 1; % lin
        elseif (v(nt) > v_threshold) && (D(nt) > D_threshold), c = 2; % fast
        elseif (v(nt) <= v_threshold) && (D(nt) <= D_threshold), c = 3; % conf
        elseif (v(nt) <= v_threshold) && (D(nt) > D_threshold), c = 4; % Brownian
        end
        plot(tab_j(ind, nt), tab_i(ind, nt), 'color', color(c, :)) % j i c
    elseif strcmp(codage, 'SCI')
        trc = detect_reconnex_to_trc(tab_param(:, nt));
        plotwithblink(trc, 'b')
        
        [v, ~, ~, ~, ~, corr_lin] = compute_speed_from_SCI(tab_param(:, nt), 0, SCI_threshold);
        
        ind_lin = (corr_lin{1}{1} > SCI_threshold);
        start_lin = find(diff([0 ind_lin]) == 1);
        end_lin = find(diff([ind_lin 0]) == -1);
        Nlin = length(start_lin);
        
        %% loop over linear events
        for nl = 1:Nlin
            if v(nl) > v_min, color = 'g';
            else color = 'k';
            end
            plotwithblink(trc(start_lin(nl):end_lin(nl), :), color)
        end
    end
    
    if mod(nt, 10) == 0
        fprintf([repmat('\b', 1, 11) '%5i/%5i'], nt, Ntrc)
        drawnow expose
    end
end % for itrc = 1:ntrc

axis image % colormap(hot), colorbar cf nanopict
fprintf([repmat('\b', 1, 11) '%5i/%5i\r'], Ntrc, Ntrc)


%% plotwithblink %%

function plotwithblink(trc, clr)
dt = diff(trc(:, 2));
t_blink = find(dt>1);
t_blink = [0; t_blink; size(trc, 1)];
for i=1:length(t_blink)-1
    tt = t_blink(i)+1:t_blink(i+1);
    plot(trc(tt, 3), trc(tt, 4), 'color', clr)
    if i<length(t_blink)-1
        tbl = [t_blink(i+1), t_blink(i+1)+1];
        plot(trc(tbl, 3), trc(tbl, 4), 'color', clr, 'linestyle', ':')
    end
end
% % % 
% % % %% plotwithoutblink %%
% % % 
% % % function plotwithoutblink(trc, clr)
% % % dt = diff(trc(:, 2));
% % % t_blink = find(dt>1);
% % % t_blink = [0; t_blink; size(trc, 1)];
% % % for i=1:length(t_blink)-1
% % %     tt = t_blink(i)+1:t_blink(i+1);
% % %     plot(trc(tt, 3), trc(tt, 4), 'color', clr)
% % % end