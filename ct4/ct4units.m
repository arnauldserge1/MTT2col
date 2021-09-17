function [par, varname, output_values] = ct4units(data_in, time_lag, pxl_size, conf_method, dirname, output_variables, side)

% [par, varname, output_values] = ct4units(data_in, timelag, pixelsize, conf_method, dirname, output_variables, side)
%
% construit les histogrammes (ou le MSD) de différentes valeurs de fit MTT
%
% ex: % ct4(file, 36, 'rel')

global N_PARAM PARAM_ALPHA PARAM_SIG2 PARAM_I PARAM_J

pdef = MTTparams_def;

if nargin<7, side = ''; end
if nargin<6, output_variables = []; end
if nargin<5, dirname = pdef{4}; end
if nargin<4, conf_method = pdef{20}; end
if nargin<3, [pxl_size, time_lag] = get_calib3; end % if nargin<3, pxl_size = str2double(pdef{24}), fprintf('Using default pixel size: %g um\r', pxl_size), end % 0.16if nargin<2, time_lag = str2double(pdef{25}); fprintf('Using default time lag: %g s\r', time_lag), end  %   timelag = eval_timing(filename);

if nargin<1
    files = dir2('*.stk');
    if isempty(files), files = dir2('*.tif'); end
    data_in = files(1).name;
    fprintf('Using first file: %s\r', files(1).name)
end

if contains(cd, 'output'), cd .., end % si on est déjà ds output..
if any(strcmp(conf_method, {'abs', 'rel'})), Nfigs = 2; else, Nfigs = 1; end  % global + conf?
nSubImg = 8;

%% ** load data
if ischar(data_in)
    filename = data_in;
    
    if contains(filename, '*') % all files: expecting a name such as '*.tif'          strcmp(filename, 'all files    ')
        tab_param = load_all_param(filename);
    else
        filename_full = [dirname filesep filename '_tab_param.mat'];
        if ~isempty(dir(filename_full))
            load(filename_full, 'tab_param') %  tab_param = fread_params_timewindow(filename_full, 1, 1, t_last);
            fprintf('%s - ', filename)
        else
            tab_param = [];
        end
    end
else % giving directly tab_param, not filename
    filename = '';
    tab_param = data_in;
end

if ~isempty(filename) %%% && ~contains(filename, '*')
    if ~contains(filename, '*') && ~isempty(side)
        img1 = imread(filename, 1);
        middle = size(img1, 2)/2;
    else
        middle = 512; 
        if ~isempty(side), disp('!!assuming middle @ 512 for all images!!'), end % or file1 = .... img1 = imread(file1, 1);????
    end
    tab_param = split_params_left_right(tab_param, middle, side); % nota: no change if side == ''
end

Ntrc = size(tab_param, 2); % length(TrcLen);

if isempty(tab_param)
    par = [];
    varname = '';
    output_values = cell(length(output_variables),1);
    disp(['No data for ' filename ' (wrong dir???)'])
    return
end

Ntrc_max = 10^4; % hence even if relative sd~1, sem~1% (otherwise, computation time may overcome several hours..)
if Ntrc > Ntrc_max
    ind_kept = round(linspace(1, Ntrc, Ntrc_max)); % homogenous sample
    fprintf('wooo, found %i trajs!! Keeping only %i\r', Ntrc, Ntrc_max)
    tab_param = tab_param(:, ind_kept); % AS 17/12/2014
end

%% ** prepare data
tab_param(PARAM_I-1:N_PARAM:end, :) = tab_param(PARAM_I-1:N_PARAM:end, :)*pxl_size;
tab_param(PARAM_J-1:N_PARAM:end, :) = tab_param(PARAM_J-1:N_PARAM:end, :)*pxl_size;

alpha = tab_param(PARAM_ALPHA-1:N_PARAM:end, :); % intensty of fitted SM
noise = sqrt(tab_param(PARAM_SIG2-1:N_PARAM:end, :)); %tab_var(5:N_PARAM:end, :);% as 4/4/7 SIGMA2
SNR = 20*log10(alpha./noise); % puis plot Log
TrcLen = sum(alpha>0)'*time_lag; % number of images in each trace, AS 11/12/7
Nfrm = tab_param(end-N_PARAM+1,1);

%% --- data, MSD & D --- % calculate mean square displacement & diff, global, conf free
r2 = calcul_r2(tab_param);
trc = detect_reconnex_to_trc(tab_param);
[msddata, ~, msdg] = msd(trc);
D = calculDinst(msddata); % D by trc
D = D(:,2)/time_lag; % um2/s
[Da, gamma, offset] = fit_anomal2(msddata);
% % % gamma = gamma(gamma>-1 & gamma<3)

if ~contains(filename, '*')
    dataname = ['ct4_by_file' filesep filename(1:end-4)  '_D_gamma_' side '.mat'];
    save(dataname, 'Da', 'gamma');
end

offset_min = 0.01*pxl_size^2; % very small values??
accuracy = sqrt(offset(offset>offset_min))/2;

data = {TrcLen alpha(:) accuracy SNR(:) r2 D msdg gamma}; % noise(:)
varname = {'TrcLen (s)' 'alpha (a.u.)' 'accuracy (um)' 'SNR (dB)' ...'noise (a.u.)'
    'r^2 (um^2)' 'D by trc (um^2/s)' 'Dglobal (um^2/s)' 'gamma'};
output_values = cell(length(output_variables),1);
n_output = 1;

%% confinement
if any(strcmp(conf_method, {'abs', 'rel'})) % ~isempty(conf_method)
    timelag_ms = time_lag*1000;
    [~, Tc, Tf, R2c, R2f, trc_c, trc_f] = ...
        MTT_probaconf(tab_param, conf_method, timelag_ms, 0, filename); % graph=0
    fprintf('confined events: ')
    [~, ~, msdc2] = msd(trc_c, 1, 1, 0, 0, 1);
    r2c = calcul_r2(trc_c); % Dc = calculDinst(msdc);
    fprintf('free events: ')
    [~, ~, msdf2] = msd(trc_f, 1, 1, 0, 0, 1);
    r2f = calcul_r2(trc_f); % Df = calculDinst(msdf);
    data(2, 1:8) = {Tc Tf sqrt(R2c) sqrt(R2f) r2c(:) r2f(:) msdc2 msdf2};
    varname(2, 1:10) = {'Tc (s)' 'Tf (s)' 'Rc (um)' 'Rf (um)' ...
        'r^2c (um^2)' 'r^2f (um^2)' 'Dc (um^2/s)' 'Df (um^2/s)' 'Nc' 'Nf'};
end

dirtitle = [' dir=', num2str(dirname)];

if isempty(side), plottitle = {cd, filename};
else, plottitle = {cd, [filename ' ' side]};
end
if ~any(strcmp(conf_method, {'abs', 'rel'})), conf_str = ''; % ' no conf. detection';
else, conf_str = [' conf.meth=' conf_method];
end
plottitle2 = {date, [dirtitle, conf_str]};

%% graphs
par = zeros(Nfigs, nSubImg+2); % +2 pour Nb trc pk evt f/c
n_graph = 1;
for u = 1:Nfigs
    figure('WindowStyle', 'docked');
    
    for v = 1:nSubImg
        subplot(nSubImg/2, 2, v);
        
        if any(strncmp(varname{u, v}, {'Dg' 'Dc' 'Df'}, 2)) % if all([u v]==[1 7]) || all([u v]==[2 5]) || all([u v]==[2 6])
            par(u, v) = ct4plot_msd(data{u, v}, time_lag); % => compute D global, for all traces
        else
            if any(strncmp(varname{u, v}, {'Trc' 'Tc ' 'Tf ' 'r^2' 'x^2' 'D b'}, 3)) % 'gam'
                do_log = 1;
            else
                do_log = 0;
            end
            par(u, v) = ct4plot_hist(data{u, v}, varname{u, v}, do_log, Nfrm); % => compute mean
        end
        
        if v==1, title(plottitle, 'Interpreter', 'none');
        elseif v==2, title(plottitle2, 'Interpreter', 'none'); end
        
        if find(output_variables==n_graph), output_values{n_output} = data{u, v}; n_output = n_output+1; end
        n_graph = n_graph+1;
    end
end

par(1, nSubImg+1) = length(alpha(alpha>0))/Nfrm;
par(1, nSubImg+2) = Ntrc; % length(TrcLen);
varname(1, nSubImg+1:nSubImg+2) = {'Npk/frm' 'Ntrc'};

if any(strcmp(conf_method, {'abs', 'rel'}))
    par(2, nSubImg+1) = length(Tc);
    par(2, nSubImg+2) = length(Tf);
end
%%%