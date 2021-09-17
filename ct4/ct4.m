function [par, varname, output_values] = ct4(data_in, timelag, conf_method, dirname, output_variables)

% [par, varname, output_values] = ct4(data_in, timelag, conf_method, dirname, output_variables)
%
% construit les histogrammes (ou le MSD) de différentes valeurs de fit MTT 
% 
% ex: % ct4(file, 36, 'rel')

global N_PARAM PARAM_ALPHA PARAM_SIG2
pdef = MTTparams_def;

if nargin<5, output_variables = []; end
if nargin<4, dirname = pdef{4}; end
if nargin<3, conf_method = pdef{20}; end
if nargin<1
    files = dir('*.stk');
    if isempty(files), files = dir('*.tif'); end
    data_in = files(1).name;
end

if strfind(cd, 'output'), cd .., end % si on est déjà ds output..
if any(strcmp(conf_method, {'abs', 'rel'})), Nfigs = 2; else Nfigs = 1; end  % global + conf?

%% ** load data
if ischar(data_in)
    filename = data_in;
    filename_full = [dirname filesep filename '_tab_param.mat'] ;
    if ~isempty(dir(filename_full))
        load(filename_full) %  tab_param = fread_params_timewindow(filename_full, 1, 1, t_last);
        fprintf('%s - ', filename)
    else
        tab_param = [];
    end
else
    filename = '';
    tab_param = data_in;
end

if isempty(tab_param)
    par = [];
    varname = '';
    output_values = cell(length(output_variables),1);
    disp('No data. wrong dir???');
    return
end

if nargin<2
    timelag = eval_timing(filename);
    disp(['time lag = ' num2str(timelag) ' ms '])
end

%% ** prepare data
alpha = tab_param(PARAM_ALPHA-1:N_PARAM:end, :); % intensty of fitted SM
noise = sqrt(tab_param(PARAM_SIG2-1:N_PARAM:end, :)); %tab_var(5:N_PARAM:end, :);% as 4/4/7 SIGMA2
SNR = 20*log10(alpha./noise); % puis plot Log
TrcLen = sum(alpha>0)'; % number of images in each trace, AS 11/12/7
Nfrm = tab_param(end-N_PARAM+1,1);

%% --- data, MSD & D --- % calculate mean square displacement & diff, global, conf free
r2 = calcul_r2(tab_param);
nSubImg = 8;
trc = detect_reconnex_to_trc(tab_param);
[msddata, ~, msdg] = msd(trc); % ,1,1,0,0,0,0 % compute sd???
D = calculDinst(msddata);
[~, gamma, offset] = fit_anomal2(msddata); % fit_anomal3??
offset_min = 0;%.01; % very small values??
accuracy = sqrt(offset(offset>offset_min))/2;

data = {TrcLen alpha(:) accuracy SNR(:) r2 D(:,2) msdg gamma}; % noise(:) 
varname = {'TrcLen (img)' 'alpha (a.u.)' 'accuracy (pxl)' 'SNR (dB)' ...'noise (a.u.)' 
    'r^2 (pxl^2)' 'D by trc (pxl^2/img)' 'Dglobal (pxl^2/img)' 'gamma'};
output_values = cell(length(output_variables),1);
n_output = 1;

%% confinement
if any(strcmp(conf_method, {'abs', 'rel'})) % ~isempty(conf_method)
    [~, Tc, Tf, R2c, R2f, trc_c, trc_f] = ...
        MTT_probaconf(tab_param, conf_method, timelag, 0, filename); % graph=0
    fprintf('confined events: ')
    [~, ~, msdc2] = msd(trc_c, 1, 1, 0, 0, 1);
    r2c = calcul_r2(trc_c); % Dc = calculDinst(msdc);
    fprintf('free events: ')
    [~, ~, msdf2] = msd(trc_f, 1, 1, 0, 0, 1);
    r2f = calcul_r2(trc_f); % Df = calculDinst(msdf);
    data(2, 1:8) = {Tc Tf sqrt(R2c) sqrt(R2f) r2c(:) r2f(:) msdc2 msdf2};
    varname(2, 1:10) = {'Tc (img)' 'Tf (img)' 'Rc (pxl)' 'Rf (pxl)' ...
     'r^2c (pxl^2)' 'r^2f (pxl^2)' 'Dc (pxl^2/img)' 'Df (pxl^2/img)' 'Nc' 'Nf'};
end

dirtitle = [' dir=', num2str(dirname)];

plottitle = {cd, filename};
if ~any(strcmp(conf_method, {'abs', 'rel'})), conf_str = ''; % ' no conf. detection'; 
else conf_str = [' conf.meth=' conf_method]; end
plottitle2 = {date, [dirtitle, conf_str]};

%% graphs
par = zeros(Nfigs, nSubImg+2); % +2 pour Nb trc pk evt f/c
n_graph = 1;
for u = 1:Nfigs
    figure('WindowStyle', 'docked');
    
    for v = 1:nSubImg
        subplot(nSubImg/2, 2, v);

        if any(strncmp(varname{u, v}, {'Dg' 'Dc' 'Df'}, 2)) % if all([u v]==[1 7]) || all([u v]==[2 5]) || all([u v]==[2 6]) 
            par(u, v) = ct4plot_msd(data{u, v}); % => compute D
        else
            if any(strncmp(varname{u, v}, {'Trc' 'Tc ' 'Tf ' 'r^2' 'x^2' 'D b'}, 3)) % 'gam'
                do_log = 1; 
            else
                do_log = 0;
            end
            par(u, v) = ct4plot_hist(data{u, v}, varname{u, v}, do_log, Nfrm);
        end
        
        if v==1, title(plottitle, 'Interpreter', 'none');
        elseif v==2, title(plottitle2, 'Interpreter', 'none'); end
        
        if find(output_variables==n_graph), output_values{n_output} = data{u, v}; n_output = n_output+1; end
        n_graph = n_graph+1;
    end
end

par(1, nSubImg+1) = length(alpha(alpha>0));
par(1, nSubImg+2) = length(TrcLen);
varname(1, nSubImg+1:nSubImg+2) = {'Npk/frm' 'Ntrc'};

if any(strcmp(conf_method, {'abs', 'rel'}))
    par(2, nSubImg+1) = length(Tc);
    par(2, nSubImg+2) = length(Tf);
end
%%%