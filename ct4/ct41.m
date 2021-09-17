function [mean_val, var_name, tab_val] = ct41(data_in, fast)

% [mean_val, output_val, tab_val] = ct41(data_in)
%
% construit les histogrammes (ou le MSD) de différentes valeurs de fit MTT 
%  var_name = {'alpha (a.u.)' 'SNR (db)' 'r^2 (pxl^2)' 'TrcLen (img)' 'asym' ...
%         'accuracy (pxl)' 'D by trc (pxl^2/img)' 'Dglobal (pxl^2/img)' 'gamma'};
% AS 2013

global N_PARAM PARAM_ALPHA PARAM_SIG2

if nargin<1, files = dir('*.tif'); data_in = files(1).name; end
if nargin<2, fast = 1; end

if strfind(cd, 'output'), cd .., end % si on est déjà ds output..

%% ** load data
if ischar(data_in)
    filename = data_in;
    filename_full = ['output23' filesep filename '_tab_param.mat'] ;
    if ~isempty(dir(filename_full))
        tab_param = importdata(filename_full);
        fprintf('%s ', filename)
    else
        tab_param = [];
    end
else
    filename = '';
    tab_param = data_in;
end

if isempty(tab_param)
    mean_val = [];
    var_name = '';
    tab_val = cell(0);
    disp('No data. wrong dir???');
    return
end

%% ** prepare data
alpha = tab_param(PARAM_ALPHA-1:N_PARAM:end, :); % intensty of fitted SM
noise = sqrt(tab_param(PARAM_SIG2-1:N_PARAM:end, :)); %tab_var(5:N_PARAM:end, :);% as 4/4/7 SIGMA2
SNR = 20*log10(alpha./noise); % puis plot Log

TrcLen = sum(alpha>0)'; % number of images in each trace, AS 11/12/7
Nfrm = tab_param(end-N_PARAM+1,1);
r2 = calcul_r2(tab_param);
asym = radius_gyration(tab_param, 0);

%% ---  MSD & D --- % calculate mean square displacement & diff
if ~fast
    [msddata, ~, msdg] = msd(detect_reconnex_to_trc(tab_param));
    D = calculDinst(msddata); % [#trace, D, erreurD, offset]
    offset = D(:,4);
    offset_min = 0.01; % very small values??
    accuracy = sqrt(offset(offset>offset_min))/2; % offset = 2 sig^2
    [~, gamma] = fit_anomal2(msddata, offset);
    
    data = {alpha(:), SNR(:), r2, TrcLen, asym, ...
        accuracy, D(:,2), msdg, gamma};
    var_name = {'alpha (a.u.)' 'SNR (db)' 'r^2 (pxl^2)' 'TrcLen (img)' 'asym' ...
        'accuracy (pxl)' 'D by trc (pxl^2/img)' 'Dglobal (pxl^2/img)' 'gamma'};
else
    data = {alpha(:), SNR(:), r2, TrcLen, asym};
    var_name = {'alpha (a.u.)' 'SNR (db)' 'r^2 (pxl^2)' 'TrcLen (img)' 'asym'};
end

tab_val = cell(size(var_name, 2), 1);

nSubImg = size(var_name, 2);
Golden = (1+sqrt(5))/2;
nSubX = floor(sqrt(nSubImg*Golden)); nSubY = ceil(nSubImg/nSubX);

%% graphs
mean_val = zeros(1, nSubImg+2); % +2 pour Nb trc pk evt f/c
figure('WindowStyle', 'docked');

for v = 1:nSubImg
    subplot(nSubX, nSubY, v)
    
    if strncmp(var_name{v}, 'Dg', 2)
        mean_val(v) = ct4plot_msd(data{v}); % => compute D
    else
        if any(strncmp(var_name{v}, {'Trc' 'r^2' 'x^2' 'D b'}, 3))
            do_log = 1;
        else
            do_log = 0;
        end
        mean_val(v) = ct4plot_hist(data{v}, var_name{v}, do_log, Nfrm);
    end
    
    if v==1, title({cd, filename}, 'Interpreter', 'none'); end
    
    if nargout>2, tab_val{n_output} = data{v}; end
end

mean_val(nSubImg+1) = length(alpha(alpha>0));
mean_val(nSubImg+2) = length(TrcLen);
var_name(nSubImg+1:nSubImg+2) = {'Npk/frm' 'Ntrc'};

%%%