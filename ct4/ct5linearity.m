function [par, varname, output_values] = ct5linearity(data_in, time_lag, pxl_size, dirname, output_variables, new_fig, clr, text_offset, side)

% function [par, varname, output_values] = ct5linearity(data_in, time_lag, pxl_size, dirname, output_variables, new_fig, clr, text_offset, side)
%
% build histograms for various linearity descriptors:
% varname = {'full path (um)', 'direct path (um)', 'instantaneous velocity (um/s)', 'resulting velocity(um/s)', 'linearity', 'asymetry'};
% returns par = mean of each descriptor
% and also returns the full distribution if 'output_variables' are specified (1 for 'full path', etc) 
%
% default: pdef = MTTparams_def;
% ct5linearity("first tif or stk file", str2double(pdef{25}), str2double(pdef{24}), pdef{4}, [], 1, 0)
% AS 2014
% see also ct4, ct5_by_file, eval_asym, calcul_linearity

global N_PARAM PARAM_ALPHA PARAM_I PARAM_J

pdef = MTTparams_def;

if nargin<9, side = ''; end
if nargin<8, text_offset = 0; end
if nargin<7, clr = [0 0 1]; end
if nargin<6, new_fig = 1; end
if nargin<5, output_variables = []; end
if nargin<4, dirname = pdef{4}; end
if nargin<3, pxl_size = str2double(pdef{24}); fprintf('Using default pixel size: %g um\r', pxl_size), end % 0.16 um
if nargin<2, time_lag = str2double(pdef{25}); fprintf('Using default time lag: %g s\r', time_lag), end % Attention, en secondes, cf timelag = eval_timing(filename) en ms!!

if nargin<1
    files = dir2('*.stk');
    if isempty(files), files = dir2('*.tif'); end
    data_in = files(1).name;
    fprintf('Using first file: %s\r', files(1).name)
end

if contains(cd, 'output'), cd .., end % si on est déjà ds output..
nSubImg = 6;

%% ** load data
if ischar(data_in)
    filename = data_in;
    if contains(filename, '*') % if strcmp(filename, 'all files    ')
        tab_param = load_all_param;
    else
        filename_full = [dirname filesep filename '_tab_param.mat'] ;
        if ~isempty(dir(filename_full))
            load(filename_full, 'tab_param') %  tab_param = fread_params_timewindow(filename_full, 1, 1, t_last);
            fprintf('%s\r', filename)
        else
            tab_param = [];
        end
    end
else
    filename = '';
    tab_param = data_in;
end

if isempty(tab_param)
    par = [];
    varname = '';
    output_values = cell(length(output_variables),1);
    disp(['No data for ' filename ' (wrong dir???)'])
    return
end

%% ** prepare data
if ~isempty(filename)
    img1 = imread(filename, 1);
    middle = size(img1, 2)/2;
    tab_param = split_params_left_right(tab_param, middle, side);
end

tab_param(PARAM_I-1:N_PARAM:end, :) = tab_param(PARAM_I-1:N_PARAM:end, :)*pxl_size;
tab_param(PARAM_J-1:N_PARAM:end, :) = tab_param(PARAM_J-1:N_PARAM:end, :)*pxl_size;

alpha = tab_param(PARAM_ALPHA-1:N_PARAM:end, :); % intensty of fitted SM
TrcLen = sum(alpha>0)*time_lag; % number of images in each trace, AS 11/12/7

[linearity, full_path, direct_path] = calcul_linearity(tab_param);
inst_velocity = full_path./TrcLen;
resulting_velocity = direct_path./TrcLen;

asym = radius_gyration(tab_param, 0); % do_plot=0

% *********************************** compil data ****************************************
data = {full_path, direct_path, inst_velocity, resulting_velocity, linearity, asym};
varname = {'full path (um)', 'direct path (um)', 'instantaneous velocity (um/s)', 'resulting velocity (um/s)', 'linearity', 'asymetry'};
% ****************************************************************************************

amin = [1e-1 1e-3 0.1 1e-4 0 0];
amax = [1e2 1e2 10 1e1 0.1 1];

output_values = cell(length(output_variables),1);
n_output = 1;

dirtitle = [' dir=', num2str(dirname)];

plottitle = {cd, [filename side]};
plottitle2 = {date, dirtitle};

%% graphs
par = zeros(1, nSubImg);
n_graph = 1;
if new_fig, figure('WindowStyle', 'docked'); end

for v = 1:nSubImg
    subplot(nSubImg/2, 2, v);
    
    if v < 5 % any(strncmp(varname{v}, {'Trc' 'Tc ' 'Tf ' 'r^2' 'x^2' 'D b'}, 3)) % 'gam'
        do_log = 1;
    else
        do_log = 0;
    end
    %*************************************************************************
    par(v) = ct4plot_hist(data{v}, varname{v}, do_log, 1, clr, 1, text_offset); % => compute mean % ct4plot_hist(data, varname, do_log, Nfrm, clr, norm, text_offset) ct4plot_hist(data, '', 1, 1, [0 0 1], 0, 0)
    %*************************************************************************
    if v==1, title(plottitle, 'Interpreter', 'none');
    elseif v==2, title(plottitle2, 'Interpreter', 'none'); end
    
    a = axis; axis([amin(v) amax(v) a(3:4)])
    
    if find(output_variables==n_graph), output_values{n_output} = data{v}; n_output = n_output+1; end
    n_graph = n_graph+1;
end

%%%