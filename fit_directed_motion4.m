function [nb, r0, D, v, chi2, v_threshold, D_threshold] ...
    = fit_directed_motion4(data_in, do_plot, amax, pxl_size, time_lag, side, use_MSDinput, align)

% function [nb, r0, D, v, chi2, v_threshold, D_threshold] ...
%     = fit_directed_motion4(data_in, do_plot, amax, pxl_size, time_lag, side, use_MSDinput, align)
%
% default: fit_directed_motion4(dir('*.tif')(1).name, 1, 40, 0.16, 0.1, 0)
%
% fit msd(t) = 2r02 + 4Dt + v2t2 and sort traces according to D & v values,
% quasi null or relevant => immobile, directed, Brownian or fast
% mean value above D & v hist is for values > threshold
% Caution: input in pxl, outputs in um!
%
% see also fit_directed_by_file, plot_all_trajs_MU

MSD_FRACTION = 0.5;

if nargin<8, align = 0; end
if nargin<7, use_MSDinput = 0; end
if nargin<6, side = ''; end
if nargin<5, [pxl_size, time_lag] = get_calib3; end
if nargin<3, amax = 40; end % max for trc axes, in pxl => amax*pxl_size in um
if nargin<2, do_plot = 1; end
if nargin<1, files = dir2('*.tif'); data_in = files(1).name; end

if ischar(data_in)
    filename = data_in;
    filename_full = ['output23' filesep filename '_tab_param.mat'] ;
    if ~isempty(dir(filename_full))
        S = load(filename_full); tab_param = S.tab_param; % if do_test, tab_param = tab_param(:, 1:10); end
        fprintf('%s - %s', filename, side)
    else
        tab_param = [];
    end
elseif use_MSDinput
    msddata = data_in; % caution, tab_param required for plot
else
    tab_param = data_in; filename = '';
end

if ~isempty(side)
    if isempty(filename)
        disp('could not split left/right wo filename')
    else
        img1 = imread(filename, 1);
        middle = size(img1, 2)/2;
        tab_param = split_params_left_right(tab_param, middle, side);
    end
end

do_print = do_plot;
if ~use_MSDinput
    msddata = msd(detect_reconnex_to_trc(tab_param), 1, do_print); % = [n t r2 sd_r2] ntraj, nimage, msd, error
end

if isempty(msddata), ntrc = 0;
else, ntrc = msddata(end, 1);
end

r0 = zeros(ntrc, 1);
D = zeros(ntrc, 1);
v = zeros(ntrc, 1);
chi2 = zeros(ntrc, 1);
nb = zeros(1, 5);

v_threshold = sqrt(3e-7)*pxl_size/time_lag; % 8.8e-4 um/s = 0.88 nm/s
D_threshold = 3e-6*pxl_size^2/time_lag; % 7.7 e-7 um2/s = 0.77 nm2/s
min_trc_length = 5;

if isempty(msddata), return, end

%% %%%%%%%%%%%%%%%%%%%%%% equation & fit opt %%%%%%%%%%%%%%%%%%%%%%
% fit_directed = inline('(p(1) + 4*p(2)*x + p(3)*x.^2).*weight', 'p', 'x', 'weight'); % = (2r02 + 4Dt + v2t2).*weight => yfit == fit_directed./weight
fit_directed = @(p, x, weight) ((p(1) + 4*p(2)*x + p(3)*x.^2).*weight); % = (2r02 + 4Dt + v2t2).*weight => yfit == fit_directed./weight
optfit = optimset('Display', 'off');
warning off Matlab:nearlySingularMatrix, warning off Matlab:SingularMatrix
warning off MATLAB:Axes:NegativeDataInLogAxis

%% for axe limits
Dmin = 1e-15*pxl_size^2/time_lag; % 2.6e-16
Dmax = 10*pxl_size^2/time_lag; % 2.6
vmin = sqrt(1e-15)*pxl_size/time_lag; % 5e-8
vmax = 10*pxl_size/time_lag; % 16

%% loop/ trajs
if (ntrc > 1), fprintf('fitting MSD for diffusion & directed motion, traj           '), end
if do_plot && (ntrc > 1), figure('WindowStyle','docked'), end

for i = 1:ntrc
    msdi = msddata(msddata(:, 1)==i, :);
    n = round(size(msdi, 1)*MSD_FRACTION);
    if (n > min_trc_length*MSD_FRACTION)
        msdi = msdi(1:n, :);
        
        tt = msdi(:, 2);
        r2 = msdi(:, 3);
        weight = 1./msdi(:, 4); % 1/error
        if (sum(weight) == 0), weight = ones(size(r2)); end
        tt(weight == inf) = []; r2(weight == inf) = []; weight(weight == inf) = []; % si err = 0..
        
        n = length(r2);
        if (n <= 2), continue, end
        
        r0ini = r2(1) - (r2(2)-r2(1));% p0 = [r0ini 3e-6 sqrt(3e-7)]; %
        p0 = [r0ini .1 .1]; % p0 = [r0ini Dini vini]; dr, Din, dr/tau???
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p = lsqcurvefit(fit_directed, p0, tt, r2.*weight, [0 0 0], [inf inf inf], optfit, weight);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        r0(i) = sqrt(p(1)/2)*pxl_size; % [i, sqrt(p(1)/2), err(1)];
        D(i) = p(2)*pxl_size^2/time_lag;
        v(i) = sqrt(p(3))*pxl_size/time_lag; % v(i, :) = [i, sqrt(p(3)), err(3)];
        chi2(i) = sum((r2 - fit_directed(p, tt, weight)./weight).^2);
        
        %% plot MSD
        if do_plot
            if (ntrc > 1)
            if (v(i) < v_threshold && D(i) < D_threshold), subplot(223)
            elseif (v(i) < v_threshold && D(i) > D_threshold), subplot(224)
            elseif (v(i) > v_threshold && D(i) > D_threshold), subplot(222)
            else, subplot(221)
            end
            end
            plot(tt, r2), hold on, plot(tt, fit_directed(p, tt, weight)./weight, 'r:')
            title(sprintf('D=%.2g v=%.2g, r0 = %.2g', D(i), v(i), r0(i)))
            figure(gcf), pause(.01)
        end
    end % if n > 2
    if (ntrc > 1), fprintf([repmat('\b', 1, 11) '%5i/%5i'], i, ntrc), end
end % for i = 1:ntrc
if (ntrc > 1), fprintf('\r'), end

nb(5) = size(v, 1); % d f B i

if do_plot && ~isempty(v) && (ntrc > 1)
    clf%figure('WindowStyle', 'docked') 24/10/2018
    
    %% all together
    subplot(2, 4, 1)
    title_str = {'All traces', ['n = ' num2str(nb(5))]};
    
    %*************************************************************************************
    plot_all_trajs_MU(tab_param*pxl_size, title_str, 1, amax*pxl_size, 0, [0 0 0], align);
    %*************************************************************************************
    
    %% histo D
    text_offset = -1; % write mean (only above threshold) (-1: dont write on plot with ct4plot_hist)
    subplot(2, 4, 2)
    ct4plot_hist(D, 'D (µm^2/s)', 1, 1, 'k', 0, text_offset);
    a = axis; axis([Dmin Dmax a(3:4)]), axis square
    set(gca, 'XTick', [1e-12 1e-6 1])
    
    Dth = D(D >= D_threshold);
    str = sprintf('D = %.2g*/%.2g µm^2/s', exp(mean(log(Dth))), exp(std(log(Dth))));
    title(str) % , 'interpreter', 'none''HorizontalAlignment', 'right',
    hold on, plot([D_threshold D_threshold], a(3:4), 'k:')
    
    %% v
    subplot(2, 4, 5)
    ct4plot_hist(v, 'v (µm/s)', 1, 1, 'k', 0, text_offset);
    a = axis; axis([vmin vmax a(3:4)]), axis square
    set(gca, 'XTick', [1e-6 1e-3 1])
    
    vth = v(v >= v_threshold);
    str = sprintf('v = %.2g*/%.2g µm/s', exp(mean(log(vth))), exp(std(log(vth))));
    title(str) % , 'interpreter', 'none''HorizontalAlignment', 'right',
    hold on, plot([v_threshold v_threshold], a(3:4), 'k:')
    view(-90, 90) % switch x, y!
    
    name = {'directed' 'fast' 'Brownian' 'immobile'};
    if strcmp(side, 'left'), color = [0 1 0; 0 0.75 0; 0 0.5 0; 0 0.25 0]; % Green (dark > light)
    elseif strcmp(side, 'right'), color = [1 0 0; 0.75 0 0; 0.5 0 0; 0.25 0 1]; % Red (dark > light)
    else, color = [0 1 1; 0 0 1; 0.5 0 1; 0 0 0]; % Cyan Blue Purple Black% % % color = [0 0 0; 0.5 0 1; 0 0 1; 0 1 1]; % Black Purple Blue Cyan  % % color = [0 0 1; 0 1 0; 1 0 0; 1 .5 0]; % B G R O
    end
end

%% sort traces
testv = (v > v_threshold);
testD = (D > D_threshold);

for i = 1:4
    switch i
        case 1, test = testv & ~testD; % lin/directed
        case 2, test = testv & testD; % fast
        case 3, test = ~testv & testD; % Brownian
        case 4, test = ~testv & ~testD; % conf/immobile % put conf last, 26/11/2015!!
    end
    nb(i) = sum(test); % size(tab_param_i, 2);
    
    if do_plot && ~isempty(v) && (ntrc > 1)
        tab_param_i = tab_param(:, test);
        
        %% correl: scatter v vs. D
        subplot(2, 4, 6)% % % %         subplot(1, 2, 1)
        if (ntrc > 1000), symbol = '.'; else, symbol = 'o'; end
        plot(D(test), v(test), symbol, 'color', color(i, :))
        hold on
        
        if (i == 4)
            plot([D_threshold D_threshold], [vmin vmax], 'k:')
            plot([Dmin Dmax], [v_threshold v_threshold], 'k:')
            axis([Dmin Dmax vmin vmax]), axis square
            set(gca, 'XTick', [1e-12 1e-6 1], 'YTick', [1e-6 1e-3 1], 'XScale', 'log', 'YScale', 'log')
            xlabel('D (�m^2/s)'), ylabel('v (�m/s)')
            title({filename, side}, 'interpreter', 'none')
        end
        
        %% plot traces
        if (i <= 2), subplot(2, 4, i+2), else, subplot(2, 4, 11-i), end
        title_str = [name{i} ' ' num2str(100*nb(i)/nb(5), 3) '%'];%         title_str = {name{i}, ['n = ' num2str(nb(i)) ' (' num2str(100*nb(i)/nb(5), 3) '%)']};
        
        %**************************************************************************
        plot_all_trajs_MU(tab_param_i*pxl_size, title_str, 1, amax*pxl_size, 0, color(i, :), align); %%% cmapi) % plot_all_trajs_MU(tab_param, name, common_origin, size_max, new_fig, traj_color)
        %**************************************************************************
    end % if do_plot && ~isempty(v)
end
%%%