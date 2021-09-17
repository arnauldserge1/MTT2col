function plot_all_trajs(data)%, common_origin)

% function plot_all_trajs(data)
%
% Plot all traces from data (either filename or trace matrix),
% using a common origin (=> all starting from 0,0)
% and/or n_row by n_col (6 by 4) first traces,
% equally reparted in a 2nd figure

global N_PARAM PARAM_I PARAM_J % PARAM_ALPHA
if isempty(N_PARAM), MTTparams_def; end


if nargin<1
    data = dir('*.tif');
    if ~isempty(data), data = data(1).name;
    else disp('nada?'), return, end
end
% if nargin<2, common_origin = 0; end

n_row = 4;
n_col = 6;

%% -- Load traces --
if ischar(data)
    filename = data;
    xls_file = [filename(1:end-4) '_filt1.xls'];
    
    if ~isempty(dir(xls_file))
        trcdata = read_xls_cell_tracks(xls_file, {'xy_sub_pxl'});
        x = trcdata{1}(:,1:2:end);
        y = trcdata{1}(:,2:2:end); %        alpha = x; %...
    else
%         n_img_max = 1000
%         n_trc_max = 1000
        trcdata = fread_params_timewindow(filename);%,1,1,n_img_max,1,n_trc_max);
        y = trcdata(PARAM_I-1:N_PARAM:end,:); % xy ji!
        x = trcdata(PARAM_J-1:N_PARAM:end,:); %        alpha = trcdata(PARAM_ALPHA-1:N_PARAM:end,:);
    end
    % else %...simul, ou run avec tab_param au lieu de filename
    %     trcdata = data;
    %     y = trcdata(PARAM_I:N_PARAM:end,:); % xy ji!
    %     x = trcdata(PARAM_J:N_PARAM:end,:); %    alpha = trcdata(PARAM_ALPHA:N_PARAM:end,:);
end
clear trcdata

ntrc = size(x, 2);
if ntrc==0, return, end


%% prepare traces: find first point t0 and set first position to 0,0
t0 = zeros(1, ntrc);

for i = 1:ntrc
    t0(i) = find(x(:, i)>0, 1); % first point defined (not NaN nor 0)
    
    x(t0(i):end, i) = x(t0(i):end, i) - x(t0(i), i);
    y(t0(i):end, i) = y(t0(i):end, i) - y(t0(i), i);
end

size_max = 1.1*max(abs([x(:); y(:)])); % 1.1 for 10% margin (or 512 pxl by default??)
cmap = jet(ntrc);

for nfig = 1:2
    figure('WindowStyle', 'docked')
    title([cd, ' \ ', filename],'interpreter','none')
    axis equal %off
    hold on
    
    %%  ** plot traces **
    for i = 1:ntrc
        xi = x(t0(i):end, i);
        yi = y(t0(i):end, i);
        
        if nfig == 1 % common_origin
            plot(xi, yi, 'color', cmap(i,:))
        else % traces dispatched, as with subplot
            [ny, nx] = integer_division(i-1, n_col);
            plot(xi + nx*size_max, yi - ny*size_max)
            plot(xi(1) + nx*size_max, yi(1) - ny*size_max, 'ko')
            if (i == n_row*n_col), return, end
            drawnow % pause
        end
    end % for i = 1:NoTrace
    %%%%%% scale_bar(calib_length, calib_text, pixelsize, W, H)
end
%%%