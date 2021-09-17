function plot_all_trajs_MU(tab_param, name, common_origin, size_max, new_fig, traj_color, align, final_position)

% function plot_all_trajs_MU(tab_param, name, common_origin, size_max, new_fig, traj_color, align, final_position)
%
% Plot all traces from data (either filename or trace matrix),
% using a common origin (=> all starting from 0,0)
% or first traces equally reparted over n_col by n_row (16 by 9)
% default: plot_all_trajs_MU(first_local_tif_file_param, '', 1, 0, 1, [], 1)
%
% AS 2013
%
% see also fit_directed_motion4

global N_PARAM PARAM_I PARAM_J PARAM_ALPHA
if isempty(N_PARAM), MTTparams_def; end

if nargin<8, final_position = 0; end
if nargin<7, align = 1; end
if nargin<6, traj_color = []; end
if nargin<5, new_fig = 1; end
if nargin<4, size_max = 0; end
if nargin<3, common_origin = 1; end
if nargin<2, name = ''; end
if nargin<1, files = dir('*.tif'); name = files(1).name; tab_param = importdata(['output23' filesep name '_tab_param.mat']); end

if ischar(tab_param), name = tab_param; tab_param = importdata(['output23' filesep name '_tab_param.mat']); end

n_col = 16;
n_row = 1;%9;

%% prepare traces: find first point t0 and set first position to 0,0~
y = tab_param(PARAM_I-1:N_PARAM:end,:); % xy ji!
x = tab_param(PARAM_J-1:N_PARAM:end,:);
alpha = tab_param(PARAM_ALPHA-1:N_PARAM:end,:);

[nfrm, ntrc] = size(x); % if ntrc==0, return, end

ncolor = size(traj_color,1);

t0 = zeros(1, ntrc);
for i = 1:ntrc
    xi = x(:,i);
    t0(i) = find((xi~=0), 1); % first point defined
    
    x(t0(i):end, i) = x(t0(i):end, i) - x(t0(i), i);
    y(t0(i):end, i) = y(t0(i):end, i) - y(t0(i), i);
end

x(alpha==0) = nan; % all_tab_param may contain nul values...
y(alpha==0) = nan;

if align
    [~, phi] = radius_gyration(tab_param, 0);
    phi = repmat(phi, nfrm, 1);
    x2 = x.*cos(phi) + y.*sin(phi);
    y = x.*sin(phi) - y.*cos(phi); % x = x2;
    for i = 1:ntrc, if nanmean(x2(:,i)) > 0, x(:,i) = x2(:,i); else x(:,i) = -x2(:,i); y(:,i) = -y(:,i); end, end
end

if size_max==0, size_max = max(abs([x(:); y(:)])); end
cmap = lines(ntrc);

if new_fig, figure('WindowStyle', 'docked'), end
if common_origin==0, axis off equal, else axis square, end
title(name, 'interpreter', 'none')
hold on

%%  ** plot traces **
for i = 1:ntrc
    xi = x(t0(i):end, i);
    yi = y(t0(i):end, i);
    xi = xi(~isnan(xi));
    yi = yi(~isnan(yi));
    
    imod = mod(i-1,ncolor)+1;
    
    if common_origin
        if ncolor==0, clr = cmap(i,:);
        elseif ncolor>1, clr = traj_color(imod,:); 
        else clr = traj_color;
        end
        
        plot(xi, yi, 'color', clr)
        if final_position
            plot(xi(end), yi(end), '.', 'markersize', 12, 'color', clr)
        end
        
    else % traces dispatched, as with subplot
        [ny, nx] = integer_division(i-1, n_col);
        xii = xi + nx*size_max;
        yii = yi - ny*size_max;
        
        if ncolor==0, clr = [0 0 0]; % if isempty(traj_color), clr = [0 0 0]; else clr = traj_color; end
        elseif ncolor>1, clr = traj_color(imod,:); 
        else clr = traj_color;
        end
        
        plot(xii, yii, 'color', clr) % plot(xii(1), yii(1), 'o', 'color', clr)
        if (i == n_row*n_col), return, end
        drawnow % pause
    end %    mean(x2(:,i)), pause
end % for i = 1:NoTrace

if common_origin, axis([-1 1 -1 1]*size_max), end % if align, a = axis; axis([0 2*a(2) a(3:4)]), end
%%%%%% scale_bar(calib_length, calib_text, pixelsize, W, H)

%%%