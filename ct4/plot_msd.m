function plot_msd(filename)

global N_PARAM PARAM_ALPHA
MTTparams_def;

%% ** load data
tab_param = fread_params_timewindow(filename);
tlag = 0.036;
pxlsz = 0.156/1.5;
min_trc_length = 2;

%% ** prepare data
alpha = tab_param(PARAM_ALPHA-1:N_PARAM:end, :); % intensty of fitted SM
TrcLen = sum(alpha>0)'; % number of images in each trace, AS 11/12/7
fprintf('discarding %i traces with less than %i step(s) \n', sum(TrcLen<min_trc_length), min_trc_length)
tab_param = tab_param(:, TrcLen>=min_trc_length); % one shot discarded, no interest!! 19/2/2010!!

%% --- data, MSD & D --- % calculate mean square displacement & diff
trc = detect_reconnex_to_trc(tab_param);
[~, ~, MSDcompil] = msd(trc);

%% graphs
figure('WindowStyle', 'docked');
% ct4plot_msd(MSDcompil); % => D
r2 = MSDcompil(MSDcompil(:,2)>0,:)*pxlsz^2;
tmax = size(r2,1);
nmax = round(tmax/2); % MSD_FRACTION
tt = (1:nmax)*tlag;

plot(tt,r2(1:nmax,2),'k') %tt,r2(1:nmax,2)+r2(1:nmax,3),'k:',tt,r2(1:nmax,2)-r2(1:nmax,3),'k:');

xlabel('\Deltat (s)'), ylabel('MSD (\mum^2)')
title(filename,'interpreter','none')
hold on

% calcul D
D = calculDinst(r2); % D = D(2);
s1 = sprintf('D_{fit} = %4.2g um^2/s',D(2)/tlag);
text(0.05,0.9,s1,'Units','normalized')
if D(2)>0, plot([1 5]*tlag, 4*D(2)*[1 5]+D(4), 'r:', 'linewidth',2), end
%%%