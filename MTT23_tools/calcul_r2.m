function [r2, r2mean] = calcul_r2(tab_param, remove_blink)

% function [r2, r2mean] = calcul_r2(tab_param)
% r2 = di.^2+dj.^2;
% r2mean = nanmean(r2);
% distance between each step of each traj.


global N_PARAM PARAM_I PARAM_J PARAM_ALPHA
if isempty(N_PARAM), MTTparams_def; end


if nargin < 2, remove_blink = 1; end


alpha = tab_param(PARAM_ALPHA-1:N_PARAM:end,:);
di = diff(tab_param(PARAM_I-1:N_PARAM:end,:));
dj = diff(tab_param(PARAM_J-1:N_PARAM:end,:));

r2 = di.^2 + dj.^2;

if remove_blink
    % blink => r2 = NaN
    r2(alpha(1:end-1,:)==0) = NaN;
    r2(alpha(2:end,:)==0) = NaN;
end

% mean (with blink)
r2mean = nanmean(r2);
if isnan(r2mean)
    r2mean = 0;
end