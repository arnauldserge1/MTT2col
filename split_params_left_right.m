function [tab_param1, tab_param2] = split_params_left_right(tab_param, middle, side, verbose)

% function [tab_param1, tab_param2] = split_params_left_right(tab_param, middle, side, verbose)
%
% side = 'left' or 'right' (or '' to get both)
% middle = 512 by default
% right: j = j-512, repositioned at [0 512]
% can return either 1 (if 'side' is specified) or 2 tables


if (nargin < 2), middle = 512, disp('Using 512 as default middle value'), end %%% size(im, 2)/2; % 512, a priori
if (nargin < 3), side = ''; end
if (nargin < 4), verbose = 1; end

if isempty(side) && (nargout == 1), tab_param1 = tab_param; return, end

global N_PARAM PARAM_J PARAM_ALPHA
if isempty(N_PARAM), MTTparams_def; end

% ii = tab_param(PARAM_I-1:N_PARAM:end, :); % i=y % ii(alpha == 0) = nan;
jj = tab_param(PARAM_J-1:N_PARAM:end, :); % j=x
alpha = tab_param(PARAM_ALPHA-1:N_PARAM:end, :);
jj(alpha == 0) = nan;

min_j = min(jj,[],1);
max_j = max(jj,[],1);

% mean_j = mean(jj(:));
% sd_j = std(jj(:));
% if (nargin < 2) && (mean_j-sd_j < 256)
%     disp('caution, using 512 as default middle value but data may be unadapted..')
% end

%% split traces
tab_param_left = tab_param(:, max_j < middle); % "green" == left == max(j) < 512
tab_param_right = tab_param(:, min_j > middle); % "red" == right == min(j) > 512
tab_j = tab_param_right(3:N_PARAM:end, :);
tab_j(tab_j > 0) = tab_j(tab_j > 0) - middle; % right: j = j-512, repositioned at [0 512] (0 left @ 0, for "not yet started traces")
tab_param_right(3:N_PARAM:end, :) = tab_j;

N_ini = size(tab_param, 2);
N_left = size(tab_param_left, 2);
N_right = size(tab_param_right, 2);
perc_left = 100*N_left/N_ini;
perc_right = 100*N_right/N_ini;
N_disc = N_ini - N_left - N_right;

if verbose
    fprintf('Splitting at %i: found %i traces below (left, %.1f%%), %i traces above (right, %.1f%%) and %i discarded (crossing middle line)\r', ...
        middle, N_left, perc_left, N_right, perc_right, N_disc)
end
if (nargout == 1)
    if strcmp(side, 'left'), tab_param1 = tab_param_left;
    elseif strcmp(side, 'right'), tab_param1 = tab_param_right;
    end
else
    tab_param1 = tab_param_left;
    tab_param2 = tab_param_right;
end