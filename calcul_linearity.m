function [linearity, full_path, direct_path] = calcul_linearity(tab_param, do_plot)%, lin_exponent)

% function [linearity, full_path, direct_path] = calcul_linearity(tab_param, do_plot)%, lin_exponent)
%
% linearity vs tortuosity:
% L = (direct_path / full_path) (^lin_exponent???)
%   = r(1_end) / sum(r(i_i+1))
%   = 1 for perfect line
% 0 < L < 1 for Brownian
% ~ 0 for confined

global N_PARAM PARAM_I PARAM_J PARAM_ALPHA
if isempty(N_PARAM), MTTparams_def; end

if nargin<2, do_plot = 0; end
% if nargin<3, lin_exponent = 1/2, end

alpha = tab_param(PARAM_ALPHA-1:N_PARAM:end,:);
tab_i = tab_param(PARAM_I-1:N_PARAM:end,:);
tab_j = tab_param(PARAM_J-1:N_PARAM:end,:);

Ntraj = size(tab_param, 2);
full_path = zeros(1, Ntraj);
direct_path = zeros(1, Ntraj);
linearity = zeros(1, Ntraj);

if do_plot, figure('WindowStyle', 'docked'), end

for itraj = 1:Ntraj
    ok = alpha(:, itraj) > 0;
    ii = tab_i(ok, itraj);
    jj = tab_j(ok, itraj);
    
    if ~isempty(ii)
        di = diff(ii);
        dj = diff(jj);
        r2 = di.^2 + dj.^2; %% blink => r2 = NaN, cf calcul_r2?????
        rr = sqrt(r2);
        
        full_path(itraj) = sum(rr);
        direct_path(itraj) = sqrt((ii(end)-ii(1))^2 + (jj(end)-jj(1))^2);
        linearity(itraj) = direct_path(itraj) / full_path(itraj);
        
        if do_plot
            clf
            plot(ii, jj, '-o')
            axis equal
            hold on
            plot([ii(1) ii(end)], [jj(1) jj(end)], '-ok', 'linewidth', 3)
            title(sprintf('full path=%.2g, direct path=%.2g, tortuosity=%.2g', full_path(itraj), direct_path(itraj), linearity(itraj)))
            pause
        end
    end
end

% hist(log10(tortuosity),2*sqrt(length(tortuosity)))
%%%