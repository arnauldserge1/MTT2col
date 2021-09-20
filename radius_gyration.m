function [asym, phi] = radius_gyration(tab_param, do_plot)

% function [asym, phi] = radius_gyration(tab_param, do_plot)
%
% asym = (R12 - R22)^2 / (R12 + R22)^2, from small & large axis of a fitted ellipse
% phi = atan(2*Txy / (Txx - Tyy)) / 2; if Txx < Tyy, phi = phi + pi/2;
% asym = 0 for 'perfectly symetric' Brownian (fitted by circle), and 1 for perfect line
%
% cf. Saxton 1993

if nargin<2, do_plot = 0; end

y_tab = tab_param(2:8:end, :);
x_tab = tab_param(3:8:end, :);
alpha_tab = tab_param(4:8:end, :);

Ntrc = size(x_tab, 2);
asym = zeros(1, Ntrc);
phi = zeros(1, Ntrc);

for ntrc = 1:Ntrc
    ind = alpha_tab(:,ntrc)>0;
    x = x_tab(ind,ntrc);
    y = y_tab(ind,ntrc);
    
    if ~isempty(x)
        % Tensor
        Txx = mean(x.^2) - mean(x)^2;
        Txy = mean(x.*y) - mean(x)*mean(y);
        Tyy = mean(y.^2) - mean(y)^2;
        
        % Principal raddii of gyration = Eigen values
        R12 = (Txx + Tyy - sqrt((Txx - Tyy)^2 + 4*Txy^2)) / 2;
        R22 = (Txx + Tyy + sqrt((Txx - Tyy)^2 + 4*Txy^2)) / 2;
        
        % Radius of gyration
        % Rg2 = R12 + R22;
        
        % asym = R22/R12
        asym(ntrc) = (R12 - R22)^2 / (R12 + R22)^2;
        
        % Main angle
        phi(ntrc) = atan(2*Txy / (Txx - Tyy)) / 2;
        if Txx < Tyy, phi(ntrc) = phi(ntrc) + pi/2; end
% %         x2 = x*cos(phi(ntrc)) + y*sin(phi(ntrc));
% %         if nanmean(x2) < 0, phi(ntrc) = phi(ntrc) + pi/2; end
        
        if do_plot
            clf
            plot(x, y, 'o-')
            hold on
            t = 0:0.01:2*pi;
            Xe = sqrt(R22) * cos(t);
            Ye = sqrt(R12) * sin(t);
            xe = Xe*cos(phi(ntrc)) - Ye*sin(phi(ntrc)) + mean(x);
            ye = Ye*cos(phi(ntrc)) + Xe*sin(phi(ntrc)) + mean(y);
            plot(xe, ye, 'r');
            plot([xe(1) xe(315)], [ye(1) ye(315)],'r:')
            title(sprintf('trc %i, asymetry = %g, major angle = %g°', ntrc, asym(ntrc), phi(ntrc)*180/pi))
            axis equal
            if ntrc>1, pause(), end
        end
    end
end

%%%

% TEST:
%
% tab_par_simul = simul_tab_param_directed(100, 100, 1000, 1, 1, 1);
% asym = radius_gyration(tab_par_simul);
% hist(asym, 20)
%
% tab_par_simul = simul_tab_param_directed(100, 100, 1000, 1, 0, 1);
% asym = radius_gyration(tab_par_simul);
% hist(asym, 20)