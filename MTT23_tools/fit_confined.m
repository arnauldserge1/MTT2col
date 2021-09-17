function [D, a, offset, R2, chi2] = fit_confined(msddata, do_plot)

% function [D, r, offset] = fit_confined(msddata, do_plot)
%
% fit pour MSD = 4/3a^2 - 128/pi^4 * a^2 * exp(-pi^2 * Dt/4a^2)+ offset 
% => Kusumi, 93
% avec offset = 2r0^2
%
% see also 
%   fit_anomal2, with linear fit on loglog: log(MSD-offset)=gamma*log(t) + log(4D),
%   fit_anomal: same code, but with error on D, gamma, offset 

if isempty(msddata)
    D = -1; a = -1; offset = -1; R2 = -1; return
end

if nargin < 2, do_plot = 0; end

MSD_FRACTION = 0.5 ;

if size(msddata, 2)==1 % vals => err==1 (...)
    msddata = [ones(size(msddata, 1), 1), (1:size(msddata, 1))', msddata, ones(size(msddata, 1), 1)];
elseif size(msddata, 2)==2 % vals & erreurs
    msddata = [ones(size(msddata, 1), 1), (1:size(msddata, 1))', msddata];
elseif size(msddata, 2)==3 % n, vals & erreurs
    msddata = [ones(size(msddata, 1), 1), msddata];
end

ntrc = msddata(end, 1);
D = zeros(1, ntrc);
a = zeros(1, ntrc);
offset = zeros(1, ntrc);
R2 = zeros(1,ntrc);

if ntrc > 1, fprintf('fitting MSD for confined diffusion, traj           '), end
if do_plot && ntrc > 1, figure('WindowStyle', 'docked'), pause(0.1), end
warning off Matlab:nearlySingularMatrix, warning off Matlab:SingularMatrix

optfit = optimset('Display', 'off');
lower_limits = [0 0 0]; 
upper_limits = [inf inf inf]; 
D0 = 1e-2/0.16^2*0.1; % pxl2/frm 
a0 = 0.1/0.16; %pxl

for i = 1:ntrc
    msdi = msddata(msddata(:, 1)==i, :);
    plot_length = round(size(msdi, 1)*MSD_FRACTION);
    if plot_length > 2
        msdi = msdi(1:plot_length, :);
        msdi(msdi(:, 4)==0, :) = []; % si err = 0..
        
        t = msdi(:, 2);
        r2 = msdi(:, 3);
        dr2 = msdi(:, 4);
        weight = 1./dr2.^2; % 1/error2
        if sum(weight) == 0, weight = ones(size(r2)); end
        
        plot_length = length(r2);
        if plot_length <= 2, continue, end
        
        f = @(p, t, weight) (((4/3)*p(2)^2 - 128/(pi^4)* p(2)^2 * exp((-pi^2*p(1)*t)/(4*p(2)^2))+p(3)).*weight); %msd modèle des mouvements confinés (carrés (2D))

        offset0 = r2(1)-(r2(2)-r2(1)); % linear interp from r2(1) & r2(2)
        p0 = [D0 a0 offset0]; % D, a, offset ini
        
        [p,chi2(i)] = lsqcurvefit(f, p0, t, r2.*weight, lower_limits, upper_limits, optfit, weight); % [0 -inf 0], [inf inf inf], optfit, weight);
        D(i) = p(1);
        a(i) = p(2);
        offset(i) = p(3);
        
        fit = (4/3)*p(2)^2 - 128/(pi^4) * p(2)^2 * exp((-pi^2*p(1)*t)/(4*p(2)^2)) +p(3);
        R2(i) = calcul_R2_fit(r2, fit);
        
        if do_plot
            r0 = sqrt(p(3)/2);
            clf, subplot(121)
            loglog(t, r2, '-', t, fit, 'r', t, (r2+dr2), 'b:', t, (r2-dr2), 'b:') 
            xlabel('log(t)'), ylabel('log(<r²>)') % ylabel('log(<r²>/t)')
            title('MSD and his fit')
            
            subplot(122)
            plot(t, r2, '.-', t, fit, 'r', t, r2+dr2, 'b:', t, r2-dr2, 'b:') % plot(tt, r2, '-', tt, fit, 'r:', tt, r2.*(1+1.*msdi(:, 4)), 'b:')
            xlabel('t'), ylabel('<r²>')
            title({['D = ' num2str(D(i), 2) ' pxl2/frm, a = ' num2str(a(i), 2) ' pxl, r_0 = ' num2str(r0, 2) ' pxl'];[' R2 = ' num2str(R2(i),2) ' pxl, chi2 = ' num2str(chi2(i),2) ' pxl']})
            drawnow
            if ntrc > 1, pause, end
        end
    end % if n > 2
    if ntrc > 1, fprintf([repmat('\b', 1, 11) '%5i/%5i'], i, ntrc), end
end % for i = 1:ntrc
if ntrc > 1, fprintf('\r'), end
