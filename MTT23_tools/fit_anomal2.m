function [D, gamma, offset] = fit_anomal2(msddata, offset, do_plot)

% function [D, gamma, offset] = fit_anomal2(msddata, offset, do_plot)
% fit MSD = 4Dt^gamma + offset => anomal coef
% offset = 2dr^2, with dr (= dx = dy) lateral accuracy 
% (error or noise on x and y)

if isempty(msddata)
    D = -1; gamma = -1; offset = -1; return
end

if nargin<3, do_plot = 0; end

MSD_FRACTION = 0.5;

if size(msddata,2)==1 % vals => err==1 (...)
    msddata = [ones(size(msddata,1),1), (1:size(msddata,1))', msddata, ones(size(msddata,1),1)];
elseif size(msddata,2)==2 % vals & erreurs
    msddata = [ones(size(msddata,1),1), (1:size(msddata,1))', msddata];
elseif size(msddata,2)==3 % n, vals & erreurs
    msddata = [ones(size(msddata,1),1), msddata];
end

ntrc = msddata(end,1);
D = zeros(1, ntrc);
gamma = zeros(1, ntrc);
if nargin<2, offset = nan(1, ntrc); end
if isempty(offset), offset = nan(1, ntrc); end % AS 2018
back = '';
if ntrc > 1
    fprintf('fitting MSD for anomalous diffusion, traj             ')
    back = repmat('\b', 1, 11);
end
if do_plot, figure('WindowStyle','docked'), pause(.1), back = '\b'; end

for i = 1:ntrc
    msdi = msddata(msddata(:,1)==i, 2:4);
    plot_length = round(size(msdi,1)*MSD_FRACTION);
    msdi = msdi(1:plot_length, :);
    msdi(msdi(:,3)==0, :) = []; % si err = 0..
    if size(msdi, 1) < 3, continue, end
    
    if isnan(offset(i))
        %% fit first 3 points to eval offset
        Results1 = weightedfit(msdi(1:3, :), 0);
        offset(i) = Results1.Intercept;
% % %         r2o = msdi(:,2);
% % %         offset(i) = r2o(1)-(r2o(2)-r2o(1));
    end
    
    %% fit log(r2) = gamma log(t) + log(4D) (with r2 = msd-offset in fact)
    tt = msdi(:,1);
    r2 = msdi(:,2) - offset(i);
    dr2 = msdi(:,3);
    dr2_log = log(1+dr2./r2); % msd or msd-offset???
    msd_log = [log([tt, r2]), dr2_log];
    msd_log = msd_log(r2 > 0, :);
    if size(msd_log, 1) < 2, continue, end % 2019
    if do_plot, clf, subplot(121), end
    
    %% weightedfit, log
    Results2 = weightedfit(msd_log, do_plot);
    D(i) = exp(Results2.Intercept)/4;
    gamma(i) = Results2.slope;
    
    if ntrc > 1, fprintf([back '%5i/%5i'], i, ntrc), end
    
    %% plot fit, log & lin
    if do_plot
        subplot(121)
        xlabel('t (log10)'), ylabel('<r²> (log10)')

        subplot(122)
        fit = 4*D(i)*tt.^gamma(i) + offset(i);
        % % %         err_pos = r2.*(1+1./weight); err_neg = max(r2.*(1-1./weight), 0); plot(tt, r2, '-', tt, err_pos, 'b:', tt, err_neg, 'b:', tt, fit, 'r:') % data +/- error & fit
        errorbar(tt, msdi(:,2), msdi(:,3), 'bs', 'MarkerFaceColor', 'b')
        hold on, plot(tt, fit, 'r.--', 'linewidth', 2)
        xlabel('t'), ylabel('<r²>')
        r0 = sqrt(offset(i)/2);
        title(['<r²> = 4Dt^{\gamma} + 2r_0², D = ' num2str(D(i),2)...
            ' \gamma = ' num2str(gamma(i),2) ' r_0 = ' num2str(r0,2)])
                
        drawnow
%         if ntrc > 1, pause(.1), end
%         if isnan(D(i)) || isinf(D(i))
%             pause
%         end
    end
end % for i=1:ntrc

if ntrc > 1, fprintf([back '%5i/%5i\r'], ntrc, ntrc), end
% D(isnan(gamma)) = 0;% gamma(isnan(gamma)) = 0;
% hist(gamma,100),hist(offset,100),hist(log(D(D>0)),100)
%%%