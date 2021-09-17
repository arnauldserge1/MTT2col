function [D, gamma, offset] = fit_anomal2(msddata, offset, do_plot)

% function [D, gamma, offset] = fit_anomal2(msddata, offset, do_plot)
% fit MSD = 4Dt^gamma + offset => anomal coef
% offset = 2dr^2, with dr (= dx = dy) lateral accuracy 
% (error or noise on x and y)

if isempty(msddata)
    D = -1; gamma = -1; return
end

if nargin<3, do_plot = 0; end

MSD_FRACTION = 0.5 ;

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

if ntrc>1, fprintf('fitting MSD for anomalous diffusion, traj             '), end
if do_plot && ntrc>1, figure('WindowStyle','docked'), back = '\b'; else, back = repmat('\b', 1, 11); end

for i=1:ntrc
    msdi = msddata(msddata(:, 1)==i, 2:4);
    n = round(size(msdi,1)*MSD_FRACTION);
    msdi = msdi(1:n, :);
    msdi(msdi(:,3)==0, :) = []; % si err = 0..
    if size(msdi, 1)<3, continue, end
    
    %% fit log(r2) = gamma log(t) + log(4D) (with r2 = msd-offset in fact)
    tt = msdi(:,1);
    r2 = msdi(:,2);
    dr2 = msdi(:,3);
    if do_plot, clf, subplot(121), end
    
    %% weightedfit, log
    Results2 = weightedfit(msd, do_plot);
    D(i) = Results2.Intercept/4;
    gamma(i) = Results2.slope;
    
    %% plot fit, log & lin
    if do_plot
        subplot(121)
        xlabel('t'), ylabel('<r²> (log10)')

        subplot(122)
        fit = 4*D(i)*tt.^gamma(i) + offset(i);
        % % %         err_pos = r2.*(1+1./weight); err_neg = max(r2.*(1-1./weight), 0); plot(tt, r2, '-', tt, err_pos, 'b:', tt, err_neg, 'b:', tt, fit, 'r:') % data +/- error & fit
        errorbar(tt, msdi(:,2), msdi(:,3), 'bs', 'MarkerFaceColor', 'b');
        hold on; plot(tt, fit, 'r.--', 'linewidth', 2);
        xlabel('t'), ylabel('<r²>')
        r0 = sqrt(offset(i)/2);
        title(['<r²> = 4Dt^{\gamma} + 2r_0², D = ' num2str(D(i),2)...
            ' \gamma = ' num2str(gamma(i),2) ' r_0 = ' num2str(r0,2)])
        drawnow expose
        figure(gcf)
        if ntrc>1, pause, end
    end
    if ntrc>1, fprintf([back '%5i/%5i'], i, ntrc), end
end % for i=1:ntrc

if ntrc>1, fprintf('\r'), end
% D(isnan(gamma)) = 0;% gamma(isnan(gamma)) = 0;
% hist(gamma,100),hist(offset,100),hist(log(D(D>0)),100)
%%%