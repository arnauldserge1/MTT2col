function [D, gamma, offset] = fit_anomal_log(msddata, do_plot)

% function [D, gamma, offset] = fit_anomal_log(msddata, do_plot)
%
% fit pour MSD = 4Dt^gamma + offset => coef anomal (<1 en général)
% soit log(MSD-offset) = gamma*log(t) + 4D
% D = [n_trc D err_D]
% offset = 2r0^2
% gamma = [n_trc gamma err_gamma]
% offset [n_trc offset err_offset]

if isempty(msddata)
    D = [1 -1 -1]; gamma = [1 -1 -1]; offset = [1 -1 -1]; return
end

if nargin<2, do_plot = 0; end
% if do_plot, figure('WindowStyle','docked'), end

MSD_FRACTION = 0.5 ;

if size(msddata,2)==1 % vals => err==1 (...)
    msddata = [ones(size(msddata,1),1), (1:size(msddata,1))', msddata, ones(size(msddata,1),1)];
elseif size(msddata,2)==2 % vals & erreurs
    msddata = [ones(size(msddata,1),1), (1:size(msddata,1))', msddata];
elseif size(msddata,2)==3 % n, vals & erreurs
    msddata = [ones(size(msddata,1),1), msddata];
end

ntrc = msddata(end,1);
D = zeros(ntrc,3);
gamma = zeros(ntrc,3);
offset = zeros(ntrc,3);
if ntrc>1, fprintf('fitting MSD for anomalous diffusion, traj           '), end
warning off Matlab:nearlySingularMatrix, warning off Matlab:SingularMatrix

for i=1:ntrc
    msdi = msddata(msddata(:,1)==i,:);
    n = round(size(msdi,1)*MSD_FRACTION);
    if n>2
        msdi = msdi(1:n,:);
        msdi(msdi(:,4)==0) = []; % si err = 0..
        
% % %         if isnan(offset(i, 2))
            %% fit first 3 points to eval offset
            Results1 = weightedfit(msdi(1:3, :), 0);
            offset(i) = Results1.Intercept;
% % %         end
        
        tt = msdi(:,2);
        r2 = msdi(:,3);
        dr2 = msdi(:,4);
        weight = 1./dr2.^2; % 1/error2
        if sum(weight)==0, weight = ones(size(r2)); end
        
        n = length(r2);
        if n<=2, continue, end
        
        f = inline('(4*p(1)+x*p(2)).*weight','p','x','weight');
        optfit = optimset('Display','off');
        offset0 = r2(1)-(r2(2)-r2(1))
        p0 = [0.2 1 offset0]
        
        [p,chi2,~,~,~,~,jacobian] =...
            lsqcurvefit(f,p0,log(tt),log(r2.*weight),[0 -inf 0],[inf inf inf],optfit,weight);
        
        Cov = inv(jacobian'*jacobian); % Covariance
        err = sqrt(max(chi2*diag(Cov),0));
        D(i,:) = [i, p(1), err(1)];
        gamma(i,:) = [i, p(2), err(2)];
        offset(i,:) = [i, p(3), err(3)];
        
        if do_plot
            clf, subplot(121)
            fit = 4*p(1)*tt.^p(2)+p(3);
            r0 = sqrt(p(3)/2);
            plot(tt,r2,'-',tt,fit,'r:',tt,r2+dr2,'b:',tt,r2-dr2,'b:')% plot(tt,r2,'-',tt,fit,'r:',tt,r2.*(1+1.*msdi(:,4)),'b:')
            xlabel('t'), ylabel('<r²>')
            title(['<r²> = 4Dt^{\gamma} + 2r_0², D = ' num2str(D(i,2),2)...
                ' \gamma = ' num2str(gamma(i,2),2) ' r_0 = ' num2str(r0,2)])
            subplot(122)
            loglog(tt,r2./tt,'-',tt,fit./tt,'r:',tt,(r2+dr2)./tt,'b:',tt,(r2-dr2)./tt,'b:') % log(r²/t) = log r²-log t = log 4Dt^(g-1)
            xlabel('log(t)'), ylabel('log(<r²>/t)')
            figure(gcf)
            if ntrc > 1, pause, end
        end
    end % if n>2
    if ntrc>1, fprintf([repmat('\b',1,11) '%5i/%5i'],i,ntrc), end
end % for i=1:ntrc
if ntrc>1, fprintf('\r'), end

% g_pos = gamma(:,2)>0;
% D = D(g_pos,:);
% gamma = gamma(g_pos,:);
% offset = offset(g_pos,:);

if isempty(gamma) %??
    D = [1 -1 -1]; gamma = [1 -1 -1]; offset = [1 -1 -1];
end


% % % msd_simul = 1 + randn(100,1) + (1:100)' + (1:100)'.^2;