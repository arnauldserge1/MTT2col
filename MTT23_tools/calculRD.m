function RD = calculRD(msddata,do_plot,discard_empty_trc)

% function RD = calculRD(msddata,do_plot)
%
% calcule RD (relative deviation from Brownian diff) pour chaque traj:
% rapport MSD(n)/(4Dn) au point t = n tau
% D = D_1_5 (diff_window = 5)
% et n = 10
% msddata = r2 ou [n t r2 dr2], 1 ou 4 colonnes
% cf. Kusumi 93
%
% AS 22/2/5 
% v1.2 12/4/2013
% v1.3 24/9/2013 correct: msdi = msddata(indi,2:4), not msdi = msddata(indi,:)

if nargin<3, discard_empty_trc = 1; end
if nargin<2, do_plot = 0; end

if isempty(msddata), RD = -1; return, end
if size(msddata,2)==1
    msddata = [ones(size(msddata,1),1), (1:size(msddata,1))', msddata, ones(size(msddata,1),1)];
end

n = 10; % time point used to compute RD
RD = -ones(1,msddata(end,1));

for i=1:msddata(end,1) % #trc
    indi = find(msddata(:,1)==i);
    if length(indi)>=n
        msdi = msddata(indi,2:4);
        D = calculDinst(msdi); % en pxl/lag!
        RD(i) = msdi(n,2)/(4*D(2)*n+D(4)); % (4Dn+sig_b)
        
        if do_plot
            N = min(12,length(indi));
            tt = 1:N;%msdi(:,2);
            r2 = msdi(1:N,2);
            dr2 = msdi(1:N,3);
            fit = 4*D(2)*tt+D(4);
            plot(tt,r2,'o-',tt,fit,'r:',tt,r2+dr2,'b^',[n n],[0 msdi(n,2)],'b:')
            xlabel('t'),ylabel('<r^2>')
            title(['<r^2> = 4Dt + r_0^2, D = ' num2str(D(2),2)...
                ' r_0 = ' num2str(sqrt(D(4)),2) ' RD = ' num2str(RD(i),2)])
             figure(gcf), pause(1)
        end
    end
end

if discard_empty_trc, RD(RD==-1) = []; end % trc too short