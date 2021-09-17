function CumulSquareDisp(files)
% CSD: cf Nature de Mellman

if nargin<1, files = '*.stk'; end
files = dir(files);
Nfile = length(files);
figure('WindowStyle','docked')

for ifile = 1:Nfile
    filename = files(ifile).name;
    tab_param = fread_all_params(filename);
    r2 = calcul_r2(tab_param); % MTT_param, x = (tab_param(PARAM_J-1:N_PARAM:end,:));
    Ntrc = size(tab_param,2);
    
    for jtrc = 1:Ntrc
        r2j = r2(:,jtrc);
        first = find(~isnan(r2j),1, 'first');
        last = find(~isnan(r2j), 1, 'last');
        if last>first+10
            clf
            plot(r2j(first:last),'o'), hold on
            r2j(isnan(r2j)) = 0;
            CSD = cumsum(r2j(first:last));%/length(r2j);
            CSD = CSD*max(r2j)/CSD(end); % norm
            plot(CSD)
            xlabel('time (frames)'), ylabel('SD & CSD (pxl^2 & a.u.)') % pxl^2/frame
            title('Square Disp. & Normalized Cumul. Sq. Disp.')
            fprintf('\b\b\b\b\b\b\b\b\b\b\b %4i/%4i ',jtrc,Ntrc)
            figure(gcf), pause(1)
        end
    end
end