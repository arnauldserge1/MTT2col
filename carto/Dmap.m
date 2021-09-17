function Dmap(tab_param) % [imean jmean Dmean] = 
% map coef D = <r2>/4 mean of 1 step dpl, over 5 points

Nav = 5;
Npar = 7;

r2min = 0.1^2; % 6nm/160nm = 0.0375???

% tab_param = fread_all_params(file);
tab_i = tab_param(2:Npar:end,:);
tab_j = tab_param(3:Npar:end,:);
tab_r2 = calcul_r2(tab_param);

Nframes = size(tab_i,1);
Ntrc = size(tab_i,2);

imean = zeros(floor(Nframes/Nav),Ntrc);
jmean = zeros(floor(Nframes/Nav),Ntrc);
Dmean = zeros(floor(Nframes/Nav),Ntrc);

% loop over traces
for itrc = 1:Ntrc
    r2i = tab_r2(:,itrc); % r2 values for current trace
    first_step = find(~isnan(r2i),1); % r2 = nan when step is not defined
    next_step = first_step;
    k = 1;
    while next_step+Nav<Nframes
        if ~isempty(next_step)
            ind = (next_step:next_step+Nav-1);
            indok = ind(r2i(ind)>r2min);
            if ~any(isnan(r2i(ind))) && ~isempty(indok) % no blink, but below r2min tolerated?????
                imean(k,itrc) = mean(tab_i(indok,itrc));
                jmean(k,itrc) = mean(tab_j(indok,itrc));
                Dmean(k,itrc) = mean(r2i(indok))/4;
                k = k+1;
            end
        end
        next_step = next_step+Nav-1+find(~isnan(r2i(next_step+Nav:end)),1);
    end
end

imean = imean(imean>0);
jmean = jmean(jmean>0);
Dmean = Dmean(Dmean>0);

Dvalues = [0.01 0.03 0.1];

figure('windowstyle','docked')
plot(imean(Dmean<=1e-2),jmean(Dmean<=1e-2),'k.')
% pause
hold on
plot(imean(Dmean>Dvalues(1)),jmean(Dmean>Dvalues(1)),'r.')
plot(imean(Dmean>Dvalues(2)),jmean(Dmean>Dvalues(2)),'g.')
plot(imean(Dmean>Dvalues(3)),jmean(Dmean>Dvalues(3)),'b.')
axis ij image

% S = contour_map(imean, jmean, Dmean);%filename)
% contourf(S)%,levels_val)
% axis ij off image % colormap colorbar, xy?