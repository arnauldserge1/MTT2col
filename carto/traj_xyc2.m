function traj_xyc2(filename, codage, dirname, min_length_ratio)

% function traj_xyc2(filename, codage, dirname, min_length_ratio)
%
% plot all traces on first or trans. image
%
% if codage=='speed' % codage movement (==speed dr/dt)
% elseif codage=='int' % codage intensity (cluster...)
% elseif codage=='time' % codage time
%
% V1.0 AS 3/5/2006
% V1.1 oct 2006 ajout du codage conf
% V1.2 dec 2006 ajout codage var
% V2.0 2013!
% V3 2014
% see also cartobyf, MTT23i, rainbow_pict


global N_PARAM PARAM_T PARAM_I PARAM_J PARAM_ALPHA
if isempty(N_PARAM), MTTparams_def; end


if nargin<4, min_length_ratio = 0; end
if nargin<3, params_def = MTTparams_def; dirname = params_def{4}; end
if nargin<2, codage = 'speed'; end
if nargin<1, files = dir('*.tif'); filename = files(1).name; end
if isempty(filename), disp('No data... Check dir & filename !'), return, end


filename_full = [dirname filesep filename '_tab_param.mat'];
tab_param = importdata(filename_full);
if isempty(tab_param), disp('No data... Check dir & filename !'), return, end

tab_t = tab_param(PARAM_T-1:N_PARAM:end, :);
tab_i = tab_param(PARAM_I-1:N_PARAM:end, :);
tab_j = tab_param(PARAM_J-1:N_PARAM:end, :);
tab_alpha = tab_param(PARAM_ALPHA-1:N_PARAM:end, :);

if strcmp(codage, 'rainbow2')
    pict = ones(512)*3;
    figure('WindowStyle', 'docked')
    h1 = imagesc(pict); axis image ij off, hold on
    title(filename, 'interpreter', 'none')
    colormap('gray')
else
    [pict, h1] = DIC_image(filename);
    pict = double(pict);
end

%% prep data
Tmax = size(tab_param, 1)/N_PARAM;
Ntrc = size(tab_param, 2);

min_trc_length = round(Tmax*min_length_ratio);
Imax = 20000;%...max - 1%??
par_def = MTTparams_def; Dmax = str2double(par_def{8}); sig_free = 2*sqrt(Dmax); Boule_free = str2double(par_def{14});
max_dr = sig_free*Boule_free; % = 1,98 pxl (or 3.17 um/s) for sJB %%% r2 = calcul_r2(tab_param); sqrt(max(r2(:)));
m = 64; % 64-elements is each colormap

%% *** codage couleur: max val ***
if strcmp(codage, 'speed') % codage par déplcmt (vitesse)
    tab_r2 = calcul_r2(tab_param);
    tab_r2b = zeros([Tmax-2, Ntrc, 2]);
    tab_r2b(:, :, 1) = tab_r2(1: end-1, :);
    tab_r2b(:, :, 2) = tab_r2(2: end, :);
    tab_r2 = [tab_r2(1, :); nanmean(tab_r2b, 3); tab_r2(end, :)]; % 2 by 2 average to get T values, not T-1...

    val_max = max_dr; % 1.98 pxl/frm ~ 3.2 um/s
    cmap = hot(1.5*m);
    cmap = cmap(round(m/8)+1:round(m/8)+m, :);
elseif strcmp(codage, 'int') % codage par intensité (cluster...)
    val_max = Imax;
    cmap = copper(m);
elseif strcmp(codage, 'time')
    val_max = Tmax;
    cmap = jet(m);
elseif strcmp(codage, 'rainbow2')
    val_max = Tmax;
    cmap = hsv(m);
else
    val_max = 1;
    cmap = zeros(m, 3); cmap(:, 3) = 1;
end

%% double colormap for pict & traces
colormap([gray(m); cmap]) %, colorbar %cf nanopict
cmin = min(pict(:));
cmax = max(pict(:));
if strcmp(codage, 'rainbow2'), cmin = 1; cmax = 4; end
C1 = min(m, round((m-1) * (pict-cmin) / (cmax-cmin)) + 1);
C2 = m + C1;
set(h1, 'CData', C1);% set(h(2),'CData',C2);
if strcmp(codage, 'speed'), set(h1, 'CData', max(C1(:)) - C1); end
if strcmp(codage, 'rainbow2'), caxis([1 2*m]) % max(C2(:))])
else caxis([min(C1(:)) max(C2(:))])
end

%% --- go through traces ---
disp('traj :            ')

for nt = 1:Ntrc
    ok = (tab_alpha(:, nt) > 0); % # des points de la traj i
    
    if sum(ok) > min_trc_length
        %% for rainbow: add larger trace below, green if complet, red if fraction of recording time
        if strcmp(codage, 'rainbow2')
            if sum(ok) < Tmax, Markerclr = [1 1 1]*.3;% unstable, short = dark gray 'g';%m+m/3;, for hsv map % green??
            else Markerclr = [1 1 1]*.9;%'w'; % stable = white %m+1; % red??
            end
%             patch('XData', [nan; tab_j(ok, nt); nan], 'YData', [nan; tab_i(ok, nt); nan], 'LineWidth', 5, 'FaceColor', 'none', 'EdgeColor', clr);
        end
        
        %%          *** codage couleur ***
        if strcmp(codage, 'speed') % codage par déplcmt (vitesse)
            val = sqrt(tab_r2(:, nt));%  if val >= max_dr, val=0, end % do_log = 1 % ???
        elseif strcmp(codage, 'int') % codage par intensité (cluster...)
            val = tab_alpha(:, nt);
        elseif strcmp(codage, 'time')
            val = tab_t(:, nt);
        elseif strcmp(codage, 'rainbow2') % hue = (n-1)/(Nimg-1); % color: from 0 to 1 for full rainbow ...
            val = (val_max - tab_t(:, nt))*0.8; %... 0 to 0.8 to skip violet & reversed (0.8 to 0) to go from blue to red
%             linew = 2???
        else
            val = ones(Ntrc, 1);
        end
        
        val = min(val, val_max); % 'ecretage', for dpl during blink, intensity 'out of range'..
        val = m+1+val*m/val_max; % scaled to [m+1 2m] % m+1, AS 5/1/2015
        if strcmp(codage, 'rainbow2')
            %             patch('XData', [nan; tab_j(ok, nt); nan], 'YData', [nan; tab_i(ok, nt); nan], 'CData', [nan; val(ok); nan], 'FaceColor', 'interp', 'EdgeColor', 'flat', 'LineWidth', 3)
            patch('XData', [nan; tab_j(ok, nt); nan], 'YData', [nan; tab_i(ok, nt); nan], 'CData', [nan; val(ok); nan], ...
                'Marker', 'o', 'MarkerFaceColor', 'flat', 'MarkerEdgeColor', Markerclr, 'MarkerSize', 8, 'LineWidth', 1);
        else
            patch('XData', [nan; tab_j(ok, nt); nan], 'YData', [nan; tab_i(ok, nt); nan], 'CData', [nan; val(ok); nan], ...
                'FaceColor', 'interp', 'EdgeColor', 'interp')
        end
    end % if length(Ni)>0
    
    if mod(nt, 10)==0
        fprintf([repmat('\b', 1, 11) '%5i/%5i'], nt, Ntrc)
        drawnow expose
    end
end % for itrc = 1:NoTrace
fprintf([repmat('\b', 1, 11) '%5i/%5i\r'], Ntrc, Ntrc)

%%%