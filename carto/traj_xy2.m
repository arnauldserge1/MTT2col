function traj_xy2(filename, codage, dirname, min_length_ratio)

% function traj_xy2(filename, codage, dirname, min_length_ratio)
%
% plot all traces on first or trans. image
%
% if codage=='dens' % codage density, n/R2, cf. Douglass
% elseif codage=='speed' % codage movement (==speed dr/dt)
% elseif codage=='int' % codage intensity (cluster...)
% elseif codage=='time' % codage time
% cf. cartobyf

% V1.0 AS 3/5/2006
% V1.1 oct 2006 ajout du codage conf
% V1.2 dec 2006 ajout codage var
% V2.0 2013!

if nargin<4, min_length_ratio = 0; end
if nargin<3, params_def = MTTparams_def; dirname = params_def{4}; end
if nargin<2, codage = 'speed'; end
if nargin<1, files = dir('*.tif'); filename = files(1).name; end
if isempty(filename), disp('No data... Check dir & filename !'), return, end

surf = (150^2)/(160^2); % surf élem, (150 nm)^2, cf. Douglass, Cell 2005, avec 1pxl = 160nm
n_max = 24; % densité "max", en pk/pxl/img

figure('WindowStyle','docked')

Ncol = 36;
cmap = colormap(hot(Ncol));
% cmap = colormap(hsv(3*Ncol)); cmap = cmap(1:Ncol, :)
Imax = 2000*Ncol;
par_def = MTTparams_def; Dmax = str2double(par_def{8}); sig_free = 2*sqrt(Dmax); Boule_free = str2double(par_def{14});
max_dr = sig_free*Boule_free; % = 1,98 pxl for sJB %%% r2 = calcul_r2(tab_param); sqrt(max(r2(:)));

%% boucle / file
filename_full = [dirname filesep filename '_tab_param.mat'] ;
tab_param = importdata(filename_full);

if strcmp(codage, 'int'), [trcdata, pkdata] = detect_reconnex_to_trc(tab_param);
else trcdata = detect_reconnex_to_trc(tab_param);
end
if isempty(trcdata), return, end

DIC_name = dicname(filename);
if ~isempty(dir(DIC_name))
    DIC = imread(DIC_name);
    sat = .002 ; % saturation 0.2% min-max du contraste
    DIC_sat = imadjust(DIC, stretchlim(DIC, [sat 1-sat]), [0 1]);
    H = fspecial('average');
    pict = imfilter(DIC_sat,H,'replicate');
else
    pict = imread(filename,1);
    pict = max(pict(:)) - pict; % invert
end

%% prep data
Tmax = max(trcdata(:,2));
ntrc = trcdata(end,1);
min_trc_length = round(Tmax*min_length_ratio);

%% *** codage couleur: max val ***
if strcmp(codage, 'dens') % codage par densité
    val_max = n_max;
elseif strcmp(codage, 'speed') % codage par déplcmt (vitesse)
    val_max = max_dr; % LOG ???
elseif strcmp(codage, 'int') % codage par intensité (cluster...)
    val_max = Imax;
elseif strcmp(codage, 'time')
    val_max = Tmax;
else val_max = inf;
end
% % % % % a=linspace(0,val_max,Ncol);imagesc(a), axis equal off,colormap(hot),colorbar,colorbar('title','speed (pxl/frm)')


%% *** met l'image de la cellule "au plancher" ***
pict = double(pict);
imagesc(pict)
colormap('gray')
axis ij image off
set(gca, 'ticklength', [0 0])
title({cd; filename ; [' codage: ' codage]}, 'interpreter', 'none')
hold on
pause(.1)

%% --- go through traces ---
disp('traj :            ')

for itrc = 1:ntrc
    
    if mod(itrc,10)==0
        fprintf([repmat('\b',1,11) '%5i/%5i'],itrc,ntrc)
        drawnow expose
    end
    
    Ni = find(trcdata(:,1)==itrc); % # des points de la traj i
    
    if length(Ni)>min_trc_length
        trci = trcdata(trcdata(:,1)==itrc,:);
        n_trci = trci(:,2)'; %  # d'image, la 2e colonne, de la trace i
        dtrci = diff(n_trci); % pour détecter les blinks
        
%%  *** intensités traji ***
        if strcmp(codage, 'int')
            pki = pkdata(trcdata(:,1)==itrc,:);
            inti = pki(:,5); % intensité des points de la traj i
        end
        
        for istep = 1:length(n_trci)-1
            
%%          coordonnées t,x,y du pas en cours
            ti = trcdata(Ni(istep),2);
            xi = trcdata(Ni(istep),3);
            xii = trcdata(Ni(istep+1),3);
            yi = trcdata(Ni(istep),4);
            yii = trcdata(Ni(istep+1),4);
            
%%          *** codage couleur ***
            if strcmp(codage, 'dens') % codage par densité
                disti = sqrt((trcdata(Ni,3)-xi).^2 + (trcdata(Ni,4)-yi).^2);
                val = sum(disti<=sqrt(surf)); % densité = nb de points ds voisinage
                
            elseif strcmp(codage, 'speed') % codage par déplcmt (vitesse)
                val = sqrt((xii-xi)^2+(yii-yi)^2); % LOG ???
                
            elseif strcmp(codage, 'int') % codage par intensité (cluster...)
                val = inti(istep);
                
            elseif strcmp(codage, 'time')
                val = ti;
                
            else val = 0; % AS 2014
            end
            
            if val>val_max, val = val_max; end % 'ecretage', for dpl during blink, intensity 'out of range'..
            val2 = round(val*(Ncol-1)/val_max)+1; % val 0 to 1, val2 1 to Ncol (2 to Ncol-1, avoid limits of cmap??)
            
%%          ** plot du segment istep **
            if dtrci(istep)==1 % pas blink
                line([xi xii], [yi yii], 'color', cmap(val2,:)) % hsv2rgb([val,1,1]))
            else % blink => LineStyle :
                line([xi xii], [yi yii], 'color', cmap(val2,:), 'LineStyle', ':')
            end
        end % for istep = 1:length(Ni)-1
    end % if length(Ni)>0
end % for itrc = 1:NoTrace

axis image % colormap(hot), colorbar cf nanopict
fprintf([repmat('\b',1,11) '%5i/%5i\r'],ntrc,ntrc)

%%%

% % % for n=1:100
% % %     in = i(:,n); jn = j(:,n); tn = t(:,n);
% % %     patch('XData', [nan;in(tn>0);nan], 'YData', [nan;jn(tn>0);nan], 'CData', [nan;tn(tn>0);nan], 'FaceColor', 'interp', 'EdgeColor', 'interp');
% % % end
