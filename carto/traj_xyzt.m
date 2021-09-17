function ok = traj_xyzt(file, color_code, min_trc_length, dirname)

% function ok = traj_xyzt(file, color_code, min_trc_length, dirname)
% traces en 3D sur image dic (ou 1e image)
% code couleur pour z, r2 ou t (ou rien)


global N_PARAM PARAM_T PARAM_I PARAM_J PARAM_K PARAM_ALPHA gamma_z zr Boule_free
if isempty(N_PARAM), MTTparams_def; end

% if nargin<1, files = dir('*.stk'); else files = dir(filename); end
% if isempty(files), disp('Hello? No data... Check dir & filename !'), return, end
if nargin<2, color_code = 'z'; end
if nargin<3, min_trc_length = 3; end
pdef = MTTparams_def; 
if nargin<4, dirname = pdef{4}; end;

pxl_size = str2double(pdef{21})/1000; % 1 pxl = 0.16 µm

% %% boucle / file % for i = 1:length(files) file = files(i).name;

ok = 0;
disp('loading data...')
filename_full = [dirname filesep file '_tab_param.dat'] ;
tab_param = fread_all_params(filename_full);
if isempty(tab_param), return, end

%% some values
ntrc = size(tab_param,2);
Tmax = size(tab_param,1)/N_PARAM; % min_trc_length = round(Tmax*min_length_ratio);

t_tab = tab_param(PARAM_T-1:N_PARAM:end,:);
alpha_tab = tab_param(PARAM_ALPHA-1:N_PARAM:end,:);
xu_tab = tab_param(PARAM_J-1:N_PARAM:end,:)*pxl_size; % x==j..
yu_tab = tab_param(PARAM_I-1:N_PARAM:end,:)*pxl_size;
zu_tab = tab_param(PARAM_K-1:N_PARAM:end,:)/1000;
% zu_tab = zu_tab - mean(zu_tab(:)); % on ramène la moyenne au niveau de la mer, pour gérer 'à l'arrache' le drift en focus!!
ref_file = dir('carto\*.stk_rect_ref.txt');
if ~isempty(ref_file)
    rect_ref = load(['carto' filesep ref_file(1).name]);
    rect_ref = rect_ref*pxl_size;
    xmin = rect_ref(1);
    ymin = rect_ref(2);
    xmax = rect_ref(1) + rect_ref(3);
    ymax = rect_ref(2) + rect_ref(4);
    
    ref = (xu_tab>xmin) & (xu_tab<xmax) & (yu_tab>ymin) & (yu_tab<ymax);
    z_ref = zu_tab(ref);
    z_ref = mean(z_ref(:));
    fprintf('z_ref = %g µm\r', z_ref)
    zu_tab = zu_tab - z_ref;
else
    disp('warning, no ''rect_ref.txt'' file found in carto!')
end

clear tab_param

zmin = 0;%(gamma_z - Boule_free*zr)/1000; % min(zu_tab(alpha_tab>0));
zmax = 2*(gamma_z + Boule_free*zr)/1000; % max(zu_tab(alpha_tab>0)); % zr = 0.2µ => max = 1.4µ

% caxis([zmin zmax]), colorbar, ylabel('z (µm)'), set(gca,'YAxisLocation','right')

% r2min = 0; r2max = 0.2;
log_r2min = -3; log_r2max = 0;

Ncol = 100;
if strcmp(color_code,'t'), traj_colors = hsv(Ncol);
elseif strcmp(color_code,'z'), traj_colors = jet(Ncol);
elseif strcmp(color_code,'r2'), traj_colors = hot(Ncol);
end

DIC_name = dicname(file);
DIC_image(file, DIC_name, 1, zmin, zmax) % triD=1

xlabel('x (µm)'), ylabel('y (µm)'), zlabel('z (µm)')
title({cd; file}, 'interpreter', 'none')

pause(.1)

%% --- go through traces ---
disp('traj :            ')

for itrc = 1:ntrc
    
    if mod(itrc,50)==0
        fprintf([repmat('\b',1,11) '%5i/%5i'],itrc,ntrc)
        pause(0.01)%figure(gcf)%drawnow
    end
    
    alpha = alpha_tab(:, itrc);
    Ni = find(alpha>0); % 'valid steps'
    
    if isempty(color_code) % traces at once: much faster!!
        plot3(xu_tab(Ni, itrc), yu_tab(Ni, itrc), zu_tab(Ni, itrc));
    else
        t = t_tab(:, itrc);
        t = t(Ni);
        dt = diff(t); % dt>1 => blink, gap in trace
        
        if length(Ni)>min_trc_length
            for istep = 1:length(Ni)-1
                %%  coordonnées x,y,z du pas en cours
                xi = xu_tab(Ni(istep), itrc);
                xii = xu_tab(Ni(istep+1), itrc);
                yi = yu_tab(Ni(istep), itrc);
                yii = yu_tab(Ni(istep+1), itrc);
                zi = zu_tab(Ni(istep), itrc);
                zii = zu_tab(Ni(istep+1), itrc);
                
                if strcmp(color_code, 't')
                    val = t(istep)/Tmax; % codage == temps % icolor = min(icolor,Tmax); icolor = max(icolor,1);
                elseif strcmp(color_code, 'z')
                    val = (zi-zmin)/(zmax-zmin);
                elseif strcmp(color_code, 'r2')
%                     r2 = (xii-xi)^2 + (yii-yi)^2 + (zii-zi)^2;
%                     val = (r2-r2min)/(r2max-r2min);
                    log_r2 = log10((xii-xi)^2 + (yii-yi)^2 + (zii-zi)^2);
                    val = (log_r2-log_r2min)/(log_r2max-log_r2min);
                end
                val = max(val,eps); val = min(val,1);
                icolor = ceil(val*Ncol);
                line_color = traj_colors(icolor,:);
                
                %%  ** plot du segment istep **
                if dt(istep)==1 % pas blink
                    plot3([xi xii], [yi yii], [zi zii], 'color', line_color)
                else % blink
                    plot3([xi xii], [yi yii], [zi zii], 'color', line_color, 'linestyle', ':')
                end
            end % for istep = 1:length(Ni)-1
        end % if length(Ni)>min_trc_length
    end % if ~color_code
end % for itrc = 1:NoTrace
fprintf([repmat('\b',1,11) '%5i/%5i\r'],ntrc,ntrc)
ok = 1;
% end % for i = 1:length(files)
%%%