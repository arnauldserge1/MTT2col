function traj_xyt_basic(filename, codage, triD)

% function traj_xyt_basic(filename, codage, triD)
% traces en 3D (ou 2D) sur image dic (ou 1e image)
% couleur selon codage: z, t (...??)


global N_PARAM PARAM_T PARAM_I PARAM_J PARAM_Z PARAM_ALPHA
if isempty(N_PARAM), MTTparams_def; end

if nargin<2, codage = 'z'; end
if nargin<3, triD = 1; end

if nargin==0, files = dir('*.stk'); else files = dir(filename); end

if isempty(files), disp('No data... Check dir & filename !'), return, end

Ncol = 100;
traj_colors = jet(Ncol);

%% boucle / file
for i = 1:length(files)
    file = files(i).name;
    disp('loading data...')
    tab_param = fread_all_params(file);
    if isempty(tab_param), return, end
    
    ntrc = size(tab_param,2);
    Tmax = size(tab_param,1)/N_PARAM; % min_trc_length = round(Tmax*min_length_ratio);
    z_tab = tab_param(PARAM_Z-1:N_PARAM:end,:);
    zmin = min(z_tab(:));
    zmax = max(z_tab(:));

    DIC_name = dicname(file);
    DIC_image(file, DIC_name, triD, Tmax)
    
    title({cd; filename ; [' codage: ' codage]}, 'interpreter', 'none')
    pause(.1)
    
    %% --- go through traces ---
    disp('traj :            ')
    
    for itrc = 1:ntrc
        
        if mod(itrc,50)==0
            fprintf([repmat('\b',1,11) '%5i/%5i'],itrc,ntrc)
            figure(gcf)%drawnow
        end
        
        trc_par = tab_param(:,itrc); % current column of params, for current trace
        alpha = trc_par(PARAM_ALPHA-1:N_PARAM:end);
        t = trc_par(PARAM_T-1:N_PARAM:end);
        Ni = find(alpha>0);
        t = t(Ni);
        dt = diff(t);
        
        if isempty(codage)
            x = trc_par(PARAM_J-1:N_PARAM:end);
            y = trc_par(PARAM_I-1:N_PARAM:end);
            x = x(Ni);
            y = y(Ni);
            plot3(x,y,t)
        else
%         if Ni>min_trc_length
            for istep = 1:length(Ni)-1
                %%          coordonnées t,x,y du pas en cours
                ni = (Ni(istep)-1)*N_PARAM - 1; % ligne des params pour Ni
                nii = (Ni(istep+1)-1)*N_PARAM - 1; % = ni + N_PARAM (sauf si blink), ligne des params pour le suivant
                ti = trc_par(ni+PARAM_T);
                tii = trc_par(nii+PARAM_T);
                xi = trc_par(ni+PARAM_J);
                xii = trc_par(nii+PARAM_J);
                yi = trc_par(ni+PARAM_I);
                yii = trc_par(nii+PARAM_I);
                
                if dt(istep)==1 % pas blink
                    
                    if strcmp(codage, 'z')
                        val = trc_par(ni+PARAM_Z);
                        val = (val-zmax)/(zmin-zmax);
                    elseif strcmp(codage, 't')
                        val = ti/Tmax;
                    end
                    
                    icolor = ceil(val*Ncol); % 0 à 1 => 1 à N
                    icolor = min(icolor,Ncol);
                    icolor = max(icolor,1);
                    
                    line_type = traj_colors(icolor);
                    
                else % blink
                    line_type = ':';
                end
                
                %%  ** plot du segment istep **
                if triD
                    if strcmp(codage, 'z')
                        plot3 ([xi xii], [yi yii], [ti tii], line_type)
                    elseif strcmp(codage, 't')
                        plot3 ([xi xii], [yi yii], [ti tii], line_type)
                    end
                elseif dtab_par_i(istep)==1 || ~no_blink
                    line ([xi xii], [yi yii], 'color', hsv2rgb([val,1,1]))
                end
            end % for istep = 1:length(Ni)-1
%         end % if length(Ni)>0
        end
    end % for itrc = 1:NoTrace
    
    axis image
    fprintf([repmat('\b',1,11) '%5i/%5i\r'],ntrc,ntrc)
    
    if isempty(dir('MCTfigs')), mkdir('MCTfigs'), end
    figname = ['MCTfigs' filesep filename 'carto3D.png'];
    saveas(gcf, figname, 'png')
    if length(files)>1, close(gcf), end
end % for i = 1:length(files)
%%%