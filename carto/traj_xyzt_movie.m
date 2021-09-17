function traj_xyzt_movie(file)
% trace les trajs (xyz) image après image (t)

global N_PARAM PARAM_BLINK PARAM_I PARAM_J PARAM_K gamma_z zr Boule_free T_off
params_def = MTTparams_def;
dirname = params_def{4};
pxl_size = str2double(params_def{21})/1000; % 1 pxl = 0.16 µm

if nargin<1, files = dir('*.stk'); file = files(1).name; end
% cd('C:\Matlab_data\Vincent Qd Lcyl') file = 'cell1d_dessus.stk';%'cell1a.stk';

zmin = (gamma_z - Boule_free*zr)/1000;
zmax = (gamma_z + Boule_free*zr)/1000;

Ncol = 100; traj_colors = jet(Ncol);

fid = fopen([dirname filesep file '_tab_param.dat'], 'rt','native') ;
value =  fscanf(fid, '%d', 2) ;
fclose(fid);
nb_t = value(2) ; % % % nb_part_max = value(1) ;

bloc_size = 50; Nbloc = floor(nb_t/bloc_size);

% DIC_name = dicname(file); %** load masque ????
DIC_image(file, dicname(file), 1, zmin, zmax) % triD=1
axis off, figure(gcf)

for bloc = 1:Nbloc % on tronconne, pour limite memoire
    t_first = (bloc-1)*bloc_size+1;
    t_last = bloc*bloc_size+1;
    
    %% trc load
    % % %     trc = detect_reconnex_to_trc(file, 1, '', max(1, t_first-Sm), min(t_last+Sm, nb_t)); % trc = [ones(8,1) (1:8)' rand(8,2)];
    disp('loading data...')
    % % %     filename_full = [dirname '\' file '_tab_param.dat'] ;
    tab_param = fread_params_timewindow(file,1,t_first, min(t_last, nb_t));
    if isempty(tab_param), return, end
    
    x_tab = tab_param(PARAM_J-1:N_PARAM:end,:)*pxl_size; % x==j..
    y_tab = tab_param(PARAM_I-1:N_PARAM:end,:)*pxl_size;
    z_tab = tab_param(PARAM_K-1:N_PARAM:end,:)/1000;
    blink_tab = tab_param(PARAM_BLINK-1:N_PARAM:end,:);
    
    %% let's go
    if isempty(dir([file(1:end-4) '_mov_xyzt'])), mkdir([file(1:end-4) '_mov_xyzt']), end
    
    for t = (t_first+1):min(t_last, nb_t) % 2:(bloc_size+1) % 
        t_modulo = t-(bloc-1)*bloc_size;
        num_parts_t = find(blink_tab(t_modulo,:)>0);
        for p = num_parts_t % 1:nb_part_max % Npeaks = sum(alpha_tab(t,:)>0);
            blk = blink_tab(t_modulo,p);
            if (blk==0) || (blk==300) || (blk<T_off) %% no trace yet, new trace or long blink
                continue
            elseif blk>0 % ok
                ls = '-';
            else %% blink
                ls = ':';
            end
            tt = [t_modulo-1, t_modulo]; % rem: x,y,z @(t-1) contient 'dernière position connue' avant blink
            xx = x_tab(tt,p);
            yy = y_tab(tt,p);
            zz = z_tab(tt,p);
            
            val = (mean(zz)-zmin)/(zmax-zmin);
            val = max(val,eps); val = min(val,1);
            lc = traj_colors(ceil(val*Ncol),:);
            
            plot3(xx, yy, zz, 'LineStyle', ls, 'Color', lc); 
        end
        
        drawnow
        
        outfile = [file(1:end-4) '_mov_xyzt/' file(1:end-4) '_xyzt' num2str(t,'%04i') '.png' ];
        print('-dpng', '-r300', outfile)
        fprintf('%5i/%5i\r', t, nb_t)
    end
end
%%%