%% MTT: Multiple Target Tracing algorithms version 2.3
%% detect_reconnex_23
%%
%%
%%
%% This programme is compatible with Matlab and Octave software
%%
%%
%% ALGORITHM AUTHORS:
%%
%% Copyright A. SERGE(1,2,3), N. BERTAUX(4,5,6)
%%
%%
%% AFFILIATIONS :
%%
%% (1) CIML, University of Marseille, F13009, FRANCE
%%
%% (2) INSERM, UMR 631, Marseille, F13009, FRANCE
%%
%% (3) CNRS, UMR 6102, Marseille, F13009, FRANCE
%%
%% (4) Fresnel Institut - PhyTI Team - MARSEILLE - F13397 - FRANCE
%%
%% (5) CNRS, UMR 6133, Marseille, F13397, FRANCE
%%
%% (6) Ecole Centrale de Marseille - France
%%
%%
%% 15/06/2010
%%
%% ==============================
%% see MTT_param.m for parameters
%% ==============================
%%%
%%%
%%% Estimation/Reconnexion
%%% See MTT_params for list & numbers of saved parameters
%%%
%%% SAVED_PARAMS=[(1) 2  3  4  5      6      7       8      9     ]
%%% tab_param = [num, t, i, j, alpha, rayon, offset, blink, sig2_b]
%%%
%%% after estimation:
%%% liste_est = [num, i, j, alpha, sigb^2, rayon, ok, offset]
%%% (~= tab_param!!)
%%%
%%%
%%% EN/ Pre-detection
%%%
%%% if blink = 0 then the particle has not been detected; otherwise
%%% blink equals the number of  consecutive presence of the particle
%%% modif from 090307
%%% notification of the info 'part is full ON'
%%% blink equals the number of consecutive presence * Nb_STK
%%% plus the number of consecutive presence (*1) full ON
%%% any new particle always  starts at Nb_STK
%%%
%%% the table of variance gives the std
%%% estimated in the past. These std allow to
%%% estimate & detect the particle which corresponds
%%% to the trajectory at time t
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FR/ Pre-detection
%%%
%%% si blink vaut 0 c'est que la particule n'a pas ete detectee
%%% sinon blink vaut le nombre de presence consecutive de la particule
%%% modif du 090307
%%% notification de l'info 'part allumee au max'
%%% blink vaut le nombre de presence consecutive * Nb_STK
%%% plus le nombre de presence consecutive (*1) allumee aux max
%%% tout nouvelle particule commence toujours a Nb_STK
%%%
%%% le tableau de variance donne les ecartypes
%%% estime dans le passe. Ces ecartypes permettent
%%% d'estimer et detecter la particule qui correspond
%%% a la trajectoire ? l'instant t
%%%
%%%
%% =============================================
%% add in path utils_SPT/ subroutines repertoire
%% =============================================



%% DO NOT MODIFY NEXT LINES
%% DO NOT MODIFY NEXT LINES
%% DO NOT MODIFY NEXT LINES
%% DO NOT MODIFY NEXT LINES
%% DO NOT MODIFY NEXT LINES


% function detect_reconnex_23(name_stk, waitbar_params.....)

include_global

%%%%%%%%%%%%%global nb_defl % AS 31/1/8 pour test ds l'article

%%%%%if nargin==0, disp('please run MTT23'), return, end

%% load parameters
% MTT_param ; %%% now done in MTT23

%% Compat Matlab/octave
if strcmp(VERSION, 'MATLAB')
    stderr = 1 ;
else
    clear stderr ;
end %if

%% define other paramter
cmd_output = ['outfile = sprintf(''%s/out_%.4d.', FORMAT_IM, ''', output_dir, t) ;' ];
if SAVE_VAR_MOY == 1
    output_file_var = [full_output_dir, '/', name_stk, '_tab_var.dat'] ;
    output_file_moy = [full_output_dir, '/', name_stk, '_tab_moy.dat'] ;
end

demi_wn = ceil(wn/2) ;

warning off images:imshow:magnificationMustBeFitForDockedFigure

%%%%%%%%%%%%%%%%%%%%%%
%% PREMIERE ITERATION
%%%%%%%%%%%%%%%%%%%%%%
last = clock;
%% gestion premiere image
%% initialisation des tableaux de sorties

if use_stk && isempty(Nb_STK)
    [im_t, Nb_STK] = tiffread(name_stk, 1) ; % test/Nb_STK, AS 9/6/8 (si defini global)
% % %     Nb_STK = min(Nb_STK, Nb_STK_read); % AS 29/3/14
% elseif use_stk && ~isempty(Nb_STK)
%     im_t = tiffread(name_stk, 1) ;
else % typically for multi-tiff
    im_t = double(imread(name_stk,1)) ;
    if isempty(Nb_STK), Nb_STK = length(imfinfo(name_stk)); end %%% ASSEZ LONG................
end

global tab_num
% % % % % % % % % if isempty(tab_num), tab_num = 1:Nb_STK ; end % cf STORM_ROI
tab_num = 1:Nb_STK ;

t = tab_num(1) ;
fprintf(stderr, 'reading of %s # %d\n', name_stk, t) ;

if CROP
    im_t = im_t(IRANGE, JRANGE) ;
end%if
[idim, jdim] = size(im_t) ;

%% nouveau fichier de sortie
fwrite_data_spt(output_file_param, 'new', name_stk) ;
if SAVE_VAR_MOY == 1
    fwrite_data_spt(output_file_var, 'new', name_stk) ;
    fwrite_data_spt(output_file_moy, 'new', name_stk) ;
end

%% Detection et Estimation sur PREMIERE image
fprintf(stderr, 'Detection & Estimation in first image\n') ;
lest = detect_et_estime_part_1vue_deflt(im_t, wn, r0, seuil_detec_1vue, nb_defl);

test_all = lest(:,7) & (lest(:,4)>seuil_alpha) & ...
    (lest(:,2)>demi_wn) & (lest(:,2)<idim-demi_wn) & ...
    (lest(:,3)>demi_wn) & (lest(:,3)<jdim-demi_wn) ;
[ind_valid, tmp] = find(test_all) ;
nb_valid = size(ind_valid,1) ; % AS 27/3/7 max(size(ind_valid)) ;
fprintf(stderr, 'nb validated particles: %d\n', nb_valid) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%pk_ok = nb_valid/size(lest,1) ;

%% t == T-T_off dans les tab_x reduits
t_red = T-T_off ;

%% on alloue les tableaux de sorties
tab_param = zeros(nb_valid, 1+N_PARAM*(t_red+1)) ;
tab_var = zeros(nb_valid, 1+N_PARAM*(t_red+1)) ;

%% on initialise les tableaux de sorties
tab_param(:,1) = (1:nb_valid)' ;
lest2tab_param((1:nb_valid), 1, t_red-1, lest, ind_valid, 'new'); % traj=(1:nb_valid), t=1, part=ind_valid
tab_param(:,2+N_PARAM*t_red) = 2*ones(nb_valid,1) ;

%% valeurs initiales (par default)
%% dans calcul_reference(traj, t, param, T)
tab_var(:,1)   = (1:nb_valid)' ;
tab_var(:,2+N_PARAM*(t_red-1):1+N_PARAM*t_red) = [ones(nb_valid,1), ...
    zeros(nb_valid,N_PARAM-3), ...
    Nb_STK*ones(nb_valid,1), ...
    zeros(nb_valid,1) ] ; %%%%lest(ind_valid,5) %% sig2_b%%%%%%%%%%%%%%%%%%%%pk_ok*ones(nb_valid,1) ] ;
tab_var(:,2+N_PARAM*t_red) = 2*ones(nb_valid,1) ;
tab_moy = tab_var ;
%% mise_a_jour_tab(t_red-1, T) ;
init_tab(t_red-1) ; %% une seul fois ! par traj



if (AFFICHAGE)
    [R,V,B] = affiche_trajectoire(im_t, t_red, max(im_t(:)), min(im_t(:)));%,liste_part,1, AFF_NUM_TRAJ) ; AS 10/10/13..
    eval(cmd_output) ;
    if strcmp(VERSION, 'MATLAB'),
        imwrite(cat(3,R,V,B)/255, outfile) ; %%% Matlab
        if (SHOW)
            imshow(cat(3,R,V,B)/255) ;
        end %if
    else
        imwrite(outfile, R, V, B, imwrite_option) ; %%% Octave
        if (SHOW)
            imshow(R/255,V/255,B/255) ;
        end %if
    end %if
    if (SHOW==inf), pause, else pause(SHOW), end % AS 10/10/13..
end %if

%% nouveau fichier de sortie
fwrite_data_spt(output_file_param, 'end', '', tab_param, t_red, t) ;
if SAVE_VAR_MOY == 1
    fwrite_data_spt(output_file_var, 'end', '', tab_var, t_red) ;
    fwrite_data_spt(output_file_moy, 'end', '', tab_moy, t_red) ;
end


%%%%%%%%%%%%%%%%%%%%%
%% BOUCLE SUR LA PILE
%%%%%%%%%%%%%%%%%%%%%
% global dt
dt = zeros(1,Nb_STK-1);
for t = tab_num(2:end)
    
    waitbar_params = time_remain(t, tab_num, waitbar_params, name_stk, MTT_version);
    
    %% lecture image courante
    fprintf(stderr, 'reading of %s # %d\n', name_stk, t) ;
    if use_stk
        im_t = tiffread(name_stk, t) ;
    else
        im_t = double(imread(name_stk,t)) ;
    end
    
    if CROP
        im_t = im_t(IRANGE, JRANGE) ;
    end%if
    [idim, jdim] = size(im_t) ;
    
    
    %% Detection et Estimation sur image courante
    fprintf(stderr, 'Detection & Estimation in current image\n') ;
    lest = detect_et_estime_part_1vue_deflt(im_t, wn, r0, seuil_premiere_detec, nb_defl);
    
    %% test sur ok et alpha>0 et position
    test_all = lest(:,7) &... %% (lest(:,4)>seuil_alpha) &\
        (lest(:,2)>demi_wn) & (lest(:,2)<idim-demi_wn) &...
        (lest(:,3)>demi_wn) & (lest(:,3)<jdim-demi_wn) ;
    [ind_valid, ~] = find(test_all) ;
    nb_valid = size(ind_valid,1) ; % max(size()) => n=1... bug corrig? le 24/4/8 AS
    fprintf(stderr, 'nb validated particles: %d\n', nb_valid) ;
    lest = lest(ind_valid, :) ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%pk_ok = nb_valid/size(lest,1) ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% BOUCLE SUR LES TRAJECTOIRES ACTIVES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %     if do_reconnect % AS 21/6/10
    fprintf(stderr, 'Reconnexion of particles in current image\n') ;
    
    %% classement par ordre decroissant de blink
    %% on commence par la plus ancienne allumee
    %% et on finit par la plus ancienne blinkee
    [tmp, part_ordre_blk] = sort(-tab_param(:,N_PARAM*(t_red-1)+PARAM_BLINK)) ;
    
    nb_traj_active = 0 ;
    nb_traj_blink = 0 ;
    for traj = part_ordre_blk'
        
        %%fprintf(stderr,'traj :%d\n',traj) ;
        %%imagesc(tab_param), figure(gcf)
        
        %% test de reconnexion
        if (tab_param(traj, N_PARAM*(t_red-1)+PARAM_BLINK) > T_off)
            part = reconnect_part(traj, t_red-1, lest, wn) ;
            nb_traj_active = nb_traj_active + 1 ;
            if part == 0
                nb_traj_blink = nb_traj_blink + 1 ;
            end%if
        else
            %% on ne test plus, traj off
            %% elle reste dans les tableaux
            part = 0;
        end %if
        
        %% mise a jour des parametres estimes
        if (part>0)
            lest2tab_param(traj, t, t_red, lest, part, 'update') ;
            %% modif + new
            if (tab_param(traj, N_PARAM*(t_red-1)+PARAM_BLINK) > 0)
                [LV_traj_part, flag_full_alpha] = rapport_detection(traj, t_red-1, lest, part, wn) ;
                tab_param(traj, N_PARAM*(t_red)+PARAM_BLINK) = tab_param(traj, N_PARAM*(t_red-1)+PARAM_BLINK) + Nb_STK ; %blk+Nb_STK
                if (flag_full_alpha==1)
                    %%fprintf(stderr, '+') ;
                    tab_param(traj, N_PARAM*(t_red)+PARAM_BLINK) = tab_param(traj, N_PARAM*(t_red)+PARAM_BLINK) + 1 ; %blk+1
                else
                    %%fprintf(stderr, '0') ;
                    tab_param(traj, N_PARAM*(t_red)+PARAM_BLINK) = ...
                        tab_param(traj, N_PARAM*(t_red)+PARAM_BLINK) - mod(tab_param(traj, N_PARAM*(t_red)+PARAM_BLINK), Nb_STK);
                end%if
                %% fin new
            else
                tab_param(traj, N_PARAM*(t_red)+PARAM_BLINK) = Nb_STK ;
            end %if
        else
            tab_param(traj, N_PARAM*(t_red)+PARAM_I) = tab_param(traj, N_PARAM*(t_red-1)+PARAM_I) ;
            tab_param(traj, N_PARAM*(t_red)+PARAM_J) = tab_param(traj, N_PARAM*(t_red-1)+PARAM_J) ;
            tab_param(traj, N_PARAM*(t_red)+PARAM_ALPHA) = 0 ;
            tab_param(traj, N_PARAM*(t_red)+PARAM_RADIUS_I) = 0 ;
            tab_param(traj, N_PARAM*(t_red)+PARAM_SIG2) = 0 ;
            
            if (tab_param(traj, N_PARAM*(t_red-1)+PARAM_BLINK) < 0)
                tab_param(traj, N_PARAM*(t_red)+PARAM_BLINK) = tab_param(traj, N_PARAM*(t_red-1)+PARAM_BLINK)-1; %blink -1
            else
                tab_param(traj, N_PARAM*(t_red)+PARAM_BLINK) = -1 ;
            end %if
            
            %% 11/07/06
            %% test trajectoire 'ephemere' : mise a off
            if (t>3)
                if (tab_param(traj, N_PARAM*(t_red-2)+PARAM_BLINK) == 0)
                    tab_param(traj, N_PARAM*(t_red)+PARAM_BLINK) = T_off ;
%                     fprintf(stderr, '--> turning OFF traj %d\n', traj) ; % commented 3/12/2013
                end %if
            end %if
            
        end %if
        
        
        %% on enleve la particule reconnectee
        %% en lui affectant des coord <0 (*-1)
        if (part>0)
            lest(part, 2) = -lest(part, 2) ;
            lest(part, 3) = -lest(part, 3) ;
        end %if
        
    end %for traj = part_ordre_blk'
    %     end % do_reconnect
    
    nb_traj_avant_new = size(tab_param, 1) ;
    
    %% boucle sur les particules restante
    %% nouvelles trajectoires
    nb_non_aff_detect = 0 ;
    dim_part = -1;
    % [dim_part, dim_tps] = size(tab_param) ;
    % tab_param = [tab_param; zeros(nb_valid, dim_tps)];
    for p=1:nb_valid
        if (lest(p, 2) > 0)
            glrt_1vue = rapport_detection(0, 0, lest, p, wn) ;
            if ((glrt_1vue > seuil_detec_1vue) && (lest(p,4) > seuil_alpha))
                nb_non_aff_detect = nb_non_aff_detect+1 ;
                %                 if do_reconnect % AS 21/6/10
                [dim_part, dim_tps] = size(tab_param) ;
                tab_param = [tab_param; dim_part+1, zeros(1,dim_tps-1)] ;
                tab_var = [tab_var; dim_part+1, zeros(1,dim_tps-1)] ; %#ok
                tab_moy = [tab_moy; dim_part+1, zeros(1,dim_tps-1)] ; %#ok
                %                 else
                %                     dim_part = dim_part+1;
                %                     [dim_var, dim_tps] = size(tab_var) ;
                %                     if dim_part>dim_var-1
                %                         tab_var = [tab_var; dim_part+1, zeros(1,dim_tps-1)] ; %#ok
                %                         tab_moy = [tab_moy; dim_part+1, zeros(1,dim_tps-1)] ; %#ok
                %                     end
                %                 end
                
                lest2tab_param(dim_part+1, t, t_red, lest, p, 'new_traj') ;
                
            end %if
        end %if
    end %for
    fprintf(stderr, 'nb validated new particles: %d\n', nb_non_aff_detect) ; % 9/1/9
    
    %% Compat Matlab/octave
    if strcmp(VERSION, 'OCTAVE')
        fflush(stdout); %%% Octave
    end %if
    
    %% decalage vers la gauche de tab_x
    new_nb_traj = size(tab_param, 1) ;
    tab_param = [tab_param(:,1), ...
        tab_param(:,2+N_PARAM*(1):1+N_PARAM*(t_red+1)), ...
        (t+1)*ones(new_nb_traj,1), zeros(new_nb_traj,(N_PARAM-1))] ;%correction arnauld
    
    %% pour tab_var
    tab_var = [tab_var(:,1), ...
        tab_var(:,2+N_PARAM*(1):1+N_PARAM*(t_red+1)), ...
        (t+1)*ones(new_nb_traj,1), zeros(new_nb_traj,(N_PARAM-1))] ;
    %% pour tab_moy
    tab_moy = [tab_moy(:,1), ...
        tab_moy(:,2+N_PARAM*(1):1+N_PARAM*(t_red+1)), ...
        (t+1)*ones(new_nb_traj,1), zeros(new_nb_traj,(N_PARAM-1))] ;
    
    init_tab(t_red-1, (1+nb_traj_avant_new):new_nb_traj) ;
    mise_a_jour_tab(t_red-1, T) ;
    
    %% enregistrement des nouveaux fichiers de sortie
    fwrite_data_spt(output_file_param, 'end', '', tab_param, t_red, t) ;
    if SAVE_VAR_MOY == 1
        fwrite_data_spt(output_file_var, 'end', '', tab_var, t_red) ;
        fwrite_data_spt(output_file_moy, 'end', '', tab_moy, t_red) ;
    end
    
    if (AFFICHAGE)
        [R,V,B] = affiche_trajectoire(im_t, t_red, max(im_t(:)), min(im_t(:)));%, liste_part, t, AFF_NUM_TRAJ) ; AS 10/10/13
        eval(cmd_output) ;
        if strcmp(VERSION, 'MATLAB'),
            imwrite(cat(3,R,V,B)/255, outfile) ; %%% Matlab
            if (SHOW)
                imshow(cat(3,R,V,B)/255) ;
            end %if
        else
            imwrite(outfile, R, V, B, imwrite_option) ; %%% Octave
            if (SHOW)
                imshow(R/255,V/255,B/255) ;
            end %if
        end %if
        if (SHOW==inf), pause, else pause(SHOW), end % AS 10/10/13..
    end %if
    
    
end %%for t=tab_num

%% enregistrement des nouveaux fichiers de sortie: mise ? jour de l'entete
fwrite_data_spt(output_file_param, 'start', '', tab_param, tab_num(end), t) ; % t_red remplac? par t = tab_num(end)...
if SAVE_VAR_MOY == 1
    fwrite_data_spt(output_file_var, 'start', '', tab_var, tab_num(end)) ;
    fwrite_data_spt(output_file_moy, 'start', '', tab_moy, tab_num(end)) ;
end

clear global tab_num % cf STORM_ROI

%%%