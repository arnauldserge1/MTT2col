% function [tab_i,tab_j,tab_alpha,tab_ray,tab_7,tab_blk] 
%   = fread_all_data_spt(filename)
%
% EN/ extracts from the input file, the info concerning
% the parameters for all trajectories
%
% Format of the output file:
% lines (modulo 8) correspond to time/image number
% columns correspond to particles,
% the number of particles increasing with time,
% at the beginning of each line, we indicate
% the current number of particles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FR/ extrait du fichier d'entree, les infos concernant
% les parametres pour toutes les trajectoires
%
% Format du fichier de sortie :
% les lignes (modulo N_PARAM+1) correspondent au temps/numero d_image
% les colonnes correspondent aux particules,
% le nombre de particules augmentant au cours du temps,
% en debut de chaque ligne on renvoie le nombre de
% particules en cours


function [tab_i,tab_j,tab_alpha,tab_ray,tab_7,tab_blk,tab_sig2]= fread_all_data_spt(filename)

%%% Estimation/Reconnexion
%%%              1    2  3  4  5      6      7   8      9
%%% tab_param = [num, t, i, j, alpha, rayon, m0, blink, sig2] 

  global N_PARAM PARAM_I PARAM_J PARAM_ALPHA PARAM_RADIUS_I PARAM_OFFSET PARAM_BLINK PARAM_SIG2 % PARAM_K PARAM_RADIUS_J
  
  fid = fopen(filename, 'rt','native') ;
  value =  fscanf(fid, '%d', 2) ;
  nb_part_max = value(1) ;
  nb_t = value(2) ;

  %% declaration des tableaux
  tab_i = zeros(nb_part_max, nb_t) ;
  tab_j = tab_i ;
  tab_alpha = tab_i ;
  tab_ray = tab_i ;
  tab_7 = tab_i ;
  tab_blk = tab_i ;
  tab_sig2 = tab_i ;

  %% on se cale sur le debut des donnees
  ligne = fgets(fid) ;
  while (strcmp(ligne(1:14), '# NEW_DATA_SPT') == 0)
    ligne = fgets(fid) ;
  end%while

  for t=1:nb_t
    for p=2:N_PARAM+1
      value =  fscanf(fid, '%d:', 2) ;
      nb_part = value(1) ;
      param_lu = value(2) ;
      
      ligne = fgets(fid) ;
      [values, cnt] =  sscanf(ligne(2:end), '%f', nb_part) ; %#ok

      if (param_lu == PARAM_I)
 	tab_i(1:nb_part, t) = values(:);
      end%if
      if (param_lu == PARAM_J)
 	tab_j(1:nb_part, t) = values(:);
      end%if
      if (param_lu == PARAM_ALPHA)
 	tab_alpha(1:nb_part, t) = values(:);
      end%if
      if (param_lu == PARAM_RADIUS_I)
 	tab_ray(1:nb_part, t) = values(:);
      end%if
      if (param_lu == PARAM_OFFSET)
 	tab_7(1:nb_part, t) = values(:);
      end%if
      if (param_lu == PARAM_BLINK)
 	tab_blk(1:nb_part, t) = values(:);
      end%if
      if (param_lu == PARAM_SIG2)
 	tab_sig2(1:nb_part, t) = values(:);
      end%if

    end%for
  %% saut de ligne  # NEW_DATA_SPT
  [ligne, cnt] = fgets(fid) ; %#ok

  end%for
  fclose(fid) ;

end %function
