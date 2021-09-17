% function out = fread_data_spt(filename, param)
%
% EN/ extracts from the input file, the info concerning
% a given parameter, for all trajectories
%
% Format of the output file:
% lines (modulo N_PARAM+1) correspond to time/image number
% columns correspond to particles,
% the number of particles increasing with time,
% at the beginning of each line, we indicate
% the current number of particles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FR/ extrait du fichier d'entree, les infos concernant
% un parametre donne, pour toutes les trajectoires
%
% Format du fichier de sortie :
% les lignes (modulo N_PARAM+1) correspondent au temps/numero d_image
% les colonnes correspondent aux particules,
% le nombre de particules augmentant au cours du temps,
% en debut de chaque ligne on renvoie le nombre de
% particules en cours


function tab_out = fread_data_spt(filename, param)

%%% Estimation/Reconnexion
%%%              1    2  3  4  5      6      7   8      9
%%% tab_param = [num, t, i, j, alpha, rayon, m0, blink, sig2] 

  global N_PARAM
  if isempty(N_PARAM), MTTparams_def; end


  fid = fopen(filename, 'rt','native') ;
  value =  fscanf(fid, '%d', 2) ;
  nb_part_max = value(1) ;
  nb_t = value(2) ;

  tab_out = zeros(nb_part_max, nb_t) ;

  %% on se cale sur le debut des donnees
  ligne = fgets(fid) ;
  while (strcmp(ligne(1:14), '# NEW_DATA_SPT') == 0)
    ligne = fgets(fid) ;
  end%while

  disp(['reading from ', filename, '...', repmat(' ',1,9)]) % if print_out % AS 11/3/9
    
  for t=1:nb_t
      if mod(t,20)==0 %|| t==nb_t%&& print_out
          fprintf([repmat('\b',1,9) '%3.0f%% done'], 100*t/nb_t)
      end
    
    for p=2:(N_PARAM+1)
      value =  fscanf(fid, '%d:', 2) ;
      nb_part = value(1) ;
      param_lu = value(2) ;
      
      ligne = fgets(fid) ;

      if (param_lu == param)
	[values, cnt] =  sscanf(ligne(2:end), '%f', nb_part) ; %#ok
 	tab_out(1:nb_part, t) = values(:);
      end%if
    end%for
  %% saut de ligne  # NEW_DATA_SPT
  [ligne, cnt] = fgets(fid) ; %#ok

  end%for
  fclose(fid) ;
  fprintf([repmat('\b',1,9) '100%% done\r']) %, if print_out, end
end %function
