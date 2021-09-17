% function fwrite_data_spt(filename, mode, text, tab_x, t)
%
% EN/ writes in ASCII format the info of a matrix tab_x
% at the end of a file, corresponding to time t (column number)
%
% openning modes:
% mode 'new': creation of the file
% mode 'end': writing of the data at the end of the file
% mode 'start': writing of the data at the beginning of the file
%
% Format of the output file:
% lines (modulo N_PARAM+1) correspond to time/image number
% columns correspond to particles,
% the number of particles increasing with time,
% at the beginning of each line, we indicate
% the current number of particles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FR/ ecrit au format ASCII les infos, d une matrice tab_x
% a la fin du fichier, correspondant au temps t (num colonne)
%
% modes d ouverture :
% mode 'new' : creation du fichier
% mode 'end' : ecriture des donnees en fin de fichier
% mode 'start' : ecriture des infos en debut de fichier
%
% Format du fichier de sortie :
% les lignes (modulo N_PARAM+1) correspondent au temps/numero d_image
% les colonnes correspondent aux particules,
% le nombre de particules augmente au cours du temps
% en debut de chaque ligne on renvoie le nombre de
% particules en cours


function fwrite_data_spt(filename, mode, text, tab_x, t_red, t)

%%% Estimation/Reconnexion
%%% SAVED_PARAMS=[1   2  3     4       5          6          7      8      9   ]
%%% tab_param = [num, t, i,    j,      alpha,     rayon,     m0,   ,blink, sig2]

global N_PARAM SAVED_PARAMS
global USE_MAT

if nargin<5,  tab_x = []; t_red = -1; end

if USE_MAT
    if nargin<6, t = t_red; end
    fwrite_data_spt_MAT(filename, mode, text, tab_x, t_red, t)
    return
end

t = t_red;




fid = -1; 

%% mode 'new' => write name, date, hour
if (strcmp(mode, 'new'))
  while fid==-1, 
      fid = fopen(filename, 'wt','native');
      pause (.01)
  end
  fprintf(fid, '############# :  (nb maxi particles, nb snapshots)\n') ;
  fprintf(fid, '# DATA_SPT : %s\n', text) ;
  fprintf(fid, '# DATA_SPT : %s ', date) ;
  clk = clock() ;
  fprintf(fid, '%.2dh%.2dm%.2ds\n', clk(4), clk(5), round(clk(6))) ;
  fclose(fid) ;
end%if

%% mode 'end' => write particule values at end of file
if (strcmp(mode, 'end'))
  while fid==-1, 
      fid = fopen(filename, 'at','native') ; 
      if fid==-1, pause (.01), end
  end
  nb_part = size(tab_x, 1) ;
  fprintf(fid, '# NEW_DATA_SPT : %d particles : %s\n', nb_part, text) ;
  for i_par = 1:N_PARAM ; %% t,i,j,alpha,...,blink, sig2
    param = SAVED_PARAMS(i_par);
    data = num2str(real(tab_x(:, N_PARAM*(t-1)+param))' , 5) ; % real ajouté 22/10/8
    fprintf(fid, '%d:%d: ', nb_part, param) ;
    fprintf(fid, '%s\n', data) ;
  end %for
  fclose(fid) ;
end%if

%% mode 'start' => write nb of particles and nb of frames at first line
if (strcmp(mode, 'start'))
  while fid==-1, 
      fid = fopen(filename, 'r+t','native') ;
      pause (.01)
      disp(['waiting for opening ' filename '... Keep on tryin...'])
  end
  nb_part = size(tab_x, 1) ;
  fprintf(fid, '%.6d %.6d', nb_part, t) ;
  fclose(fid) ;
end%if

end %function
