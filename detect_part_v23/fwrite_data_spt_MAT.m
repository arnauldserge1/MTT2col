% function fwrite_data_spt_MAT(filename, mode, text, tab_x, t_red, t)
%
% EN/ writes in Matlab .mat format the info of a matrix tab_x
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
% FR/ ecrit au format Matlab .mat les infos, d une matrice tab_x
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


function fwrite_data_spt_MAT(filename, mode, text, tab_x, t_red, t)

%%% Estimation/Reconnexion
%%% SAVED_PARAMS=[1   2  3  4  5     6      7   8      9   ]
%%% tab_param = [num, t, i, j,alpha, rayon, m0, blink, sig2]

global N_PARAM SAVED_PARAMS

%% mode 'new'
if (strcmp(mode, 'new'))
    s1 = sprintf('############# :  (nb maxi particles, nb snapshots)') ; %#ok
    s2 = sprintf('# DATA_SPT : %s', text) ; %#ok
    s3 = sprintf('# DATA_SPT : %s ', date) ; %#ok
    clk = clock() ;
    s4 = sprintf('%.2dh%.2dm%.2ds', clk(4), clk(5), round(clk(6))) ; %#ok
    save(filename, 's1', 's2', 's3', 's4')
end%if

%% mode 'end'
if (strcmp(mode, 'end'))
    nb_part = size(tab_x, 1) ;
%     s1 = sprintf('# NEW_DATA_SPT : %d particles : %s\n', nb_part, text) ;
    
    liste_param = SAVED_PARAMS; % cf. MTT_params_v23 for list of saved parameters
    data = [nb_part*ones(size(liste_param')) liste_param' real(tab_x(:, N_PARAM*(t_red-1)+liste_param))'] ; %#ok
    data_name = sprintf('data%05i', t);%     frame_number = data(1,1);
    eval([data_name ' = data ;'])
    
    save(filename, data_name, '-append') % , 's1'
end%if

%% mode 'start'
if (strcmp(mode, 'start'))
    nb_part = size(tab_x, 1) ;
    s1 = sprintf('%.6d %.6d :  (nb maxi particles, nb snapshots)', nb_part, t) ; %#ok
    save(filename, 's1', '-append') % écrase et remplace s1
end%if

end %function