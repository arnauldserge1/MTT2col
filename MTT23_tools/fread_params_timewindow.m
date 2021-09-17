% function out = fread_params_timewindow(filename, print_out, t_first, t_last, n_first, n_last)
%
% extract du fichier d'entree, les infos concernant
% tous les parametres pour toutes les trajectoires
%
% Format du fichier de sortie :
% les lignes (modulo N_PARAM) correspond au temps/numero d_image
% les colonnes correspondent aux particules,
% le nombre de particules augmente au cours du temps
% en debut de chaque ligne on renvoie le nombre de
% particule en cours
% Adaptation de fread_data_spt ? tous les params. AS

function tab_out = fread_params_timewindow(filename, print_out, t_first, t_last, n_first, n_last)

%%% Estimation/Reconnexion (num non renvoy?!)
%%% SAVED_PARAMS-1(0)  1  2     3       4          5          6      7      8
%%% tab_param = (num) [t, i,    j,      alpha,     rayon,     m0,   ,blink, sig2]

global N_PARAM
global USE_MAT


if nargin<2, print_out = 1 ; end % affiche % effectu?
if nargin<3, t_first = 1 ; end
if nargin<4, t_last = inf ; end
if nargin<5, n_first = 1 ; end
if nargin<6, n_last = inf ; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if USE_MAT, tab_out = fread_params_timewindow_MAT(filename, print_out, t_first, t_last, n_first, n_last); return, end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~strcmp(filename(end-3:end),'.dat') % si "filename short", .tif, .stk.. AS 25/3/9
    %    select_output_folder_version
    params_default = MTTparams_def; dirname = params_default{4};
    % % % % %     if ~isdir(dirname), dirname = 'output22'; disp('caution, using good-old output22!!!'), end
    filename = [dirname filesep filename '_tab_param.dat'] ;
end

fid = fopen(filename, 'rt','native') ;
if fid==-1, disp('gloups...no data'), tab_out = [] ; return, end

value =  fscanf(fid, '%d', 2) ; % [nb part (=n_last), nb images(=t_last)
if isempty(value), disp('fit incomplet...'), tab_out = [] ; return, end

str = '';
if print_out
    if t_first>1 || t_last<inf, str = [', image ' num2str(t_first) ' to ' num2str(t_last)]; end
    if n_first>1 || n_last<inf, str = [str ', traces ' num2str(n_first) ' to ' num2str(n_last)]; end
    disp(['reading from ', filename, str, '............'])
end

n_last = min(n_last, value(1)) ; % n_last = inf par d?faut, pour tout lire
nb_part_tot = n_last-n_first+1;
t_last = min(t_last, value(2)) ; % t_last = inf par d?faut, pour tout lire
nb_t = t_last-t_first+1;

if nb_t*N_PARAM*nb_part_tot>=memavailable/8
    tab_out = [];
    fprintf('file too large: %g frames & %g particles...\r', nb_t, nb_part_tot)
    return
end

tab_out = zeros(nb_t*N_PARAM, nb_part_tot) ; %t_last corrig? 9/1/9!!!

%% on se cale sur le debut des donnees
ligne = fgets(fid) ;
while (strcmp(ligne(1:14), '# NEW_DATA_SPT') == 0)
    ligne = fgets(fid) ;
end%while

current_line = 1;

%% saute d?but
% N_PARAM lignes + une ligne de comments (#)
% for t=1:(t_first-1)*(N_PARAM+1) % NB: skipped if t_first=0 ou 1...
%     ligne = fgets(fid) ; %#ok
% end
for t=1:t_first-1
    for p=2:(N_PARAM+1)
        value =  fscanf(fid, '%d:', 2) ;
        nb_part = value(1) ;
        fscanf(fid, '%f', nb_part) ;
    end%for
    %% saut de ligne  # NEW_DATA_SPT
    [ligne, cnt] = fgets(fid) ; %#ok
    [ligne, cnt] = fgets(fid) ; %#ok
end%for

%% read
for t=t_first:t_last
    if mod(t,10)==0 && print_out
        fprintf([repmat('\b',1,9) '%3.0f%% done'], 100*(t+1-t_first)/nb_t)
    end
    
    for p=2:(N_PARAM+1)
        value =  fscanf(fid, '%d:', 2) ;
        nb_part = value(1) ;        %param_lu = value(2) ;
        
        values =  fscanf(fid, '%f', nb_part) ;
        
        if nb_part < n_last
            tab_out(current_line, 1:nb_part-n_first+1) = values(n_first:nb_part) ;
        else
            tab_out(current_line, 1:n_last-n_first+1) = values(n_first:n_last) ;
        end
        current_line = current_line+1;
    end%for
    %% saut de ligne  # NEW_DATA_SPT
    fgets(fid) ;
    fgets(fid) ;
    
end%for
fclose(fid) ;

if print_out, fprintf([repmat('\b',1,9) '100%% done\r']), end

end %function