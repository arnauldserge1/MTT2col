% function tab_out = fread_all_params(filename, print_out)
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
% Adaptation de fread_data_spt à tous les params. AS
%%% Estimation/Reconnexion (num non renvoyé!)
%%%              (0)  1  2  3  4      5      6   7      8
%%% tab_param = (num)[t, i, j, alpha, rayon, m0, blink, sig2]

function tab_out = fread_all_params(filename, print_out)

global N_PARAM
global USE_MAT
if isempty(N_PARAM), MTTparams_def; end

if nargin<2, print_out = 1 ; end% affiche % effectué

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if USE_MAT, tab_out = fread_params_timewindow_MAT(filename, print_out); return, end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~strcmp(filename(end-3:end),'.dat') % si "filename short", .tif, .stk.. AS 25/3/9
    params_def = MTTparams_def; dirname = params_def{4};
    filename = [dirname filesep filename '_tab_param.dat'] ;
end

tab_out = [] ;

fid = fopen(filename, 'rt','native') ;
if fid==-1, disp('gloups...no data, check dir, MTT version...'), return, end

value =  fscanf(fid, '%d', 2) ;
if isempty(value), disp('fit incomplet...'), return, end

nb_part_max = value(1) ;
nb_t = value(2) ;

if nb_t*N_PARAM*nb_part_max>=memavailable/8
    fprintf('file too large: %g frames & %g particles...\r', nb_t, nb_part_max)
    return
end

tab_out = zeros(nb_t*N_PARAM, nb_part_max) ; %%(nb_part_max*N_PARAM, nb_t) ;

%% on se cale sur le debut des donnees
ligne = fgets(fid) ;
while (strcmp(ligne(1:14), '# NEW_DATA_SPT') == 0)
    ligne = fgets(fid) ;
end%while
if print_out
    disp(['reading from ', filename, '...', repmat(' ',1,9)])
end

current_line = 0;

for t=1:nb_t
    if print_out && mod(t,20)==0
        fprintf([repmat('\b',1,9) '%3.0f%% done'], 100*t/nb_t)
    end
    
    for p=2:(N_PARAM+1)
        value = fscanf(fid, '%d:', 2) ;

        nb_part = value(1) ;

        values = fscanf(fid,'%f', nb_part);
        
        current_line = current_line+1;

        tab_out(current_line, 1:nb_part) = values(:) ;
        %end%if
    end%for
    %% saut de ligne  # NEW_DATA_SPT
    fgets(fid) ; % [ligne, cnt] = fgets(fid) ;
    fgets(fid) ;

end%for
fclose(fid) ;

if print_out, fprintf([repmat('\b',1,9) '100%% done\r']), end

end %function