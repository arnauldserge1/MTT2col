% function out = fread_params_timewindow_MAT(filename, print_out, t_first, t_last, n_first, n_last)
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

function tab_out = fread_params_timewindow_MAT(filename, print_out, t_first, t_last, n_first, n_last)

%%% Estimation/Reconnexion (num non renvoyé!)
%%%              (0)  1  2     3       4          5          6      7
%%% tab_param = (num)[t, i,    j,      alpha,     rayon,     m0,   ,blink]

global N_PARAM ;

if ~strcmp(filename(end-3:end),'.mat') % si "filename short", .tif, .stk.. AS 25/3/9
    filename = ['output23\' filename '_tab_param.mat'] ;
end

if nargin<2, print_out = 1 ; end % affiche % effectué
if nargin<3, t_first = 1 ; end
if nargin<4, t_last = inf ; end
if nargin<5, n_first = 1 ; end
if nargin<6, n_last = inf ; end

if isempty(dir(filename)), disp('gloups...no data'), tab_out = [] ; return, end

load('-mat', filename, 's1')
if strcmp(s1(1),'#'), disp('fit incomplet...'), tab_out = [] ; return, end

if print_out
    if t_first>1 || t_last<inf, str = [', image ' num2str(t_first) ' to ' num2str(t_last)]; else str = ''; end
    if n_first>1 || n_last<inf, str = [str ', traces ' num2str(n_first) ' to ' num2str(n_last)]; end
    disp(['reading from ', filename, str, '............'])
end

% n_files = whos('-file', filename);
% t_tot = n_files - 4; % car 4 fichiers text s1..4

n_tot = str2double(s1(1:6));
t_tot = str2double(s1(7:7+6));

n_last = min(n_last, n_tot) ; % n_last = inf par défaut, pour tout lire
n_max = n_last-n_first+1;
t_last = min(t_last, t_tot) ; % t_last = inf par défaut, pour tout lire
t_max = t_last-t_first+1;

if t_max*N_PARAM*n_max>=memavailable/8
    tab_out = [];
    fprintf('file too large: %g frames & %g particles...\r', t_max, n_max)
    return
end

tab_out = zeros(t_max*N_PARAM, n_max) ;
data_name = cell(1,t_max);
list_name = cell(1,t_max);

for t=1:t_first-1
    list_name{t} = '';
end
for t=t_first:t_last
    data_name{t} = sprintf('data%05i',t);
    list_name{t} = [', ''' data_name{t} ''''];
end

list_names = cell2mat(list_name);
eval(['load(''-mat'', filename' list_names ')'])

for t=t_first:t_last
    
    if exist(data_name{t},'var')
        eval(['data_t = ' data_name{t} '(:,3:end) ;'])
    else
        disp(['Hey, ' data_name{t} ' not found!!'])
        continue
    end
    
    ind = (t-1)*N_PARAM+1:t*N_PARAM;
    nb_part = size(data_t,2); %#ok nb of particles in current frame t
    
    if nb_part<n_last % current frame has less particles
        tab_out(ind, 1:nb_part-n_first+1) = data_t ;
    else % current frame has more particles (=> when collecting only a subset, n_last<inf)
        tab_out(ind, 1:n_last-n_first+1) = data_t(:,n_first:n_last) ;
    end
    
end%for
%%%