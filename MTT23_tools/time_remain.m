function waitbar_params = time_remain(t, tab_num, waitbar_params, name_stk, MTT_version)

% function waitbar_params = time_remain(t, tab_num, waitbar_params, name_stk, MTT_version)
% Calcul le temps restant pour traiter n images et n fichiers (supposés équivalents!)
% par extrapolation du temps déjà passé
% Met à jour la waitbar en conséquence

if nargin<5, MTT_version = ''; end

Nb_img = length(tab_num);
ifile = waitbar_params{1};
Nfiles = waitbar_params{2};
waitbar_handle = waitbar_params{3};
t_start = waitbar_params{4};
N_img_by_stk = waitbar_params{5};


%% *** calcul temps restant ***
if t==tab_num(2)
    N_img_by_stk(ifile) = Nb_img;
    waitbar_params{5} = N_img_by_stk;
end

t_elapsed = datenum(clock-t_start); % temps déjà écoulé ; unité datenum: 1=1jour
N_img_done = sum(N_img_by_stk(1:ifile-1))+t-tab_num(1); % images traitées dans ce fichier et les précédents
mean_t_by_img = t_elapsed/N_img_done; % temps moyen / image

N_img_remaining = tab_num(end)-t+1 + mean(N_img_by_stk(N_img_by_stk>0))*(Nfiles-ifile); % extrapolation for other files
t_remain_tot = N_img_remaining*mean_t_by_img; % ajoute les fichiers restants (supposés de taille équivalente)

global TEST_TIMER, TEST_TIMER(end+1)=t_remain_tot; % pour debug >> global TEST_TIMER, TEST_TIMER=[], MTT23i({'*.tif'})

%% *** affiche temps restant ***
hours_remain = num2str(floor(t_remain_tot*24));
t_remain_str = [hours_remain ':' datestr(t_remain_tot, 'MM:SS')];
waitbar_text = sprintf('Image%4d/%4d, file%3d/%3d, time remaining %s',...
    t, tab_num(end), ifile, Nfiles, t_remain_str); % suppose moins d'un jour... (t_remain_tot<1 )
disp(waitbar_text)

if ishandle(waitbar_handle) % check si fermée par utilisateur
    Nimag_done = (ifile-1)*Nb_img + t-tab_num(1);
    waitbar(Nimag_done/(Nb_img*Nfiles), waitbar_handle, waitbar_text);
    if t==tab_num(2), set(waitbar_handle,'name',['fitting ' name_stk ', MTT v' MTT_version]), end
end