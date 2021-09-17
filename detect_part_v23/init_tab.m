% function init_tab(t, new_traj)
%
% EN/ initialisation of tables of values
% mean of parameters and variances (std)
%
% if input new_traj non null
% then only this traj is init
% otherwise all trajectories are (at the beginning)
%
%
% FR/ initialisation des tableaux de valeurs
% moyenne des parametres et des variances (std)
%
% si entree new_traj non nulle
% alors seule cette traj est init
% sinon toutes les trajectoires le sont (au debut)


function init_tab(t, new_traj) %% (t, T, new_traj)

global tab_param ;
global tab_moy ;
global tab_var ;
global sig_free ;

global N_PARAM PARAM_I PARAM_J PARAM_ALPHA PARAM_RADIUS_I PARAM_BLINK % PARAM_T PARAM_OFFSET PARAM_SIG2 PARAM_K PARAM_RADIUS_J

 coef = 0.5; % Ratio_def std/mean 1/7/8

if (nargin < 2)
  nb_traj = size(tab_param, 1) ;
  tab_traj = 1:nb_traj ;
  new = 0 ;
else
  tab_traj = new_traj ;
  new = 1 ;
end%if

%% boucle sur les particules
for traj = tab_traj

%% alpha
param = PARAM_ALPHA ;
local_param = tab_param(traj, N_PARAM*t+param) ; % alpha
tab_moy(traj, N_PARAM*t+param) = local_param + sqrt(-1)*local_param ;% moyenne,max
tab_var(traj, N_PARAM*t+param) = coef*local_param ; % std

%% r
param = PARAM_RADIUS_I ;
local_param = tab_param(traj, N_PARAM*t+param) ;
tab_moy(traj, N_PARAM*t+param) = local_param ;
tab_var(traj, N_PARAM*t+param) = coef*local_param ; 
  
%% i,j
param = PARAM_I ;
tab_var(traj, N_PARAM*t+param) = sig_free ; 
param = PARAM_J ;
tab_var(traj, N_PARAM*t+param) = sig_free ; 

%% blink pour info
param = PARAM_BLINK ;
tab_moy(traj, N_PARAM*t+param) = tab_param(traj, N_PARAM*t+param) ;
tab_var(traj, N_PARAM*t+param) = tab_param(traj, N_PARAM*t+param) ;

%% affection identique a t-1
%% pour compat avec mise_a_jour_tab
if (new)

  %% alpha
  param = PARAM_ALPHA ;
  local_param = tab_param(traj, N_PARAM*t+param) ; % alpha
  tab_moy(traj, N_PARAM*(t-1)+param) = local_param + sqrt(-1)*local_param ;% moyenne,max
  tab_var(traj, N_PARAM*(t-1)+param) = coef*local_param ; % std
  
  %% r
  param = PARAM_RADIUS_I ;
  local_param = tab_param(traj, N_PARAM*t+param) ;
  tab_moy(traj, N_PARAM*(t-1)+param) = local_param ;
  tab_var(traj, N_PARAM*(t-1)+param) = coef*local_param ; 

  %% i,j
  param = PARAM_I ;
  tab_var(traj, N_PARAM*(t-1)+param) = sig_free ; 
  param = PARAM_J ;
  tab_var(traj, N_PARAM*(t-1)+param) = sig_free ; 
  
  %% blink pour info
  param = PARAM_BLINK ;
  tab_moy(traj, N_PARAM*(t-1)+param) = tab_param(traj, N_PARAM*t+param) ;
  tab_var(traj, N_PARAM*(t-1)+param) = tab_param(traj, N_PARAM*t+param) ;
  
end %if


end %for

end %function

