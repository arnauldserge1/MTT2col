% function part = reconnect_part(traj, t, lest, wn)
%
%
% EN/ function that search/detect the particle
% among those pre-detected, which best corresponds
% to the given trajectory
% returns the number of the particle in lest
%
%
% FR/ fonction qui cherche/detecte la particule
% parmi celles predetectee qui correspond a la
% trajectoire donnee
% renvoie le numero de la particule dans lest


function part = reconnect_part(traj, t, lest, wn)

%% evite les allocations sur des constantes
% global tab_param ;
% global tab_var ;
global Nb_combi ;
global stderr % AS 12/12/7

%%% Pre-detection
%%% liste_param = [num, i, j, alpha, sig^2, rayon, ok]

%%% Estimation/Reconnexion
%%%              1    2  3     4       5          6          7      8
%%% tab_param = [num, t, i,    j,      alpha,     rayon,     m0,   ,blink] 
%%% tab_var =   [num, t, sig_i,sig_jj, sig_alpha, sig_rayon, sig_b ,blink] 

%% reconnexion de la particule num part

%% recherche des particules dans la boule
%% de diffusion libre de la traj
ind_boule = liste_part_boule([], traj, lest, t) ;
nb_part_boule = size(ind_boule, 1) ;


if (nb_part_boule == 0)
  %% la trajectoirn a blinkee
  part = 0 ;
else
  %% plusieurs particule candidate
  %% on prend celle la plus probable

  %% on recherche si d'autre trajectoires sont
  %% en competition sur les particules de la boule
  %% cad si leur boule de recherche intersecte
  %% les particules en test
  vec_traj_inter = liste_traj_inter(ind_boule, traj, lest, t) ;
  vec_traj = [traj; vec_traj_inter] ;
  nb_traj = size(vec_traj,1) ;

  %% limitation du nb de trajectoire
  %% en competition
  if (nb_traj > Nb_combi)
    fprintf(stderr, '--> limitation combi for traj : traj %d\n', traj) ;
    vec_traj = limite_combi_traj_blk(vec_traj, t) ;
    %%vec_traj = limite_combi_traj_dst(vec_traj, t) ;
    nb_traj = Nb_combi ;
  end %if

  %% prise en compte des particules
  %% qui seraient dans les boules de recherche
  %% des trajectoires en competition
  %% "On travail ici a l'ordre 1"
  ind_boule_O1 = ind_boule ;
  for ntraj = 2:nb_traj
    traj_inter = vec_traj(ntraj) ;
    ind_boule_O1 =  liste_part_boule(ind_boule_O1, traj_inter, lest, t) ;
  end %for
  nb_part_boule_O1 = size(ind_boule_O1, 1) ;

  %% passage a l_ordre 1
  ind_boule = ind_boule_O1 ;
  nb_part_boule = nb_part_boule_O1 ;

  %% limite du nb de particules
  %% en competition
  if (nb_part_boule > Nb_combi)
    fprintf(stderr, '--> limitation combi for part : traj %d\n', traj) ;
     ind_boule = limite_combi_part_dst(traj,ind_boule,lest,t) ;
     nb_part_boule = Nb_combi ;
  end %if

  %% recherche de la meilleur config par 
  %% maximum de vraisemblance
  if nb_traj <= nb_part_boule
    %% pas le blink
    %% mais possiblilite de nouvelle particule (<)
    vec_part = ind_boule ;
  else
    %% nb_traj > nb_part_boule
    %% des particules ont blinkees
    vec_part = [ind_boule; zeros(nb_traj-nb_part_boule,1)] ;
  end %if

  part = best_config_reconnex(vec_traj, vec_part, t, lest, wn);


end %if

end %function



%% renvoie les indices des N point les plus proche 
%% de C(ic,jc)
%% dans l'ordre d'eloignement

function [indice, dist2] = N_plus_proche(ic, jc, liste_i, liste_j, N)
  dim_liste = size(liste_i(:),1) ;
  %% distances au point C
  sq_dist = (liste_i-ic).^2 + (liste_j-jc).^2 ;
  %% classement
  [sq_dist_classe, ind_classe] = sort(sq_dist) ;
  if (dim_liste > N)
    indice = ind_classe(1:N) ;
    dist2 = sq_dist_classe(1:N) ;
  else
    indice = ind_classe ;
    dist2 = sq_dist_classe ;
  end%if
end%function

%% limite le nombre de part
%% a celle les plus proche de la ref traj
%% au nombre Nb_combi (global)
function vec_part_out = limite_combi_part_dst(traj, vec_part_in, lest, t)
  global Nb_combi ;
  global tab_param ;
  global N_PARAM PARAM_I PARAM_J;

  ic = tab_param(traj, N_PARAM*t+PARAM_I) ;
  jc = tab_param(traj, N_PARAM*t+PARAM_J) ;
  tabi = lest(vec_part_in, 2) ;
  tabj = lest(vec_part_in, 3) ;
  indice = N_plus_proche(ic, jc, tabi, tabj, Nb_combi) ;
  vec_part_out = vec_part_in(indice) ;
end%function

%% limite le nombre de traj 
%% a celle les plus proche de la ref (vec_traj(1))
%% au nombre Nb_combi (global)
% function vec_traj_out = limite_combi_traj_dst(vec_traj_in, t)
%   global Nb_combi ;
%   global tab_param ;
%   tabi = tab_param(vec_traj_in, N_PARAM*t+PARAM_I) ;
%   tabj = tab_param(vec_traj_in, N_PARAM*t+PARAM_J) ;
%   indice = N_plus_proche(tabi(1), tabj(1), tabi(2:end), tabj(2:end), Nb_combi-1) ;
%   vec_traj_out = [vec_traj_in(1); vec_traj_in(indice+1)] ;
% end%function

%% limite le nombre de traj 
%% a celle les plus anciennes (en blink)
%% au nombre Nb_combi (global)
%% la premiere reste : ref
function vec_traj_out = limite_combi_traj_blk(vec_traj_in, t)
  global Nb_combi ;
  global tab_param ;
  global N_PARAM PARAM_BLINK

  tab_blk = tab_param(vec_traj_in(2:end), N_PARAM*t+PARAM_BLINK) ;
  [blk,indice] = sort(tab_blk, 'descend') ; 
  vec_traj_out = [vec_traj_in(1); vec_traj_in(indice(1:(Nb_combi-1))+1)] ;
end%function



%% function qui renvoie les trajectoires (non reconnectees)
%% dont les boules (espace libre) de recherche intersecte
%% au moins une des particules de la liste

function liste_traj = liste_traj_inter(liste_part, traj_ref, lest, t)

global tab_param ;
% global sig_free ;
global Boule_free ; 
% global T_off ; 
global N_PARAM PARAM_I PARAM_J PARAM_BLINK

%% init des trajectoires a tester
traj_boule = zeros(size(tab_param,1), 1) ;

nb_traj = size(tab_param, 1) ;
boules = Boule_free * sigij_free_blink(1:nb_traj, t) ;


for p = liste_part'

  ic = lest(p,2) ;
  jc = lest(p,3) ;

  %% fenetre de recherche des trajectoires
  %% boule de recherche en ecartype
  icm = max(ic - boules, 0) ; 
  icM = ic + boules ; 
  jcm = max(jc - boules, 0) ; 
  jcM = jc + boules ;
  
  traj_boule = traj_boule | ...
      ((tab_param(:, N_PARAM*t+PARAM_I) < icM) & (tab_param(:, N_PARAM*t+PARAM_I) > icm) & ...
      (tab_param(:, N_PARAM*t+PARAM_J) < jcM) & (tab_param(:, N_PARAM*t+PARAM_J) > jcm)) ;

 
end %for

%% trajectoires a tester
%% celles non encore reconnectees
traj_boule = traj_boule & (tab_param(:, N_PARAM*(t+1)+PARAM_BLINK) == 0)  ;

%% test 22 11 07 : celles qui ne sont pas blinkee
traj_boule = traj_boule & (tab_param(:, N_PARAM*t+PARAM_BLINK) > 0)  ;

%% on enleve la ref
traj_boule(traj_ref) = 0 ;

%% les trajectoires en question
liste_traj = find(traj_boule) ;

end %function



%% function qui renvoie les particules
%% appartenant a la boule de la traj
%% et qui ne sont pas deja dans la liste de reference
%%
%% si liste_par_ref == [] alors pas de ref

function liste_part = liste_part_boule(liste_part_ref, traj, lest, t) 

%% evite les allocations sur des constantes
global tab_param ;
% global tab_var ;
% global sig_free ;
global Boule_free ;
global N_PARAM PARAM_I PARAM_J

nb_part = size(lest, 1) ;
dim_liste = size(liste_part_ref,1) ;

if (dim_liste ~= 0)
  %% generation d_un masque
  masque_part_ref = ones(nb_part, 1) ;
  masque_part_ref(liste_part_ref) = 0 ;
end %if


%% reconnexion de la particule num part
%% centrage trajectoire
ic = tab_param(traj, N_PARAM*t+PARAM_I) ;
jc = tab_param(traj, N_PARAM*t+PARAM_J) ;
%% fenetre de recherche des particules candidates
%% boule de recherche pour reconnexion en ecartype
boule = Boule_free * sigij_free_blink(traj, t) ; 
icm = max(ic - boule, 0) ; % en ecart type sur i
icM = ic + boule ; % en ecart type
jcm = max(jc - boule, 0) ; % en ecart type sur j
jcM = jc + boule ; % en ecart type

%% liste locale a la boule (carree!!) <<<====================================!!!
part_boule = ...
    (lest(:,2) < icM) & (lest(:,2) > icm) & ...
    (lest(:,3) < jcM) & (lest(:,3) > jcm) ;

%% on enleve les part deja referencees
if (dim_liste ~= 0)
  part_boule = part_boule & masque_part_ref ;
end %if

liste_new_part = find(part_boule) ;
liste_part = [liste_part_ref; liste_new_part] ;


end %function
