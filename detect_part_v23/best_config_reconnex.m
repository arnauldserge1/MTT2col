% function part = best_config_reconnex(vec_traj, vec_part, t, lest, wn)
%
% EN/ This function looks for the best config
% according to the ML and sends back 
% the most probable particle (O if blinked)
% The refering trace IS the first one
% of the list in vect_traj 
%
% Caution with combinatory tests beyond nb_part = 4
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FR/ cette fonction regarde la meilleur config
% au sens du MV et renvoie la particule la
% plus probable (0 si blinkee)
% la trajectoire de reference EST la premiere
% de la liste dans vec_traj
%
% attention aux tests combinatoire au dela de nb_part = 4
%
% =============================================
% 
% Exemple de calcul de Log-vraisemblance
% de reconnexion a 3 trajectoires/Particules
% reprenant les notations de la figure 3
% 
% 
% - T_domain = P_domain
% 3 trajectoires (a,k,b) pour 3 particules (1,2,3)
% on cherche la meilleure combinaison en comparant leur vraisemblance
% respective :
% L(a,1)  +  L(k,2)  +  L(b,3)  <>
% L(a,1)  +  L(k,3)  +  L(b,2)  <>
% L(a,2)  +  L(k,1)  +  L(b,3)  <>
%   etc.
% 
% - T_domain > P_domain
% 3 trajectoires (a,k,b) pour 2 particules : On en cree une 3eme OFF
% (1,2,OFF).
% Alors les tests deviennent :
% L(a,1) + L(k,2)   + L(b,OFF)  <>
% L(a,1) + L(k,OFF) + L(b,2)    <>
% L(a,2) + L(k,1)   + L(b,OFF)  <>
%   etc.
% Dans ce cas, seulement 2 termes de vraisemblance sont
% presents car le troisieme L(T,OFF), probabilite de disparition d'une 
% trajectoire, est identique pour toute trajectoires T.
% 
% - T_domain < P_domain
% 2 trajectoires (a,k) pour 3 particules (1,2,3) : On cree une nouvelle
% trajectoire (a,k,NEW)
% Alors les tests deviennent :
% L(a,1) + L(k,2)  + L(NEW,3)  <>
% L(a,1) + L(k,3)  + L(NEW,2)  <>
% L(a,2) + L(k,1)  + L(NEW,3)  <>
%   etc.
% La encore, dans ce cas, seulement 2 termes de vraisemblance sont
% presents car le troisieme terme L(NEW,P), probabilite d'apparition
% d'une nouvelle trajectoire dans l'image, est suppose constant.
% Elle ne depend ni de la position, ni de l'intensite.
% 
% ==============================================


function part = best_config_reconnex(vec_traj, vec_part, t, lest, wn)%,T)

  %% permutation circulaire sur vec_part
  %% limite a la taille nb_traj
  % nb_part = size(vec_part(:), 1) ;
  nb_traj = size(vec_traj(:), 1) ;
  best_part = vec_part(1) ;
  vec_part_ref = vec_part(1:nb_traj) ;
  vrais =  vrais_config_reconnex(vec_traj, vec_part_ref, t, lest, wn);%,T)


  tab_perms_vec_part = perms(vec_part)' ;
  nb_perms = size(tab_perms_vec_part, 2) ;

  for p=2:nb_perms
    vec_part_perm = tab_perms_vec_part(:,p) ;
    %% si plus de particules que de traj
    %% cas de nouvelle traj
    %% on test les nb_traj premiere
    vec_part_perm = vec_part_perm(1:nb_traj) ;
    vrais_tmp =  vrais_config_reconnex(vec_traj, vec_part_perm, t, ...
				       lest, wn);%, T) ;

    if (vrais_tmp > vrais)
      vrais = vrais_tmp ;
      best_part = vec_part_perm(1) ;
    end %if
  end %for

part = best_part ;


end %function





%% function qui renvoie la vraisemblance
%% d'une configuration de reconnexion/blink
%% pour plusieurs trajectoires et particules
%%
%% les cardinaux des deux vecteurs d'entrees 
%% sont forcement egaux

function vrais = vrais_config_reconnex(vec_traj, vec_part, t, lest, wn)%,T)

  %% vec_part peut avoir des valeurs nulles
  %% correspondant aux particules blinkees
  %% la proba de blink est sans a priori
  %% et c'est la meme pour tout traj
  %% toute comparaison ce fera a nombre de
  %% blink egal, donc on fait rien

  vrais = 0 ;
  nb_part = size(vec_part(:), 1) ;
  for p = 1:nb_part
    part = vec_part(p) ;
    traj = vec_traj(p) ;
    if part ~= 0
      vrais_p = rapport_detection(traj, t, lest, part, wn) ;
      vrais = vrais + vrais_p ;
    end %if
  end %for

end %function
