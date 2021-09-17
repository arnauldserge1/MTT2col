% function glrt = rapport_detection(traj, t, lest, num_test, wn)
%
%
% EN/ returns the likelihood for reconnection between 
% the particle num_test (in list lest) and the trajectory traj
%
% flag_full_alpha indicates if the particle
% is ON (FULL, belongs to the Gaussian)
%
%
% FR/ renvoie la vraisemblance de reconnexion d un
% particule num_test (dans la liste lest) a une 
% trajectoire traj
%
% flag_full_alpha indique si la particule
% est allumee a max (appartient a la gaussienne)


function [out, flag_full_alpha] = rapport_detection(traj, t, lest, part, wn)%, T)


%% evite les allocations sur des constantes
global tab_param ;
global tab_moy ;
global tab_var ;
global im_t ;
% global sig_free ; %% la diffusion libre
global T_off ;
% global Nb_STK ;
global Poids_melange_aplha ;
global Poids_melange_diff ;
global N_PARAM PARAM_I PARAM_J PARAM_ALPHA PARAM_RADIUS_I PARAM_BLINK

N = wn*wn ;% carree

%% ==============================================
%% test glrt avec parametre estime pour new traj
%% ==============================================
%%
%% glrt_1vue : le meme que dans carte_H0H1
%% mais ici calcule pour les parametres estimes
%% pour les particules nouvelles

if (traj<=0)

  %% P(x|H1)
  %% deja calcule lors de l'estimation 1 vue
  sig2_H1 = lest(part, 5 ) ;
  LxH1 = -N/2*log(sig2_H1) ; % -N/2

  %% P(x|H0)
  %% probleme lors des deflations!!!
  Pi = round(lest(part, 2)) ;
  Pj = round(lest(part, 3)) ;
  di = (1:wn)+Pi-floor(wn/2) ;
  dj = (1:wn)+Pj-floor(wn/2) ;
  im_part = im_t(di, dj) ;
  sig2_H0 = var(im_part(:)) ;
  LxH0 = -N/2*log(sig2_H0) ; % -N/2

  % glrt = -2*( LxH0 - LxH1 ) ;
  out = -2*( LxH0 - LxH1 ) ;
  return ;

end%if 


%% ==============================================
%% vraixemblance de la reconnexion traj <-> part
%% ==============================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% proba de reapparition pendant un blink (blinch)
%% gaussienne > 0
%% nb_blink
if (tab_param(traj, N_PARAM*t+PARAM_BLINK) < 0)
  nb_blink = -tab_param(traj, N_PARAM*t+PARAM_BLINK)  ;
  sig_blink = -T_off/3 ;
  Pblink = 2*inv(sqrt(2*pi)*sig_blink) * exp(-inv(2*sig_blink^2)*nb_blink^2) ;
  Lblink = log(Pblink) ;
else
  % nb_blink = 0 ;
  Lblink = 0 ;
end %if

%%%%%%%%%%%%%%%%%%%%%
%% intensite des pics
%% P(alpha|H1)
%% c'est un melange de loi uniforme et gaussienne
%% on travail a tres faible nombre d'echantillons
%% estimer la proportion uni/gauss est tres difficile

%% ancienne version _old
%% on fait plutot un test de vraisemblance uni/gauss
%% on est don, soit gaussienne, soit uniforme

%% nouvelle version : la loi est une combinaison
%% d'une loi uniforme et gaussienne
%% en effet, l'estimation conjointe des parametres
%% avec melange n'est possible qu'a grand nombre d echantillons
%% ici : a nombre d echantillons reduits on test si 
%% le comportement et plutot gaussien ou uniforme (MV)
%% si uniforme, on garde les anciens parametres
%% si gaussien, on met a jour moyenne et variance
%%
%% Attention, les para metres des lois (les stats m, var)
%% sont mis a jour par la fonction mise_a_jour_tab.m
%% dans les tableaux tab_moy et tab_var
%%
%% si modif verifer mise_a_jour_tab.m

alpha =  lest(part, 4 ) ;
alpha_moy = real(tab_moy(traj, N_PARAM*t+PARAM_ALPHA)) ; %% sur tout l'echantillon
sig_alpha = tab_var(traj, N_PARAM*t+PARAM_ALPHA) ; %% sur tout l'echantillon
%%alpha_max = imag(tab_moy(traj, N_PARAM*t+PARAM_ALPHA)) ;
alpha_max = alpha_moy ;

%% gaussienne (1)
if sig_alpha==0, Palpha_gaus = inf;
else Palpha_gaus = inv(sqrt(2*pi)*sig_alpha) * exp(-inv(2*sig_alpha^2)*(alpha-alpha_moy)^2) ;
end
%% uniforme (2)
if (alpha < alpha_max)
  Palpha_univ = 1/alpha_max ;
else
  Palpha_univ = 0.0 ;
end%if

poids = Poids_melange_aplha ;
Lalpha = log(poids*Palpha_gaus + (1-poids)*Palpha_univ) ;

if (Palpha_gaus > Palpha_univ)
  flag_full_alpha = 1 ;
else
  flag_full_alpha = 0 ;
end%if


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% vraisemblance de traj a intensitee full
%% nb_alpha_full = mod(tab_param(traj,N_PARAM*t+8), Nb_STK);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% vraisemblance position / mvt brownien (libre/confine)
%% P(n0|H1)
i0 = lest(part, 2) ;
j0 = lest(part, 3) ;
ic = tab_param(traj, N_PARAM*t+PARAM_I) ; % derniere position
jc = tab_param(traj, N_PARAM*t+PARAM_J) ;
% sig_ij_ref = tab_var(traj, N_PARAM*t+3);
sig_ij_ref = sigij_blink(traj, t) ;%% avec blk
sig_free_blk = sigij_free_blink(traj, t) ;%% avec blk  !!! 

poids = Poids_melange_diff ; %% entre gaussienne ref et gaussienne libre
if sig_ij_ref==0, Pn_ref = inf;
else Pn_ref = inv(2*pi*sig_ij_ref^2) * exp(- inv(2*sig_ij_ref^2) * ((i0-ic)^2 + (j0-jc)^2)) ;
end
if sig_free_blk==0, Pn_free = inf;
else Pn_free = inv(2*pi*sig_free_blk^2) * exp(- inv(2*sig_free_blk^2) * ((i0-ic)^2 + (j0-jc)^2)) ;
end
Ln0 = log(poids*Pn_ref + (1-poids)*Pn_free);

%%%%%%%%%%
%% P(r|H1)
r =  lest(part, 6) ;
r_ref = tab_moy(traj, N_PARAM*t+PARAM_RADIUS_I) ;
sig_r_ref = tab_var(traj, N_PARAM*t+PARAM_RADIUS_I) ;
if (sig_r_ref ~= 0)
  Lr = -0.5*log(sig_r_ref) - inv(2*sig_r_ref^2) * (r-r_ref)^2;
else
  Lr = 0 ;
end%if

%% vraisemblance
out = Lalpha + Ln0 + Lblink + 0*Lr;


end %function



