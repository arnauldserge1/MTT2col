%% function [test, glrt, test_acc] = detect_conf_lent(r2,Na,Nb,seuil,blkr2)
%% 
%% 
%% Programme de detection de confinements
%% dans un mouvement brownien
%% Hypothese de ralentissement du mouvement
%% correspondant au cas de confinements forts.
%%
%% Base sur les travaux de Mastere d'Arnaud Woiselle
%% Detection d'un ralentissement dans le vecteur des
%% deplacement carrees (r^2), a partir d'un test
%% sur une fenetre glissante decomposee en deux zones
%%
%%          Na     Nb
%%     |---------||--|
%%
%% avec Na de l'ordre de 10 a 15 echantillons et
%% Nb de 3 a 5 (ralentissement de courte duree)
%% Les deplacements carres (r2) d'un mouvement brownien
%% correspondent a des va gamma. On cherche ainsi
%% des ruptures de valeurs moyennes dans les r2
%%
%% 
%%
%% valeurs des seuils pour des pfa fixees (+ pd pour info)
%% pour un contraste de 10 (sur les r2) ce qui correspond
%% a un ralentissement detectable a "l'oeil"
%%
%% NA = 10, NB = 3, C=10
%% pfa     10^-3    10^-2    10^-1
%% pd      0.2      0.63     0.95
%% seuil   5.7      3.5      1.43
%%
%%
%% NA = 10, NB = 4, C=10
%% pfa     10^-3    10^-2    10^-1
%% pd      0.42     0.81     0.985
%% seuil   5.64     3.47     1.42
%%
%%
%% NA = 10, NB = 5, C=10
%% pfa     10^-3    10^-2    10^-1
%% pd      0.62     0.9      0.94
%% seuil   5.64     3.44     1.41
%%
%%
%% NA = 15, NB = 3, C=10
%% pfa     10^-3    10^-2    10^-1
%% pd      0.23     0.68     0.976
%% seuil   5.68     3.5      1.43
%%
%%
%% NA = 15, NB = 4, C=10
%% pfa     10^-3    10^-2    10^-1
%% pd      0.46     0.86     0.994
%% seuil   5.63     3.46     1.414
%%
%%
%% NA = 15, NB = 5, C=10
%% pfa     10^-3    10^-2    10^-1
%% pd      0.68     0.94     0.998
%% seuil   5.6      3.43     1.4
%%
%%
%% Parametres par defauts
%%
%% Na=10, Nb=3, seuil=4.5 
%% vecteur de blkr2 : pas de blink
%%
%% Parametre de sorties
%%
%% le test binaire de detection et
%% pour info (genre indice de confinement) le glrt.
%% Par exemple, test.*glrt donne les zones detectees
%% comme zones de ralentissment plus l'info sur la force
%% du confimement (contraste de vitesse)
%%
%%
%%
%% NOTE : Cette fonction gere le phenomene de blink
%% si on lui donne le vecteur temps pour retrouver les zones de blink
%% correspondant a la trajectoire de la particule testee
%%
%% version 1.0
%% NB 
%% Institut Fresnel (PhyTI) / CIML
%% Marseille - France

%% Procedure de test (exemple)
%% > r2 = [10*gamma_rnd(1,1,1,150), gamma_rnd(1,1,1,4), 10*gamma_rnd(1,1,1,150)];
%% > test=detect_conf_lent(r2);
%% > plot([r2; test-2; glrt-8]');


function [test,glrt,test_acc] = detect_conf_lent(r2,Na,Nb,seuil,blkr2)

  r2 = r2(:)' ;
  stack = size(r2, 2) ; % stack-1 en fait...

  %% parametres par defauts 
  if (nargin < 5)
    blkr2 = zeros(1,stack) ; %% pas de blink
  end%if
  if (nargin < 4)
    seuil = 4.5 ; 
  end%if
  if (nargin < 3)
    Nb = 3 ;
  end%if
  if (nargin < 2)
    Na = 10 ;
  end%if
  
  warning off MATLAB:log:LogOfZero;

  blkr2 = blkr2(:)';

  %% fenetres glissantes
  N=Na+Nb;
  Ta=[ones(1,Na)/Na,zeros(1,Nb)];
  Tb=[zeros(1,Na),ones(1,Nb)/Nb];
%   T=ones(1,N)/N;
  

  if stack<N % r2 trop court, < Na+Nb... AS 11/12/2009
      test = []; glrt = []; test_acc = [];
      return
  end
  
  
  %% correlations : moyennes
  TTA = [Ta, zeros(1, stack-N)] ;
  TTB = [Tb, zeros(1, stack-N)] ;
  ma=real(ifft( conj(fft(TTA)) .* fft(r2) )); %%%max(size(r2))
  mb=real(ifft( conj(fft(TTB)) .* fft(r2) ));
  m=(Na*ma+Nb*mb)/(Na+Nb);
  
  %% decallage de Na pour se replacer sur le bord (en Na+1)
  m=[m(end-Na+1:end),m(1:end-Na)]; % 9 (Na-1) derniers mis en tete
  ma=[ma(end-Na+1:end),ma(1:end-Na)];
  mb=[mb(end-Na+1:end),mb(1:end-Na)];
	
  %% calcul du rapport de vraisemblance
  glrt = log( (m.^N) ./ ((ma.^Na) .* (mb.^Nb)) ) ;
  glrt(blkr2==1) = 0 ;
  
  %% calcul du test 
  test = (glrt>seuil) .* (mb<ma) ; %%(que les ralentissements)
  test_acc = (glrt>seuil) .* (mb>ma) ; %% accélération??

  %% gestion du blink
  %% le Test n'est pas valide a l'instant t si 
  %% la particule est blinkee au moins une fois entre t et t+Nb-1
  %% On elargie les zones invalides seulement par rapport a la
  %% taille de la fenetre Nb. En effet, un blink dans la fenetre B
  %% conduira a une baisse de la moyenne et donc a des fausses alarme.
  %% La presence de blink dans la fenetre A fait baisser la Pd
  %% D'ou la zone de validite suivante :
  blink = conv([blkr2,zeros(1,Nb+1)],ones(1,Nb+1))>0;
  blink = blink(Nb+1:end-Nb-1);
  test = test .* ~blink ;


end%function
