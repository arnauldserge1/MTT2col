%% MTT: Multiple Target Tracing algorithms
%% PARAMETERS FILE
%% 
%%
%%
%% This programme is compatible with Matlab and Octave software
%%
%%
%% ALGORITHM AUTHORS:
%% 
%% Copyright A. SERGE(1,2,3), N. BERTAUX(4,5,6)
%%
%%
%% AFFILIATIONS :
%%
%% (1) CIML, University of Marseille, F13009, FRANCE
%%
%% (2) INSERM, UMR 631, Marseille, F13009, FRANCE
%%
%% (3) CNRS, UMR 6102, Marseille, F13009, FRANCE
%%
%% (4) Fresnel Institut - PhyTI Team - MARSEILLE - F13397 - FRANCE
%%
%% (5) CNRS, UMR 6133, Marseille, F13397, FRANCE
%% 
%% (6) Ecole Centrale de Marseille - France
%%
%%
%% last modifications 03/12/07



%% =============================================
%% add in path utils_SPT/ subroutines repertoire
%% =============================================


%% =====================================================
%% =====================================================
%% FILE PARAMETERS
%% =====================================================
%% =====================================================

%% =====================================================
%% positions of parameters
%% =====================================================

include_global
% global N_PARAM PARAM_T PARAM_I PARAM_J PARAM_ALPHA PARAM_RADIUS_I ...
%     PARAM_OFFSET PARAM_BLINK PARAM_SIG2 PARAM_K PARAM_RADIUS_J
% global SAVED_PARAMS SAVE_VAR_MOY

% % % % % N_PARAM = 8 (up to 11) ; % total number of output parameters

%%% PARAM_PARTICLE NUMBER = 1 ; not save in tab_param.dat
%%% thus param numbers start at 2

PARAM_T = 2 ; % = frame number
PARAM_I = 3 ; % position along lines, sub-pixel
PARAM_J = 4 ; % position along columns, sub-pixel
PARAM_ALPHA = 5 ; % signal, see definition of the Gaussian
PARAM_RADIUS_I = 6 ; % Gaussian sd along i axis
PARAM_OFFSET = 7 ;
PARAM_BLINK = 8 ;

%%% new params, from version MTT2.3 %%%
PARAM_SIG2 = 9 ; % power of noise (in tab_param, no longer in tab_var)
PARAM_K = 10 ; % or z, along optical axis, if evaluated %%% placed in position 9 for compatibility!!
PARAM_RADIUS_J = 11 ; % ri & rj different when using astigmatic, cylindrical lens to evaluate k

% PARAM_LAMBDA = 12 ; % for multicolor measure (i.e. 2 populations, red & green) (default: 0)
%%% add if required...

SAVED_PARAMS = [PARAM_T PARAM_I PARAM_J PARAM_ALPHA PARAM_RADIUS_I ...
    PARAM_OFFSET PARAM_BLINK PARAM_SIG2];% PARAM_K PARAM_RADIUS_J];
% SAVED_PARAMS_v22 = [PARAM_T PARAM_I PARAM_J PARAM_ALPHA PARAM_RADIUS_I ...
%     PARAM_OFFSET PARAM_BLINK];
N_PARAM = length(SAVED_PARAMS); % => 8

SAVE_VAR_MOY = 0 ; % save or not tab_var and tab_moy (default from version MTT2.3: save only tab_param.dat)


%% version program
%% comment corresponding ligne
%%VERSION = 'OCTAVE' ;
VERSION = 'MATLAB' ;

MTT_version = '2.3';

OPTIM_R = 0; % 29 2 2013

%% =============================================
%% DISPLAY PARAMETERS
%% =============================================

%% output as images, format ppm
AFFICHAGE = 0 ; 

%% output file extension (ppm, jpg, gif, etc.)
FORMAT_IM = 'png' ;

%% Display particles number
AFF_NUM_TRAJ = 1 ; %% 0/1 

SHOW = 0 ; %% 0/1 inline dislay output images

%% image size output option
%%imwrite_option = ''; 
imwrite_option = '-resize 320x240' ;

%% display of the tracking of a limited number of particles
%%
%% all particules
liste_part = 0 ; 
%%
%% only particle number (>0)
%% liste_part = [1 3 8 9 12 13 15 16 17] ;
%%
%% all without particle number (<0)
%% liste_part = -[18 19 5 20 21 34] ;



%% ===============================================
%% DATA INPUT
%% ===============================================

%% directory
repertoire = '' ;

%% filename
%%%stack = 'EGFR-Qd605.stk' ; AS 12/12/7
filename = '*.tif';

%% number of images
%%%Nb_STK = 300 ; AS 12/12/7
Nb_STK = [];

%% limitation of the zone of interest
CROP = 0 ; %% boolean 0/1
IRANGE = 180 + (1:120) ; 
JRANGE = 15  + (1:80) ; 


%% ===============================================
%% OUTPUT
%% ===============================================
output_dir = 'output23' ; %'../data/output23' ;
full_output_dir = ['./', output_dir] ;



%% ===============================================
%% Tracking parameters
%% ===============================================

%% ====================
%% detection parameters
%% ====================

%% Pre-Detection threshold
seuil_premiere_detec = 24 ; %% 10^-6

%% Final detection threshold
seuil_detec_1vue = 28 ;  %% 10^-7

%% Gaussian radius in pixel (r0)
r0 = 1.1 ; 

%% Size Windows en pixel (Ws)
wn = floor(7*r0) ; % 7 % AS 17/9/8

%% =======================
%% Reconnection parameters
%% =======================

%% temporal sliding window (Wt)
T = 5 ; %% taille fenetre glissante temporelle

%% number of deflation loop
nb_defl = 1 ;

%% disappearance prob of blinking (tau_off)
T_off = -15 ;

%% Maximum diffusion coef Dmax
Din = 0.08 ; % for JamB at 100ms, 16/4/2013
sig_free = 2*sqrt(Din) ; %% in pixel

%% Reference diametre for research set particles/trajectories
%% diameter = Boule_free*sig_free 
Boule_free = 3.5 ;
 
%% Limitation of combinations number 
%% Nb_combi define maximum number of particle/trajectories
%% becarefull the complexity change as fact(Nb_combi) (4!=24 , 5!=120)
Nb_combi = 4 ;

%% validation of pre-detected particles
%% & new ones
seuil_alpha = 0; % 4000 ; %% environ 20 dB // ramené à 0 le 25/4/8 pour pouvoir traiter data à 1ms (env. 20dB justement)
% 4000 reduisait faux departs sur particules tres faibles

%% weight of likelihood of alpha
%% between uniform and gaussian law
Poids_melange_aplha = 0.5 ;

%% weight of likelihood between
%% maximum and local diffusion
Poids_melange_diff = 0.9 ;

%% min. trace length kept for analysis
Tmin = 5;

%% min. diffusion coef kept for analysis (lower assumed as immobile, discarded)
Dmin = 0.001;


%% =====================================================
%% =====================================================
%% END PARAMETERS FILE
%% =====================================================
%% =====================================================


