function lest2tab_param(traj, t, t_red, lest, part, mode)

% function [tab_param, tab_var] = lest2tab_param(tab_param, traj, t_red, lest, part, mode)
% Reattribute estimated parameters from lest (fit results) to tab_param,
% for current trajectory & particle
% using global variables for each parameter position in the output table
% (see MTT_param)
% 7/12/2009

global tab_param ;
global Nb_STK ;

global N_PARAM PARAM_T PARAM_I PARAM_J PARAM_ALPHA PARAM_RADIUS_I ...
    PARAM_OFFSET PARAM_BLINK PARAM_SIG2 % PARAM_K PARAM_RADIUS_J
global SAVED_PARAMS

% PARAM_T = 2 ; % = frame number
% PARAM_I = 3 ; % position along lines, sub-pixel
% PARAM_J = 4 ; % position along columns, sub-pixel
% PARAM_ALPHA = 5 ; % signal, see definition of the Gaussian
% PARAM_RADIUS_I = 6 ; % Gaussian sd along i axis
% PARAM_OFFSET = 7 ;
% PARAM_BLINK = 8 ;
%
% %%% new params, from version MTT2.3 %%%
% PARAM_SIG2 = 9 ; % power of noise (in tab_param, no longer in tab_var)
% PARAM_K = 10 ; % or z, along optical axis, if evaluated %%% placed in position 9 for compatibility!!
% PARAM_RADIUS_J = 11 ; % ri & rj different when using astigmatic, cylindrical lens to evaluate k

for i_par = 1:N_PARAM
    switch SAVED_PARAMS(i_par)
        case  PARAM_T
            tab_param(traj, N_PARAM*(t_red)+PARAM_T) = t ; %t
        case PARAM_I
            tab_param(traj, N_PARAM*(t_red)+PARAM_I) = lest(part, 2) ; %i
        case PARAM_J
            tab_param(traj, N_PARAM*(t_red)+PARAM_J) = lest(part, 3) ; %j
        case PARAM_ALPHA
            tab_param(traj, N_PARAM*(t_red)+PARAM_ALPHA) = lest(part, 4) ; %alpha
        case PARAM_RADIUS_I
            tab_param(traj, N_PARAM*(t_red)+PARAM_RADIUS_I) = lest(part, 6) ; %r
        case PARAM_OFFSET
            tab_param(traj, N_PARAM*(t_red)+PARAM_OFFSET) = lest(part, 8) ; %offset (m0)
        case PARAM_BLINK
            if (strcmp(mode, 'new')) || (strcmp(mode, 'new_traj')) % but not for 'update'
                tab_param(traj, N_PARAM*(t_red)+PARAM_BLINK) = Nb_STK ;                
            end
        case PARAM_SIG2
            tab_param(traj, N_PARAM*(t_red)+PARAM_SIG2) = lest(part, 5) ; %sig2_b
%             tab_var(traj, N_PARAM*(t_red)+7) = lest(part, 5) ; %sig2_b
    end
end