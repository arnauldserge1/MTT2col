function [left_peaks right_peaks] = selectDVpeaks(tab_param, Nj, bord)

global N_PARAM
if isempty(N_PARAM), MTTparams_def; end

if nargin<3
    MTT_param % => wn = 7
    bord = ceil(wn/2) ;
end

%%% centrageDV???
tab_j = tab_param(3:N_PARAM:end,:) ;

left_peaks = tab_param(tab_j < Nj/2 - bord);
right_peaks = tab_param(tab_j > Nj/2 + bord);