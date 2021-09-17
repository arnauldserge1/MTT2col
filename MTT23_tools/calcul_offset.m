%% calcul de la moyenne locale m0 (plancher/offset de la gaussienne)    AS 29/3/7
% function offset = calcul_offset(lest, input_deflt, wn)
% AS 29/3/7

function offset = calcul_offset(lest, input_deflt, wn)

demi_wn = ceil(wn/2) ; % => 4!
range = 1-demi_wn:demi_wn-1 ;

nb_part = size(lest,1) ;
offset = zeros(nb_part,1);

for part = 1:nb_part % nb_valid
    irange = round(lest(part,2))+range ;
    jrange = round(lest(part,3))+range ;
    irange = irange(irange>0 & irange<=size(input_deflt,1)) ; % écarte éventuelles val foireuses
    jrange = jrange(jrange>0 & jrange<=size(input_deflt,2)) ;
    fond = input_deflt(irange, jrange) ; % fond local: input_deflt = im(wn x wn) - gauss(i,j)
    if ~isempty(fond)
        offset(part) = mean(fond(:)) ; 
    end
end

%tab_param(traj, 7*(t_red)+7) = mean(fond(:)) ; %offset : moyenne locale (bruit, sans signal, après deflation)