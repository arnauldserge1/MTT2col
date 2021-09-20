%function [liste_param, reste]
%   = estim_param_part_GN(im, wn, liste_info_param, r0, bornes_ijr)
%
% EN/ sub-pixel estimate of the peak position by Gauss Newton regression
% liste_info_param is a line of the matrix liste_detect
% wn odd
%
%
% FR/ estimation sub-pixel de la position du pic par Gauss Newton
% liste_info_param est une ligne de la matrice liste_detect
% wn impaire


%%% version 11 07 05

function [liste_param, reste] = estim_param_part_GN(im, wn, liste_info_param, r0, bornes_ijr)

scale_factor = r0/1.1 ; % bornes prop. à r0

if (nargin < 5)
     bornes_ijr(1) = -1.5 ;
     bornes_ijr(2) = 1.5 ;
     bornes_ijr(3) = -1.5 ;
     bornes_ijr(4) = 1.5 ;
     bornes_ijr(5) = 0.3 ;
     bornes_ijr(6) = 3.0 ; 
     bornes_ijr = bornes_ijr * scale_factor ; % AS 16/9/8
end%if

Pi = liste_info_param(2) ;
Pj = liste_info_param(3) ;
ii = (1:wn)+Pi-floor(wn/2) ;
jj = (1:wn)+Pj-floor(wn/2) ;
im_part = im(ii, jj) ;


r = r0 ;
i = 0.0 ;
j = 0.0 ; 
dr = scale_factor ;
di = scale_factor ;
dj = scale_factor ;
fin = 0.01 ;
sig2 = inf ;
cpt = 0 ;
test = 1 ;
ITER_MAX = 50 ;
while (test)
   %%[r, i, j, dr, di, dj, alpha, sig2] = deplt_GN_estimation (r, i, j, im_part) ;
   [r, i, j, dr, di, dj, alpha, sig2, reste] = deplt_GN_estimation (r, i, j, im_part, sig2, dr, di, dj) ;
   cpt = cpt + 1 ;
   test = max([abs(di), abs(dj), abs(dr)]) > fin ;
   if (cpt > ITER_MAX) 
     test = 0 ;
   end%if

   %% on stop si l_on sort des bornes  
   result_ok = ~((i < bornes_ijr(1)) || (i > bornes_ijr(2)) || ...
		 (j < bornes_ijr(3)) || (j > bornes_ijr(4)) || ...
		 (r < bornes_ijr(5)) || (r > bornes_ijr(6)) ) ;
   test = test & result_ok ;

end%while


% liste_info_param = [num, i, j, alpha, sig^2]
% liste_param = [num, i, j, alpha, sig^2, rayon, ok]

liste_param = [liste_info_param(1), ...
               Pi+i , ...
               Pj+j , ...
               alpha , ...
               sig2 , ...
	           r , ...		      
               result_ok ];

global SHOW_PLOT_FIT
if SHOW_PLOT_FIT, plot_fit(im_part, reste, wn, i, j, r, alpha, sig2, result_ok), end

end%fonction


function plot_fit(im_part, reste, wn, i, j, r, alpha, sig2, result_ok)
% AS 22/6/10

colormap(gray), clf

subplot(2,3,1)
imagesc(im_part), axis equal off
hold on
if result_ok, plot(j+wn/2+0.5, i+wn/2+0.5, 'g+')
else plot(j+wn/2+0.5, i+wn/2+0.5, 'r+'), end
title('input')

subplot(2,3,4)
surf(im_part), axis tight

subplot(2,3,2)
imagesc(alpha*gausswin2(r, wn, wn, i, j)), axis equal off
title('fit')

subplot(2,3,5)
surf(alpha*gausswin2(r, wn, wn, i, j)), axis tight
title(sprintf('i=%.2g, j=%.2g, r=%.3g', i, j, r))

subplot(2,3,3)
imagesc(reste), axis equal off
title('reste')

subplot(2,3,6)
surf(reste), axis tight
title(sprintf('sig^2=%.2g, fit ok=%i', sig2, result_ok))

pause%(1)

end%fonction

%%%