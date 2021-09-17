function z = cyl_focus(sig_x, sig_y, sig_0, z_r, gamma_z)
% function z = cyl_focus(sig_x, sig_y, sig_0, z_r, gamma_z)
% Compute axial/focal position z according to sig_x and sig_y as eval
% by using a cylindrical lens inducing an axial astigmatism gamma_z, 
% for diffraction limited peaks of width sig_0 (1.22 lambda / 2 NA),
% with a focal depth z_r (3 sig_0 by def.)
% all parameters in pixel (160 nm)
% In our setup, cylindric axis along y, thus z>0 <=> sig_x>sig_y

% cf. Holtzer, Meckel & Schmidt, APL 2007
% AS 25/11/8

if nargin<3, sig_0 = 1.22*605/(2*1.4)/160 ; end 
if nargin<4, z_r = 3*sig_0 ; end % ~500 nm?? ~1.5 lambda??
if nargin<5, gamma_z = 1.25 ; end % ~200 nm??

% equiv radius
% sig_r = sqrt(sig_x*sig_y) ;

% ellipticity
% ellipt = sqrt(sig_y/sig_x) ;

if sig_x > sig_y  % ellipt<1 % z<0, peak is below focus
    z = z_r * sqrt(sig_x^2/sig_0^2 - 1) + gamma_z ;
else % z>0, peak is above
    z = - z_r * sqrt(sig_y^2/sig_0^2 - 1) + gamma_z ;
end
