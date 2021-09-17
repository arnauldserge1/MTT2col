function [pict, h] = DIC_image(filename, DIC_name, triD, zmin, zmax, viewDV, verbose, new_fig)

% function pict = DIC_image(filename, DIC_name, triD, zmin, zmax, viewDV, verbose, new_fig)
% plot DIC (or phase, or first fluo image of the stack, if no pict in
% 'dic' folder), in 2D or 3D axe

if nargin < 8, new_fig = 1; end
if nargin < 7, verbose = 1; end
if nargin < 6, viewDV = ''; end
if nargin < 3, triD = 0; zmin = 0; zmax = 0; end % donc 2D!
if nargin < 2, DIC_name = dicname(filename, verbose); end

% pxl_size = 0.16; % 1 pxl = 0.16 um

if new_fig, figure('WindowStyle', 'docked'), end

if ~isempty(DIC_name)
    if ischar(DIC_name), DIC = imread(DIC_name); end   %DIC = uint16(tiffread...
    sat = 0.002 ; % saturation 0.2% min-max du contraste
    DIC_sat = imadjust(DIC, stretchlim(DIC, [sat 1-sat]), [0 1]);
    H = fspecial('average');
    pict = imfilter(DIC_sat, H, 'replicate');
else
    pict = double(imread(filename, 1)); % 1e image du stack...
    pict = max(pict(:)) - pict; % invert
end

if (length(DIC_name) > 4) && ~strcmp(DIC_name(1:5),'trans') % DIC image double, take half. trans image already single
    DIC_width = size(pict, 2);
    if contains(viewDV, 'left'), pict = pict(:, 1:DIC_width/2); % 1/2 image de gauche !
    elseif contains(viewDV, 'right'), pict = pict(:, DIC_width/2+1:end); end % strcmp(viewDV, '_right')...
end

if ~triD % 2D donc...
    h = imagesc(pict); % brighten(+0.3)??
else
%% *** met l'image de la cellule "au plancher" ***
% % % % % % % %     pict(:, :, 2) = 0; % juste pour creer une 3D...
    [ny, nx] = size(pict);
    [x, y, z] = meshgrid((1:nx)*pxl_size, (1:ny)*pxl_size, [zmin zmax]);
    v = zeros(ny, nx, 2);
    v(:, :, 1) = pict;
    h = slice(x, y, z, v, [], [], zmin); % coupe a zmin
    set(h, 'EdgeColor', 'none')
% % %     zlabel('frames')
    axis([0 nx*pxl_size 0 ny*pxl_size zmin zmax])
%     view(65, 55) % azimut, hauteur
end

axis image ij off %([1 size(DIC_sat, 1) 1 size(DIC_sat, 2)]) % else axis equal
title({cd; filename}, 'interpreter', 'none')
colormap('gray')
hold on
if new_fig, drawnow, end

if nargout < 1, pict = []; end
%%%