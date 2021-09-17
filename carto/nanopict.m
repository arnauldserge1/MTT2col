function nanopict(file, ROI, check_roi_only, codage, dr_in, pixelsize, Nimg_max, Nsub_img)

% function nanopict(file, ROI, check_roi_only, codage, dr_in, pixelsize, Nimg_max, Nsub_img)
%
% Génère une image résolue 'à l'échelle moléculaire', en exploitant les data
% de type MTT : position et précision de localisation, traduit en rayon
% des disques de localisation individuel (cf. principe PALM, STORM, GSDIM...)
% nota : l'info de reconnexion (éventuelle) des traces n'est pas exploitée ici.
%
% dr_in def = 0.1 pxl (selon SNR et/ou Dmin)
%
% Image restreinte à une région d'intérêt : ROI = [x y w h] (coin et taille)
% (déf : taille 1/8e, centré)
% utiliser check_roi_only = 1 pour juste tester la position de la ROI
%
% codage : 'alpha' (intensité, défaut) ou 't' (time) (ou speed dr/dt??)
%
% AS 15/12/8

MTT_param %global N_PARAMif isempty(N_PARAM), , end

if nargin==0
%     files = dir('*.stk');%     if isempty(files), files = dir('*.tif'); end%     file = files(1).name;
    [file, file_path] = uigetfile({'*.tif'; '*.stk'});
    cd(file_path)
    user_input = 1;
else
    user_input = 0;
end
if nargin<2, ROI = []; end
if nargin<3, check_roi_only = 0; end
if nargin<4, codage = 'alpha'; end
if user_input
    answer = inputdlg('Codage? (alpha or t)', '', 1, {codage}); % length(imfinfo(file))?
    codage = answer{1};
end
if nargin<5, dr_in = 0.1; end % dr_in def = 0.1 pxl (selon SNR et/ou Dmin)
if nargin<6, magnification = 100*1.5; pixelsize = 16/magnification; end
if nargin<7, Nimg_max = length(imfinfo(file)); end 
if nargin<8, Nsub_img = 200; end % charge data par partie
% if user_input
%     answer = inputdlg('Nombre d''images à traiter? (all => inf)', '', 1, {num2str(Nsub_img)}); % length(imfinfo(file))?
%     Nimg_max = str2double(answer{1});
%     if (Nimg_max == inf), Nimg_max = length(imfinfo(file)); end
% end

n_loop = ceil(Nimg_max/Nsub_img);

image1 = imread(file);
[H_pict_ini, W_pict_ini] = size(image1);
if isempty(ROI)
    ROI = ceil([H_pict_ini*7/16, W_pict_ini*7/16, H_pict_ini/8, W_pict_ini/8]); % ceil([max(all_j)/2, max(all_i)/2, max(all_j)/8, max(all_i)/8]); % x=j & y=i!!
end


%% Max image + ROI
fprintf('computing image of max (using first 10%% frames) ....')
Nimg_for_max = Nimg_max/10;
img_max = image_max(file, 1:Nimg_for_max, 1); % figure('WindowStyle','docked') % imagesc(img_max)

% gouts & couleurs...
if strcmp(codage,'alpha'), cmap = 'hot';
elseif strcmp(codage,'t'), cmap = rainbow; cmap(1,:) = [0 0 0];
end
colormap(cmap)
axis image

if user_input % isempty(ROI) && nargin>0
    fprintf('\n  select ROI please...')
    ROI = getrect(gca);
    ROI = round(ROI);
    fprintf('x0=%g, y0=%g, width=%g, height=%g\n', ROI)
end

subplot(3,4,1)
imagesc(img_max) % full pict max
axis image, hold on
rectangle('Position', ROI, 'EdgeColor', 'w')

subplot(3,4,5)
imagesc(img_max) % ROI pict max
axis image
axis([ROI(1) ROI(1)+ROI(3) ROI(2) ROI(2)+ROI(4)])
drawnow

%% kernel = PSF gaussienne (suréchantillonée, amplitude 1)
sampling = 10;
W_nanopict = ROI(3)*sampling;
H_nanopict = ROI(4)*sampling;
sig = dr_in*sampling;
wn = round(2*3.5*sig);

if wn<min(W_nanopict, H_nanopict),
    g = sqrt(pi)*sig*gausswin2(sig, wn); % max(g)=1
    g = expand_w(g, H_nanopict, W_nanopict) ;
else
    g = sqrt(pi)*sig*gausswin2(sig, H_nanopict, W_nanopict); % max(g)=1
end

if strcmp(codage,'t'), g = g>exp(-1/2) ; end % gaussienne => disque binaire

%% loop to load data by part
Image_out = [];

for i_loop = 1:n_loop
    t_first = (i_loop-1)*Nsub_img+1;
    t_last = i_loop*Nsub_img;
    %% data
    data = fread_params_timewindow(file, 1, t_first, t_last); % print_out = 1;
    if isempty(data), return, end
    pos_i = data(PARAM_I-1:N_PARAM:end,:);
    pos_j = data(PARAM_J-1:N_PARAM:end,:);
    blink = data(PARAM_BLINK-1:N_PARAM:end,:); %     if dr_in==0, di = data_in{2}(3:N_PARAM:end,:); end
    
    if strcmp(codage,'alpha'), val = data(PARAM_ALPHA-1:N_PARAM:end,:);
    elseif strcmp(codage,'t'), val = data(PARAM_T-1:N_PARAM:end,:);
    else disp(['method ' codage ' not known...']), return
    end
    
    %% full pict, with enlighted pixels at each peak
    all_i = ceil(pos_i(blink>0)); % ceil => 1 à 512
    all_j = ceil(pos_j(blink>0));
    
    full_pict_i = uint16(zeros(H_pict_ini, W_pict_ini));%max(all_i), max(all_j)));
    ind = sub2ind(size(full_pict_i), all_i, all_j);
    full_pict_i(ind) = val(blink>0);
    
    if i_loop == 1
        full_pict = full_pict_i;
    else
        full_pict = full_pict + full_pict_i;
    end
    % portrait ou paysage?    %     if W_pict_ini<H_pict_ini, subplot(2,1,1), else subplot(1,2,1), end
    subplot(3,4,9) % full nanopict
    imagesc(full_pict)
    axis image, hold on
    colormap(cmap)
    
    %% ROI def & scale bar
    rectangle('Position', ROI, 'EdgeColor', 'w')
    scale_bar(10, '10 µm', pixelsize, max(all_i), max(all_j))%     drawnow
    
    if check_roi_only
        disp('nice place for a ROI, isnt it??')
        fprintf('  ROI:%i %i %i %i',ROI)
        return
    end
    
    %% Inside & Non blinking only!
    ok = (pos_i>ROI(2) & pos_i<ROI(2)+ROI(4) ...
        & pos_j>ROI(1) & pos_j<ROI(1)+ROI(3) & blink>0);
    
    pos_i = pos_i(ok); % nota: linéarise la matrice
    pos_j = pos_j(ok);
    val = val(ok);
    
    %% adding positions (for all frames at once!)
    Image_i = zeros(H_nanopict, W_nanopict); % image (sur)échantillonée
    ROI_i = ceil((pos_i-ROI(2))*sampling); % positions, ramenées à ROI suréch => donc précision 0.1 pixel (1 pixel dans image suréchantillonée)
    ROI_j = ceil((pos_j-ROI(1))*sampling);
    ind = sub2ind(size(Image_i), ROI_i, ROI_j);
    Image_i(ind) = val;
    
    %% correlation by Gaussian & plot nanopict
    Image_i = real(fftshift(ifft2(fft2(Image_i) .* fft2(g)))) ;
    if strcmp(codage,'t'),  Image_i(Image_i>max(val)) = max(val); end
    
    if i_loop == 1
        Image_out = Image_i;
    else
        Image_out = Image_out + Image_i;
    end
    
    %     if W_pict_ini<H_pict_ini, subplot(2,1,2), else subplot(1,2,2), end
    subplot(3,4,[2 3 4 6 7 8 10 11 12])
    %     imagesc(log((Image_out-min(Image_out(:)))))
    Im = (Image_out-min(Image_out(:)))/(max(Image_out(:))-min(Image_out(:)));
    gamma1 = 0.5;
    imagesc(Im.^gamma1) % ROI nanopict
    axis image, hold on
    scale_bar(1*sampling, '1 µm', pixelsize, W_nanopict, H_nanopict)
    title(file, 'interpreter', 'none')
    if i_loop == 1, colorbar('location','SouthOutside'), end %,'xlabel',codage
    if (i_loop == 1), c = caxis; caxis([0 c(2)/2]); end % sat...
    if user_input && (i_loop == 1 || i_loop == n_loop)
        ok = 0;
        while ~ok
            answer = inputdlg('change gamma? (keep current value to continue)', '', 1, {num2str(gamma1)});
            gamma2 = str2double(answer{1});
            imagesc(Im.^gamma2)
%             answer = inputdlg('stop (1) ou encore (0)??', '', 1, {'1'}); ok = str2double(answer{1});
            if (gamma2 == gamma1), ok = 1; else gamma1 = gamma2 ; end
        end
    end
end % i_loop

%%%

function scale_bar(calib_length, calib_text, pixelsize, W, H)
% draw a line in bottom left corner, with corresponding length & text,
% in current pict, of size W by H

calib = calib_length/pixelsize;
if min(W, H)>300, offset = min(W, H)/30; else offset = 10; end% offset from pict border
line([offset offset+calib], [H-offset/2 H-offset/2], 'LineWidth', 2, 'color', 'w')
text(offset+calib/2, H-offset, calib_text, 'color', 'w', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')

%%%

function cmap = rainbow(n)
% rainbow color map, from magenta to red (cf. HSV, but reversed and stopped at red)

if nargin==0, n = 64; end

p = ceil(n/5);
zz = zeros(1,p);
oo = ones(1,p);
ramp = 0:1/(p-1):1;

red = [1-ramp zz zz ramp oo]';
green = [zz ramp oo oo 1-ramp]';
blue = [oo oo 1-ramp zz zz ]';

cmap = [red green blue];
cmap = cmap(1:n,:);
%%%