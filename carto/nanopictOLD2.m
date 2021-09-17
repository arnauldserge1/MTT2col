function Image_out = nanopict(data_in, ROI, check_roi_only, codage, dr_in)

% function Image_out = nanopict(data_in, ROI, check_roi_only, codage, dr_in)
%
% Génère une image résolue à l'échelle moléculaire, en exploitant les data
% de type MTT : position et précision de localisation, traduit en rayon
% des disques de localisation individuel (cf. principe PALM, STORM, GSDIM...)
% nota : l'info de reconnexion (éventuelle) des traces n'est pas exploitée ici.
%
% data_in = filename ('toto.stk') ou {tab_param, tab_var} 
% (donner tab_var si dr_in non donné)
%
% Image restreinte à une région d'intérêt : ROI = [x y w h] (coin et taille)
% (déf : taille 1/8e, centré)
% codage : 'alpha' (intensité) ou 't' (time)
%
% AS 15/12/8

if nargin<2, ROI = []; end
if nargin<3, check_roi_only = 0; end
if nargin<4, codage = 'alpha'; end
if nargin<4, dr_in = 0; end

dirname = 'output22'; chg_dir = 0;

%% data
if ischar(data_in) % nom du fichier
    if isempty(strfind(cd, dirname)), cd(dirname), chg_dir = 1; end
    
    disp('loading i...')
    pos_i = fread_data_spt([data_in '_tab_param.dat'], 3);
    disp('loading j...')
    pos_j = fread_data_spt([data_in '_tab_param.dat'], 4);
    if dr_in==0
        disp('loading di...')
        di = fread_data_spt([data_in '_tab_var.dat'], 3); % (nota: di=dj)
    end
    disp('loading blink...')
    blink = fread_data_spt([data_in '_tab_param.dat'], 8); % (7 => offset)
    
    disp(['loading ' codage ' (as intensity values)...'])
    if strcmp(codage,'alpha')
        val = fread_data_spt([data_in '_tab_param.dat'], 5);
    elseif strcmp(codage,'t')
        val = fread_data_spt([data_in '_tab_param.dat'], 2);
    end
    
    if chg_dir, cd .., end
else % tableau de data, directement
    if ~iscell(data_in),data_in = {data_in}; end %...
    pos_i = data_in{1}(2:7:end,:);
    pos_j = data_in{1}(3:7:end,:);
    if dr_in==0, di = data_in{2}(3:7:end,:); end
    blink = data_in{1}(7:7:end,:);
    
    if strcmp(codage,'alpha'), val = data_in{1}(4:7:end,:);
    elseif strcmp(codage,'t'), val = data_in{1}(1:7:end,:);
    else disp(['method ''' codage ''' not known...']), return
    end
    
    clear data_in;
end

%% full pict, with enlighted pixels at each peak
all_i = ceil(pos_i(blink>0)); % ceil => 1 à 512
all_j = ceil(pos_j(blink>0));

full_pict = uint16(zeros(max(all_i), max(all_j)));
ind = sub2ind(size(full_pict), all_i, all_j);
full_pict(ind) = val(blink>0);

figure('WindowStyle','docked')
% portrait ou paysage?
if max(all_i)<max(all_j), subplot(2,1,1), else subplot(1,2,1), end
% gouts & couleurs...
if strcmp(codage,'alpha'), cmap = 'hot';
elseif strcmp(codage,'t'), cmap = rainbow; cmap(1,:) = [0 0 0];
end

imagesc(full_pict)
axis image, hold on
colormap(cmap)

%% ROI def & scale bar
if isempty(ROI)
    ROI = ceil([max(all_j)/2, max(all_i)/2, max(all_j)/8, max(all_i)/8]); % x=j & y=i!!
end

rectangle('Position', ROI, 'EdgeColor', 'w')

pixelsize = 0.16;
stamp_calib_bar(10, '10 µm', pixelsize, max(all_i), max(all_j))
drawnow

if check_roi_only, disp('nice place for a ROI, isnt it??'), return, end

%% Inside & Non blinking only!
ok = (pos_i>ROI(2) & pos_i<ROI(2)+ROI(4) ...
    & pos_j>ROI(1) & pos_j<ROI(1)+ROI(3) & blink>0);

pos_i = pos_i(ok); % nota: linéarise la matrice
pos_j = pos_j(ok);
val = val(ok);

if dr_in==0
    di = di(ok); % dj = dj(ok);
    di_def = 0.7; % default value, non relevant..
    dr_avg = mean(di(di~=di_def));
else
    dr_avg=dr_in;
end
disp(['average position accuracy = ' num2str(dr_avg) ' pixels = ' ...
    num2str(dr_avg*pixelsize*1000) ' nm'])

%% kernel = PSF gaussienne (suréchantillonée, amplitude 1)
sampling = 10;
W = ROI(3)*sampling;
H = ROI(4)*sampling;
sig = dr_avg*sampling;
wn = round(2*3.5*sig);
if wn<min(W, H),
    g = sqrt(pi)*sig*gausswin2(sig, wn); % max(g)=1
    g = expand_w(g, H, W) ;
else
    g = sqrt(pi)*sig*gausswin2(sig, H, W); % max(g)=1
end

if strcmp(codage,'t'), g = g>exp(-1/2) ; end % gaussienne => disque binaire

%% adding positions (for all frames at once!)
Image_out = zeros(H, W); % image (sur)échantillonée
ROI_i = ceil((pos_i-ROI(2))*sampling); % positions, ramenées à ROI suréch
ROI_j = ceil((pos_j-ROI(1))*sampling);
ind = sub2ind(size(Image_out), ROI_i, ROI_j);
Image_out(ind) = val;

%% correlation by Gaussian & plot nanopict
Image_out = real(fftshift(ifft2(fft2(Image_out) .* fft2(g)))) ;
if strcmp(codage,'t'),  Image_out(Image_out>max(val)) = max(val); end

if max(all_i)<max(all_j), subplot(2,1,2), else subplot(1,2,2), end
imagesc(Image_out)
axis image, hold on
stamp_calib_bar(1*sampling, '1 µm', pixelsize, W, H)
colorbar('location','SouthOutside') %,'xlabel',codage

%%%

function stamp_calib_bar(calib_length, calib_text, pixelsize, W, H)
% draw a line in bottom left corner, with corresponding length & text,
% in current pict, of size W by H

calib = calib_length/pixelsize;
if min(W, H)>300, offset = min(W, H)/30; else offset = 10; end% offset from pict border
line([offset offset+calib], [H-offset/2 H-offset/2], 'LineWidth', 2, 'color', 'w')
text(offset+calib/2, H-offset, calib_text, 'color', 'w', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')

%%%

function cmap = rainbow(n)
% rainbow color map, from magenta to red (cf. HSV, but reversed and stoped at red)

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