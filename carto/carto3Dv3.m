function S = carto3Dv3(trc, filename, do_subplot)
% function S = carto3Dv3(trc, filename, do_subplot)
% calcule la carte S(xyz) des niveaux de confinement avec z=Log10(L)
% avec lissage gaussien
% (+ lignes de niveau et dot plot si do_subplot)
% AS, 28/2/2007

if nargin<3, do_subplot = 0; end

%% ** ROI **
%     if ~ischar(filename), ROI = [1 128 1 128]; % simul nico
%     elseif timing==36, ROI = [1 512 1 512];
%     else ROI = [floor(min(x)) ceil(max(x)) floor(min(y)) ceil(max(y))];
%     end
 [lg,ht] = stk_size(filename); ROI = [1 lg 1 ht];

%% *** conf ***
timing = eval_timing(filename);
Lc = probaconf(trc, 0, timing); % 'Nico'
%%%%%%%% REDONDANT???????
% Lmax = 370;
% Lc(Lc>Lmax) = 0; % on supprime les vals trop élevées (immobiles, non pertinentes)
trc(:,5) = log10(Lc); % PASSAGE AU LOG: z=Log(Lconf)
trc(trc(:,5)<0,:) = [];% supprime points nuls ou de log neg
if isempty(trc), S = zeros(lg,ht); disp('gloups, y''a rien..'), return, end

%% moyenne temporelle sur troncons de trajs, longs de wT (tlast-tfirst+1)
% if t_last == inf % une seule image, dc pas de moyenne ??
%     x = trc(:,3);%     y = trc(:,4);%     z = trc(:,5);
% else
    Ntrc = trc(end,1);
    x = zeros(Ntrc,1); y = zeros(Ntrc,1);z = zeros(Ntrc,1);
    for i=1:Ntrc
        indi = find(trc(:,1)==i);
        if ~isempty(indi)
            indm = round(mean(indi)); % point médian du tronçon
            x(i) = trc(indm,3);
            y(i) = trc(indm,4);
            zi = trc(indi,5); % vals pour trci, à moyenner
            z(i) = mean(zi);
        end
    end
% end
clear trc
%z = max(z,0);
x(z==0) = []; y(z==0) = []; z(z==0) = [];

%% ** sélection pics dans masque **
DIC_name = dicname(filename);
maskdir = ['DIC\masques en position basse\Masque IN\MIB.' DIC_name(5:end)];% remove 'dic\'
if ~isempty(dir(maskdir))
    maskIN = imread(maskdir);
    Rdil = 5;
    maskINdil = imdilate(maskIN,strel('disk',Rdil)); % on dilate d'une marge de Rdil pxl

    if do_subplot, figure('WindowStyle','docked'), plot(y,x,'b.'), hold on, end % total
    [y x ind] = select_pk_in_mask(y, x, maskINdil); % NB y x == i j
    z = z(ind>0);
    if do_subplot, plot(y,x,'r.'), end % in mask

%% contour
    Rdil2 = 10;
    maskINdil2 = imdilate(maskINdil,strel('disk',Rdil2)); % on dilate d'une marge de Rdil pxl
    maskCONTdil = bwperim(maskINdil2);% imdilate(maskCONT,strel('disk',Rdil));
    [im jm] = find(maskCONTdil);
    xm = jm(1:15:end); ym = im(1:15:end); % on échantillonne à 1/15 => 200 pts IMPAIR DE PREF!!
    x(end+1:end+size(xm)) = xm;
    y(end+1:end+size(xm)) = ym;
    z(end+1:end+size(xm)) = zeros(size(xm)); % contour = min = 0 (log(1)!!), comme bords
    if do_subplot, plot(ym,xm,'k.'), end
else
    disp('no mask found..') 
end

%% ajout des bords
xb = ROI(1); xe = ROI(2); yb = ROI(3); ye = ROI(4);
xx = (xb:20:xe); yy = (yb:20:ye); N = 2*length(xx)+2*length(yy);
x(end+1:end+N) = [xx, ones(size(yy))*xb, xx, ones(size(yy))*xe]';
y(end+1:end+N) = [ones(size(xx))*yb, yy, ones(size(xx))*ye, yy]';
z(end+1:end+N) = zeros(N,1); % bords = min, plancher

%% ****************** compute surf ***********************
%if t_last == inf, delta=1/3; else delta=1; end
x1 = ROI(1):ROI(2); y1 = ROI(3):ROI(4); %:delta:
[Xi,Yi] = meshgrid(x1,y1);
S = griddata(x,y,z,Xi,Yi,'cubic');%'v4');%,{'Qt','Qbb','Qc','Qz'});%
S = max(S,0); % supprime NaN!
%%%