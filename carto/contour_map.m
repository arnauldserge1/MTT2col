function S = contour_map(x,y,z,filename)
% calcule la carte S(xyz) des niveaux de confinement
% avec lissage gaussien
% (+ lignes de niveau et dot plot si do_subplot)
% AS, 28/2/2007

if nargin<4, filename = ''; end

%% ** ROI **
if ~isempty(filename)
    [lg,ht] = stk_size(filename);
    ROI = [1 lg 1 ht];
else
    ROI = [floor(min(x)) ceil(max(x)) floor(min(y)) ceil(max(y))];
end

%% ** sélection pics dans masque **
if ~isempty(filename)
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
end

%% ****************** compute surf ***********************
x1 = ROI(1):ROI(2); y1 = ROI(3):ROI(4);
[Xi,Yi] = meshgrid(x1,y1);
S = griddata(x,y,z,Xi,Yi,'cubic');%'v4');%,{'Qt','Qbb','Qc','Qz'});%
S = max(S,0); % supprime NaN!
%%%