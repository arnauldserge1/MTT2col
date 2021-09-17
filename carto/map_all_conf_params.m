function map_all_conf_params(filename, tab_param)
% maps duration 'T'& size 'R' for 'free' & 'conf' events,
% as defined by MTT_proba_conf
% by interpolation through (ijk) points, with
% i,j barycentre of each event and k associated par, T or R


if nargin==0, files = dir('*.stk'); filename = files(f).name; end
if nargin<2, tab_param = fread_all_params(filename); end

[Lconf Tc Tf R2c R2f trc_c trc_f] = MTT_probaconf(tab_param, 'rel');

% parameters = ;type =
AFMmap = load_AFMmap;
n_levels = size(AFMmap,1);
n_half = round(n_levels/2);
conf_map = AFMmap(1:n_half,:);%#ok
free_map = AFMmap(n_half+1:end,:);
free_map(1,:) = [0 0 .3];%#ok % plancher bleu marine

% i_plot = 0;
for i_par = 1:2
    switch i_par
        case 1
            par = 'T';
        case 2
            par = 'R';
    end
    
    for i_type = 1:2
        switch i_type
            case 1
                trc_events = trc_f;
                event_type = 'free';
            case 2
                trc_events = trc_c;
                event_type = 'conf';
        end
        
        S = map_conf_params(trc_events, par, filename);
        
%         i_plot = i_plot+1;%         subplot(2,2,i_plot)
        % levels_val = [0 log10(zmax)/256:0.4:log10(zmax)];
        figure('windowstyle','docked')
        contourf(log10(S)) % ,levels_val
        title([filename ' ' par, event_type])
        axis ij off image
        
%         switch event_type
%             case 'conf'%                 colormap(conf_map)
%             case 'free'%                 colormap(free_map)
%         end
        eval(['colormap(' event_type '_map)'])
        colorbar
    end
end

%% ** surface **
function S = map_conf_params(trc_events, par, filename)
% maps either duration 'T', frequency 'f' (=1/T) or size 'R'
% for 'free' or 'conf' event_type, as defined by MTT_proba_conf
% by interpolation through (ijk) points, with
% i,j barycentre of each event and k associated par, T, f or R

global N_PARAM
if isempty(N_PARAM), MTTparams_def; end

if nargin<2, par = 'T'; end

N_events = size(trc_events,2);
i_val = zeros(N_events,1);
j_val = zeros(N_events,1);
k_val = zeros(N_events,1);

for i_event = 1:N_events
    event = trc_events(:,i_event);
    i = event(2:N_PARAM:end);
    j = event(3:N_PARAM:end);
    alpha = event(4:N_PARAM:end);
    i = i(alpha>0);
    j = j(alpha>0);
    i_val(i_event) = mean(i);
    j_val(i_event) = mean(j);
    
    switch par
        case 'T'
            k_val(i_event) = numel(i);
        case 'f'
            k_val(i_event) = 1/numel(i);
        case 'R'
            k_val(i_event) = sqrt(calcul_size(event));
    end
end % for i_event = 1:N_events

%% ** masque **
DIC_name = dicname(filename);
maskdir = ['DIC\masques en position basse\Masque IN\MIB.' DIC_name(5:end)];
if ~isempty(dir(maskdir))
    mask = imread(maskdir);
    Rdil = 5;
    mask = imdilate(mask,strel('disk',Rdil)); % on dilate d'une marge de Rdil pxl
    [i_val, j_val, ind] = select_pk_in_mask(i_val, j_val, mask);
    k_val = k_val(ind>0);
    ROI = [1 size(mask,1) 1 size(mask,2)];
else
    ROI = [floor(min(i_val)) ceil(max(i_val)) floor(min(j_val)) ceil(max(j_val))];
    disp('no mask found..')
end

% hist(log10(k_val))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ** sélection pics dans masque **

if ~isempty(dir(maskdir))
%% contour
    Rdil2 = 10;
    maskINdil2 = imdilate(mask,strel('disk',Rdil2)); % on dilate d'une marge de Rdil pxl
    maskCONTdil = bwperim(maskINdil2);% imdilate(maskCONT,strel('disk',Rdil));
    [im jm] = find(maskCONTdil);
    im = im(1:15:end); jm = jm(1:15:end); % on échantillonne à 1/15 => 200 pts IMPAIR DE PREF!!
    i_val(end+1:end+size(im)) = im;
    j_val(end+1:end+size(jm)) = jm;
    k_val(end+1:end+size(im)) = ones(size(im))*min(k_val); 
end

%% ajout des bords
xb = ROI(1); xe = ROI(2); yb = ROI(3); ye = ROI(4);
xx = (xb:20:xe); yy = (yb:20:ye); N = 2*length(xx)+2*length(yy);
i_val(end+1:end+N) = [xx, ones(size(yy))*xb, xx, ones(size(yy))*xe]';
j_val(end+1:end+N) = [ones(size(xx))*yb, yy, ones(size(xx))*ye, yy]';
k_val(end+1:end+N) = ones(N,1)*min(k_val); % bords = min, plancher

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% % % cf S = carto3Dv3(event)
x1 = ROI(1):ROI(2); y1 = ROI(3):ROI(4); %:delta:
[Xi,Yi] = meshgrid(x1,y1);
S = griddata(i_val,j_val,k_val,Xi,Yi,'cubic');%'v4');%,{'Qt','Qbb','Qc','Qz'});%
S = max(S,min(k_val));
%%%