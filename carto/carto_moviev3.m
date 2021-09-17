function carto_moviev3(filename, do_contour, nstart, extraframes, dT)%!!! dT~=1 et N_av~=1 non garanti.....
% sauve images pour video avec carto3Dv3
% ou une seule image globale, si nstart=0

if nargin<2, do_contour = 1; end
if nargin<3, nstart = 1; end
if nargin<4, extraframes = 3; end
if nargin<5, dT = 1; end %if ~do_contour, dT = 1; else dT = 110; wT = wT+0, end

zmax = 370; % log10(zmax)=2.57
do_subplot = 0;
if strfind(cd,'output'), cd .., end % si on est déjà ds output..
[lg, ht, Ntiff] = stk_size(filename);
params_def = MTTparams_def; dirname = params_def{4}; % dirname = 'output23';
codage = 'rel';
dir_carto = ['carto_v3_' codage '_' dirname];

if nstart>0 % image ou video
    s = which(['MTT_conf_index_' codage]); param_dir = fileparts(s);
    params = load([param_dir filesep 'proba_conf_params_' codage '.dat']); % params = load('\\Amon\Picsl\HETMLAB\DMLAB\Matlab-scripts\Arnauld\carto\proba_conf_params_varM.dat')
    wT = params(3)+extraframes; % wconf(=4) +3 = 7 (NB ok pour calculer Dlocal vs Dmin)
    dir_carto = [dir_carto filesep filename(1:end-4) '_movie'];
    Nframes = floor((Ntiff-wT)/dT); % Ntiff-N_av-wT+1
    kmax = 3;
else % pour une seule image, au lieu d'un movie
    nstart = 1; % ...
    wT = inf;
    Nframes = 1;
    kmax = 1;
end


%%  *** colormap ***
n = 256; sat = 4/3;
r = [(0:n*(1-1/2/sat)-1)'/(n*(1-1/2/sat)); ones(n/2/sat,1)];
g = [zeros(n*(1-3/4/sat),1); (0:n/2/sat-1)'/(n/2/sat); ones(n/4/sat,1)];
b = [zeros(n*(1-1/2/sat),1); (0:n/2/sat-1)'/(n/2/sat)];
AFMmap = [r g b]; AFMmap(1,:) = [0 0 .3]; % plancher bleu marine

levels_val = [0 log10(zmax)/256:0.4:log10(zmax)]; % 0:0.4:log10(zmax); % pour fixer les niveaux ad vitam, 1e au min AS 28/8/8

%% ini
if isempty(dir(dir_carto)), mkdir(dir_carto), end
disp(['& c parti pour ' num2str(Nframes-nstart+1) ' image(s)'])
Surf3Dk = zeros(lg,ht,kmax);

t_first = (nstart-1)*dT+1;
t_last = (nstart-1)*dT+wT+kmax-2;
trc = detect_reconnex_to_trc(filename,1,dirname,t_first,t_last);
trc_mobiles = discard_fix(trc, filename); % ON EJECTE LES QD FIXES, 2/9/8
tt = trc_mobiles(:,2);

for k = 1:kmax-1 % calcule kmax-1, pour 1e image
    t_first = (nstart-1)*dT+1;
    t_last = (nstart-1)*dT+wT;
    trck = trc_mobiles(tt>=t_first+k-1 & tt<=t_last+k-1,:);
    Surf3Dk(:,:,k) = carto3Dv3(trck, filename, do_subplot);
end

%% boucle sur frames
for n_image = nstart:Nframes % 1:995
    t_first = (n_image-1)*dT+1;
    t_last = (n_image-1)*dT+wT;
    fprintf('frame %g / %g ',n_image,Nframes)
    
    trc = detect_reconnex_to_trc(filename,1,dirname,t_first+kmax-1,t_last+kmax-1);
    trc_mobiles = discard_fix(trc, filename); % ON EJECTE LES QD FIXES, 2/9/8
    tt = trc_mobiles(:,2);
    
    trcn = trc_mobiles(tt>=t_first+kmax-1 & tt<=t_last+kmax-1,:);
    Surf3Dk(:,:,kmax) = carto3Dv3(trcn, filename, do_subplot);
    Surf3D = mean(Surf3Dk,3); % moy. temp. sur kmax frames
    Surf3Dk(:,:,1:kmax-1) = Surf3Dk(:,:,2:kmax); % permute
    
    %% *** contour, dot plot, save fig ***
    if n_image==nstart, figure('WindowStyle','docked'), end % 1e image
    if ~do_contour
        surf(Surf3D,'linestyle','none','facecolor','interp')
        view(-20,80), axis ij off tight % colorbar
        a = axis; axis([a(1:4) 0 log10(zmax)])
    else %do_contour
        contourf(Surf3D,levels_val)
        axis ij off image % colorbar, xy?
    end
    colormap(AFMmap), caxis([0 log10(zmax)]), drawnow
    
    %     save(['video_conf\frame' num2str(n_image,'%04.0f') '.txt'], 'Surf3D', '-ascii'), 0+0 %-append
    
    if do_contour, outfile = [dir_carto filesep filename(1:end-4) '_Cont' num2str(n_image) '.png'];
    else outfile = [dir_carto filesep filename(1:end-4) '_3D' num2str(n_image) '.png']; end
    %if ~isempty(dir(outfile)), disp ([outfile ' already done']), continue, end
    print('-dpng', outfile) % '-r300',    % saveas(gcf,outfile,'png') %tif??
    pause(.1)
end % for n_image=1:Nframes %% boucle sur frames
%%
crop_png(filename, do_contour)
%%%