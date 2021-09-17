function traj_xyt_movie(file)
% trace les trajs (xy) image après image (t)

if nargin==0
    cd('M:\Arnauld Serge\2008-09-23 CO+LatA')
    file = 'ctl cell5 1000frames.stkkk';
    %     DIC_name = dicname('ctl cell5.stk');
end

params_def = MTTparams_def; dirname = params_def{4};
color_dif = [.5 0 0]; % bordeaux
color_conf = [1 .5 0]; % orange
timing = eval_timing(file);
bloc_size = 50;

codage = 'rel';
s = which(['MTT_conf_index_' codage]); param_dir = fileparts(s);
params_conf = load([param_dir filesep 'proba_conf_params_' codage '.dat']); % params = load('\\Amon\Picsl\HETMLAB\DMLAB\Matlab-scripts\Arnauld\carto\proba_conf_params_varM.dat')
Lmin = params_conf(1);
Sm = round(params_conf(3)*36/timing);

DIC_name = dicname(file);

DIC_image(file)
% % % figure('WindowStyle','docked');
% % % %% DIC
% % % if ~isempty(DIC_name)
% % %     DIC = imread(DIC_name);
% % %     sat = .002 ; % saturation 0.2% min-max du contraste
% % %     DIC_sat = imadjust(DIC,stretchlim(DIC, [sat 1-sat]),[0 1]);
% % %     H = fspecial('average');
% % %     DIC_sat_f = imfilter(DIC_sat,H,'replicate');
% % %     imagesc(DIC_sat_f)    %brighten(+.3) % 24/4/7
% % %     axis tight image %([1 size(DIC_sat,1) 1 size(DIC_sat,2)])
% % % else axis equal
% % % end
% % % axis ij, colormap('gray'), hold on
% % % drawnow

%% ** load masque **
maskIN = [];
maskdir = ['DIC\masques en position basse\Masque IN\MIB.' DIC_name(5:end)]; % remove 'dic\'
if ~isempty(dir(maskdir)) % && usemask
    maskIN = imread(maskdir);
    disp(['using Masque IN\MIB.' DIC_name(5:end)]);
end

fid = fopen([dirname filesep file '_tab_param.dat'], 'rt','native') ;
value =  fscanf(fid, '%d', 2) ;
fclose(fid);
nb_t = value(2) ; % nb_part_max = value(1) ;
Nbloc = floor(nb_t/bloc_size);

for bloc=1:Nbloc % on tronconne, pour limite memoire
    t_first = 1+(bloc-1)*bloc_size;
    t_last = bloc*bloc_size;
    
    %% trc: load, select in mask, proba conf
    trc = detect_reconnex_to_trc(file, 1, '', max(1, t_first-Sm), min(t_last+Sm, nb_t)); % trc = [ones(8,1) (1:8)' rand(8,2)];
    
    if ~isempty(trc) && ~isempty(maskIN) %&& usemask
        [itoto jtata ind] = select_pk_in_mask(trc(:,4), trc(:,3), maskIN); % NB y x == i j
        trc = trc(ind>0,:);
    end
    
    Lconf = probaconf(trc, 0, timing, codage);
    trc(:,5) = Lconf;
    
    %% blink
    dn = diff(trc(:,1)); % # trc => diff=1 (ou +..) pour chaque new trc
    trc(:,1) = [0; dn]; % !!! on remplace #trc par diff !!!
    tt = trc(:,2); % # frame (temps)
    dt = diff(tt);
    trc(:,6) = [1; dt-1]; % diff(tt) = 1 normal, diff(tt)~=1 = blink ou new trc (diff(nn)=1)
    %     Nframes = max(tt);
    
    %% let's go
    if isempty(dir([file(1:end-4) '_mov_xyt'])), mkdir([file(1:end-4) '_mov_xyt']), end
    
    for t=t_first:min(t_last, nb_t) %1:Nframes
        indi = find(tt==t); % # des pics de l'image t
        pki = trc(indi,:);
        pk_previous = pki;
        pk_previous(indi>1,:) = trc(indi(indi>1)-1,:); % pics img precedente (on conserve pki si c le début)
        
        for j=1:size(pki,1)
            if pki(j,1)==0 % dn>0 => new trc, segment à ne pas tracer
                if pki(j,6)~=0, ls = ':'; else ls = '-'; end
                if pki(j,5)>Lmin, lc = color_conf; else lc = color_dif; end
                xx = [pk_previous(j,3) pki(j,3)];
                yy = [pk_previous(j,4) pki(j,4)];
                plot(xx, yy, 'LineStyle', ls, 'Color', lc);
            end
        end
        
        outfile = [file(1:end-4) '_mov_xyt' filesep file(1:end-4) '_xyt' num2str(t,'%04i') '.png' ];
         print('-dpng', '-r300', outfile) % saveas(gcf,outfile,'png')
        %         f = getframe(gca); % rect(roi)??
        %         fig_i = im2uint8(f.cdata);
        %         imwrite(fig_i, outfile, 'png')
        %         imwrite(fig_i, outfile, 'tif', 'WriteMode', 'append', 'compression', 'none')
        fprintf('%5i/%5i\r',t,nb_t) %fprintf([repmat('\b',1,11) '%5i/%5i'],itrc,ntrc)
    end
end