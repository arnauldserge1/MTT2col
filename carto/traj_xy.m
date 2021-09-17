function traj_xy(filename, DIC_name, codage, dirname)
% function traj_xy(filename, DIC_name, codage, dirname)
%
% traces sur image dic
% if codage=='conf' % codage par index de confinement, Dt/R2, cf. Saxton
% elseif codage=='varM' % codage par index de confinement, Dt/var(R2) cf. Meilhac
% ou enfin, en exclusivite, codage='Nico'!!
% see also cartobyf

% V1.0 AS 3/5/2006
% V1.1 oct 2006 ajout du codage conf
% V1.2 dec 2006 ajout codage var
% V2 fev 2007
% v3 26/3/9


if nargin<2, DIC_name = dicname(filename); end
if nargin<3, codage = ''; end % rel??
if ~isempty(codage), disp(['method used: ' codage]); end
if nargin<4, params_def = MTTparams_def; dirname = params_def{4}; end

if isempty(filename), disp('No data... Check dir & filename !'), return, end

figure('WindowStyle','docked')

%% *** met l'image de la cellule "au plancher" ***
if ~isempty(DIC_name)
    DIC = imread(DIC_name);    %DIC = uint16(tiffread...
    sat = .002 ; % saturation 0.2% min-max du contraste
    DIC_sat = imadjust(DIC,stretchlim(DIC, [sat 1-sat]),[0 1]);
    H = fspecial('average');
    pict = imfilter(DIC_sat,H,'replicate');
else
    pict = imread(filename,1); % ...
    pict = max(pict(:)) - pict; % invert
end

imagesc(pict)    %brighten(+.3) % 24/4/7
axis tight image ij %([1 size(DIC_sat,1) 1 size(DIC_sat,2)]) % else axis equal
colormap('gray'), hold on
drawnow

%% ** load masque **
maskdir = ['DIC\masques en position basse\Masque IN\MIB.' DIC_name(5:end)]; % remove 'dic\'
if ~isempty(dir(maskdir))
    maskIN = imread(maskdir);
    str_mask = [', using Masque IN\MIB.' DIC_name(5:end)];
else
    maskIN = [];
    str_mask = ''; %   str_mask = ', no mask';%  disp(['no mask found for ' maskdir])
end

%% -- Load traces --
if ischar(filename)
    title([cd, ' \ ', filename, ', codage: ', codage, str_mask], 'interpreter', 'none')
    full_filename = [dirname filesep filename '_tab_param.mat'];
    load(full_filename) % tab_param = fread_all_params(full_filename);
else %...simul, ou run avec tab_param au lieu de filename
    title([cd, ', codage: ', codage, str_mask],'interpreter','none')%%, ' \ ', 'simul data.mat', ', codage: ', codage, str_mask],'interpreter','none')
    tab_param = filename; % cf. f.tab_param = ...; trc=detect_reconnex_to_trc(f);
end

if ~isempty(tab_param)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_traces(tab_param, maskIN, codage, filename)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else %(out of mem..)
    tab_param = 'toto'; isub = 1;
    while ~isempty(tab_param)
        full_filename = [dirname filesep filename '_tab_param.dat'];
        tab_param = fread_params_timewindow(full_filename,1,1,inf,100*isub-99,100*isub);
        if isempty(tab_param) && isub==1, disp('no traces...'), return, end
        isub = isub+1;
        
        plot_traces(tab_param, maskIN, codage, filename)
    end % while ~isempty(trcdata)
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_traces(tab_param, maskIN, codage, filename)
% Plot traces for a given trcdata matrix, checking for blink, contour mask
% of cell, confinment domains

global N_PARAM PARAM_I PARAM_J PARAM_ALPHA
if isempty(N_PARAM), MTTparams_def; end

ntrc = size(tab_param,2);

if ~isempty(codage)
    color_dif = [.5 0 0]; % bordeaux ([0 0 1] Marek)
    color_conf = [1 .5 0]; % orange
else
    color_dif = [0 0 1]; % bleu
end

timing = eval_timing(filename);

disp('traj :            ')

%% --- go through traces ---
for itrc = 1:ntrc
    if mod(itrc,200)==0
        fprintf([repmat('\b',1,11) '%5i/%5i'],itrc,ntrc)
        drawnow
    end
    
    trci = tab_param(:,itrc); % reshape(trcdata(:,itrc)',7,[])';
    
    %%  *** test hors masque de la cellule (ou immobile/traj (HIV) correcte) ***
    if ~isempty(trci) && ~isempty(maskIN)
        ii = trci(PARAM_I-1:N_PARAM:end);
        jj = trci(PARAM_J-1:N_PARAM:end);
        alpha  = trci(PARAM_ALPHA-1:N_PARAM:end);
        
        [~, ~, ind] = select_pk_in_mask(ii, jj, maskIN, alpha); % NB y x == i j
        
        trci(PARAM_ALPHA-1:N_PARAM:end) = alpha.*(ind>0); % trci(2:N_PARAM:end) = ii.*(ind>0); trci(3:N_PARAM:end) = jj.*(ind>0);
        %             plot(jj.*(ind>0),ii.*(ind>0),'.'), pause
    end
    
    %% ****************** calcul index Lconf *********************
    if ~isempty(codage)
        [~ , ~, ~ , ~ , ~, trc_conf, trc_free] = ...
            MTT_probaconf(trci, codage, timing, 0, filename); % graph=0
    else
        trc_conf = []; trc_free = trci;
    end
    
    if isempty(trc_conf) && isempty(trc_free), continue, end
    clear Lconf Tc Tf R2c R2f
    
    %%  ** plot traces **
    if isempty(trc_conf) %&& ~isempty(trc_free)
        plotwithblink(trc_free, color_dif)
    elseif isempty(trc_free) %&& ~isempty(trc_conf)
        plotwithblink(trc_conf, color_conf)
    else % NOTA: Trajs tout conf ou tout free écartées par methode de Nico
        Nevent = size(trc_free,2);%{1}
        
        for ievent = 1:Nevent
            trcfreei = trc_free(:,ievent);
            plotwithblink(trcfreei, color_dif)
        end
        
        Nevent = size(trc_conf,2);
        
        for ievent = 1:Nevent % trace conf apres, par dessus!
            trcconfi = trc_conf(:,ievent);
            plotwithblink(trcconfi, color_conf)
        end
    end
    %     disp(['traj ' num2str(itrc)]), pause
end % for itrc = 1:NoTrace
fprintf([repmat('\b',1,11) '%5i/%5i'],ntrc,ntrc)
fprintf('\r')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotwithblink(trc, clr)

global N_PARAM PARAM_T PARAM_I PARAM_J PARAM_ALPHA
if isempty(N_PARAM), MTTparams_def; end

plot_blink = 1; % 0 for Marek...

t = trc(PARAM_T-1:N_PARAM:end);
y = trc(PARAM_I-1:N_PARAM:end);
x = trc(PARAM_J-1:N_PARAM:end);
alpha = trc(PARAM_ALPHA-1:N_PARAM:end); % xy ji!
t = t(alpha>0); x = x(alpha>0); y = y(alpha>0);
%      disp([t x y]), pause

dt = diff(t);
t_blink = find(dt>1);
t_blink = [0; t_blink; length(t)];

for i=1:length(t_blink)-1
    tt = t_blink(i)+1:t_blink(i+1);
    plot(x(tt),y(tt),'Color',clr,'LineWidth',1)
    if (i<length(t_blink)-1) && (plot_blink)
        tbl = [t_blink(i+1), t_blink(i+1)+1];
        plot(x(tbl),y(tbl),'Color', clr, 'LineStyle', ':')
    end
end
%%%