function traj_xyt(filename, pict, codage, bin, seuil, triD, viewDV, timing, min_length_ratio)
% function traj_xyt(filename, pict, codage, bin, seuil, triD, viewDV, timing, min_length_ratio)
%
% traces en 3D sur image dic
% if codage=='dens' % codage par densité, n/R2, cf. Douglass
% elseif codage=='speed' % codage par déplcmt (vitesse)
% elseif codage=='int' % codage par intensité (cluster...)
% elseif codage=='conf' % codage par index de confinement, Dt/R2, cf. Saxton
% elseif codage=='var' % codage par index de confinement, Dt/var(R2) cf. Meilhac
% elseif codage=='time' % codage par t
% bin: binaire, 2 couleurs, délim. par seuil relatif (def 1/2)
% arc en ciel, sinon: 1, magenta, slow... 6, red, fast (defaut)
% cf. cartobyf (et carto)

% V1.0 AS 3/5/2006
% V1.1 oct 2006 ajout du codage conf
% V1.2 dec 2006 ajout codage var

%%%%%%%%%%% EN CHANTIER!!! %%%%%%%%%%%%%%%%%%
%%%%%%%%%%% EN CHANTIER!!! %%%%%%%%%%%%%%%%%%
%%%%%%%%%%% EN CHANTIER!!! %%%%%%%%%%%%%%%%%%

global N_PARAM PARAM_I PARAM_J PARAM_ALPHA
if isempty(N_PARAM), MTTparams_def; end

if nargin<2, pict = []; end
if nargin<3, codage = 'rel'; end
if nargin<4, if strcmp(codage,'rel'), bin = 1; else bin = 0; end, end
if nargin<5, seuil = 1/2; end
if nargin<6, triD = 1; end
if nargin<7, viewDV = []; end
if nargin<8, timing = 36; end % ms timing = eval_timing(filename)
if nargin<9, min_length_ratio = 0; end % utile pour cell_track

if strcmp(codage, 'time'), no_blink = 1; else no_blink = 0; end

% if strcmp(filename(end-4:end),'_left'), filename(end-4:end) = []; viewDV = '_left'; end
% if strcmp(filename(end-5:end),'_right'), filename(end-5:end) = []; viewDV = '_right'; end

if nargin==0, files = dir('*.stk'); else files = dir(filename); end

if isempty(files), disp('No data... Check dir & filename !'), return, end

surf = (150^2)/(160^2); % surf élem, (150 nm)^2, cf. Douglass, Cell 2005, avec 1pxl = 160nm
n_max = 24/(timing/36); % densité "max", en pk/pxl/img. Empiriquement, 24 ok pour 36 ms
% si dt(=timing) diminue, dr diminue, donc n_max augmente...

Sm = 12; graph = 0; %tc = 6*timing/36; %Lc = 3; % Pour proba_conf, cf. Saxton

traj_colors = 'mbcgyr'; % rainbow order
Ncol = length(traj_colors); % 6...
Imax = 1300*2*Ncol*timing/36; % Imax deteremine par ct2 aussi??? (cf msdturbo)
if ~isempty(findstr(cd,'HTHLAB')) || strfind(filename,'_filt'), Imax = 255; end

%% boucle / file
for i = 1:length(files)
    file = files(i).name;
    disp('loading data...')
    tab_param = fread_all_param(file);
    if isempty(tab_param), return, end
    
    ntrc = size(tab_param,2); % sum(diff(trcdata(:,1))>0)+1;
    Tmax = size(tab_param,1)/N_PARAM;
    min_trc_length = round(Tmax*min_length_ratio);
    
    nf = strfind(file,'_filt');
    if ~isempty(nf) % enleve '_filt.tif'
        DIC_name = dicname([file(1:nf-1) '.tif']);
    else
        DIC_name = dicname(file);
    end
    
    DIC_image(filename, DIC_name, triD)
%     figure('WindowStyle','docked')
% 
%     if isempty(pict) && ~isempty(dir(DIC_name))
%         DIC = imread(DIC_name);
%         sat = .002 ; % saturation 0.2% min-max du contraste
%         DIC_sat = imadjust(DIC,stretchlim(DIC, [sat 1-sat]),[0 1]);
%         H = fspecial('average');
%         pict = imfilter(DIC_sat,H,'replicate');
%     else
%         % xy_max = max([pkdata(:,3); pkdata(:,4)]); % pict = zeros(ceil(xy_max));
%         if strfind(filename,'.stk')
%             pict = tiffread(filename,1);
%         else
%             pict = imread(filename,1);
%         end
%         pict = max(pict(:)) - pict; % invert
%     end
%     
%     %% *** met l'image de la cellule "au plancher" ***
%     %         pict_width = size(pict,2);
%     %         if strcmp(viewDV,'_left')%             pict = pict(:,1:pict_width/2); % 1/2 image de gauche !
%     %         elseif strcmp(viewDV,'_right')%             pict = pict(:,pict_width/2+1:end);%         end
%     pict = double(pict);
%     
%     if ~triD % 2D donc...
%         imagesc(pict)
%     else
%         pict(:,:,2) = 0; % juste pour créer une 3D...
%         h = slice(pict,[],[],1); % coupe à z( soit t)=1
%         set(h,'EdgeColor','none')
%         zlabel('frames')
%         axis([1 size(pict,1) 1 size(pict,2) 0 Tmax])
%         view(65,55) % azimut, hauteur
%     end
%     colormap('gray'), hold on
%     axis ij image % axis tight, axis equal


    title({cd; filename ; [' codage: ' codage]}, 'interpreter', 'none')
    pause(.1)
    
    %% calcul de r2 (si utile, selon codage)
    if strcmp(codage, 'speed')
        r2 = calcul_r2(tab_param); max_dr = sqrt(max(r2));
    end
    
    %% --- go through traces ---
    disp('traj :            ')
    
    for itrc = 1:ntrc
        
        if mod(itrc,50)==0
            fprintf([repmat('\b',1,11) '%5i/%5i'],itrc,ntrc)
            figure(gcf)%drawnow
        end
        
        trci = tab_param(:,itrc);
        alpha = trci(PARAM_ALPHA:N_PARAM:end);
        
        if length(alpha)>min_trc_length
            
            %%  *** calcul index Lconf de Saxton ***
            if (strcmp(codage, 'conf') || strncmp(codage, 'var', 3)) && length(n_trci)>Sm
                Lconf = MTT_probaconf(trci, codage, timing, graph);
            end
            
            for istep = 1:length(n_trci)-1
                
                if (strcmp(codage, 'conf') || strncmp(codage, 'var', 3)) && (istep<=Sm/2 || istep>length(n_trci)-Sm/2), continue, end
                % début ou fin de trace, index Lconf non défini
                
                %%          coordonnées t,x,y du pas en cours
                ti = trcdata(Ni(istep),2); % trci(istep,2)??
                tii = trcdata(Ni(istep+1),2);
                xi = trcdata(Ni(istep),3);
                xii = trcdata(Ni(istep+1),3);
                yi = trcdata(Ni(istep),4);
                yii = trcdata(Ni(istep+1),4);
                
                if dtrci(istep)==1 % pas blink
                    
                    %%          *** calcul densité locale ***
                    if strcmp(codage, 'dens')
                        disti = sqrt((trcdata(Ni,3)-xi).^2 + (trcdata(Ni,4)-yi).^2);
                        ndensi = sum(disti<=sqrt(surf)); % densité = nb de points ds voisinage
                    end
                    %%          *** calcul des pas élém. (vitesses) ***
                    if strcmp(codage, 'speed')
                        dr = sqrt((xii-xi)^2+(yii-yi)^2);
                        dr = max(dr,eps); % if dr==0...
                    end
                    %%          *** codage couleur ***
                    if strcmp(codage, 'dens') % codage par densité, 1, magenta, slow... 6, red, fast
                        val = min(ndensi/n_max,1);
                        
                    elseif strcmp(codage, 'speed') % codage par déplcmt (vitesse)
                        val = max((1-dr/max_dr),eps);
                        
                    elseif strcmp(codage, 'int') % codage par intensité (cluster...)
                        val = min(inti(istep)/Imax,1);
                        
                    elseif (strcmp(codage, 'conf') || strncmp(codage, 'var', 3)) % codage par index de conf/variance, cf. Saxton ou Meilhac
                        val = Lconf(istep-Sm/2);
                        
                    elseif strcmp(codage, 'time') % codage par intensité (cluster...)
                        val = ti/Tmax;
                    end
                    
                    %%               ** codage binaire **
                    if bin  % codage par densité, binaire : 4 vert, slow ou 6, rouge, fast
                        if val<seuil, icolor = 4; else icolor = 6; end
                    else
                        icolor = ceil(val*Ncol); % 0 à 1 => 1 à 6
                        icolor = min(icolor,Ncol);
                        icolor = max(icolor,1);
                    end
                    
                    line_type = traj_colors(icolor);
                    
                else % blink
                    line_type = ':k'; icolor = 4;
                end
                
                %%              ** plot du segment istep **
                if triD
                    plot3 ([xi xii], [yi yii], [ti tii], line_type)
                    %                 else
                elseif dtrci(istep)==1 || ~no_blink
                    sat = (Ncol-icolor)/Ncol; % saturation de la couleur, 1/3 ou 0 si bin; entre 0 et 1 sinon
                    if bin && icolor==4, sat = 1/2; end %...
                    
                    if strcmp(viewDV,'_left') % en vert, contrairement aux conventions maritimes... vif=conf, pale=diff
                        line ([xi xii], [yi yii], 'color', [0 1/2 0]+sat)
                    elseif strcmp(viewDV,'_right')% idem, en rouge
                        line ([xi xii], [yi yii], 'color', [1/2 0 0]+sat)
                    else % bleu
                        if strcmp(codage, 'time')
                            line ([xi xii], [yi yii], 'color', hsv2rgb([val,1,1]))
                        else
                            line ([xi xii], [yi yii], 'color', [0 0 1/2]+sat)
                        end
                        %                         plot ([xi xii], [yi yii], line_type)
                    end
                end
            end % for istep = 1:length(Ni)-1
        end % if length(Ni)>0
    end % for itrc = 1:NoTrace
    
    axis image
    fprintf([repmat('\b',1,11) '%5i/%5i\r'],ntrc,ntrc)
    
    if isempty(dir('MCTfigs')), mkdir('MCTfigs'), end
    figname = ['MCTfigs' filesep filename 'carto3D.png'];
    saveas(gcf, figname, 'png')
    if length(files)>1, close(gcf), end
end % for i = 1:length(files)
