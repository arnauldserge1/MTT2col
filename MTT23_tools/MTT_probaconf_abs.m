function [L lenConf lenFree freqConf sizeConf trc_conf trc_free beginConf endConf sizeFree Dlocal] = probaconf(trc, graph, method, timelag, params)

% function [L lenConf lenFree freqConf sizeConf trc_conf trc_free begin_conf end_conf sizeFree Dlocal] 
%           = probaconf(trc, graph, method, timelag, params)
%
% Disclaimer: Developped for Matlab, not checked for Octave!
%
% EN/ Allows to determine if trajectories display confinement along time. 
% Using a sliding window along the trajectory, we calculate the function 
% of Saxton & Simson / Meilhac (according to the method used). 
% Then we threshold this function, which leads to
% confined events. This calculation is done for all trajs long enough
% and we build the associated statistics (length of confinement...)
%
% Input:   requires the matrix trc 
%               graph: plot or not data
%               Method = 'conf' (Saxton), 'dens' (Douglass) or 'abs' (Meilhac, def)
%               timelag in ms (def 36)
%               param = [Lc tc Sm], threshold, min. time and calculation window
%       def: params = load('../carto/proba_conf_params_varM.dat');
%
% Output:   3 graphs, trajectory with zone(s) of confinement in orange, 
%               the diffusion coefficient vs time,
%               and the function of confinement vs time, with the threshold
%               used
%
% cf Simson 95 et Saxton 93, Douglass 2005, Meilhac 2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FR/ Elle permet de déterminer s’il y a au cours du temps des confinements des
% trajectoires. A l’aide d’une fenêtre glissante qui parcours la trajectoire
% on calcul la fonction de Saxton & Simson/Meilhac (selon method). 
% Puis on seuil cette fonction ce qui nous
% donne les épisodes confinés. On effectue ce calcule sur toutes les
% trajectoires assez longues et on calule les statistiques associées
% (longueur de ces confinements...)
%
% En entrée :   nécessite la matrice .trc 
%               graph: plot ou non data
%               Method = 'conf' (Saxton), 'dens' (Douglass) ou 'abs' (Meilhac, def)
%               timelag en ms (def 36)
%               param = [Lc tc Sm], seuil, temps min. et fenetre de calcul
%       def: params = load('../carto/proba_conf_params_abs.dat');
%
% En sortie :   On a 3 graphiques, la trajectoire avec la/les zone(s) de
%               confinement en orange, le coefficient de diffusion en
%               fonction du temps, et la fonction de confinement au cours du
%               temps et le seuil utilisé.
%
% cf Simson 95 et Saxton 93, Douglass 2005, Meilhac 2006


% V1.2 dec 2006 ajout codage var
% V1.3 oct 2007 pour var(R), R2=(x(i)-x(Sm/2))2+(y y), sur fenetre de Sm (wn) points (Sm+1 initialmt)

if nargin<1
    output_dir = 'output23'; files = dir('*.tif'); filename = files(1).name;
    trc = detect_reconnex_to_trc([output_dir filesep filename]);%(['..' filesep 'data' filesep 'EGFR-Qd605.stk']); 
end
if nargin<2, graph = 1; end
if nargin<3, method = 'abs'; end
if nargin<4, timelag = 36; end % ms..
if nargin<5
    s = which(['MTT_conf_index_' method]);
    param_dir = fileparts(s); % s.path;
    params = load([param_dir filesep 'proba_conf_params_' method '.dat']);
end
Lc = params(1);
tc = round(params(2)*36/timelag); % temps conf min, en # image (=69 chez simson)
Sm = round(params(3)*36/timelag); 

if isempty(trc), Ntraj = 0; else Ntraj = trc(end,1); end
global TRAJN
% TRAJN=121;
if isempty(TRAJN), TRAJN = 1:Ntraj; end % n° of trajs to test
if length(TRAJN)==Ntraj, newfigs = 0; else newfigs = 1; end % si sélection restreinte, figs indep 

Tmax = size(trc,1);
L = zeros(Tmax,1);%cell(1,Ntraj);
freqConf = zeros(1,Ntraj);
lenConf = cell(1,Ntraj);
lenFree = cell(1,Ntraj);
sizeConf = cell(1,Ntraj);
sizeFree = cell(1,Ntraj);
beginConf = cell(1,Ntraj); % 31/5/7
endConf = cell(1,Ntraj);
trc_conf = [];
trc_free = [];
Nconf = zeros(1,Ntraj);
Dlocal = zeros(Tmax,1);%cell(1,Ntraj);

if isempty(trc), 
    lenConf = []; lenFree = []; sizeConf = []; beginConf = []; endConf = []; sizeFree = [];
    if graph, disp('no traces...'), end
    return
end

Dglobal = 0.25*timelag/36; % pxl/lag, rescalé!!!timelag AS 5/4/7 %%%%%%% (0.22*0.16^2)/0.036 ; %en µm²/

ifree = 1; iconf = 1;
length_prev_trc = 0;

color_dif = [.5 0 0]; % bordeaux
color_conf = [1 .5 0]; % orange
    
if Ntraj>trc(1,1) % => plusieurs trajs
    disp('detect conf dans traj:             '), 
end

for itraj = TRAJN %1:Ntraj%
    
    indi = find(trc(:,1)==itraj);
    trci = trc(indi,:);
    nfinal = size(trci,1);
    
    if mod(itraj,20)==0 && Ntraj>trc(1,1)
        fprintf(repmat('\b',1,13))
        fprintf('%6i/%6i', itraj, Ntraj)
    end

%% ****** index de confinement Li ******
    if nfinal>Sm
        Li = zeros(nfinal-Sm+1,1);
        Dlocali = zeros(nfinal-Sm+1,1);      

        for n = 1:nfinal-Sm+1 % NB fenetre glissante #n
            trcn = trci(n:n+Sm-1,:); % bout de trace #n
            xn = trcn(:,3); yn = trcn(:,4);
            R = sqrt((xn-xn(ceil(Sm/2))).^2 + (yn-yn(ceil(Sm/2))).^2); % ref: 1e point (cf Simson) (milieu du tronçon, Sm/2??) 

            if graph || nargout>10
                msdn = msd(trcn,1,0);
                DD = calculDinst(msdn,Sm-1); % [n D sD]
                Dlocali(n) = max(DD(2),0); % Diff instantanée pour ce troncon (<0 foireux...)
            end

%%  **** fct de Saxton ou variantes de Meilhac / Douglass ****
            switch method
%                 case 'conf'
%                     Li(n) = 2.5117*(Dglobal)*Sm / R2mean(n) - 0.2048 - 1;
%                     %Dtrc(itraj)*.16^2/(timelag/1000)
                case 'abs'
                warning off MATLAB:divideByZero;
                    Li(n) = Dglobal*Sm / var(R); % Meilhac
%                 case 'dens'
%                     distn = sqrt((xi-xi(n)).^2 + (yi-yi(n)).^2); 
%                     surf = (150^2)/(160^2); % surf élem, (150 nm)^2, cf. Douglass, Cell 2005, avec 1pxl = 160nm
%                     Li(n) = sum(distn<=sqrt(surf)); % densité = nb de points ds voisinage
                otherwise
                    disp('Wrong method name!')
            end
        end % for n = 1:nfinal-Sm % NB fenetre glissante #n

%% détermination des épisodes confinés %%
        Nconf(itraj) = sum(diff(Li>Lc)==1); % nombre d'épisodes conf. (transitions Li<=Lc à Li>Lc)
        begin_conf = find(diff(Li>Lc)==1)+1; % index (vectoriel) des moments où Li>Lc, cad qu'on devient conf (image suivante, dc +1)
        end_conf = find(diff(Li>Lc)==-1); % plus précisément, index des moments où l'on redevient non conf, sur le critère Li<=Lc

        if Li(1)>Lc, begin_conf = [1; begin_conf]; end %else end_conf = [1; end_conf]; end % on commence confiné (ou non!)
        if Li(end)>Lc, end_conf = [end_conf; length(Li)]; end %else begin_conf = [begin_conf; length(Li)]; end % on termine...

        length_conf_events = end_conf - begin_conf; % durée des confinements, par définition... (R), y compris too short, à ce niveau

%% cut off the too short events
        ind_too_short = (length_conf_events < tc);
        begin_conf (ind_too_short) = [];
        end_conf (ind_too_short) = [];

        length_conf_events = end_conf - begin_conf;
        length_free_events = begin_conf(2:end) - end_conf(1:end-1);
        Nconf(itraj) = length(end_conf);

        if ~isempty(begin_conf) %% gestion 1e et dernier
            if Li(1)<=Lc, length_free_events = [begin_conf(1); length_free_events]; end
            if Li(end)<=Lc, length_free_events = [length_free_events; length(Li) - end_conf(end)]; end
        else 
            if Li(1)<=Lc, length_free_events = length(Li); end % free tout du long sur cette traj
        end
        if isempty(end_conf) && Li(1)>Lc, length_conf_events = length(Li); end        % conf " " " (rem: cas déjà pris en compte + haut??..)

%% calcul taille, durée & fréq.
        if Nconf(itraj)>0

            freq = zeros(1,Nconf(itraj)-1);
            size_conf_events = zeros(1,Nconf(itraj));
            size_free_events = zeros(1,Nconf(itraj)+1);
            
            %% gestion 1e et dernier
            if Li(1)<=Lc, size_free_events(1) = calcul_size(trci(1:begin_conf(1),:)); end
            if Li(end)<=Lc, size_free_events(end) = calcul_size(trci(end_conf(end):length(Li),:)); end

            for paulo = 1:Nconf(itraj) % boucle sur confs de itraj
                if paulo<Nconf(itraj)
                    freq(paulo) = 1/(begin_conf(paulo+1)-begin_conf(paulo)); 
                    size_free_events(paulo+1) = calcul_size(trci(end_conf(paulo):begin_conf(paulo+1),:));
                end                    
                size_conf_events(paulo) = calcul_size(trci(begin_conf(paulo):end_conf(paulo),:)); %% taille (R2 max) du ss dom confiné #itraj  end_conf(paulo)-1??
            end
            
            size_free_events = size_free_events(size_free_events>0);
            if Nconf(itraj)>1, freqConf(itraj) = mean(freq); end
            if ~isempty(length_conf_events), lenConf{itraj} = length_conf_events'; end
            if ~isempty(length_free_events), lenFree{itraj} = length_free_events'; end
            if ~isempty(size_conf_events), sizeConf{itraj} = size_conf_events; end
            if ~isempty(size_free_events), sizeFree{itraj} = size_free_events; end
            
            if ~isempty(begin_conf), beginConf{itraj} = begin_conf+length_prev_trc+ceil((Sm-1)/2); end
            if ~isempty(end_conf), endConf{itraj} = end_conf+length_prev_trc+ceil((Sm-1)/2); end
            
%% matrices trajs free/conf
            if begin_conf(1)>1 % on commence non confiné, bug corr. 14/2/7 (StV!)
                index_free = 1:begin_conf(1);
                trc_free = [trc_free; [ifree*ones(length(index_free),1) trci(index_free,2:end)]];
                ifree = ifree+1;
            end

            for paulo = 1:Nconf(itraj)

                if paulo<Nconf(itraj)
                    index_free = end_conf(paulo):begin_conf(paulo+1);
                    trc_free = [trc_free; [ifree*ones(length(index_free),1) trci(index_free,2:end)]];
                    ifree = ifree+1;
                end

                index_conf = begin_conf(paulo):end_conf(paulo);
                %RDpaulo = calculRD(msd(trci(index_conf,2:end,1,0)));
                %if RDpaulo<RDmax % écarte trajs linéaires (>m+3std)
                    trc_conf = [trc_conf; [iconf*ones(length(index_conf),1) trci(index_conf,2:end)]];
                    iconf = iconf+1;
                %end
            end

            if end_conf(end)<length(Li) % on termine non confiné, ajouté 14/2/7 (StV!)
                index_free = end_conf(end):length(Li);
                trc_free = [trc_free; [ifree*ones(length(index_free),1) trci(index_free,2:end)]];
                ifree = ifree+1;
            end

        else % pas de changement d'état: tout free ou ...
            if Li(1)<=Lc
                trc_free = [trc_free; [ifree*ones(nfinal,1) trci(:,2:end)]];
                ifree = ifree+1;
                lenFree{itraj} = length(Li);
                sizeFree{itraj} = calcul_size(trci);
            else %... tout conf
                trc_conf = [trc_conf; [iconf*ones(nfinal,1) trci(:,2:end)]];
                iconf = iconf+1;
                lenConf{itraj} = length(Li);
                sizeConf{itraj} = calcul_size(trci);
             end
        end % if Nconf(itraj)>0
        
%% ***************** PLOTS **************
        Tc = sum(lenConf{itraj}); Tf = sum(lenFree{itraj});
        alpha = Tc/(Tc+Tf);
        if graph %&& Nconf(itraj)>2 && Nconf(itraj)<10 && alpha>.3 && alpha<.7 && length(Li)>100
%% traj
            if newfigs || itraj==TRAJN(1), figure('WindowStyle','docked'), else clf, end
            subplot(1,2,1)
            x = trci(ceil(Sm/2):end-ceil(Sm/2),3)*.16; y = trci(ceil(Sm/2):end-ceil(Sm/2),4)*.16; % current trc, #i
            plot(x-x(1),y-y(1),'Color',color_dif), hold on, plot(0,0,'ob')
            axis ij % hold on
            
%% ajout des ss domaines confs
            if ~isempty(Lc)
                xconf = cell(1,Nconf(itraj)); yconf = cell(1,Nconf(itraj));
                for nconf=1:Nconf(itraj)
                    xconf{nconf} = x(begin_conf(nconf):end_conf(nconf)-1);
                    yconf{nconf} = y(begin_conf(nconf):end_conf(nconf)-1);
                    plot(xconf{nconf}-x(1),yconf{nconf}-y(1),'Color',color_conf)
                end
            end
            hold off
            xlabel('x (µm)'), ylabel('y (µm)'), title(['traj ' num2str(itraj)])
            axis image
            
%% index Li
            subplot(2,2,2)
            tt = (1:length(Li))*(timelag/1000);
            semilogy(tt,Li,'k') % Li(t)
            hold on
            if ~isempty(Lc), semilogy(tt,Lc*ones(size(Li)),'Color',color_conf,'LineStyle','--'), end % seuil
            a = axis; axis([0 tt(end) a(3:4)]) % nfinal-Sm
            xlabel('time (s)'), ylabel('conf. index')
            title(['alpha_{conf} = ' num2str(alpha)])

            semilogy([0 tt(end)],[a(4) a(4)],'Color',color_dif,'LineWidth',5) % moments dif.
            if ~isempty(Lc)
                for nconf=1:Nconf(itraj)
                    tn = [begin_conf(nconf) end_conf(nconf)]*(timelag/1000);
                    semilogy(tn,[a(4) a(4)],'Color',color_conf,'LineWidth',7) % moments conf.
                    %a = axis; axis([0 tn(2) a(3:4)])
                end
            end
            hold off

%% Diff
            subplot(2,2,4)
            semilogy(tt(Dlocali>0),Dlocali(Dlocali>0)*.16^2/(timelag/1000),'k')
            xlabel('time (s)'), ylabel('D_{inst} (µm²/s)')
            a = axis; axis([0 tt(end) a(3:4)])
            
            if ~newfigs, fprintf('\r strike any key...'), pause, end
        end % if graph
        
        L(indi) = [zeros(ceil((Sm-1)/2),1); Li; zeros(floor((Sm-1)/2),1)];%L{itraj} = ;
        Dlocal(indi) = [zeros(ceil((Sm-1)/2),1); Dlocali; zeros(floor((Sm-1)/2),1)];%Dlocal{itraj} = ;
        % points début et fin complétés par 0 pour garder meme taille
    end % if nfinal>Sm
    length_prev_trc = length_prev_trc+nfinal;
end % for itraj=1%:Ntraj

if Ntraj>trc(1,1)
    fprintf(repmat('\b',1,13))
    fprintf('%6i/%6i done \r', Ntraj, Ntraj) % affiche final
end

global USECELL
if isempty(USECELL), USECELL = 0; end

if ~USECELL
    lenConf = cell2mat(lenConf);
    lenFree = cell2mat(lenFree);
    sizeConf = cell2mat(sizeConf);
    sizeFree = cell2mat(sizeFree);
    beginConf = cell2mat(beginConf');
    endConf = cell2mat(endConf');
end

clear global TRAJN

function trc_size = calcul_size(trc)
x = trc(:,3);
y = trc(:,4);
R2 = (x-mean(x)).^2+(y-mean(y)).^2;
trc_size = mean(R2);
%%%