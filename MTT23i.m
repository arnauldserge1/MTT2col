%% MTT23i(params): Multiple Target Tracing algorithms version 2.3i
% (i: with user interface)
% This programme is compatible with Matlab and Octave software
%
% ALGORITHM AUTHORS:
% N. BERTAUX, A. SERGE
%
% lance detect_reconnex_23, puis carto et ct4 (histos D, I, Tfree/conf...)
% sur tous les fichiers .stk du répertoire local
% filename = '*.stk' par défaut
% refitdata = 0 par défaut, n'écrase pas data déjà fittées
% use_stk = 1 par défaut, fichiers .stk (0: fichiers .tif)
% par défaut:
% MTT23i({'*.stk','0',cd,'output23','24','28','7','0.1225','5','',...
% '1.1','1','-15','3.5','4','0','0.5','0.9','0','rel','1'})
%
% paramètres prédéfinis (D, ro...) pour:
% params = 1 => MCT
% params = 2 => STORMTT
%
% AS 18/12/7
% interface pour input params 12/6/8
% v2.3 ne sauve plus var et moy 3/12/9

function ok = MTT23i(params)

% % % path2version('MTT23')

%% define global parameters
include_global

global Nb_STK ;
global T ;
global T_off ;
global sig_free ;
global Boule_free ;
global Nb_combi ;
global Poids_melange_aplha ;
global Poids_melange_diff ;

if nargin==0, params = ''; end
% % if isnumeric(params), opt = params; else opt = -1; end% AS 2014, default chooice for future ana

ok = 1; % 2/5/17

[filename refitdata output_dir seuil_premiere_detec seuil_detec_1vue wn sig_free T Nb_STK r0 ... %%%% _def
    nb_defl T_off Boule_free Nb_combi seuil_alpha Poids_melange_aplha Poids_melange_diff SHOW conf_method ...
    OPTIM_R Tmin Dmin pxl_size time_lag use_gui name_param1 name_param2 params] = MTT23_dialog_box(params); %#ok
if isempty(params), ok = 0; return, end

%%--------------------------------------------------
%% data = .stk or .tiff??
use_stk = 1; %#ok
if ischar(filename), if strcmp(filename, '*.tif'), if isempty(dir('*.tif')), filename = '*.stk'; end, end, end % shift default from tif to stk, if no tif found 1/10/2013
if ischar(filename), if strcmp(filename(end-3:end),'.tif'), use_stk = 0; end, end %#ok

%% ini boucle sur fichiers stk, AS 26/3/7
% if use_stk
files = dir(filename);
Nfiles = length(files);
% else % liste de TIFF: toto_001.tif... AS 5/6/7
%     Nfiles = 1;
% end

if Nfiles==0, disp('ya qqun? (check dir ou testez ''clear all global'' ??)'), return, end % AS 21/5/7

if isempty(Nb_STK) && use_gui % isempty(params)
    disp('Did you/would you create mask (using masque.m) for the cells ?? (then press ctl+C)')
    %%%pause(2)
end

%% waitbar initialisation
waitbar_handle = waitbar(0,'Please wait, ...','Name', 'starting fit algo MTT SuperReknctor v2.3'); %if ~AFFICHAGE, end
waitbar_params{4} = clock; % t_start: temps de départ, pour calcul temps restant
waitbar_params{5} = zeros(1,Nfiles); % N_img_by_stk
 
if (SHOW > 0), AFFICHAGE = 1; figure('windowstyle','docked'), end %#ok

if isempty(dir(output_dir)), mkdir(output_dir), end
full_output_dir = ['.', filesep, output_dir] ;

save_params([name_param1; name_param2] , params, [params{4}, filesep, 'MTTparams.txt']'') % 24/2/9

%%%%%%%%%%%%%%%%%%%%%%
%% BOUCLE SUR FICHIERS
%%%%%%%%%%%%%%%%%%%%%%
for ifile = 1:Nfiles
    
    name_stk = files(ifile).name ; % filei ds v1.0...
    if strcmp(name_stk(1), '.'), continue, end
    
    if USE_MAT
        output_file_param = [full_output_dir, filesep, name_stk, '_tab_param.mat'] ;
    else
        output_file_param = [full_output_dir, filesep, name_stk, '_tab_param.dat'] ;
    end
    
    %% Test si fichier déjà complété, AS 27/3/7
    fid = fopen(output_file_param, 'rt','native') ;
    if fid>0 
        if ~refitdata
            value = fscanf(fid, '%d', 2) ;
            if ~isempty(value)
                disp([name_stk ' already done, skipped'])
                fclose(fid);
                continue
            end
        end
        fclose(fid);
    end
    
    waitbar_params(1:3) = {ifile, Nfiles, waitbar_handle};
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    detect_reconnex_23 %(name_stk, waitbar_params.....)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Nb_STK = []; %%%%%if ~use_gui, end % RAZ, AS 18/6/10

end % for ifile = 1:Nfiles

if ishandle(waitbar_handle), delete(waitbar_handle); pause(0.1), end

%% select traces long enough, fast enough, and save as .mat
sort_traces(filename, output_dir, Tmin, Dmin);

%% adding suppl. analyses AS 27/3/7

if ~strcmp(output_dir(end-2:end), 'out') && ~isempty(conf_method) % isempty(strfind(output_dir, 'DinDout')) && (Tmin > 1)% if refitdata==0 & 1 %&& Nfiles>1 && files_done<Nfiles

    disp('fit terminé, on passe aux cartes & autres histos')
    
    if strcmp(conf_method, 'coloc') || contains(conf_method, '2colors'), Nrun = 2; % 2 channels => 2 runs
    else, Nrun = 1;
    end
    
% %     %% fit Dv % ! requires Optimization Toolbox !
% %     amax = 40;
% %     for nc = 1:Nrun
% %         if (Nrun == 2) && (nc == 1), side = 'left';
% %         elseif (Nrun == 2) && (nc == 2), side = 'right';
% %         else, side = '';
% %         end
% %         disp('Sorting traces for D & v'), fit_directed_by_file(filename, amax, pxl_size, time_lag, side); % traces, sorted by v & D
% %     end
    
    %% histos
    for nc = 1:Nrun
        if (Nrun == 2) && (nc == 1), side = 'left';
        elseif (Nrun == 2) && (nc == 2), side = 'right';
        else, side = '';
        end
        
        disp('histos')
        ct4_by_file(filename, output_dir, conf_method, [], pxl_size, time_lag, side); % output_variables = []
        ct5_by_file(filename, output_dir, [], pxl_size, time_lag, side);
    end

    %% carto
    if ~contains(cd, 'split') % if isempty(strfind(cd, 'split')) % pb recalage origine xy si crop... COMPATIBILITY ISSUE??
        if strcmp(conf_method, 'nano')
            disp('nanopicts'), nanopictbyf(filename)
        elseif strcmp(conf_method, 'SCI')
            compil_SCI(filename)
        else
            disp('cartes 2D'), cartobyf(conf_method, filename, output_dir) % disp('cartes 3D'), carto3Dv3(filename)
        end
    end

    miniatures % for raw videos

end
%%%