function MTT_2_colors(params, JAMC_side)

%% function MTT_2_colors(params, JAMC_side)
%
% run MTT,
% detect contour of contact from trajs of given side: 'left' or 'right'
% and save associated trajs for the other side
% (default, 'right', assuming JAM-C in red, hence selecting JAM-B trajs, in 'lef't)
%
% params: -1 for default, empty for user interface, or detailled
%
%% See article:
%
%% Single molecule analysis unveils the dynamics of leukemic stem/bone-marrow stromal cell contacts
% Oksana Gorshkova1, Jessica Cappaï1 , Lorianer Maillot2 and Arnauld Sergé1,2*
% 1 Centre de Recherche en Cancérologie de Marseille, Institut Paoli-Calmettes, Inserm U1068, CNRS UMR7258, Aix-Marseille Université UM105, France.
% 2 Laboratoire Adhésion et Inflammation, Inserm U1067 CNRS UMR 7333, Aix-Marseille Université, France.
% * Correspondence should be addressed to A.S. (arnauld.serge@univ-amu.fr).
%
%% (c) Arnauld Sergé 2020
%
% see also MTT23i, detect_cell_contact


global coloc_dist_max coloc_time_min
if isempty(coloc_dist_max), coloc_dist_max = 2; end
if isempty(coloc_time_min), coloc_time_min = 2; end

if nargin < 1, params = ''; end
if nargin < 2, JAMC_side = 'right'; end % default: JAM-C on KG1 in red, so right image, and JAM-B on MS5 in green, left image

% % % if (params == 0), params = ''; end % to be removed??
% % % if (params == 0)
% % %     params = MTTparams_def;
% % %     params{20} = 'coloc';
% % %     disp('using coloc for method, not default for one color(?)')
% % % end 

if isempty(params)
    [filename, refitdata, output_dir, seuil_premiere_detec, seuil_detec_1vue, wn, sig_free, T, Nb_STK, r0, ... %%%% _def
        nb_defl, T_off, Boule_free, Nb_combi, seuil_alpha, Poids_melange_aplha, Poids_melange_diff, SHOW, ~, ...
        OPTIM_R, Tmin, Dmin, pxl_size, time_lag, use_gui, name_param1, name_param2, params] = MTT23_dialog_box(params); %#ok
end

if isempty(params), return, end % cancelled by user

params_default = MTTparams_def;
if length(params)<length(params_default)
    params = [params(:); params_default(length(params)+1:end)];
end

%% prepare_coloc % nota bene: requires Matlab R2014 or later
if ~strcmp(JAMC_side, 'right') && ~strcmp(JAMC_side, 'left'), disp('whats that side??'), return, end

curdir = cd;
if strcmp(curdir(end-1:end), [filesep 'z']), disp('caution, running from z folder, going up to main one!!'), cd .., end

%% %%%%%% MTT %%%%%%%
if iscell(params)
% % %     params{20} = 'coloc'; % "conf detec" method
    %%%%%%%%%%%%%%%%%%%%
    MTT23i(params); %% detect_reconnex_23, carto & histos
    %%%%%%%%%%%%%%%%%%%%
elseif (params == -1)
    params = MTTparams_def;
    params{20} = 'coloc';
    disp('using coloc for method, not default for one color(?)')
    %%%%%%%%%%%%%%%%%%%%
    ok = MTT23i(params); % coloc => ct4, ct5, fit_directed_by_file for left & right + carto_coloc, with Red-Green distances & coloc_N_duration
    %%%%%%%%%%%%%%%%%%%%
    if ok == 0, return, end
end
%%%%%%%%%%%%%%%%%%%

merge_colors % & max

if isfolder('dic - Copy')
    detect_cell_contact2 % detect for contact def with ImageJ black circles
else
    detect_cell_contact(JAMC_side); % stats (ct4, ct5, fit_directed_by_file) only at cell contact
end

build_mosaic('carto')
if strcmp(params{20}, 'Dv_2colors'), build_mosaic('carto_Dv_2colors'), end
build_mosaic('max')

%% z
if isfolder('z'), zdir = 'z';
elseif isfolder('Z'), zdir = 'Z';
else, zdir = '';
end

if ~isempty(zdir)
    cd(zdir)
    params = MTTparams_def; % 3/5/17
    params{20} = ''; % conf detec method...
    params{22} = '3'; % "Tmin" => Zmin
    params{23} = '0'; % Dmin
    
    %%%%%%%%%%%%%%%%%%%
    MTT23i(params); % MTT_z ???
    %%%%%%%%%%%%%%%%%%%
    gradient_z2; %(JAMC_side)?
    merge_colors % & max
    
    cd ..
end
%%%