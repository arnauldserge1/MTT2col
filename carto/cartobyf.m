function cartobyf(codage, filename, dirname)%, show_colorbar)%useDV)%, useDIC)

% cartobyf(codage, filename, dirname, show_colorbar)
%
% lance trajxyt.m sur chaque stk (selon filename, '*.stk' par d?f.)
%
% !! la syntaxe du tif correspondant peut etre a ajuster !!
% convention usuelle pour nom de l'image DIC (ou trans)
% blabla cell1a.stk => DIC\blabla cell1.tif
%
% see also MTT23i, dicname

% 29/6/6


global redo
if isempty(redo), redo = 0; end

%% criteria for coloc
global coloc_dist_max coloc_time_min

if nargin<1, codage = 'coloc_magenta'; end %%%'SCI' 'Dv' 'rel' 'abs' 'speed' 'int' ''
if nargin<2, filename = '*.tif'; end
if nargin<3, params_default = MTTparams_def; dirname = params_default{4}; end% if nargin<4, show_colorbar = 0; end

if contains(codage, 'coloc')
    if isempty(coloc_dist_max), coloc_dist_max = 2; end
    if isempty(coloc_time_min), coloc_time_min = 2; end
end

DIC_folder = 'DIC';
if isempty(dir(DIC_folder)), DIC_folder = 'trans'; end % 17/2/21

%% registration of the 2 channels
searchType = 'file';
if ~isempty(dir(DIC_folder)) && exist('imregtform', searchType)%exist('imregtform', 'builtin')????? 10/5/2017
    string_reg = 'reg';
    cd(DIC_folder), reg_2_colors; cd .. % mean reg from all files
    % % elseif contains(cd, 'EQ_Charafe_Ginestier') || contains(cd, 'Violette')% %     string_reg = '';
else
%     disp(['no ' DIC_folder ' folder or no imregtform function, no registration possible, so no 2-color cartography..'])%     string_reg = '';
    disp(['Caution : no ' DIC_folder ' folder or no imregtform function, no registration...'])%     string_reg = '';
    string_reg = 'no reg';
%     return % 21/2/18 % disabled 1/3/2021, for using single camera, sequentially
end

if ~any(strcmp(codage, {'', 'SCI', 'Dv', 'rel', 'abs', 'speed', 'int', 'time', 'coloc', 'coloc_magenta', 'rainbow2', 'Dv_2colors', 'anomal_2colors'}))
    disp(['unknokn codage ' codage ' for carto'])
    return
end

if contains(cd,'output'), cd .., end % si on est deja ds output..

%% get files
files = dir2(filename);
Nfiles = length(files);
if isempty(files), disp('oh?, y''a qqchose??'), return, end

%% output dir
dir_carto = 'carto';
if ~isempty(codage), dir_carto = [dir_carto '_' codage]; end
if strcmp(codage, 'coloc'), dir_carto = [dir_carto '_' num2str(coloc_dist_max) string_reg]; end
if strcmp(codage, 'coloc_magenta'), dir_carto = [dir_carto '_' num2str(coloc_dist_max) string_reg]; end
% if strcmp(codage, 'Dv_2colors'), dir_carto = [dir_carto '_' num2str(coloc_dist_max) string_reg]; end
if isempty(dir(dir_carto)), mkdir(dir_carto), end

% % % figure

%% loop over files
for nf = 1:Nfiles
    file = files(nf).name;
    
    outfile = [dir_carto filesep file(1:end-4) '.png'];
    if ~isempty(dir(outfile)) && ~redo, disp ([outfile ' already done']), continue, end
    
    %% get D & v or D & gamma values if required
    D = []; v = []; gamma = [];
    if contains(codage, 'Dv')% strcmp(codage, 'Dv')
        dataname = ['Brownian_vs_directed' filesep file(1:end-4) ' trc sort data'];
        if ~isempty(dir([dataname '.mat']))
            s = load([dataname '.mat']);
            D = s.D; v = s.v;
        else % two sides, left, right
            sl = load([dataname 'left.mat']);
            sr = load([dataname 'right.mat']);
            D = [sl.D; sr.D]; v = [sl.v; sr.v];
        end
    elseif contains(codage, 'anomal')
        dataname = ['ct4_by_file' filesep file(1:end-4)  '_D_gamma'];
        if ~isempty(dir([dataname '.mat']))
            s = load([dataname '.mat'], 'Da', 'gamma');
            D = s.Da; gamma = s.gamma;
        elseif ~isempty(dir([dataname '_left.mat']))% two sides, left, right
            sl = load([dataname '_left.mat'], 'Da', 'gamma');
            sr = load([dataname '_right.mat'], 'Da', 'gamma');
            D = [sl.Da(:); sr.Da(:)]; gamma = [sl.gamma(:); sr.gamma(:)];
%         else [D, gamma] = fit_anomal2...
        end
    end
    
    ok = 1;
% % %     clf % figure
    
    %% *************** carto, with selected option *********************
    if strcmp(codage, 'Dv') || strcmp(codage, 'SCI')
        ok = traj_xyc(file, codage, dirname, D, v);
    elseif contains(codage, 'coloc') % meaning 2 colors SPT: red, green & coloc, color coded for time
        [N_coloc, coloc_duration] = traj_xy_coloc(file, dirname, string_reg); % codage = diff py def
    elseif strcmp(codage, 'Dv_2colors') % meaning 2 colors SPT: red, green & coloc, with color code for Dv values
        [N_coloc, coloc_duration] = traj_xy_Dv_2colors(file, D, v);%, dirname, string_reg);
    elseif strcmp(codage, 'anomal_2colors') % meaning 2 colors SPT: red, green & coloc, with color code for D & gamma values (anomal diff)
        [N_coloc, coloc_duration] = traj_xy_anomal_2colors(file, D, gamma);
    else
        traj_xyc2(file, codage, dirname) % coded for speed, int, time or rainbow
    end
    %% ********************************************************************
    
    %% save coloc duration
    if contains(codage, 'coloc') || strcmp(codage, 'Dv_2colors') || strcmp(codage, 'anomal_2colors')
        if (nf == 1) && redo, wmode = 'w'; else, wmode = 'a'; end
        fid = fopen([dir_carto filesep 'coloc_N_duration.txt'], wmode); % ouvre et ferme a chaque fois => moins de soucis si interrompu
        if (nf == 1) %% *** create txt file & print entete ***
            fprintf(fid, '%s\t\r\n', date);  % coloc_dist???
            fprintf(fid, '%s\t%s\t%s\t\r\n', 'file', 'N coloc', 'duration of each coloc (frames)');
        end
        if (fid == -1), disp('Couldn''t open file'); continue, end % error('Couldn''t open file')???
        
        fprintf(fid, '%s\t%g\t', file, N_coloc);
        fprintf(fid, '%g\t', coloc_duration); % 0, 1 or more, according to N_coloc
        fprintf(fid, '\r\n'); % retour ligne
        fclose(fid);
    end
    
    %% save carto image
    if ok
        % %         a = axis;% %         w = a(2) - a(1);% %         h = a(4) - a(3); % w = h usually, 2h for two cameras
        % %         if show_colorbar, colorbar, w2 = 2; else w2 = 0; end
        % %         h_out = 20; % setting h = 20cm
        % %         set(gcf,'PaperSize',[w*h_out/h+w2 h_out+1],'PaperPosition',[0 0 w*h_out/h+w2 h_out+1],'units','centimeters') % +1cm for title
        % %         set(gca,'units','centimeters','position',[0 0 w*h_out/h h_out])
        disp(['saving ' outfile])
        if strcmp(codage, 'rainbow2')
            saveas(gcf, outfile, 'png')
        else
            title('')
            set(gcf,'Renderer','OpenGL')
            im = print('-RGBImage', '-r600', '-painters'); % print(gcf, '-dpng', '-r600', outfile)  %      no_save = 1
            s = sum(im,3)/255/3;
            l = find(min(s,[],2) < 1); % lines for image, w/o border to crop only around axes
            c = find(min(s,[],1) < 1);
            im2 = im(l(1):l(end), c(1):c(end), :);
            imwrite(im2, outfile)
        end
%         if (Nfiles > 1), close gcf, end
    end
end % ifile
close gcf
%%%