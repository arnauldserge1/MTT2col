function carto_xyzt_byf(filename, color_code, dirname)

% function carto_xyzt_byf(filename, color_code, dirname)
% lance traj_xyzt.m sur chaque stk
% code couleur pour z, r2 ou t (ou rien)

global redo
if isempty(redo), redo = 0; end

if nargin<1, filename = '*.stk'; end
if nargin<2, color_code = 'z'; end
if nargin<3, dirname = 'output3D'; end

if strfind(cd,'output'), cd .., end % si on est déjà ds output..
files = dir(filename);
Nfiles = length(files);
if isempty(files), disp('ohé, y''a qqchose??'), return, end

dir_carto = 'carto';
% % % if ~strcmp(dirname, 'output3D'), dir_carto = ['carto_' dirname(10:end)], end;

min_trc_length = 3;

for ifile=1:Nfiles
    filei = files(ifile).name;
    
    outfile3D = [dir_carto '/' filei '_' color_code '_3D.png'];
    outfile2D = [dir_carto '/' filei '_' color_code '_2D.png'];
    if ~isempty(dir(outfile3D)) && ~redo, disp ([outfile3D ' already done']), continue, end
    
    %% ************************************************
    ok = traj_xyzt(filei, color_code, min_trc_length, dirname);
    %% ************************************************
    
    if ~ok, continue, end
    
    if isempty(dir(dir_carto)), mkdir(dir_carto), end
    disp(['saving ' outfile3D])
    saveas(gcf, outfile3D, 'png') % saveas(gcf,[outfile3D '.fig'],'fig') trop lourdingue...
    
    view(0,89.9), drawnow % pb avec view(2).....
    saveas(gcf, outfile2D, 'png')
    
    if ifile<Nfiles, close gcf, else view(3), end % Nfiles>1 pour économ. mémoire, si plusieurs fichiers (on garde la derniere par défaut, pour ROI éventuelle...)
end % ifile
%%%