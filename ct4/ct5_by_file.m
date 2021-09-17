function output_values = ct5_by_file(filename, dirname, output_variables, pxl_size, time_lag, side)

% function output_values = ct5_by_file(filename, dirname, output_variables, pixel_size, time_lag, side)
%
% lance ct5 sur chaque fichier vidéo,
% sauve fig en png et valeurs moyennes en txt

global redo
if isempty(redo), redo = 0; end

if nargin<1, filename = '*.stk'; end
if isempty(dir(filename)), filename = '*.tif'; end
if nargin<2, pdef = MTTparams_def; dirname = pdef{4}; end
if nargin<3, output_variables = []; end
if nargin<5, [pxl_size, time_lag] = get_calib3; end
if nargin<6, side = ''; end

if isempty(side), side_str = side; else, side_str = ['_' side]; end

outputdir = 'ct5linearity_by_file';
% % % % % if ~strcmp(dirname, 'output3D'), outputdir = ['ct4_by_file' dirname(10:end)]; end;

curdir = cd; findslash = strfind(curdir,filesep);
curdirname = curdir(findslash(end)+1:end);
savefile = [outputdir filesep curdirname '_ct5', dirname, '.txt'];

files = dir2(filename);
files = sort_nat({files.name});
Nfiles = length(files);

if isempty(files), disp('Hello! check dir, data, filename, weather forecast...'), return, end
if isempty(dir(outputdir)), mkdir(outputdir), end

output_values = cell(Nfiles,1);

%% ************ file by file ******
for nfile=1:Nfiles+1
    if nfile == Nfiles+1 % all together now!
        if isempty(side)
            filei = filename; % 'all files    '; % 4 spaces instead of '.stk'...            %         data = load_all_param;
            disp('all together now')
        else
            continue
        end
    else
        filei = files{nfile};
    end

    figname = [outputdir filesep filei(1:end-4) side_str '_ct5.png']; % set(gcf,'name',figname)
    if ~isempty(dir(figname)) && ~redo, disp([figname ' already done']), continue, end
    
    if nfile==1 && redo, fid = fopen(savefile,'w');
    else, fid = fopen(savefile,'a'); % ouvre et ferme à chaque fois => moins de soucis si interrompu
    end
    if (fid==-1), error('Couldn''t open file'); end
    
    %%%************%%%************%%%************%%%************%%%************%%%************%%%******%%%******%%%******
    [par, varname, output_values{nfile}] = ct5linearity(filei, time_lag, pxl_size, dirname, output_variables, side);
    %%%************%%%************%%%************%%%************%%%************%%%************%%%******%%%******%%%******
    
    if nfile==1
        %% *** print entêtes ***
        if ~isempty(files)
            fprintf(fid, '%s\t', date, dirname);
            fprintf(fid, '\r\n\t');
            varname = varname'; varname = varname(:)'; % remet les lignes bout à bout....
            for ii = 1:numel(varname), fprintf(fid, '%s\t', varname{ii}); end
            fprintf(fid, '\r\n');
        end
    end
    
    if ~isempty(par)
        fprintf(fid, '%s\t', [filei side_str]);
        par = par'; par = par(:)';
        for ii = 1:numel(par)
            fprintf(fid,'%g\t',par(ii));
        end
        fprintf(fid,'\r\n'); % retour ligne
        
        saveas(gcf, figname, 'png')
        delete(gcf) % pour éviter 36 fenetres!
    end
    
    fclose(fid);
end
%%%