function output_values = ct4_by_file(filename, dirname, conf_method, output_variables, pxl_size, time_lag, side)

% function output_values = ct4_by_file(filename, dirname, conf_method, output_variables, pixel_size, time_lag, side)
% default: output_values = ct4_by_file('*.stk', 'output23', '', [], 0.1375, 0.14, '')
%
% lance ct4 sur chaque fichier vidéo,
% sauve fig en png et valeurs moyennes en txt

global redo
if isempty(redo), redo = 0; end

pdef = MTTparams_def;

if nargin<1, filename = '*.stk'; end
if isempty(dir(filename)), filename = '*.tif'; end

if nargin<2, dirname = pdef{4}; end

if nargin<3, conf_method = ''; end
if strcmp(conf_method, 'nano'), conf_method = ''; end % ou histos dédiés??
if any(strcmp(conf_method, {'abs', 'rel'})), Nfigs = 2; else, Nfigs = 1; end % global & conf

if nargin<4, output_variables = []; end
if nargin<6, [pxl_size, time_lag] = get_calib3; end % if nargin<5, pxl_size = str2double(pdef{24}), fprintf('Using default pixel size: %g um\r', pxl_size), end % 0.16 if nargin<6, time_lag = str2double(pdef{25}); fprintf('Using default time lag: %g s\r', time_lag), end  %   timelag = eval_timing(filename);
if nargin<7, side = ''; end

if isempty(side), side_str = side; else, side_str = ['_' side]; end

outputdir = 'ct4_by_file';
if ~strcmp(dirname, 'output3D'), outputdir = ['ct4_by_file' dirname(10:end)]; end

curdir = cd; findslash = strfind(curdir,filesep);
curdirname = curdir(findslash(end)+1:end);
savefile = [outputdir filesep curdirname '_Ct4', dirname, '.txt'];

files = dir2(filename);
files = sort_nat({files.name});
Nfiles = length(files);

if isempty(files), disp('Hello! check dir, data, filename, weather forecast...'), return, end
if isempty(dir(outputdir)), mkdir(outputdir), end

figname = cell(Nfigs,1);
output_values = cell(Nfiles,1);

%% ************ file by file ******

for nfile = 1:Nfiles+1
    
    if nfile == Nfiles+1 % all together now!
% % %         if isempty(side) % removed 19/5/2017
            filei = filename; % 'all files    '; % 4 spaces instead of '.stk'...   % data = load_all_param;
            disp('all together now')
% % %         else
% % %             continue
% % %         end
    else
        filei = files{nfile};
    end
    % % %     timelag_ms = eval_timing(filei);
    
    for ff = Nfigs:-1:1 % if Nfigs == 1, figname{1} = [outputdir filesep filei(1:end-4) side '_ct4.png'];???
        figname{ff} = [outputdir filesep filei(1:end-4) side_str '_ct4_' num2str(ff) '.png'];
        figname{ff}(strfind(figname{ff}, '*')) = '@'; % all files: expecting a name such as '*.tif' => '@.tif_ct4_ff.png'
    end
    if ~isempty(dir(figname{Nfigs})) && ~redo, disp([figname{Nfigs} ' already done']), continue, end
    
    if (nfile == 1) && redo, fid = fopen(savefile,'w');
    else, fid = fopen(savefile,'a'); % ouvre et ferme à chaque fois => moins de soucis si interrompu
    end
    if (fid == -1), error('Couldn''t open file'); end
    
    %%%************%%%************%%%************%%%************%%%************%%%************%%%******%%%******%%%******
    [par, varname, output_values{nfile}] = ct4units(filei, time_lag, pxl_size, conf_method, dirname, output_variables, side);
    %%%************%%%************%%%************%%%************%%%************%%%************%%%******%%%******%%%******
    
    if (nfile == 1)
        %% *** print entêtes ***
        if ~isempty(files)
            fprintf(fid, '%s\t', date, dirname);
            fprintf(fid, '\r\nfile\t');
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
        
        for ff = Nfigs:-1:1
            saveas(gcf, figname{ff}, 'png')
            delete(gcf) % pour éviter 36 fenetres!
        end
    end
    
    fclose(fid);
end
%%%