function [nb, r0, D, v, chi2] = fit_directed_by_file(filename, amax, pxl_size, time_lag, side)

% function [nb, r0, D, v, chi2] = fit_directed_by_file(filename, amax, pxl_size, time_lag, side)
%
% default: fit_directed_by_file('*.tif', 40, 0.16, 0.1, '') (see MTTparams_def)
% run fit_directed_motion4 over files:
% fit msd(t) = 2r02 + 4Dt + v2t2 and sort traces according to D & v values,
% quasi null or relevant => immobile, directed, Brownian or fast
%
% AS 2014
% see also fit_directed_motion4, MTT23i

global redo
if isempty(redo), redo = 0; end

if nargin<1, params_default = MTTparams_def; filename = params_default{1}; end % *.tif
if nargin<2, amax = 40; end % axis size, in pxl 
if nargin<3, params_default = MTTparams_def; pxl_size = eval(params_default{24}); end % if nargin<5, [pxl_size, time_lag] = get_calib3; end ????
if nargin<4, params_default = MTTparams_def; time_lag = eval(params_default{25}); end
if nargin<5, side = ''; end

do_plot = 0;

files = dir2(filename);
if isempty(files), files = dir2('*.stk'); end
if strcmp(filename(end-2:end), 'stk'), amax = 20; end % => Claire, PCM1

Nfiles = length(files);

dir_out = 'Brownian_vs_directed';
if isempty(dir(dir_out)), mkdir(dir_out), end

curdir = cd; findslash = strfind(curdir,filesep);
curdirname = curdir(findslash(end)+1:end);
savefile = [dir_out filesep curdirname '_trc_sorting.txt'];

for i = 1:Nfiles
    filei = files(i).name; % & all files??
    
    outfile = [dir_out filesep filei(1:end-4) ' trc sort data' side '.mat']; % % % [dir_out filesep filei side '_trc.png'];
% % %     aaa = 'erase!!'
    if ~isempty(dir(outfile)) && ~redo
        disp ([outfile ' already done'])
        continue
    end
    
    %% ******************************************************************
    [nb, r0, D, v, chi2, v_threshold, D_threshold]...
        = fit_directed_motion4(filei, do_plot, amax, pxl_size, time_lag, side);
    %% ******************************************************************
    %% OR D, gamma = fit_anomal() ?????????? %%
    %% ******************************************************************

    %% save vm Dm %
    fid = fopen(savefile,'a'); % ouvre et ferme � chaque fois => moins de soucis si interrompu
    if i==1
        fprintf(fid, '%s\t', 'file', '% directed', '% fast D,v', '% immobile', '% Brownian', 'N total', 'D mean (um2/s)', 'v mean(um/s)');
        fprintf(fid,'\r\n'); % retour ligne
    end
    
    fprintf(fid, '%s\t', [filei side]);
    Dm = mean(D > D_threshold)*pxl_size^2/time_lag;
    vm = sqrt(mean(v > v_threshold))*pxl_size/time_lag;
    perc = nb(1:4) ./ nb(5);
    par = [perc, nb(5), Dm, vm];
    for ii = 1:numel(par)
        fprintf(fid,'%g\t',par(ii));
    end
    fprintf(fid,'\r\n'); % retour ligne
    fclose(fid);
    
    disp(['saving ' outfile])
    save(outfile, 'nb', 'r0', 'D', 'v', 'chi2') %  save([dir_out filesep filei(1:end-4) ' trc sort data' side '.mat'], 'nb', 'r0', 'D', 'v', 'chi2')  

% % %     if do_plot, print(gcf, '-dpng', outfile), close gcf, end % print(gcf, '-dpng', '-r600', outfile) % saveas(gcf, outfile) % def resol 150
end % ifile

%%%