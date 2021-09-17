function CompilConf(files)

if nargin<1, files = '*.stk'; end
files = dir(files);
Nfile = length(files);
conf_method = 'rel';
graph = 0;
Tc = cell(Nfile,1);
Tf = cell(Nfile,1);

for ifile = 1:Nfile
    
    filename = files(ifile).name;
    timelag = eval_timing(filename);
    
    t_last = check_mem(filename);
%     if t_last<inf, sprintf('memory limit: keeping only up to %i frames', t_last), end
    tab_param = fread_params_timewindow(filename, 1, 1, t_last);
    
    %% confinement
    [Lc Tc{ifile} Tf{ifile}] = ... % R2c R2f trc_c trc_f
        MTT_probaconf(tab_param, conf_method, timelag, graph, filename);
    
end

Tcs = cell2mat(Tc)*timelag/1000;
Tfs = cell2mat(Tf)*timelag/1000;

figure('WindowStyle','docked')
do_log = 1;
subplot(2,1,1)
ct4plot_hist(Tcs, 'T_c (s)', do_log);
subplot(2,1,2)
ct4plot_hist(Tfs, 'T_f (s)', do_log);
%%%