% % % function ct4_by_dir

% function ct4_by_dir
%
% run ct41 over all subdirs
% (ie pooling data for each experimental condition)
%
% AS 17/07/2013

% cd('/Users/serge/Documents/labo/Geneve/Alexandre')
if ispc, cd('V:\ARNAULD'), else cd('/Volumes/IMAGERIE/Arnauld'), end

dirs = dir;
ok = zeros(size(dirs));

for ndir = 3:length(dirs)
    ok(ndir) = dirs(ndir).isdir;
end

dirs = dirs(logical(ok));
Ndir = length(dirs);
    
var_name = {'alpha (a.u.)' 'SNR (db)' 'r^2 (pxl^2)' 'TrcLen (img)' 'asym' 'Npk/frm' 'Ntrc'};
% var_name = {'alpha (a.u.)' 'SNR (db)' 'r^2 (pxl^2)' 'TrcLen (img)' 'asym' ...
%     'accuracy (pxl)' 'D by trc (pxl^2/img)' 'Dglobal (pxl^2/img)' 'gamma'  'Npk/frm' 'Ntrc'};
Nvar = size(var_name, 2);
mean_val = zeros(Ndir, Nvar);

for ndir = 1:Ndir
    dirname = dirs(ndir).name;
    cd(dirname)
    
    all_tab_param = load_all_param;
    if ~isempty(all_tab_param)
        mean_val(ndir, :) = ct41(all_tab_param);
    end
    
    cd ..
end

nSubImg = 9%size(var_name, 2);
Golden = (1+sqrt(5))/2;
nSubX = 3%floor(sqrt(nSubImg*Golden)); 
nSubY = ceil(nSubImg/nSubX);

%% histos
figure('WindowStyle', 'docked');

for nplot = 1:nSubImg
    subplot(nSubX, nSubY, nplot)
    barh(mean_val(:, nplot))
    title(var_name{nplot})
    axis tight
    if mod(nplot,nSubY)==1, set(gca,'YTick',1:31,'YTickLabel',{dirs.name,'FontSize',6}), end
end

c = cell(size(dirs));
i = cell(size(dirs));
for ndir = 1:31
    if strfind((dirs(ndir).name),'confl')
        c{ndir} = mean_val(ndir, :);
    elseif strfind((dirs(ndir).name),'isol')
        i{ndir} = mean_val(ndir, :);
    end
end
c = cell2mat(c);
i = cell2mat(i);
ci = [c;i];
for nplot = 1:nSubImg
    subplot(nSubX, nSubY, nplot)
    bar(ci(:, nplot))
    title(var_name{nplot})
end
%%%