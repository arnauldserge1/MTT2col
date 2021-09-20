function all_tab_param = load_all_param(file_specification, crop2max, MTT_dirname)

% function all_tab_param = load_all_param(file_specification, crop2max, MTT_dirname)
%
% load all MTT params from local dir, concatenated in a single table
% Caution: if different frame #, crop to max frames by default #!
%
% AS 17/07/2013

if nargin < 1, file_specification = '*.tif'; end
if nargin < 2, crop2max = 1; end
if nargin < 3, MTT_dirname = 'output23'; end % dirname = MTTparams_def{4};

% files = dir(file_specification);
files = dir([MTT_dirname filesep file_specification '_tab_param.mat']); % 3/5/17
Nfile = length(files);
all_tab_param = cell(1, Nfile);
n_frm = zeros(1, Nfile);

for nfile = 1:Nfile
%     filename = files(nfile).name;
%     filename_full = [MTT_dirname filesep filename '_tab_param.mat'] ;
    filename_full = [MTT_dirname filesep files(nfile).name];
    
    if ~isempty(dir(filename_full))
        all_tab_param{nfile} = importdata(filename_full);
    end
    
    [n_frm(nfile), nt_rc] = size(all_tab_param{nfile}); n_frm(nfile) = n_frm(nfile)/8;
    fprintf('%i frames & %i traces for file %s\r', n_frm(nfile), nt_rc, filename_full)
end

if crop2max
    t_max = max(n_frm);
    if (min(n_frm) < max(n_frm)), fprintf('Caution, padding with zeros to max = %i frames for concatenation!\r', t_max), end
    for nfile = 1:Nfile
        [l, c] = size(all_tab_param{nfile});
        all_tab_param{nfile} = [all_tab_param{nfile}; zeros(t_max*8-l, c)];
    end
else
    t_min = min(n_frm);
    if (min(n_frm) < max(n_frm)), fprintf(' Caution, cropping to min = %i frames for concatenation!\r', t_min), end
    for nfile = 1:Nfile
        all_tab_param{nfile} = all_tab_param{nfile}(1:t_min, :);
        tab_alpha = all_tab_param{nfile}(4:8:end, :);
        all_tab_param{nfile} = all_tab_param{nfile}(:, sum(tab_alpha)>0);
    end
end

all_tab_param = cell2mat(all_tab_param);

% imagesc(log(abs(all_tab_param)))
%%%