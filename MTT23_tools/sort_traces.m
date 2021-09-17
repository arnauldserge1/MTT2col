function [n0, n1] = sort_traces(filename, dirname, Tmin, Dmin, two_colors, redo_sort)

% function [n0, n1] = sort_traces(filename, dirname, Tmin, Dmin, two_colors, redo_sort)
% remove traces below Tmin (def = 5) or Dmin (def = 0.001)
% and save as mat file (for fast loading)
% n0/n1 = number of peaks before/after sort, resp.
% AS 2013

global N_PARAM PARAM_ALPHA %%% redo

params_def = MTTparams_def;
if nargin < 1, filename = params_def{1}; end
if nargin < 2, dirname = params_def{4}; end
if nargin < 3, Tmin = str2double(params_def{22}); end % cf diff_window = 5; for Dinst
if nargin < 4, Dmin = str2double(params_def{23}); end % = 0.001 pxl2/img = 2.56 um2/s log10(Dmin) = -3.6
if nargin < 5, two_colors = 0; end
if nargin < 6, redo_sort = 0; end %%%if isempty(redo), redo = 0; end

files = dir2(filename);
Nf = length(files);

if two_colors == 0
    Nc = 1;
    side = {''};
else % for Violette
    Nc = 2;
    side = {'left' 'right'};
end
n0 = zeros(Nf, Nc);
n1 = zeros(Nf, Nc);

for nf = 1:Nf
    filei = files(nf).name;
    filename_out = [dirname filesep filei '_tab_param.mat'];
    
    if ~isempty(dir(filename_out)) && ~redo_sort %%%%% && ~two_colors ?????
        disp([filename_out ' already sorted out & saved'])
    else
        for nc = 1:Nc
            tab_param = fread_all_params([dirname filesep filei '_tab_param.dat']); % tab_param_in
        
            im1 = imread(filei,1); middle = size(im1,2)/2;
            tab_param = split_params_left_right(tab_param, middle, side{nc}); % for Violette
            %% Tmin
            if (Tmin < 1) % ratio, not abs value, cf cell_track
                min_length_ratio = Tmin;
                Tmax = size(tab_param, 1)/N_PARAM;
                Tmin = round(Tmax*min_length_ratio);
            end
            
            alpha = tab_param(PARAM_ALPHA-1:N_PARAM:end, :); % intensity of fitted SM
            TrcLen = sum(alpha > 0, 1)'; % number of images in each trace, AS 11/12/7
            tab_param = tab_param(:, TrcLen >= Tmin); % tab_param_out
            fprintf('keeping %i traces among %i, with at least %i step(s) \n', sum(TrcLen >= Tmin), length(TrcLen), Tmin)
            
            %% Dmin
            if Dmin > 0
                trc = detect_reconnex_to_trc(tab_param);
                msd_by_file = msd(trc);
                D = calculDinst(msd_by_file);
                tab_param = tab_param(:, D(:,2) >= Dmin); % tab_param_out
                fprintf('keeping %i traces with D > %3g pxl2/img \n', size(tab_param,2), Dmin)
            end
            
            if ~two_colors, save(filename_out, 'tab_param'), end
            
            n0(nf, nc) = length(TrcLen);
            n1(nf, nc) = size(tab_param,2);
        end
    end
end