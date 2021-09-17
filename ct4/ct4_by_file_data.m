function par = ct4_by_file_data(filename, dirname, conf_method, output_variables, pxl_size, time_lag)

%  par = ct4_by_file_data(filename, dirname, conf_method, output_variables, pixel_size, time_lag)
%
% lance ct4 sur chaque fichier vidéo, get params


pdef = MTTparams_def;

if nargin<1, filename = '*.tif'; end
if nargin<2, dirname = pdef{4}; end
if nargin<3, conf_method = 'coloc'; end
if nargin<4, output_variables = []; end
if nargin<6, [pxl_size, time_lag] = get_calib3; end % if nargin<5, pxl_size = str2double(pdef{24}), fprintf('Using default pixel size: %g um\r', pxl_size), end % 0.16 if nargin<6, time_lag = str2double(pdef{25}); fprintf('Using default time lag: %g s\r', time_lag), end  %   timelag = eval_timing(filename);

files = dir(filename);
files = sort_nat({files.name});
Nfiles = length(files);

par = zeros(Nfiles, 10, 2);

for ns = 1:2
    if ns == 1, side = 'left'; else, side = 'right'; end
    
    for nf = 1:Nfiles
        filei = files{nf};
        results = ct4units(filei, time_lag, pxl_size, conf_method, dirname, output_variables, side);
        if ~isempty(results), par(nf, :, ns) = results; end
        close(gcf)
    end
end

%%%