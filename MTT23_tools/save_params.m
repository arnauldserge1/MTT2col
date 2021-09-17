function save_params(name_parameters, params, name)

if nargin <3, name = 'params.txt'; end

%   while fid==-1,
fid = fopen(name, 'wt'); %, 'wt','native');
%       pause (.01)
%   end

fprintf(fid, '%s @ %s ', name, date) ;
clk = clock() ;
fprintf(fid, '%.2dh%.2dm%.2ds\n\n', clk(4), clk(5), round(clk(6))) ;

for i=1:length(params)
    if iscell(params{i})
        params{i} = cell2mat(params{i});
    end
    if ischar(params{i})
        template = '%s : %s\n';
    else
        template = '%s : %g\n';
    end
    fprintf(fid, template, name_parameters{i}, params{i}) ;
end

fclose(fid) ;
