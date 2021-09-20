function data_mat = cell2mat_extended(data_cell, use_nan)

% function data_mat = cell2mat_extended(data_cell, use_nan)
%
% Convert from cell to mat:
% create a matrix and asign each element of data_cell to a column of data_mat
% output size is thus max(numel(data_cell)) by length(data_cell)
% if elements are of different size, matrix is completed by zeros (or NaN)
% AS

if nargin < 2, use_nan = 0; end

data_cell = data_cell(:);

nc = length(data_cell);
sz = zeros(size(data_cell));

for i=1:nc
    sz(i) = numel(data_cell{i});
end

if use_nan, data_mat = nan(max(sz), nc);
else, data_mat = zeros(max(sz), nc);
end

for i=1:nc
    data_mat(1:sz(i),i) = data_cell{i};
end
%%%