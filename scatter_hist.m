function [mean_data, sd_data, data_wo_outliers] = scatter_hist(data, x_tags, y_tag, data_pos, do_log, use_sem, color, x_offset, remove_outliers, show_n, color_mean_sd, symbol, plot_bar)

% function [mean_data, sd_data, data_wo_outliers] = scatter_hist(data, x_tags, y_tag, data_pos, do_log, use_sem, color, x_offset, remove_outliers, show_n, color_mean_sd, symbol, plot_bar)
% data sets ranged in columns, with associated names (and common y_tag)
% i.e scatter_hist(data(:, [22 23 7]), {'D_{free}', 'D_{conf}', 'D (\mum^2/s)'})
%
% default: scatter_hist(data, {1:N}, '', 0, 0, 1, -1, 0, 0, 0, 'k', {'o', 's', '*', 'd', '+', 'p', 'h', '.', 'x', '^', 'v', '>', '<'}, 0)
% def colors: gradual_blue
% AS 31/3/9

global outlier_method legend_keys
if isempty(outlier_method), outlier_method = 'med'; end

if nargin < 2, x_tags = num2cell(1:size(data, 2)); end
if nargin < 3, y_tag = ''; end
if nargin < 4, data_pos = 0; end
if nargin < 5, do_log = 0; end % 0: 21/6/2016
if nargin < 6, use_sem = 1; end
if nargin < 7, color = -1; end
if nargin < 8, x_offset = 0; end
if nargin < 9, remove_outliers = 0; end
if nargin < 10, show_n = 0; end
if nargin < 11, color_mean_sd = 'b'; end
if nargin < 12, symbol = {'o', 's', '*', 'd', '+', 'p', 'h', '.', 'x', '^', 'v', '>', '<'}; end
if nargin < 13, plot_bar = 0; end

if iscell(data)
    data_wo_outliers = cell(size(data));
    use_nan = 1;
    data = cell2mat_extended(data, use_nan);
else
    data_wo_outliers = zeros(size(data));
end % nan: 21/6/2016
if ~iscell(symbol), symbol = {symbol}; end

N_data_sets = size(data, 2);
% % % if data_pos% % %     for n = 1:N_data_sets% % %         ind = find(data(:, n)>0, 1);% % %         if isempty(ind)% % % % % %             data(:, n) = []; % 10/12/13% % %             N_data_sets = N_data_sets-1;% % %         end% % %     end% % % end

scatter_width = 0.1; % or .3 for rand??
% cla reset
hold on
set(gca, 'TickDir', 'out')
if do_log, set(gca, 'YScale', 'log'), else, set(gca, 'YScale', 'lin'), end
ylabel(y_tag)%, 'interpreter' ,'tex')%, 'interpreter' ,'none')

ind = 1:size(data, 1);
yy = cell(N_data_sets,1);

mean_data = nan(N_data_sets, 1);
sd_data = nan(N_data_sets, 1);
legend_keys = zeros(N_data_sets, 1);

for n = 1:N_data_sets
    if data_pos, ind = find(data(:, n)>0); end
    yy{n} = data(ind, n);
    
    if (color == -1)
        blue_level = n/N_data_sets*0.8; c = [0 0 blue_level];
    else
        if isnumeric(color) && (size(color, 2) > 1)
            c = color(n,:); 
        elseif ischar(color) && (length(color) > 1)%size(color, 1) == 1 && size(color, 2) > 1
            c = color(n); 
        else
            c = color;
        end % one or several colors for input?
        if nargin < 12, symbol = {'.'}; end%{'o'};
    end
    if size(color_mean_sd, 1) > 1, cc_mean_sd = color_mean_sd(n,:); else, cc_mean_sd = color_mean_sd; end
    
    n2 = mod(n-1, length(symbol))+1;
    
    if remove_outliers
        TF = isoutlier(yy{n}, outlier_method);%, 'gesd');%, 'grubbs');% 'mean');
        if remove_outliers < 2
            yy_outliers = yy{n}(TF);
            xx_outliers = n+(rand(length(yy_outliers), 1)-0.5)*scatter_width+x_offset;
            plot(xx_outliers, yy_outliers, 'x', 'color', 'r')%symbol{n2}, 'color', 'r')
        end
        yy{n} = yy{n}(~TF);
        
        if iscell(data_wo_outliers)
            data_wo_outliers{n} = yy{n}; % data(ind, n);
        else
            data_wo_outliers(1:length(yy{n}), n) = yy{n}; % data(ind, n);
        end
    end
    randxx = randn(length(yy{n}), 1);
    randxx = (mod(randxx-2, 4))-2; % rnormal distrib, but restricted to -2 2
    xx = n+randxx*scatter_width+x_offset; % xx = n+(rand(length(yy{n}), 1)-0.5)*scatter_width+x_offset; % xx = n+randn(length(yy{n}), 1)*scatter_width+x_offset; put rand, not randn, 4/4/18
    
    plot(xx, yy{n}, symbol{n2}, 'color', c)
    
    % default color for mean & sd: black, gray, red, same as for dots??
    %     if (clr == -1), c_mean_sd = c+0.2
    %     else
    % if size(clr, 1) > 1, c_mean_sd =.........???
    %         c_mean_sd = 'k';%[1 1 1]*0.3;%if size(clr, 1) > 1, c2 = clr(n,:); else, c2 = clr; end % one or several colors for input?
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [mean_data(n), sd_data(n), legend_keys(n)] = plot_mean_err(n+x_offset, yy{n}, cc_mean_sd, scatter_width, do_log, use_sem, plot_bar);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

a = axis;
if do_log, tag_offset = a(3) * (a(3)/a(4))^(1/10);
else, tag_offset = a(3) + (a(3)-a(4))/12;
end

if show_n
    for n = 1:N_data_sets
        % % %         if remove_outliers        % % %             n_tot = sum(~isnan(yy));%length(TF);        % % %             n_left = sum(~TF | ~isnan(data(ind, n)));%length(yy);        % % %             text(n, -20, sprintf('n=%i/%i', n_left, n_tot))        % % %         else
        n_tot = sum(~isnan(yy{n}));
% % %         xlabel({x_tags{n}, n_tot})
%         text(n+0.25, tag_offset, sprintf('n=%i', n_tot))
        text(n, tag_offset, sprintf('n=%i', n_tot), 'HorizontalAlignment', 'center')
        % % %         end
    end
end
if N_data_sets > 0, set(gca, 'XLim', [0.5 N_data_sets+0.5], 'XTick', 1:N_data_sets, 'XTickLabel', x_tags), end
if ~verLessThan('matlab','8.1')
    set(gca, 'TickLabelInterpreter', 'none')
end
% % % labels = cellfun(@(x) strrep(x,'um','\n'), x_tags,'UniformOutput',false);% % % a = gca;% % % a.XTickLabel = labels;% % % % % % if ~do_log, a = axis; if (a(4) > 0), axis([a(1) a(2) 0 a(4)]), end, end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% add stats (for 2 data sets)
result = stat_test(yy, do_log);
title(result)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if N_data_sets == 3
    add_stats({data(:,1) data(:,2) data(:,3)})
end
if N_data_sets == 4
    add_stats({data(:,1) data(:,2) data(:,3) data(:,4)})
end% legend(leg,legend_text)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [m, s, legend_key] = plot_mean_err(n, yy, c, scatter_width, do_log, use_sem, plot_bar)

linewidth_mean = 3; 
linewidth_sd = 1; 

yy(abs(yy)==inf) = [];

if ~isempty(yy)
    if do_log
        m = exp(nanmean(log(yy)));
        s = exp(nanstd(log(yy)));
        if use_sem, s = s^(1/sqrt(length(yy))); end
        
        legend_key = line(n+scatter_width*[-1 1]*2.5, m*[1 1], 'color', c, 'linewidth', linewidth_mean); % m
        if plot_bar
            line(n-scatter_width*[1 1]*2.5, m*[0 1], 'color', c, 'linewidth', linewidth_mean);
            line(n+scatter_width*[1 1]*2.5, m*[0 1], 'color', c, 'linewidth', linewidth_mean);
        end
        line(n+scatter_width*[-1 1], m*s*[1 1], 'color', c, 'linewidth', linewidth_sd) % sd, log
        line(n+scatter_width*[-1 1], m/s*[1 1], 'color', c, 'linewidth', linewidth_sd) % sd
        line(n*[1 1], [m/s m*s], 'color', c, 'linewidth', linewidth_sd)
    else
        m = nanmean(yy);
        s = nanstd(yy);
        if use_sem, s = s/sqrt(length(yy)); end
        
        legend_key = line(n+scatter_width*[-1 1]*2.5, m*[1 1], 'color', c, 'linewidth', linewidth_mean);
        if plot_bar
            line(n-scatter_width*[1 1]*2.5, m*[0 1], 'color', c, 'linewidth', linewidth_mean);
            line(n+scatter_width*[1 1]*2.5, m*[0 1], 'color', c, 'linewidth', linewidth_mean);
        end
        line(n+scatter_width*[-1 1], (m+s)*[1 1], 'color', c, 'linewidth', linewidth_sd)
        line(n+scatter_width*[-1 1], (m-s)*[1 1], 'color', c, 'linewidth', linewidth_sd)
        line(n*[1 1], m+s*[-1 1], 'color', c, 'linewidth', linewidth_sd)
    end
else
    m = NaN; s = NaN; legend_key = 0; % mean([])=NaN;
end
%%%