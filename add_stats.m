function add_stats(data, do_log)
%% function add_stats(data, do_log)
% compute statistics by checking normality and choosing adequate test
% and plot a line with stars and p value
% see also scatter_hist, stat_test

if nargin < 2, do_log = 0; end

if size(data, 2) == 3
    [~,~,stars12] = stat_test({data{1} data{2}}, do_log);
    [~,~,stars23] = stat_test({data{2} data{3}}, do_log);
    [~,~,stars13] = stat_test({data{1} data{3}}, do_log);
    
    a = axis; ymax = a(4)*0.99;
    if do_log, yy = [1 1]*ymax; % ??
    else, yy = [1 1]*ymax; end
    
    plot([1.05 1.95], yy, 'k', 'linewidth', 1)
    plot([2.05 2.95], yy, 'k', 'linewidth', 1)
    plot([1.05 2.95], yy*1.1, 'k', 'linewidth', 1)
    
    text(1.5, ymax*0.95, stars12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')%, 'units', 'normalized')
    text(2.5, ymax*0.95, stars23, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
    text(2, ymax*1.05, stars13, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
    
elseif  size(data, 2) == 4 % assuming paired: 1 & 2 / 3 & 4
    [~,~,stars12] = stat_test({data{1} data{2}}, do_log);
    [~,~,stars34] = stat_test({data{3} data{4}}, do_log);
    
    a = axis; ymax = a(4)*0.93;
    if do_log, yy = [1 1]*ymax;
    else, yy = [1 1]*ymax; end
    
    plot([1.05 1.95], yy, 'k', 'linewidth', 1)
    plot([3.05 3.95], yy, 'k', 'linewidth', 1)
    
    if strcmp(stars12, 'ns'), text_pos = 1; else, text_pos = 0.95; end
    text(1.5, ymax*text_pos, stars12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
    if strcmp(stars34, 'ns'), text_pos = 1; else, text_pos = 0.95; end
    text(3.5, ymax*text_pos, stars34, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
end