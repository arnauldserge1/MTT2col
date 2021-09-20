function [txt_result, p, stars] = stat_test(data, do_log, do_plot, selected_test)
% function [txt_result, p, stars] = stat_test(data, do_log, do_plot, selected_test)
% test normality (Single sample Kolmogorov-Smirnov goodness-of-fit hypothesis test)
% & check p value for d1 & d2 (t test or rank)
% assuming data = {d1 d2}
%
% Default:[txt_result, p, stars] = stat_test(data, 0, 0, 1)
% selected_test = 1 => ttest2, use Student by def, if both norm. distributions, else ranksum
% selected_test = 2 => kstest2, KS
% selected_test = 3 => vartest2, Var (needs norm. distributions in principle)
%
% see also scatter_hist

if (nargin < 4), selected_test = 1; end % => use Student by def
if (nargin < 3), do_plot = 0; end
if (nargin < 2), do_log = 0; end

if ~iscell(data)
    if size(data, 2) == 2, data = {data(:,1), data(:,2)}; % 2 columns !! caution if size(data) = (2,2)... !! (ok, then that's not a lot of data..)
    elseif size(data, 1) == 2, data = {data(1,:), data(2,:)}; % 2 rows
    end
end

if isempty(find(selected_test == 1:3, 1)), fprintf('Selected test #%i not available, soooory....\r', selected_test), end

%% check # datasets
N_data_sets = length(data);
% % % N_data_sets = size(data, 2);% % % yy = cell(N_data_sets,1);

%% check outliers
% % % ind = 1:size(data, 1);
% % % for n = 1:N_data_sets
% % %     if data_pos, ind = find(data(:, n) > 0); end
% % %     yy{n} = data(ind, n);
% % %     if remove_outliers% % %         TF = isoutlier(yy{n});%, outlier_method);%, 'gesd');%, 'grubbs');% 'mean');% % %         yy{n} = yy{n}(~TF);% % %     end
% % % end

if N_data_sets == 2
    data1 = data{1}(~isnan(data{1}) & ~isinf(data{1}));
    data2 = data{2}(~isnan(data{2}) & ~isinf(data{2}));
    
    %% normality test for each
    % data must be normalized, for test with G(0,1)!! 26/2/2019
    %     if isempty(data1) || isempty(data2)    %         h1 = -1; h2 = -1;    %     else
    if do_log
        data1 = data1(data1 > 0);
        data1l = log(data1);
        data1n = (data1l-mean(data1l))/std(data1l);
        if ~isnan(data1n), h1 = kstest(data1n); else, h1 = 0; end % 0 = norm % see adtest
        
        data2 = data2(data2 > 0);
        data2l = log(data2);
        data2n = (data2l-mean(data2l))/std(data2l);
        if ~isnan(data2n), h2 = kstest(data2n); else, h2 = 0; end
    else
        data1n = (data1-mean(data1))/std(data1);
        if ~isnan(data1n), h1 = kstest(data1n); else, h1 = 0; end % 0 = norm % see adtest
        data2n = (data2-mean(data2))/std(data2);
        if ~isnan(data2n), h2 = kstest(data2n); else, h2 = 0; end
    end
    %     end
    
    if (h1 == 0), txt_norm1 = 'data1 norm.'; else, txt_norm1 = 'data1 not norm.'; end
    if (h2 == 0), txt_norm2 = 'data2 norm.'; else, txt_norm2 = 'data2 not norm.'; end
    
    if do_plot
        figure
        [cdf1, xx1] = ecdf(data1n);
        subplot(221), plot(xx1, cdf1, 'b', xx1, normcdf(xx1, 0, 1), 'r:')
        legend('data1 normalized', 'Standard Normal CDF', 'Location', 'SE')
        title(txt_norm1)
        
        [cdf2, xx2] = ecdf(data2n);
        subplot(222), plot(xx2, cdf2, 'g', xx2, normcdf(xx2, 0, 1), 'r:')
        legend('data2 normalized', 'Standard Normal CDF', 'Location', 'SE')
        title(txt_norm2)
    end
    
    if (h1 == 0) && (h2 == 0) && length(data1) > 6 && length(data2) > 6
        %% param => T test
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if do_log, data1n = log(data1); data2n = log(data2); % normalized to log
        else, data1n = data1; data2n = data2;
        end
        
        switch selected_test
            case 1, [~, p] = ttest2(data1n, data2n);
            case 2, [~, p] = kstest2(data1n, data2n);% nota: yy is rid off outliers, eventually
            case 3, [~, p] = vartest2(data1n, data2n);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('Norm. data (n1=%i & n2=%i), using param. t-test: ', length(data1), length(data2))
        % %         test_used = 'S';%'t-test';
    else
        %% non param => Wilcoxon (1945) rank sum test = Mann-Whitney U test (1947) ref: https://fr.wikipedia.org/wiki/Test_de_Wilcoxon-Mann-Whitney
        %         if do_log, p = ranksum(log(data1), log(data2));else, ... end % ranktest not log sensitive
        %%%%%%%%%%%%%%%%%%%%%%%%%
        if isempty(data1) || isempty(data2), p = 0;
        else
            switch selected_test
                case 1, p = ranksum(data1, data2);
                case 2, [~, p] = kstest2(data1, data2); % nota: yy is rid off outliers, eventually % if do_log, [~, p] = kstest2(log(data1), log(data2));
                case 3, [~, p] = vartest2(data1n, data2n);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('%s (n1=%i), %s (n2=%i), using non-param. Wilcoxon test: ', txt_norm1, length(data1), txt_norm2, length(data2)) % % test_used = 'W';%'WMW-test';
    end
    
    stars = 'ns';
    if p < 0.05, stars = '*'; end
    if p < 0.01, stars = '**'; end
    if p < 0.001, stars = '***'; end    % %     if p < 0.0001, stars = '****'; end % ???
    txt_result = sprintf('%s p = %.2g', stars, p); % result = sprintf('%s %s p = %.2g', stars, test_used, p);
    fprintf('%s p = %.2g\r', stars, p)
    
    if do_plot
        subplot(212)
        if do_log
            [cdf1, xx1] = ecdf(data1);
            [cdf2, xx2] = ecdf(data2);
            semilogx(xx1, cdf1, 'b', xx2, cdf2, 'g')
        else
            plot(xx1, cdf1, 'b', xx2, cdf2, 'g')
        end
        legend('data1', 'data2', 'Location', 'SE')
        title(txt_result)
    end
else
    %% ANOVA....
    txt_result = '';
end
%%%

%% test on Gaussian values!!
% mu = 10; sig = 3; n = 100;
% normdist = normrnd(mu, sig, n, 1);
% subplot(121), histogram(normdist)
% kstest(normdist)
%
% normdist2 = (normdist-mu)/sig;
% subplot(122), histogram(normdist2)
% kstest(normdist2)
%%%