function [mean_data, sd_data, value_at_max] = ct4plot_hist(data, varname, do_log, Nfrm, clr, norm, text_offset, Nbins, smooth_span, cum_prob)
% function [mean_data, sd_data, value_at_max] = ct4plot_hist(data, varname, do_log, Nfrm, clr, norm, text_offset, Nbins, smooth_span, cum_prob)
% default: ct4plot_hist(data, '', 1, -1, [0 0 1], 0, 0, 0, 0, 0)
% see also ct4

if nargin<10, cum_prob = 0; end
if nargin<9, smooth_span = 0; end
if nargin<8, Nbins = 0; end
if nargin<7, text_offset = 0; end
if nargin<6, norm = 0; end
if nargin<5, clr = [0 0 1]; end
if nargin<4, Nfrm = -1; end
if nargin<3, do_log = 1; end
if nargin<2, varname = ''; end

data = data(:); % 27/4/2017
data(abs(data)==inf) = [];
data(isnan(data)) = [];
data(data==0) = []; % if ~strcmp(varname, 'z (nm)') data = data(data>0);
if do_log, data(data<0) = []; end

gamma_min = -2;
gamma_max = 4;
if strcmp(varname, 'gamma') % AS 23/6/2021: keep data only around 1 (theoretical range: [0 2]). Ousiders are put at limit, not discarded, for graph histo and for mean/sd
    data(data < gamma_min) = gamma_min;
    data(data > gamma_max) = gamma_max;
end

xmin = min(data);
xmax = max(data);
if Nbins == 0, Nbins = ceil(2*sqrt(length(data))); end % nb 'harmonieux' de barres: 2sqrt(N) % if length(data) > 400, Nbins = 80, end

if cum_prob, Ndata = length(data); nn = (1:Ndata)/Ndata; end
    
if do_log
    step = (log(xmax)-log(xmin))/Nbins;
    x = exp(log(xmin):step:log(xmax));
else
    step = (xmax-xmin)/Nbins;
    x = (xmin:step:xmax);
end

if (length(x) > 1) && (xmax > xmin)
    p = hist(data, x); %     p = p/max(p); % if norm, end
    xxx = [x(1), x, 2*x(end)-x(end-1)]; % rajoute les bords de l'histo (manquent chez Matelabe, par defaut... cf doc stairs)
    [xx, pp] = stairs(xxx, [0 p 0]);
    value_at_max = mean(x((p==max(p)))); % (sommet de la distribution de proba => valeur satistiquement la + probable)
        
    if (smooth_span > 0), pp = smooth(pp, smooth_span); end
    if norm, pp = pp/max(pp); end
    
    if do_log
            if cum_prob, semilogx(sort(data), nn, 'color', clr)
            else, semilogx(xx, pp, 'color', clr)
            end
        mean_data = exp(mean(log(data)));% moyenne geometrique (log)
        sd_data = exp(std(log(data)));
        sep_string = '*/';
    else
        if cum_prob, plot(sort(data), nn, 'color', clr)
        else, plot(xx, pp, 'color', clr)
        end
        if strcmp(varname, 'z (nm)')
            a = axis;
            ind = [1; (xx~=value_at_max)];
            if (norm == 0), axis([a(1) a(2) 0 max(pp(ind(1:end-1) & ind(2:end)))]), end % ecarte val. par defaut, gamma...
        end
        mean_data = mean(data);
        sd_data = std(data);
        sep_string = '+/-';
    end
    
    s1 = sprintf('%s = %.3g%s%.3g', varname, mean_data, sep_string, sd_data);
    if any(strncmp(varname, {'TrcLen' 'Mean Speed' 'Time Coloc' 'Tc (img)' 'Tf (img)'}, 5))
        s2 = sprintf('N = %i', length(data));
    elseif strncmp(varname, 'alpha' , 5)
        s2 = sprintf('N = %.1f pk/frm', length(data)/Nfrm);
    else
        s2 = '';
    end
%     if (text_offset == 0), txtclr = 'k'; else, txtclr = clr; end
    if (text_offset > 0), text(0.05, 0.9-text_offset, {s1 s2}, 'Units', 'normalized')%, 'color', clr);
    elseif (text_offset == 0), text(0.05, 0.8, {s1 s2}, 'Units', 'normalized'); end
%     elseif (text_offset < 0), end

    if strcmp(varname, 'gamma')
        hold on
        a = axis;
        plot([1 1], [a(3) a(4)], ':')
    end
else
    mean_data = -1; sd_data = -1; value_at_max = -1;
end

xlabel(varname)