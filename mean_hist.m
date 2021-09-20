function mean_hist(data, Np)
% function mean_hist(data, Np)
% build histogram with mean and sd

if nargin < 2, Np = 100; end

Nd = length(data);

min_x = inf;
max_x = -inf;
for nd = 1:Nd
    min_x = min(min_x, min(data{nd}(:)));
    max_x = max(max_x, max(data{nd}(:)));
end
xh = linspace(min_x, max_x, Np);
Nx = length(xh);

yh = zeros(Nd, Nx);
for nd = 1:Nd
    yh(nd, :) = hist(data{nd}(:), xh);
end

xh2 = [xh(1), xh, 2*xh(end)-xh(end-1)]; % rajoute les bords de l'histo (manquent chez Matelabe, par défaut... cf doc stairs)
[xx, yy] = stairs(xh2, [0 mean(yh) 0]);
plot(xx, yy, 'k')

hold on
bar_width = 0.5*(xh(2)-xh(1));

for nx = 1:Nx-1
    m = mean(yh(:, nx));
    s = std(yh(:, nx))/sqrt(Nd);
    xbar = mean(xh(nx: nx+1));
    plot(xbar + bar_width*[-1 1], (m+s)*[1 1], 'k')%, 'color', c, 'linewidth', 2)
    plot(xbar*[1 1], m + s*[0 1], 'k')
end