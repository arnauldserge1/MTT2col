function AFMmap = load_AFMmap

%%  *** colormap ***
n = 256; sat = 4/3;
r = [(0:n*(1-1/2/sat)-1)'/(n*(1-1/2/sat)); ones(n/2/sat,1)];
g = [zeros(n*(1-3/4/sat),1); (0:n/2/sat-1)'/(n/2/sat); ones(n/4/sat,1)];
b = [zeros(n*(1-1/2/sat),1); (0:n/2/sat-1)'/(n/2/sat)];
AFMmap = [r g b];

AFMmap(1,:) = [0 0 .3]; % plancher bleu marine
