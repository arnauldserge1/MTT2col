function [mean_conf, var_conf, n_pos, visites_map] = mean_var_conf(filename,Nframes)
% compute mean & var of Lconf & #revisit at each pixel, using carto3Dv3, as in carto_moviev3

%cd('M:\Arnauld Serge\2005-09-29 SPT specif'),[m v n vm] = mean_var_conf('cell1c.stk');
% x = 70; y = 140; w = 400; h = 340; crp = [x+1 x+w y+1 y+h]; axis(crp)
% x = 100; y = 190; w = 170; ROI = [x+1 x+w y+1 y+w]; axis(ROI) %(zoom x3)
% 
% cd('M:\Arnauld Serge\2005-06-03 LatB'), [m v n vm] = mean_var_conf('ctl4a-10nM.stk');
% x = 30; y = 190; w = 128; ROI = [x+1 x+w y+1 y+w]; axis(ROI) %(zoom x4)
% x = 250; y = 350; w = 128; ROI2 = [x+1 x+w y+1 y+w]; axis(ROI2) %(zoom x4)%%%

%%  *** colormap ***
n = 256; sat = 4/3;
r = [(0:n*(1-1/2/sat)-1)'/(n*(1-1/2/sat)); ones(n/2/sat,1)];
g = [zeros(n*(1-3/4/sat),1); (0:n/2/sat-1)'/(n/2/sat); ones(n/4/sat,1)];
b = [zeros(n*(1-1/2/sat),1); (0:n/2/sat-1)'/(n/2/sat)];
AFMmap = [r g b]; AFMmap(1,:) = [0 0 .3]; % plancher bleu marine

[lg, ht, Ntiff] = stk_size(filename);
params = load('\\Amon\Picsl\HETMLAB\DMLAB\Matlab-scripts\Arnauld\carto\proba_conf_params_varM.dat')
wT = params(3)+3; % (wconf+3)
if nargin==1, Nframes = Ntiff-wT+1; Ntiff = Nframes; end %...

%% moyenne & std
mean_conf = zeros(lg,ht); sum_square = mean_conf; n_pos = mean_conf;

for n_image=1:Nframes
    fprintf('frame %g / %g ',n_image,Nframes)
%     Surf3D = carto3Dv3(filename, 0, n_image, n_image+wT-1); % Lc=Log(DT/var(r2)) (:,:,n_image)
    trcn = detect_reconnex_to_trc(filef, 1, dirname, n_image, n_image+wT-1);
    Surf3D = carto3Dv3(trcn, filef, 0); % Lc=Log(DT/var(r2)) (:,:,n_image)
    mean_conf = mean_conf + Surf3D;
    sum_square = sum_square + Surf3D.^2;
    n_pos = n_pos + (Surf3D>0);
end 

mean_conf = mean_conf./n_pos;
var_conf = (sum_square./n_pos) - (mean_conf.^2);
mean_conf(n_pos<Nframes-10) = 0; % supprime bords
var_conf(n_pos<Nframes-10) = 0;

%% DIC
% DIC_name = dicname(filename); DIC = imread(DIC_name);
% sat = .002 ; DIC_sat = imadjust(DIC,stretchlim(DIC, [sat 1-sat]),[0 1]);
% H = fspecial('average'); DIC_sat_f = imfilter(DIC_sat,H,'replicate');
%% revisites
if Ntiff<=300
    trc = detect_reconnex_to_trc(filename);
%     trc = discard_fix(filename);
    visites_map = revisites(trc,lg,ht);
    visites_map = visites_map(:,:,1);
else
    visites_map = zeros(lg,ht);
    for i=1:350:Ntiff
        trc = detect_reconnex_to_trc(filename,1,'',i,i+349);
        visites_map_i = revisites(trc,lg,ht);
        visites_map = visites_map+visites_map_i(:,:,1);
    end
end

%% graphs 
%%%AFMmap(1,:) = [0 0 .8]; % plancher bleu vif, pour revisit!!!
figure('WindowStyle','docked'), colormap(AFMmap)

subplot(2,3,1), imagesc(mean_conf), axis image, title('mean(L{conf})'), colorbar
subplot(2,3,2), imagesc(var_conf./mean_conf), axis image, title('VMR'), colorbar
subplot(2,3,3), imagesc(visites_map), axis image, title('# revisites'), caxis([5 60]), colorbar
%subplot(2,3,4), imagesc(DIC_sat_f), axis image, title('DIC'), colorbar % colormap('gray')...

subplot(2,3,4), plot(mean_conf(1:30:end), var_conf(1:30:end)./mean_conf(1:30:end),'.')
xlabel('mean(L{conf})'), ylabel('VMR'), title(cd), axis([0 2.1 0 0.3])

mp = mean_conf(visites_map>0);
vmp = visites_map(visites_map>0);
r = 0.2*rand(numel(vmp),1)-.1; % rand pour disperser en y..
subplot(2,3,5), plot(mp(:), vmp+r,'.') 
xlabel('mean(L{conf})'), ylabel('# revisites')
if ischar(filename), title(filename), else title('raw data (simul.mat?)'), end

vmr = var_conf(visites_map>0)./mp;
subplot(2,3,6), plot(vmr(:), vmp+r,'.') 
xlabel('VMR'), ylabel('# revisites')
%%%