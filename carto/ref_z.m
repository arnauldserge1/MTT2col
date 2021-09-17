function ref_z(file)
% affiche image et pics, demande rectangle de référence (en dehors de la
% cellule, cf masque??) , calcule z moyen et sauve dans [fichier _z_ref.txt]
% => z_ref = load(['carto\' file '_z_ref.txt']);

global N_PARAM PARAM_I PARAM_J % PARAM_K PARAM_ALPHA gamma_z zr Boule_free
if isempty(N_PARAM), MTTparams_def; end

if nargin==0, [file file_dir] = uigetfile('*.stk'); cd(file_dir), end
% if nargin==0, file = ('*.stk'); end

% files = dir(file);
% for nf = 1:length(files)
%     file = files(nf).name;
im1 = imread(file,1);
figure('windowstyle','docked')
imagesc(im1)
axis image
colormap(1-gray)

tab_param = fread_all_params(file);

x_tab = tab_param(PARAM_J-1:N_PARAM:end,:); % x==j..
y_tab = tab_param(PARAM_I-1:N_PARAM:end,:);
% z_tab = tab_param(PARAM_K-1:N_PARAM:end,:);

hold on
plot(x_tab(:),y_tab(:),'.','markersize',2)

disp('select rectangle for ref, please (allez quoi!)')
rect = getrect;
fprintf('x=%g y=%g w=%g h=%g\r',rect)

if isempty(dir('carto')), mkdir('carto'), end
save(['carto\' file '_rect_ref.txt'], 'rect', '-ascii')

% xmin = rect(1);
% ymin = rect(2);
% xmax = rect(1) + rect(3);
% ymax = rect(2) + rect(4);
% 
% ref = (x_tab>xmin) & (x_tab<xmax) & (y_tab>ymin) & (y_tab<ymax);
% z_ref = z_tab(ref);
% z_ref = mean(z_ref(:));
% fprintf('z_ref = %g nm\r', z_ref)
% 
% save(['carto\' file '_z_ref.txt'], 'z_ref', '-ascii')
% end