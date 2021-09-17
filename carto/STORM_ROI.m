function ROI = STORM_ROI(name_stk, N_img_trial)


include_global 
global tab_num tab_param CROP IRANGE JRANGE

if nargin<1
    name_stk = uigetfile({'*.tif';'*.stk'});
%     user_input = 1;% else%     user_input = 0;
end

if nargin<2
    options.Interpreter = 'tex';
    str = {'Nombre d''images test? (0 => skip)'};
    answer = inputdlg(str, '', 1, {'10'}, options); 
    N_img_trial = str2double(answer{1});% if user_input, end    N_img_trial = 100;
end
if strcmp(name_stk(end-3:end),'.stk')
    [~, N_img_tot] = tiffread(name_stk, 1);
else
    N_img_tot = length(imfinfo(name_stk)); % cf Nb_STK
end
tab_num = round(N_img_tot/2) + (1:N_img_trial);

MTTparams_trial = MTTparams_def(2);
MTTparams_trial{1} = name_stk;
MTTparams_trial{4} = 'output_full_pict';
if N_img_trial>0, MTTparams_trial{2} = '1'; end

MTT23i(MTTparams_trial)

n_trc_tot = size(tab_param, 1);
x = tab_param(:,end-2*N_PARAM+PARAM_J-1) + 1; % x <=> j!
y = tab_param(:,end-2*N_PARAM+PARAM_I-1) + 1;

image_max(name_stk, tab_num); hold on
plot(x,y,'.','markersize',1)
axis equal

test = 0;
while ~test
    fprintf('\n  select ROI please... ')
    ROI = getrect(gca);
    ROI = round(ROI);
    fprintf('x0=%g, y0=%g, width=%g, height=%g\n',ROI)
    
    n_trc_ROI = sum(x>ROI(1) & x<ROI(1)+ROI(3) & y>ROI(2) & y<ROI(2)+ROI(4));
    prop = dir(['output_full_pict\' name_stk '_tab_param.dat']);
    trial_size = prop.bytes;
    
    expected_size = (N_img_tot/N_img_trial)*(n_trc_ROI/n_trc_tot)*trial_size;
    fprintf('test file (in output_full_pict) : %g Mo for %g particles in %g images\n', ...
        trial_size/2^20, n_trc_tot, N_img_trial)
    fprintf('expected size for output file tab_param.dat, using ROI: %g Mo for about %g particles in %g images\n', ...
        expected_size/2^20, n_trc_ROI, N_img_tot)
    answer = inputdlg('go??','',0.1,{'y'});
    if strcmp(answer,'y'), test = 1; end
end

%%%%
CROP = 1;
IRANGE = ROI(1):ROI(1)+ROI(3);
JRANGE = ROI(2):ROI(2)+ROI(4);
tab_num = 1:N_img_tot;

MTTparams = MTTparams_def(2);
MTTparams{1} = name_stk;

MTT23i(MTTparams)
nanopictbyf(name_stk)
%%%