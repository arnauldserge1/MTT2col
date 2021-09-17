function [filename, refitdata, output_dir, seuil_premiere_detec, seuil_detec_1vue, wn, sig_free, T, Nb_STK, r0,...
    nb_defl, T_off, Boule_free, Nb_combi, seuil_alpha, Poids_melange_aplha, Poids_melange_diff, SHOW, conf_method,...
    OPTIM_R, Tmin, Dmin, pxl_size, time_lag,  use_gui, name_param1, name_param2, params] = MTT23_dialog_box(params)


%% -------------------- dialog box ---------------------
name_param1 = {...
    'filename';...
    'refit data? (0/1)';...
    'data directory';...
    'output folder';...
    'pre-threshold PFA_peak';...
    'final threshold PFA_orphan';...
    'spatial sliding search window wn (pxl)';...
    'diff max (pxl²/lag)';...
    'temporal sliding search window wt (frame)';...
    '# of frames/stack (empty = auto detect)';...
    'peak radius r0 (std, in pxl)';...
    };
name_param2 = {...
    '# deflation loops';...
    'disappearance prob of blinking T_off (frame)';...
    'ref diametre for research set part/traj (R = n.r_max)';...
    'max # combi part/traj';...
    'min intensity';...
    'weight of Gauss (on) vs uniform (blink) intensity law';...
    'weight of local vs max diff';...
    'show detect results (0, XX seconds per frame, or inf for manual)';... AS 10/10/13..
    'confinement detection method (coloc SCI int time rel abs speed)';...
    'optimize Gaussian radius r0 (0/1)';...
    'min. trace length (frame)';...
    'min. diffusion coef. (pxl2/frame)';...
    'pixel size (um)';... % not nm!!! AS 11/9/2014
    'time lag (s)';...
    %     'do reconnect';...
    };

%% load parameters / lecture params par défaut dans script MTT_param

params_default = MTTparams_def;
params_def1 = params_default(1:length(name_param1));
params_def2 = params_default((1:length(name_param2))+length(name_param1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(params) %%% && strcmp(VERSION,'MATLAB')
    params1 = inputdlg(name_param1,'Parameters - list 1/2',1, params_def1,'on');
    if isempty(params1)
        filename = ''; refitdata = ''; output_dir = ''; seuil_premiere_detec = ''; seuil_detec_1vue = ''; wn = ''; sig_free = ''; T = ''; Nb_STK = ''; r0 = '';
        nb_defl = ''; T_off = ''; Boule_free = ''; Nb_combi = ''; seuil_alpha = ''; Poids_melange_aplha = ''; Poids_melange_diff = ''; SHOW = ''; conf_method = '';
        OPTIM_R = ''; Tmin = ''; Dmin = ''; pxl_size = ''; time_lag = '';  use_gui = ''; name_param1 = ''; name_param2 = ''; params = '';
        return
    end
    
    params2 = inputdlg(name_param2,'Parameters - list 2/2',1, params_def2,'on');
    if isempty(params2)
        filename = ''; refitdata = ''; output_dir = ''; seuil_premiere_detec = ''; seuil_detec_1vue = ''; wn = ''; sig_free = ''; T = ''; Nb_STK = ''; r0 = '';
        nb_defl = ''; T_off = ''; Boule_free = ''; Nb_combi = ''; seuil_alpha = ''; Poids_melange_aplha = ''; Poids_melange_diff = ''; SHOW = ''; conf_method = '';
        OPTIM_R = ''; Tmin = ''; Dmin = ''; pxl_size = ''; time_lag = '';  use_gui = ''; name_param1 = ''; name_param2 = ''; params = '';
        return
    end
    params = [params1; params2];
    use_gui = 1;
elseif isnumeric(params)
    params = MTTparams_def(params); % assuming params === opt !!
    use_gui = 0;
else% complete, si seulement premiers champs fournis
    if length(params)<length(params_default)
        params = [params(:); params_default(length(params)+1:end)];
    end
    use_gui = 0;
end
if isempty(params),return, end

%% ----------------- reattribute parameters --------------------
filename = params{1};
refitdata = str2double(params{2});

% % % % cd(params{3})
output_dir = params{4};

seuil_premiere_detec = str2double(params{5});
seuil_detec_1vue = str2double(params{6});
wn = str2double(params{7});
Dmax = str2double(params{8}); % pxl/lag => µ2/s?????
sig_free = 2*sqrt(Dmax);

T = str2double(params{9});
Nb_STK = str2double(params{10});
if strcmp(params{10},''), Nb_STK = []; end

r0 = str2double(params{11});
nb_defl = str2double(params{12});
T_off = str2double(params{13});
Boule_free = str2double(params{14});
Nb_combi = str2double(params{15});
seuil_alpha = str2double(params{16});
Poids_melange_aplha = str2double(params{17});
Poids_melange_diff = str2double(params{18});
SHOW = str2double(params{19});
conf_method = params{20};
OPTIM_R = str2double(params{21});
Tmin = str2double(params{22});
Dmin = str2double(params{23});
pxl_size = str2double(params{24});
time_lag = str2double(params{25});
%%%