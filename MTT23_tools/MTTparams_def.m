function params_default = MTTparams_def(opt)

global Nb_STK ;
global coloc_dist_max

% global USE_MAT

if nargin<1, opt = 0; end % opt = 1=> cell track

%% load parameters / lecture params par défaut dans script MTT_param
MTT_param
    
%% refit ou skip already fitted? 
refitdata = 0;
% do_reconnect = 1; %#ok % not for STORM???....

%% -------------------- dialog box fields ---------------------
params_default = {...
    filename;... % 1
    num2str(refitdata);... % 2
    cd;... % 3
    output_dir;... % 4
    num2str(seuil_premiere_detec);... % 5
    num2str(seuil_detec_1vue);... % 6
    num2str(wn);... % 7
    num2str(Din);... % 8
    num2str(T);... % 9
    num2str(Nb_STK);... % 10
    
    num2str(r0);... % 11
    num2str(nb_defl);... % 12
    num2str(T_off);... % 13
    num2str(Boule_free);... % 14
    num2str(Nb_combi);... % 15
    num2str(seuil_alpha);... % 16
    num2str(Poids_melange_aplha);... % 17
    num2str(Poids_melange_diff);... % 18
    num2str(SHOW);... % 19
    'coloc';... % 20 % coloc SCI int time rel abs speed Dv... %%% conf_method !!!
    num2str(OPTIM_R);... % 21
    num2str(Tmin);... % 22
    num2str(Dmin);... % 23
%     num2str(do_reconnect);... % 22
    '0.1375';% 95B=11um*2binning/160, optovar %% Cascade: '0.16';... pxl_size = str2double(params{24});  % 16 um at 100x magn
    '0.14';... 95B: 50msx2channels + delai  %% time_lag = str2double(params{25}); default 50 ms x 2 + 40 ms delay
    };

if isempty(coloc_dist_max), params_default{20} = 'speed'; end % int, time?? BUT NOT  coloc !!

%% Multi Cell Track
if opt==1
    %     p{1} = filtered_file; %%%% VOIR CELL_TRACK.M %%%%
    %         p{2} = '1'; % refit     %%%_Toff_150
    %         p{6} = '15'; % pfa%         p{7} = '20'; spfa
    %         p{8} = num2str(3*r0); % wn
    params_default{8} = '1.36'; % Din = (r/3)^2/4 (??) (r=7 taille cell detectees)
    %         p{12} = num2str(r0);
    params_default{12} = '0'; % #defl
    params_default{13} = '0'; % Toff
    params_default{20} = 'time'; % coloc SCI int time rel speed
    
    %% nanoscope - STORMTT
elseif opt==2
%     USE_MAT = 1;
    params_default{1} = '*.tif';
    params_default{8} = '0.01';  % Din = d_min^2/4 & d_min = 0.2 pxl @ 22 dB
    params_default{11} = '2';  % r0 = 650/2/1.49/(16000/150)
    params_default{20} = 'nano';
    params_default{21} = '0';
%     params_default{22} = '1'; % reconnect, par défaut 

%% Claire, PCM1
elseif opt==3  
    params_default{1} = '*.tif';
%     params_default{8} = '0.01';  % Din = d_min^2/4 & d_min = 0.2 pxl @ 22 dB
%     params_default{5} = '30'; % pfa
%     params_default{6} = '34';
    params_default{20} = 'SCI';
    params_default{24} = num2str(16/63); % pxl 16 um at 63x magn
    params_default{25} = '0.029'; % time_lag, s
    
%% Avais, Paxilin-Venus 100x!!!
elseif opt==4
%      params_default{8} = '0.01';  % Din = d_min^2/4 & d_min = 0.2 pxl @ 22 dB 100x
     params_default{8} = '0.004';  % Din = d_min^2/4 & d_min = 0.2*63/100=0.126 pxl @ 22 dB 63x
     params_default{23} = '0.001';  % D_min = d_min^2/4 & d_min = 0.2*63/100=0.126 pxl @ 22 dB 63x
%     params_default{5} = '16'; % pfa
%     params_default{6} = '20';
%     params_default{19} = '1'; % show results 1s for tests
    params_default{25} = '60'; % time_lag, s
%     params_default{20} = 'rainbow2'; % def in cell_track
    
%% Avais, 63x!!!
elseif opt==5
%      params_default{8} = '0.01';  % Din = d_min^2/4 & d_min = 0.2 pxl @ 22 dB 100x
     params_default{8} = num2str(0.04*(63/100)^2);  % Din = d_min^2/4 & d_min = 0.2*63/100=0.126 pxl @ 22 dB 63x
     params_default{23} = num2str(0.001*(63/100)^2);  % D_min = d_min^2/4 & d_min = 0.2*63/100=0.126 pxl @ 22 dB 63x
%     params_default{19} = '1'; % show results 1s for tests
    params_default{25} = '60'; % time_lag, s
end