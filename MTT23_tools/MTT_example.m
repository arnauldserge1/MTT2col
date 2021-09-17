function MTT_example(file_name)
%%% Basic examples showing how to recover MTT output results
%%% to plot each trace and to build the histogram 
%%% of fluorescence intensities
%%% AS 2010

if nargin<1
    files = dir('*.tif');
    if isempty(files), files = dir('*.stk'); end
    if isempty(files), disp('no data in current dir'), return, end
    file_name = files(1).name;
end
file_param = [file_name '_tab_param.dat'];

%% load data
cd('output23')
% [tab_i,tab_j,tab_alpha,tab_ray,tab_7,tab_blk] = fread_all_data_spt(file_param);
tab_i = fread_data_spt(file_param,3);
tab_j = fread_data_spt(file_param,4);
tab_alpha = fread_data_spt(file_param,5);
tab_blk = fread_data_spt(file_param,8);
cd ..

%% loop over traces
N_traces = size(tab_i,1);

for itrc = 1:N_traces 
    No_blk_index = tab_blk(itrc,:)>0; % non blinking steps only
    plot(tab_i(itrc, No_blk_index),tab_j(itrc, No_blk_index))
    xlabel('i (pixel)'), ylabel('j (pixel)')
    title(['trace # ' num2str(itrc)])
    pause
end

%% fluo hist
N_datapoints = sum(tab_blk(:)>0); % non blinking steps only
hist(tab_alpha(tab_blk>0),2*sqrt(N_datapoints))
xlabel('intensity (a.u.)'), ylabel('occurence')
title('histogram of particles fluorescence intensity')