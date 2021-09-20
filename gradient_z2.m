function [N1, N2] = gradient_z2(do_plot)

% function [N1, N2] = gradient_z2(do_plot)
%
% For a z-stack, compute axial gradient = number of detected molecules at
% each image, for each side (for 2 cameras, N1 at left, N2 at right)
% and mean gradient, for all files
% save graph N(z) for green & red, and values = matrix N(z,file), for green & red
%
% see also MTT_2_colors


global redo
if isempty(redo), redo = 0; end

if nargin < 1, do_plot = 1; end

outputdir = 'grad_z';
if ~isdir(outputdir), mkdir(outputdir), end
outputfile = [outputdir filesep 'gradient_z'];

if do_plot, figure('WindowStyle', 'docked'), end

files = dir2('*.tif');
Nf = length(files);
Nz_max = 0;
N1 = zeros(Nz_max, Nf);
N2 = zeros(Nz_max, Nf);
alreadydone = 0;

for nf = 1:Nf
    file = files(nf).name;
    if strcmp(file(1), '.'), continue, end
    
    %% check for multiple KG1 cells %%
    % %     if isdir('../DIC - Copy') % DIC images with black holes @ each KG1
    % %         dic_image_with_holes = imread(['../DIC - Copy' filesep file]);
    % %         BW = (dic_image_with_holes(:,1:512) > 0);
    % %         S = regionprops(~BW, 'Boundingbox');
    % %         N_cells = length(S);% count KG1 cells
    % %         %    imagesc(dic_image_with_holes(:,1:512)), title({file N_cells ' KG1 cells'})
    % %     end
    
    figname = [outputdir filesep file(1:end-4) '_grad.png'];
    if ~isempty(dir(figname)) && ~redo, disp([figname ' already done']), alreadydone = 1; continue, end
    
    %% load data, count SM, plot N(z) for left & right
    tab_param = importdata(['output23/' file '_tab_param.mat']);
    img1 = imread(file, 1); middle = size(img1, 2)/2;
    [tab_param_1, tab_param_2] = split_params_left_right(tab_param, middle, '', 0); % if strcmp(side, 'right'), tab_param_JAMB = tab_param_1; tab_param_JAMC = tab_param_2; else tab_param_JAMC = tab_param_1; tab_param_JAMB = tab_param_2; end
    tab_alpha_1 = tab_param_1(4:8:end, :);
    tab_alpha_2 = tab_param_2(4:8:end, :);
    Nz = size(tab_alpha_2, 1);
    Nz_max = max(Nz, Nz_max);
    N1(1:Nz, nf) = sum(tab_alpha_1 > 0, 2);
    N2(1:Nz, nf) = sum(tab_alpha_2 > 0, 2);
    
    %% graph
    if do_plot
        plot(N1(1:Nz, nf), 1:Nz, 'g', N2(1:Nz, nf), 1:Nz, 'r'), xlabel('number of molecules'), ylabel('z (frame)'), title(file), drawnow
        saveas(gcf, figname)
    end
end

%% graph of mean N & all N
if do_plot
    figname1 = [outputdir filesep 'mean_grad.png'];
    figname2 = [outputdir filesep 'all_grads.png'];
    
    if ~isempty(dir(figname2))
        disp([figname1 ' & ' figname2 ' already done'])
    else
        %% mean
        plot(mean(N1, 2), 1:Nz_max, 'g', mean(N2, 2), 1:Nz_max, 'r'), xlabel('mean number of molecules'), ylabel('z (frame)'), title('all files'), drawnow
        saveas(gcf, figname1)
        %% all
        plot(N1, 1:Nz_max, 'g', N2, 1:Nz_max, 'r'), xlabel('number of molecules'), ylabel('z (frame)'), title('all files'), drawnow
        saveas(gcf, figname2)
    end
    delete(gcf)
end

%% *** create txt file & print headers *** save number of green & red molecules at each z plane
if ~alreadydone
    save([outputfile '.mat'], 'N1', 'N2') % 28/3/2019
    
    if redo, wmode = 'w'; else, wmode = 'a'; end
    fid2 = fopen([outputfile '.txt'], wmode);
    if (fid2 == -1), warning('Couldn''t open file %s!!', [outputfile '.txt']), return, end 
    fprintf(fid2, '%s\t\r\n', date);
    
    %% save N left
    fprintf(fid2, '%s\t%s\t\r\n', 'file', 'N left (green channel, for ascending z)');
    for nf = 1:Nf
        fprintf(fid2, '%s\t', files(nf).name);
        fprintf(fid2, '%g\t', N1(1:Nz, nf));
        fprintf(fid2, '\r\n'); % retour ligne
    end
    %% save mean N
    fprintf(fid2, 'all_files_left\t');
    fprintf(fid2, '%g\t', mean(N1, 2));
    fprintf(fid2, '\r\n'); % retour ligne
    
    %% save N right
    fprintf(fid2, '%s\t%s\t\r\n', 'file', 'N right (red channel, for ascending z)');
    for nf = 1:Nf
        fprintf(fid2, '%s\t', files(nf).name);
        fprintf(fid2, '%g\t', N2(1:Nz, nf));
        fprintf(fid2, '\r\n'); % retour ligne
    end
    
    %% save mean N
    fprintf(fid2, 'all_files_right\t');
    fprintf(fid2, '%g\t', mean(N2, 2));
    fprintf(fid2, '\r\n'); % retour ligne
    
    fclose(fid2);
end
%%%