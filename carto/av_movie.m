%% av_movie

cd('\\Pcdm109506\20061218 Backup poste SPT\Arnauld\SPT\2005-09-29 SPT specif\carto\cell1c_3Dmovie_Gm')
im_name = 'cell1c_3D';
%mkdir('av')
n_av = 3; % ???
n_frames = 994;

im = cell(n_av);
disp(['averaging' im_name '..'])

for i=1:n_frames-n_av+1
    im_av = uint8(zeros(901,1201,3));
    for j=1:n_av
        im{j} = imread([im_name num2str(i+j-1) '.png']);
        im_av = im_av+im{j}/n_av;
    end
    imwrite(im_av, ['av' filesep im_name 'av_' num2str(i) '.png']);
    fprintf('\b\b\b%3.0i', i)
end

