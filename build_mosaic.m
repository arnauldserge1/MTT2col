function build_mosaic(image_type)
% function build_mosaic(image_type)
% build and save one image from several ones

if (nargin == 0), image_type = 'carto'; end

warning off images:imshow:magnificationMustBeFitForDockedFigure % warning off images:initSize:adjustingMag

switch image_type
%     case 'raw'
        % % dir_img = '.'; img_type = 'tif'; b = 16; crop_val = []; % 16 bits...sat = 0.002 ;     imsat = imadjust(im, stretchlim(im, [sat 1-sat]), [0 1]);
    case 'max'
        dir_img = 'max'; img_type = 'tif'; nbits = 16;
    case 'carto'
        dir_img = 'carto_coloc_magenta_2reg'; % dir_img = 'carto_coloc_2reg';
        if ~isfolder(dir_img), dir_img = 'carto_coloc_2'; end
        if ~isfolder(dir_img)
            dir_img = dir('carto_*'); 
            if isempty(dir_img), disp('carto_* not found'), return, end
            dir_img = dir_img(1).name; disp(['using ' dir_img])
        end
        img_type = 'png'; nbits = 8;
    case 'carto_Dv_2colors'
        dir_img = 'carto_Dv_2colors'; img_type = 'png'; nbits = 8;
    otherwise
        fprintf('%s unknown for mosaic...\r', image_type)
end

dir_ini = cd;
if ~isfolder(dir_img), disp([dir_img ' not found']), return, else, cd(dir_img), end
if ~isempty(dir('mosaic.png')), disp('mosaic already done'), cd(dir_ini), return, end % 23/5/2017

files = dir2(['*.' img_type]); % remove possible previous mosaic image? or keep generating 'fractal' result....
files = sort_nat({files.name});
Nf = length(files);

Nmax = 4*3;
N_mosaic = ceil(Nf/Nmax);
im_tot = cell(4, 3, N_mosaic); % & transpose!

if ~isempty(dir(sprintf('mosaic_%i.png', N_mosaic))), disp('mosaic already done'), cd(dir_ini), return, end % 2018!

a = zeros(Nmax,1); b = a;
a2 = zeros(Nmax,1); b2 = a2;
n_ok = 0;

for nf = 1:Nf
    file = files{nf};
    if contains(file, 'mosaic.png') || contains(file, '._'), continue, end
    disp(['building mosaic: ' file]) 
    
    im = imread(file);
    [a(nf), b(nf), ~] = size(im);
    n_ok = n_ok+1;
    [i,j,k] = ind2sub([4,3,N_mosaic], n_ok); % k = ceil(n_ok/Nmax);    j = ceil(n_ok-(k-1)*Nmax/4);

    if contains(image_type, 'carto') %~isempty(crop_val)
        bw = (mean(im, 3) < 127);
        max_line = max(bw);
        max_col = max(bw, [], 2);
        first_line = find(max_col, 1);
        last_line = find(max_col, 1, 'last');
        first_col = find(max_line, 1);
        last_col = find(max_line, 1, 'last');
        im_crop = im(first_line:last_line, first_col:last_col, :);
        im_size = last_col-first_col+1;
        im_tot{i,j,k} = imresize(im_crop, 512/im_size);
    else
        size2 = round(size(im)*1.02); % add marg
        im2 = uint16(ones(size2)*2^nbits);
        [A, B, ~] = size(im2);
        a0 = round(A/2-a(nf)/2);
        b0 = round(B/2-b(nf)/2);
        im2(a0+1:a0+a(nf), b0+1:b0+b(nf), :) = im;
        im_tot{i,j,k} = im2;
    end
    [a2(n_ok), b2(n_ok), ~] = size(im_tot{n_ok});
    % if n_ok >= Nmax, fprintf('caution, %s mosaic stopped at %i\r', dir_img, Nmax), break; end % 2nd mosaic??
end

for nf = 1:n_ok
    im = ones(max(a2)+10, max(b2)+10, 3)*2^nbits;
    if contains(image_type, 'carto'), im = uint8(im); else, im = uint16(im); end 
    im(round((max(a2)+10-a2(nf))/2+1:round(max(a2)+10+a2(nf))/2), round((max(b2)+10-b2(nf))/2+1:round(max(b2)+10+b2(nf))/2), :) = im_tot{nf};
    [i,j,k] = ind2sub([4,3,N_mosaic], nf);
    im_tot{i,j,k} = im;
end
% if (std(a(a>0)) > 0) && (std(b(b>0)) > 0), cd(dir_ini), disp('not all at same size, mosaic cancelled'), return, end %  => put at max size??
for nf = (n_ok+1):(N_mosaic*Nmax)
    [i,j,k] = ind2sub([4,3,N_mosaic], nf);
    im_tot{i,j,k} = max(im_tot{1}, 2^nbits); % completed with white
end

for nm = 1:N_mosaic
    im_tot2 = cell2mat(im_tot(:,:,nm)'); % imshow(im_tot2)
    if ~isempty(im_tot2), fprintf('Saving mosaic_%i.png\r', nm), imwrite(im_tot2, sprintf('mosaic_%i.png', nm)), end
end
cd(dir_ini)
%%%