cd('M:\Tom Trombik\2009-05-13 LatrA Cascade128 Fab LotA')

S=128 ; % size(img)
N=1000 ; % Nimg
P=36/2 ; % time ratio
Q=(floor(N/P)) ;
N=P*Q ;

files = dir('*.stk');

for f =1:length(files)
    file = files(f).name;%'LatrA_2ms_a.stk';
    name_out = [file(1:end-4) '_2to36ms.tif'] ;%name_out = 'seq_2to36ms.tif' ;
    if ~isempty(dir(name_out)), delete(name_out), end
    
    for q = 1:Q
        im_p =  zeros(S,S) ;
        for p = 1:P
            im_p = im_p + double(tiffread(file,(q-1)*P+p));
        end%for
        im_p = im_p / P;
        if (0)
            imagesc(im_p) ;
            axis image
            colormap('gray')
            pause ;
        end%if
        imwrite(uint16(im_p), name_out, 'WriteMode', 'append', 'Compression', 'none') ;
        disp([name_out num2str(q)]) ;
    end%for
end%for