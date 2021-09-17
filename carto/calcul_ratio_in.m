function ratio = calcul_ratio_in
% calcule le nb (fraction) de pics ds le masque correspondant
% AS & MCB 10/12/7

files = dir('*.stk');
ratio = zeros(length(files),1);

for i=1:length(files)
    filename = files(i).name;
    trc = detect_reconnex_to_trc(filename);
    
    DIC_name = dicname(filename);
    maskdir = ['DIC\masques en position basse\Masque IN\MIB.' DIC_name(5:end)];
    maskIN = imread(maskdir);
    i_in = select_pk_in_mask(trc(:,4), trc(:,3), maskIN);

    density_in = length(i_in)/numel(maskIN(maskIN>0)); % nb de pics/pixel
    density_out = (size(trc,1)-length(i_in))/numel(maskIN(maskIN==0));
    ratio(i) = density_in/(density_in+density_out);

    disp(['ratio of peaks in mask = ' num2str(ratio(i))])
    pause(.1)
end

save('ct2_by_file\ratio_in.txt', 'ratio', '-ascii');