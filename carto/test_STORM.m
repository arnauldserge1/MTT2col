function test_STORM

path2version('MTT23')

%     file = uigetfile({'*.tif';'*.stk'});
file = '*.tif';

files = dir(file);

opt = MTTparams_def(2);
wn = eval(opt{7});
r0 = eval(opt{11});
pfa = eval(opt{6});
n_deflt = 1;

n_part_defl = zeros(length(files),1);
for i_file = 1:length(files)
    
    pict = double(imread(files(i_file).name, 1));
    
    [~, ~, ~, n_part_defl(i_file)] = detect_et_estime_part_1vue_deflt(pict, wn, r0, pfa, n_deflt);

    fprintf('  %i peaks found after deflation in %s\n\n', n_part_defl(i_file), files(i_file).name)
end

disp('dense files:')
for i_file = 1:length(files)
    if n_part_defl(i_file)>10
        fprintf('  %i peaks found after deflation in %s\n', n_part_defl(i_file), files(i_file).name)
    end
end