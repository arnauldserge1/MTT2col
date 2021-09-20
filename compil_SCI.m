function compil_SCI(filename)

% cd('D:\Users\UTILISATEURS\Claire_Acquaviva\2014-09-18\Matlab')
if nargin<1, filename = '*.stk'; end

files = dir(filename);

if ~isdir('SCI_results'), mkdir('SCI_results'), end

for nf = 1: length(files)
    [v, f, d, t, N, corr_lin, csi_max] = compute_speed_from_SCI(files(nf).name);
    
%     outfile = ['SCI_results' filesep files(nf).name '_v.txt'];
%     save(outfile, 'v', '-ascii', '-tabs')
%     outfile = ['SCI_results' filesep files(nf).name '_f.txt'];
%     save(outfile, 'f', '-ascii', '-tabs')
%     outfile = ['SCI_results' filesep files(nf).name '_d.txt'];
%     save(outfile, 'd', '-ascii', '-tabs')
%     outfile = ['SCI_results' filesep files(nf).name '_t.txt'];
%     save(outfile, 't', '-ascii', '-tabs')
%     outfile = ['SCI_results' filesep files(nf).name '_N.txt'];
%     save(outfile, 'N', '-ascii', '-tabs')
    
    outfile = ['SCI_results' filesep files(nf).name '_SCI.mat'];
    save(outfile, 'v', 'f', 'd', 't', 'N', 'csi_max')
end