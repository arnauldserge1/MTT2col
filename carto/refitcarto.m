function refitcarto(filename)

if nargin==0, filename = '*.stk'; end

cd ('M:\Marie-Claire Blache\camera128')%Arnauld Serge')

folders=dir;
folders(1:2)=[]; % . et ..

params_default = MTTparams_def;
refit = 1;
params_default{3} = filename;
params_default{end} = refit;

for i=10:14%1:length(folders) % [8 10 16 17 20]% 6 7 dossiers pertinents
    if isdir(folders(i).name) % && strfind(f{i},'chol')
        cd (folders(i).name)
        if ~isempty(dir(filename))
            disp([num2str(i) ' - ' folders(i).name])
            params_default{2} = cd;
            %ct2_by_file(filename)
            MTT22i(params_default)

        end
        cd ..
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% global USE_CONDITIONS, USE_CONDITIONS = 0;
% %global redo, redo = 1;
% %dispmod = [46.2 46.7 46.1 46.6 46.8 46.4]; [50.1 %[82 47.2 47.7 47.9 47.8 47.1 47.6 47.4]
%
% cd('\\Pcdm109506\20061218 Backup poste SPT\Arnauld\SPT\2005-06-03 LatB')
% %carto_movie('ctl1a-10nM.stk'),carto_movie('ctl2a-10nM.stk'),carto_movie('ctl4a-10nM.stk')%,36,1,dispmod
% ct2_by_file('ctl*.stk')
%
% %cd ..,cd('2005-06-02 LatA')
% %carto_movie('ctl6b.stk'), BOFFF... 1nM...
% %carto_movie('ctl7a.stk'),carto_movie('ctl7b');% 6c carto_movie('cell6b-10nM.stk')
% %ct2_by_file('ctl*.stk')
%
% cd ..,cd('2005-09-29 SPT specif')
% %carto_movie('cell1c.stk'),carto_movie('cell2c.stk'),carto_movie('cell3a.stk')
% %('cell1a.stk')
% ct2_by_file('cell*.stk')
%
% cd ..,cd('2005-10-03 SPT EGF-EGFR')
% %carto_movie('EGFR-5c.stk') %('EGFR-4a.stk')
% %%%carto3D('EGFR-1a.stk')
% %ct2_by_file('EGFR*.stk')
%
% cd ..,cd('2006-07-06 cholox+latB')
% %%%carto3D('ctl2a.stk');
% ct2_by_file('ctl*.stk')
%
% % FIT CRASHED!!!!!!!!!!!!!!!!!!!!!!!!!
%
% cd('\\Pcmicrovideo\datad\Arnauld\2007-02-27 gamme temps acquis'),
% %%%%carto_movie('cell3 36ms 5000img.stk'),carto_movie('cell5 20ms 6000img 20pc fluo.stk')
% %('cell4 36ms 300img.stk')
% %%%carto3D('cell2 36ms 300img.stk');carto3D('cell3 36ms 300img.stk');
% ct2_by_file('cell*.stk')
%
% return
