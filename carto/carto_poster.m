%%% carto_poster

cd ('C:\20061218 Backup poste SPT\Arnauld\SPT\2007-02-27 gamme temps acquis')
%\\Pcmicrovideo\data\Arnauld\2007-02-27 gamme temps acquis')
filename = 'cell2 36ms 300img.stk';

%trc = detect_reconnex_to_Leiden(filename);%#ok ,1,dirname
%L = probaconf(filename,0,'var',36);
full = [1 512 1 512];
ROI = [100 200 100 200];

carto3D(filename,2,full)%,40);
carto3D(filename,1,ROI)%,4);

%cartobyf(filename)
%axis(ROI)
% function carto3D(trc, L, s, ROI, do_ech, Sm, Nvoisins, dist_max)
% >>> function carto3D(filename, lissage, ROI, Nech, do_subplot, t_first, t_last, nf)