cd \\amenophis\matlab$\Collaborateurs\Marek Cebecauer
filename = '4088 timelapse 0.5s 555 V1mChcl3 FlatBckg 071214.stk';
cartobyf('filename')
cd \\amenophis\matlab$\Collaborateurs\Marek Cebecauer\carto_v3_varM_output22\fig JCS review
title ''
axis off
axis([2 256 1 255])
saveas(gcf,'4088 v2','png')
axis([107 107+32 121 121+32])
saveas(gcf,'4088 v2 zoom','png')

% time rainbow??
cd \\amenophis\matlab$\Collaborateurs\Marek Cebecauer
traj_xyt(filename, '', 'time', 0, 0, 0)
cd \\amenophis\matlab$\Collaborateurs\Marek Cebecauer\carto_v3_varM_output22\fig JCS review
title ''
axis off
axis([2 256 1 255])
saveas(gcf,'4088 v2 time colors','png')
axis([107 107+32 121 121+32])
saveas(gcf,'4088 v2 time colors zoom','png')