cd('P:\HETMLAB\DMLAB\Arnauld\Rémi - cellules T et B\rl100603bstiff')
mkdir('doubles_cartos')

for i=5:12
    figure('WindowStyle','docked')
    plot_all_traces(['rl100603bs_P' num2str(i) '_20x_GFPL_filt1.tif'],0,'g')
    plot_all_traces(['rl100603bs_P' num2str(i) '_20x_TRTC_filt1.tif'],0,'r')
    axis image
    figname = ['doubles_cartos\rl100603bs_P' num2str(i) '_GFPL_&_TRITC_carto.png'];
    saveas(gcf, figname, 'png'), close(gcf)
end
