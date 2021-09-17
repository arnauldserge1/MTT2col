function miniatures
% plot la 1e image de chaque stk, en contraste inversé (=> pics noirs sur fond clair)
% AS 30/3/7

Nlignes = 4 ; Ncolonnes = 6 ;
files = dir('*.stk') ;
if isempty(files), files = dir('*.tif'); end
Nfiles = length(files) ;
Nfig = ceil(Nfiles/(Nlignes*Ncolonnes)) ;

for ifig=1:Nfig
    if Nfig>1, savename = ['miniatures_' num2str(ifig)];
    else savename = 'miniatures'; end
    
    if isempty(dir([savename '.png']))
        figure, colormap('gray')
        for i=1:Nlignes*Ncolonnes
            iplot = i + (ifig-1)*Nlignes*Ncolonnes ;
            if iplot<= Nfiles
                subplot(Nlignes,Ncolonnes,i)
                img = imread(files(iplot).name,1) ;
                imagesc(max(img(:))-img)
                title(files(iplot).name,'interpreter','none','FontSize',6)
                axis image off
            end
        end

        %set(gcf,'WindowStyle', 'normal','units','normalized','position',[0 0 1 1])
        saveas(gcf, savename, 'png')    %saveas(gcf, savename, 'emf')
        %saveas(gcf, savename) % .fig Matlab
        pause(.1), close
    else
        if ifig==Nfig, disp('mini c bon!'), end
    end
end % for ifig=1:Nfig
pause(.1)