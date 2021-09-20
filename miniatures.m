function miniatures
% plot the fist image of each  stk, in inverted contrast  (=> black peaks on light background)
% AS 30/3/7

Nlignes = 4 ; Ncolonnes = 6 ;
files = dir('*.stk');
if isempty(files), files = dir('*.tif'); end
files = sort_nat({files.name});
Nfiles = length(files);
Nfig = ceil(Nfiles/(Nlignes*Ncolonnes)) ;

for ifig=1:Nfig
    if Nfig>1, savename = ['miniatures_' num2str(ifig)];
    else, savename = 'miniatures'; end
    
    if isempty(dir([savename '.png']))
        figure, colormap('gray')
        for i=1:Nlignes*Ncolonnes
            iplot = i + (ifig-1)*Nlignes*Ncolonnes ;
            if iplot <= Nfiles
%                 subplot(Nlignes,Ncolonnes,i) 
                left = mod(i-1, Ncolonnes)/Ncolonnes; bottom = (i-1-mod(i-1,Ncolonnes))/Ncolonnes/Nlignes; width = 1/Ncolonnes; height = 1/Nlignes;
                subplot('Position',[left bottom width height])
                file = files{iplot};
                if strcmp(file(1), '.'), continue, end % for hidden files such as ._toto.tif...
                img = imread(file,1);
                img(:, 1) = max(img(:));
                img(1, :) = max(img(:));
                imagesc(max(img(:))-img)
%                 title(files(iplot).name,'interpreter','none','FontSize',6)
                text(0,0,file,'interpreter','none','FontSize',6,'units','normalized')
                axis image off
            end
        end

        %set(gcf,'WindowStyle', 'normal','units','normalized','position',[0 0 1 1])
        saveas(gcf, savename, 'png')    %saveas(gcf, savename, 'emf')        %saveas(gcf, savename) % .fig Matlab
        pause(.1), close(gcf)
    else
        if ifig==Nfig, disp('mini c bon!'), end
    end
end % for ifig=1:Nfig
pause(.1)