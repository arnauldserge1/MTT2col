function movie_traj(file, overlay, imagephase)
% cree une video des trajs 
% eventuellement superposees a l'image de fluo (overlay=1)
% ou de phase (overlay=2) 
% AS 23/6/4

if nargin==1
    overlay = 0;
end

trc = load(['trc\' file '.trc']);

if overlay==1
    debut = 1;
    [imagen, fin] = tiffread(file,1);
else
    debut = min(trc(:,2));
    fin = max(trc(:,2));
end

clf
hold on
axis ij, axis equal, axis tight, axis off

colormap('gray')

for n = debut:fin
%for n=fin
    if overlay==1
        clf
        if n>1, imagen = tiffread(file,n); end
        imagesc(imagen)
        hold on
        axis ij, axis equal, axis tight, axis off
    elseif overlay==2 && n==debut
        imagesc(imagephase)
    end
    
    title(n)
    
    indn = trc(trc(:,2)==n,1); % # des trajs concernant l'image n
    %for p=1:length(indn) % seulement trajs concernant image n
        %trcnp=trc(find(trc(:,1)==indn(p) & trc(:,2)<=n),:);
    for p = 1:max(indn) % toutes trajs de 1 a n
        trcnp = trc((trc(:,1)==p & trc(:,2)<=n),:);
        plot(trcnp(:,3),trcnp(:,4),'g') %,'linewidth',3)
        pause(.01)
    end
end
