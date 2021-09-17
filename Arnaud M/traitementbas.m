function masque_traite_enbas = traitementbas(image_obtenue)

%% Explication

% Cette fonction permet de recupérer un masque de la cellule en position
% basse, c'est àdire avec les filopodes...

% Elle fonctionne sur cette base : 

% On crée un masque grossier que l'on applique ensuite à l'image initiale
% afin d'enlever les bruits résultant du filtrage par fft.

% Puis par une suite de dilatation erosion, on recupère un filtre qui
% semble plus précis puisqu'il arrive à différencier certains filopodes.

%% Masque grossier

[Nbrligne_image Nbrcolonne_image] = size(image_obtenue) ;

% On supprime les points trop petits
bw1 = bwareaopen(image_obtenue,45,8); %figure, imagesc(bw1), title('bw1') ; colormap(gray)
% On rejoints
bw2 = imfill(bw1,8,'holes'); %figure, imagesc(bw2), title('bw2') ; colormap(gray)
% On dilate grossièrement
bw3 = imdilate(bw2,strel('disk',8)); %figure, imagesc(bw3), title('bw3') ; colormap(gray)
% On remplis
bw4 = imfill(bw3,4,'holes'); %figure, imagesc(bw4), title('bw4') ; colormap(gray)
% On erode pour arrondir
bw5 = imerode(bw4,strel('diamond',4)); %figure, imagesc(bw5), title('bw5') ; colormap(gray)
% On supprime les petites zones
masque_grossier = bwareaopen(bw5,2000,8); %figure, imagesc(masque_grossier), title('masque_grossier') ; colormap(gray)



%% Traitement fin

% On supprime les zones non-intéressantes
bw7 = masque_grossier.*image_obtenue; %figure, imagesc(bw7), title('bw7') ; colormap(gray)
% Création d'un masque de bord
bwprim = ones(Nbrligne_image,Nbrcolonne_image); bwprim(Nbrligne_image,1:2)=0;
bwprim(1:Nbrligne_image,Nbrcolonne_image-1:Nbrcolonne_image)=0; bwprim(1:2,1:Nbrcolonne_image)=0;
bwprim(Nbrligne_image-1:Nbrligne_image,1:Nbrcolonne_image)=0; 
% On supprime les bords
bw8 = bwprim.*bw7; %figure, imagesc(bw8),title('bw8') ; colormap(gray)
% On dilate suivant l'horizontale
bw9 = imdilate(bw8,strel('line',3,0)); %figure, imagesc(bw9), title('bw9') ; colormap(gray)
% On erode
bw10 = imerode(bw9,strel('line',3,0)); %figure, imagesc(bw10), title('bw10') ; colormap(gray)
% On remplis
bw11 = imfill(bw10,4,'holes'); %figure, imagesc(bw11), title('bw11') ; colormap(gray)
% On dilate suivant la verticale
bw12 = imdilate(bw11,strel('line',3,90)); %figure, imagesc(bw12), title('bw12') ; colormap(gray)
% On erode
bw13 = imerode(bw12,strel('line',3,90)); %figure, imagesc(bw13), title('bw13') ; colormap(gray)
% On remplis
bw14 = imfill(bw13,4,'holes'); %figure, imagesc(bw14), title('bw14') ; colormap(gray)
% On dilate suivant l'horizontale
bw9 = imdilate(bw14,strel('line',3,0)); %figure, imagesc(bw9), title('bw9') ; colormap(gray)
% On erode
bw10 = imerode(bw9,strel('line',3,0)); %figure, imagesc(bw10), title('bw10') ; colormap(gray)
% On remplis
bw11 = imfill(bw10,4,'holes'); %figure, imagesc(bw11), title('bw11') ; colormap(gray)
% On dilate suivant la verticale
bw12 = imdilate(bw11,strel('line',3,90)); %figure, imagesc(bw12), title('bw12') ; colormap(gray)
% On erode
bw13 = imerode(bw12,strel('line',3,90)); %figure, imagesc(bw13), title('bw13') ; colormap(gray)
% On remplis
bw14 = imfill(bw13,4,'holes'); %figure, imagesc(bw14), title('bw14') ; colormap(gray)
% On dilate suivant l'horizontale
bw9 = imdilate(bw14,strel('line',3,0)); %figure, imagesc(bw9), title('bw9') ; colormap(gray)
% On erode
bw10 = imerode(bw9,strel('line',3,0)); %figure, imagesc(bw10), title('bw10') ; colormap(gray)
% On remplis
bw11 = imfill(bw10,4,'holes'); %figure, imagesc(bw11), title('bw11') ; colormap(gray)
% On dilate suivant la verticale
bw12 = imdilate(bw11,strel('line',3,90)); %figure, imagesc(bw12), title('bw12') ; colormap(gray)
% On erode
bw13 = imerode(bw12,strel('line',3,90)); %figure, imagesc(bw13), title('bw13') ; colormap(gray)
% On remplis
bw14 = imfill(bw13,4,'holes'); %figure, imagesc(bw14), title('bw14') ; colormap(gray)
% On dilate
bw15 = imdilate(bw14,strel('diamond',1)); %figure, imagesc(bw15), title('bw15') ; colormap(gray)
% On remplis
bw16 = imfill(bw15,4,'holes'); %figure, imagesc(bw16), title('bw16') ; colormap(gray)
% On erode
bw17 = imerode(bw16,strel('diamond',1)); %figure, imagesc(bw17), title('bw17') ; colormap(gray)
% On supprime les points
bw18 = bwareaopen(bw17,5,4); %figure, imagesc(bw18), title('bw18') ; colormap(gray)
% On dilate
bw19 = imdilate(bw18,strel('disk',1)); %figure, imagesc(bw19),title('bw19') ; colormap(gray)
% On remplis
bw20 = imfill(bw19,4,'holes'); %figure, imagesc(bw20),title('bw20') ; colormap(gray)
% On erode
bw21 = imerode(bw20,strel('diamond',1)); %figure, imagesc(bw21), title('bw21') ; colormap(gray)
% On supprime les zones éloignées
bw21 = bwareaopen(bw21,30,4); % figure, imagesc(bw21), title('bw210') ; colormap(gray)
% On dilate
bw22 = imdilate(bw21,strel('disk',2)); %figure, imagesc(bw22),title('bw22') ; colormap(gray)
% On remplis
bw23 = imfill(bw22,4,'holes'); %figure, imagesc(bw23),title('bw23') ; colormap(gray)
% On erode
bw24 = imerode(bw23,strel('diamond',1)); %figure, imagesc(bw24),title('bw24') ; colormap(gray)
% On dilate
bw25 = imdilate(bw24,strel('disk',2)); %figure, imagesc(bw25),title('bw25') ; colormap(gray)
% On remplis
bw26 = imfill(bw25,4,'holes'); %figure, imagesc(bw26),title('bw26') ; colormap(gray)
% On erode
bw27 = imerode(bw26,strel('disk',3)); %figure, imagesc(bw27),title('bw27') ; colormap(gray)
% On supprime les grosses zones
% Resultat :
masque_traite_enbas = bwareaopen(bw27,10000,8); %figure, imagesc(masque_traite_enbas),title('masque_traite_enbas') ; colormap(gray)