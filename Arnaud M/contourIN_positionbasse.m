function [masque_IN_positionbasse , masque_contourIN_positionbasse ] = ...
    contourIN_positionbasse(image_travail,Taille_matricegauss,Taille_median,bas,haut)

%% Explication
%
% ce programme fait ressortir les détails de l'image.
% On convolue chaque pixel aux pixels proche d'une zone
% (Taille_matricegauss,Taille_matricegauss) avec les coeff d'une gaussienne...
%
% 1) On crée la matrice 2D gaussienne, puis on la centre dans l'image.
%
% 2) Pour chaque pixel, on convolue avec cette matrice
%
% 3) On calcule une matrice qui pour chaque pixel prend la médiane de ses
% Taille_median^2 plus proche voisin, puis on divise, pour chaque pixel,
% cette matrice par la précedente.
%
% 4) histogramme : On supprime le pic central correspondant au bruit de
% fond
%
% 5) On passe l'image dans le programme traitement, qui ressort le masque.
%
% 6) On affiche l'image  avec le contour du masque...

%% Paramètres et préliminaires

if (nargin < 1)
    disp('mettre une image dic en paramètre')
end
Tr_image = double(image_travail);
masque_contourIN_positionbasse = imadjust(image_travail);
[Nbrligne_image Nbrcolonne_image] = size(Tr_image);

%% 1) Création de la matrice (Taille_matricegauss,Taille_matricegauss) centrée.

matricegauss = gausswin2(1,Taille_matricegauss,Taille_matricegauss);
matricegauss = expand_w(matricegauss,Nbrligne_image,Nbrcolonne_image)./sum(matricegauss(:));

%% 2) Convolution du filtre et de l'image

tf_image = fft2(Tr_image);
tf_matricegauss = fft2(matricegauss);

tf_final = real(fftshift(ifft2(tf_image.*tf_matricegauss)));

%% 3) division de la medianne par la gaussienne

median_image = medfilt2(Tr_image, [Taille_median Taille_median]) ;

div_image = median_image./tf_final;


%% 4) Histogramme

seuil_div_image = (div_image > bas) & (div_image < haut);

seuil_div_image = ones(Nbrligne_image,Nbrcolonne_image) - seuil_div_image;

if (bwarea(seuil_div_image) > 80000)
    haut = haut + 0.0005;
    seuil_div_image = (div_image > bas) & (div_image < haut);
    seuil_div_image = ones(Nbrligne_image,Nbrcolonne_image) - seuil_div_image;
end

%% 5) Traitement

% Appel du programme de traitement
masque_IN_positionbasse = traitementbas(seuil_div_image);

% Contour du masque
BWoutline = bwperim(masque_IN_positionbasse);
    % Matrice pour tracer le contour de manière visible
    se1 = zeros(3,3);
    se1(1,2) = 1; se1(2,1) = 1; se1(2,2) = 1;
% Superposition des deux images
masque_contourIN_positionbasse(imdilate(BWoutline, se1)) = max(masque_contourIN_positionbasse(:));
% Affichage de l'image
%figure, imshow(masque_contourIN_positionbasse), title(['contour du masque de position basse de l ''image '    ' avec bas = ' num2str(bas) ' avec haut = ' num2str(haut) ]); colormap(gray)
