function masque(filename)

%% Explication
%
% Ce programme, à partir d'images DIC (filename = '*.tif', pr défaut), 
% crée les contours IN et CONT des cellules.
%
% Il fait appel au programme contourIN_positionbasse qui trouve son
% [masque,contour]. On affiche les contours pour pouvoir valider ou non les
% maques.

% contour = file de fer du masque dilaté, et non fil de fer dilaté du masque
% non valide non sauvé AS 24/8/7

%% Préliminaires

if nargin==0, filename = '*.tif'; end

dossier = cd;
dirchanged = 0;
if ~strcmp(dossier(end-2:end),'dic'), if isdir('dic'), cd dic, dirchanged = 1; end, end
%if dossier==0, return, end% Correspond à l'abandon

% Chargement des images dans le workspace
files = dir(filename); % [dossier '\*.tif']);
Nb_image = length(files);

% Si mauvaise selection
if (Nb_image==0)
    disp('aucune image DIC trouvée, le dossier n''est pas valide'), return
end

% Taille des images :
premiere_image = imread(files(1).name);%premiere_image = imread([files(1).directory, files(1).name]);
[Nbrligne_image Nbrcolonne_image] = size(premiere_image) ;

% Paramètres pour contourIN_positionbasse :
Taille_matricegauss = 3;
Taille_median = 11;
Masque_IN_enbas = ones(Nbrligne_image,Nbrcolonne_image*Nb_image);
Contours_IN_positionbasse = ones(Nbrligne_image,Nbrcolonne_image*Nb_image);
Masque_CONT_positionbasse = ones(Nbrligne_image,Nbrcolonne_image*Nb_image);
bas = 0.985;
if ~isempty(strfind(cd,'HIV')) % HeLa, en nunc, a priori
    haut = 1.0035 -.0015;
else % COS, a priori
    haut = 1.0035-.0022; % as 9/3/10
end
% im_nonvalide = zeros(Nb_image,1);
compteur_6=1; nfig = 1;
f(nfig) = figure('Position',[2 40 1275 910],'Name','images DIC 1 à 6');

%% Création des masques IN

disp('Contour des cellules en cours..   ')
% Pour 1 jusqu'au nombre d'image :
for l=1:Nb_image
    file_l = files(l).name;
    II = 1:Nbrligne_image;
    JJ = 1+Nbrcolonne_image*(l-1):Nbrcolonne_image*l;
    % On recupère le masque et le contour dessiné sur l'image à partir du
    % programme contourIN_positionbasse.
    [Masque_IN_enbas(II,JJ), Contours_IN_positionbasse(II,JJ)] = ...
        contourIN_positionbasse(imread(file_l),Taille_matricegauss,Taille_median,bas,haut);

    % Création de la figure

    if compteur_6>6
        nfig = nfig+1; % =floor(l/6)...
        f(nfig) = figure('Position',[2 40 1275 910],'Name',['images DIC ' num2str(l) ' à ' num2str(l+5)]);
        compteur_6=1;
    end

    % Affichage des images sur la figure
    subplot(2,3,compteur_6), imagesc(Contours_IN_positionbasse(II,JJ))
    title({file_l, ['image ' num2str(l) ' sur ' num2str(Nb_image)]},'interpreter','none')
    colormap(gray)
    compteur_6 = compteur_6+1;
    fprintf('\b\b\b%3i', l), drawnow
end

%% Sélection des images à retravailler

% Ouverture d'une fenêtre :
prompt=['Selectionner les images à retravailler : Rentrer dans la zone de texte le numéro des images'...
    ' sous la forme 1 2 8 24 '];
name='Selection des images à retravailler';
numlines=1;
defaultanswer={'1'};
options.Resize='on';
im_sel = (inputdlg(prompt,name,numlines,defaultanswer,options)); % fenêtre
if ~isempty(im_sel), im_sel = str2num(im_sel{1,1}); end %#ok<ST2NM>

%% Retravaillage des images

% Nombre d'images à retravailler
nbr_im_sel = length(im_sel);
% Différentes options
% str = {'Elargir le contour (+5)','Resserrer le contour (-5)','Elargir le contour grossierement (+20)',...
%     'Resserrer le contour grossierement (-20)','Elargir le contour finement (+1)','Resserrer le contour finement (-1)'};
str1 = 'Coef. de retouche du contour';
str2 = 'Elargir (>0) ou resserrer le contour (<0) ? Valeurs typiques: +/-5, grossierement +/-20, finement +/-1';

% Boucle : on parcoure les images à retravailler
for e=1:nbr_im_sel
    % Création de la figure ou apparait le contour
    h = figure;
    set(h,'Position',[5 40 Nbrcolonne_image Nbrcolonne_image],'Name','contour à retravailler');
    iptsetpref('ImshowBorder','tight')
    JJ2 = 1+Nbrcolonne_image*(im_sel(e)-1):Nbrcolonne_image*im_sel(e);
    imagesc(Contours_IN_positionbasse(II,JJ2)),...
        title(['contour de l'' image ' num2str(im_sel(e)) ' sur ' num2str(Nb_image)]); colormap(gray)
    pasencorebon = 1;

    % Tant que pasencorebon = 1, on relance la fenetre de retouchage
    while pasencorebon == 1
        % Fenetre de retouchage
        options.Resize='on';
        options.WindowStyle='normal';
        answer = inputdlg(str2,str1,1,{'+5'},options);
        
        haut = haut - str2double(answer)/1e4;
        
%         [elarg_resser,bouton1] = listdlg('name',['retouche de l image ' num2str(im_sel(e))],...
%             'ListSize',[250 82],...
%             'PromptString','Que voulez vous faire :',...
%             'SelectionMode','single',...
%             'CancelString','Abandonner',...
%             'ListString',str);

%         if bouton1 == 0  % Correspond à l'abandon
%             disp(['Le contour de l image ' num2str(im_sel(e)) ' sera considéré comme non-valide']);
%             pasencorebon = 0;
%             im_nonvalide(l) = 1; 
%             
%         elseif bouton1 == 1  % Si l'on demande une action :
%             switch elarg_resser
%                 case 1  % Elargir le contour
%                     haut = haut - 0.0005;
%                 case 2  % Resserrer le contour
%                     haut = haut + 0.0005;
%                 case 3  % Elargir le contour grossierement
%                     haut = haut - 0.002;
%                 case 4  % Resserrer le contour grossierement
%                     haut = haut + 0.002;
%                 case 5  % Elargir le contour finement
%                     haut = haut - 0.0001;
%                 case 6  % Resserrer le contour finement
%                     haut = haut + 0.0001;
%             end

               % On rappelle "contourIN_positionbasse" pour retrouver le
               % nouveau masque
            [Masque_IN_enbas(II,JJ2), Contours_IN_positionbasse(II,JJ2)] = ...
                contourIN_positionbasse(imread(files(im_sel(e)).name),...
                Taille_matricegauss,Taille_median,bas,haut);
            imagesc(Contours_IN_positionbasse(II,JJ2))
            title(['contour de l''image ' num2str(im_sel(e)) ' sur ' num2str(Nb_image)]), colormap(gray)
            % On demande si le nouveau masque affiché convient
            button2 = questdlg('Le contour convient-il ?',...
                'Continue Operation','Yes','No','Yes');
            if strcmp(button2,'Yes'), pasencorebon = 0;
            elseif strcmp(button2,'No'), pasencorebon = 1;
            end
%         end
    end % while pasencorebon == 1
    close(h);
end % for e=1:nbr_im_sel

close(f) % all % toutes les fenetres se ferment à la fin du programme

%% Création des masques contours

for l = 1:Nb_image
    Rdil = 10;
    % Je dilate de Rdil pixel le masque, paramètre fixé, 
    % sans justification préalable. Il se peut qu'il soit à modifier...
    JJ = 1+Nbrcolonne_image*(l-1):Nbrcolonne_image*l;
    Masque_IN_dilate(II,JJ) = imdilate(Masque_IN_enbas(II,JJ),strel('disk',Rdil)); %#ok
    
    % Pour chaque image, je prend le contour du masque IN.
    Masque_CONT_positionbasse(II,JJ) = bwperim(Masque_IN_dilate(II,JJ));

%     Masque_CONT_positionbasse(II,JJ) = ...
%         imdilate(Masque_CONT_positionbasse(II,JJ),strel('disk',Rdil));
end

%% Sauvegarde des images contourées et des masques

% On veut sauvegarder tous les masques de toutes les images dans le dossier :
% \\Pcmicrovideo\data\HIV\2007-01-25\DIC\masques en position basse\Masque IN, par exemple.
% Et sous le nom de fichier : MIB.name.tif . ou MCONT.name.tif 
%    M -> Masque
%    I -> Interieur   et CONT -> Contour
%    B -> Bas, masque en position basse

% On enregistre également tous les contours des masques dans le dossier :
% \\Pcmicrovideo\data\HIV\2007-01-25\DIC\masques en position
% basse\verification
% Les contours sont enregistré sous : CIB.name.tif .
%    C -> Contour
%    I -> Interieur
%    B -> Bas, masque en position basse

% Création des dossiers et des noms de masques et de contours
rajout_name_IN = 'MIB.';
rajout_name_CONT = 'MCONTB.';
rajout_name_contour = 'CIB.';

directory_name_IN = ['masques en position basse' filesep 'Masque IN' filesep] ;
directory_name_CONT = ['masques en position basse' filesep 'Masque CONT' filesep];
directory_name_contour = ['masques en position basse' filesep 'verification' filesep];
if isempty(dir('masques en position basse')),
    mkdir('masques en position basse');
    mkdir(directory_name_IN);
    mkdir(directory_name_CONT);
    mkdir(directory_name_contour);
end

% Pour 1 jusqu'au nombre d'image :
for l=1:Nb_image
    file_l = files(l).name;
%     if (im_nonvalide(l) == 0) % Si l'image était valide
        save_IN = strcat(rajout_name_IN,file_l);
        save_CONT = strcat(rajout_name_CONT,file_l);
        save_contour = strcat(rajout_name_contour,file_l);
%     end
    save_string_IN = [directory_name_IN,save_IN];
    save_string_CONT = [directory_name_CONT,save_CONT];
    save_string_contour = [directory_name_contour,save_contour];
    
    M = Masque_IN_enbas((1:Nbrligne_image),(1+Nbrcolonne_image*(l-1):Nbrcolonne_image*l));
    CONT = Masque_CONT_positionbasse((1:Nbrligne_image),(1+Nbrcolonne_image*(l-1):Nbrcolonne_image*l));
    C = Contours_IN_positionbasse((1:Nbrligne_image),(1+Nbrcolonne_image*(l-1):Nbrcolonne_image*l));
    C = uint16(C);
    
    % On sauve les images :
    imwrite(M,save_string_IN,'tif')
    imwrite(CONT,save_string_CONT,'tif')
    imwrite(C,save_string_contour,'tif')
end
%close(h), close(f)
if dirchanged, cd .., end

msgbox([' FIN DU PROGRAMME : Tous les masques ont été enregistrés dans le dossier ',dossier]);
%%%