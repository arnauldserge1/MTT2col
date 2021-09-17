function modificationmatricepositionbasse

%% Explication :

% Cette fonction appelle une fenetre qui permet de selectionner
% trois matrices dans un dossier de mnanip comme Par exemple : 
% \\Pcmicrovideo\dataD\HIV\2007-05-23\output21\ .

% A partir de la, le programme cherche le masque correspondant;

% Puis par appel de modificationmatrice_parampositionbasse va créer les
% nouvelles matrices.


%% Recupération des trois matrices à reformater en fonction du masque

a = 1;
while a == 1
    % Tant que les matrices ne sont pas valides, une fenetre s'ouvre dans
    % laquelle l'utilisateur doit sélectionner trois matrices auxquelles il
    % veut applique un masque
    [A,dossier_output] = uigetfile({'*.dat','text Files'},...
        'Selectionner les 3 fichiers textes de données (moy, param, var)','\\Pcmicrovideo\DataD\HIV\','MultiSelect', 'on');
    if dossier_output == 0
        disp('abadon')
        break 
    end
    
    if length(A) ~= 3
        h = msgbox('Il y a trois matrices à selectionner : tab_moy, tab_param et tab_var','problème');
        uiwait(h);
    elseif length(A) == 3
        a = 0;
        distance = findstr(A{1},'.stk_tab'); % numéro du caractère à partir duquel on peut lire .stk...
        for n = 1:3
        MOY = findstr(A{n},'.stk_tab_moy.dat'); % On regarde si .stk_tab_moy.dat est présent dans A(n)
                                                % , et si oui, la fonction findstr renvoie le rang correspondant
        PARAM = findstr(A{n},'.stk_tab_param.dat');
        VAR = findstr(A{n},'.stk_tab_var.dat');
            if MOY == distance % si le rang est la distance, alors
                matrice_a_modifier_moy =  A{n}; % La matrice _moy est bien celle la.
            end
            if PARAM == distance
                matrice_a_modifier_param =  A{n};
            end
            if VAR == distance
                matrice_a_modifier_var = A{n};
            end
        end
    end
end

%% Recupération du masque

% On retrouve le nom du fichier
nom_fichier = regexprep(matrice_a_modifier_moy,'.stk_tab_moy.dat', '');
% On trouve le dossier de manip
dossier = regexprep(dossier_output, [filesep 'output22' filesep], '');
% On en déduit le dossier des images DIC.
dossier_masque = [dossier,'DIC' filesep 'masques en position basse' filesep];
% On ouvre une fenetre demandant de sélectionner le masque. Pour aider, on
% fait apparaitre le nom du fichier dans l'intitutlé de la fenetre.
[masquecorrespondant,dossier_masque] = uigetfile({'*.tif','masque binaires'},...
        ['Selectionner le masque correspondant aux trois fichiers :',nom_fichier],dossier_masque);
masquecorrespondant = imread([dossier_masque,masquecorrespondant]);
    
%% Appel de modificationmatrice_masquecorespondant

% Appel de modificationmatrice_parampositionbasse pour réécrire les
% nouvelles matrices.
modificationmatrice_parampositionbasse(matrice_a_modifier_moy,matrice_a_modifier_param,matrice_a_modifier_var,masquecorrespondant,dossier_output);

disp('Les trois matrices modifiées en fonction du masque ont été enregistrées dans le dossier "\output22\matrice modifié suivant un MIB"');

end