function modificationmatrice_parampositionbasse(matrice_a_modifier_moy,matrice_a_modifier_param,matrice_a_modifier_var,masquecorrespondant,dossier)

%% Explication :

% Cette fonction permet de ressortir trois matrices dans lesquelles les pics
% étant à l'exterieur du masaque correspondant sont mis à zero...

%% Cretion des trois nouvelles matrices

mkdir(dossier,'matrice modifié suivant un MIB'); % Creation du dossier
dossier_destination = [ dossier,filesep,'matrice modifié suivant un MIB' ];

nouvelle_matrice_moy = [ dossier_destination,filesep,'MIB',matrice_a_modifier_moy];
nouvelle_matrice_param = [ dossier_destination,filesep,'MIB',matrice_a_modifier_param];
nouvelle_matrice_var = [ dossier_destination,filesep,'MIB',matrice_a_modifier_var];

ancienne_matrice_moy = [dossier,filesep,matrice_a_modifier_moy];
ancienne_matrice_param = [ dossier,filesep,matrice_a_modifier_param];
ancienne_matrice_var = [ dossier,filesep,matrice_a_modifier_var];

fid_moy = fopen(nouvelle_matrice_moy, 'wt','native') ; % On ouvre un nouveau fichier texte
fid_param = fopen(nouvelle_matrice_param, 'wt','native') ;
fid_var = fopen(nouvelle_matrice_var, 'wt','native') ;

%% Initialisation des fichiers textes

fprintf(fid_moy, '############# :  (nb maxi particules, nb snapshots)\n') ;
fprintf(fid_moy, '# DATA_SPT : %s\n', matrice_a_modifier_moy) ;
fprintf(fid_moy, '# DATA_SPT : %s', date) ;
clk = clock() ;
fprintf(fid_moy, ' %.2dh%.2dm%.2ds\n', clk(4), clk(5), round(clk(6))) ;
fclose(fid_moy) ;

fprintf(fid_param, '############# :  (nb maxi particules, nb snapshots)\n') ;
fprintf(fid_param, '# DATA_SPT : %s\n', matrice_a_modifier_param) ;
fprintf(fid_param, '# DATA_SPT : %s', date) ;
clk = clock() ;
fprintf(fid_param, ' %.2dh%.2dm%.2ds\n', clk(4), clk(5), round(clk(6))) ;
fclose(fid_param) ;

fprintf(fid_var, '############# :  (nb maxi particules, nb snapshots)\n') ;
fprintf(fid_var, '# DATA_SPT : %s\n', matrice_a_modifier_var) ;
fprintf(fid_var, '# DATA_SPT : %s', date) ;
clk = clock() ;
fprintf(fid_var, ' %.2dh%.2dm%.2ds\n', clk(4), clk(5), round(clk(6))) ;
fclose(fid_var) ;


%% Boucle

fid_old_moy = fopen(ancienne_matrice_moy, 'rt','native') ;
fid_old_param = fopen(ancienne_matrice_param, 'rt','native') ;
fid_old_var = fopen(ancienne_matrice_var, 'rt','native') ;

value =  fscanf(fid_old_param, '%d', 2) ;
nb_part_max = value(1);
nb_image_video = value(2) ;

fscanf(fid_old_moy, '%d', 2) ;
fscanf(fid_old_var, '%d', 2) ;

  %% on se cale sur le debut des donnees
ligne_old_moy = fgets(fid_old_moy) ; 
ligne_old_param = fgets(fid_old_param) ;
ligne_old_var = fgets(fid_old_var) ; 

while (strcmp(ligne_old_moy(1:14), '# NEW_DATA_SPT') == 0)
    ligne_old_moy = fgets(fid_old_moy) ;
end%while

while (strcmp(ligne_old_param(1:14), '# NEW_DATA_SPT') == 0)
    ligne_old_param = fgets(fid_old_param) ;
end%while

while (strcmp(ligne_old_var(1:14), '# NEW_DATA_SPT') == 0)
    ligne_old_var = fgets(fid_old_var) ;
end%while

for t=1:nb_image_video
    
    for p=2:8
      value =  fscanf(fid_old_param, '%d:', 2) ;
      nb_part = value(1) ;
      ligne_old_param = fgets(fid_old_param); 
 	  tab_out_old_param(p-1,(1:nb_part)) = sscanf(ligne_old_param(2:end), '%f', nb_part);  %#ok<AGROW>
      
      fscanf(fid_old_moy, '%d:', 2);
      if p == 5
          ligne_old_moy = fgets(fid_old_moy); 
 	      tab_out_old_moy_temp = sscanf(ligne_old_moy(2:end), '%f - %f i');   
          for temp = 1:nb_part
              tab_out_old_moy(p-1,temp) = tab_out_old_moy_temp(2*temp-1) - tab_out_old_moy_temp(2*temp)*1i; %#ok<AGROW>
          end
      elseif p ~= 5
          ligne_old_moy = fgets(fid_old_moy) ;
 	      tab_out_old_moy(p-1,(1:nb_part)) = sscanf(ligne_old_moy(2:end), '%f', nb_part);  %#ok<AGROW> 
      end
      
      fscanf(fid_old_var, '%d:', 2) ;
      ligne_old_var = fgets(fid_old_var) ;
 	  tab_out_old_var(p-1,(1:nb_part)) = sscanf(ligne_old_var(2:end), '%f', nb_part) ; %#ok<AGROW>
      
    end%for
    
    for numero_pic = 1:nb_part
        if masquecorrespondant(int8(tab_out_old_param(3,numero_pic)),int8(tab_out_old_param(2,numero_pic))) == 0
            tab_out_old_moy((4:7),numero_pic) = 0; %#ok<AGROW>
            tab_out_old_param((4:7),numero_pic) = 0; %#ok<AGROW>
            tab_out_old_var((4:7),numero_pic) = 0; %#ok<AGROW>
        end
    end
    
    fid_moy = fopen(nouvelle_matrice_moy, 'at','native') ;
    fprintf(fid_moy, '# NEW_DATA_SPT : %d particules : %s\n', nb_part, matrice_a_modifier_moy) ;
    fid_param = fopen(nouvelle_matrice_param, 'at','native') ;
    fprintf(fid_param, '# NEW_DATA_SPT : %d particules : %s\n', nb_part, matrice_a_modifier_param) ;
    fid_var = fopen(nouvelle_matrice_var, 'at','native') ;
    fprintf(fid_var, '# NEW_DATA_SPT : %d particules : %s\n', nb_part, matrice_a_modifier_var) ;
    
    for param = 2:8 ; %% t,i,j,alpha,...,blink
       remplir_moy = num2str( (tab_out_old_moy(param-1, :)) , 5);
       fprintf(fid_moy, '%d:%d: ', nb_part, param) ;
       fprintf(fid_moy, '%s\n', remplir_moy) ;
       remplir_param = num2str( (tab_out_old_param(param-1, :)) , 5);
       fprintf(fid_param, '%d:%d: ', nb_part, param) ;
       fprintf(fid_param, '%s\n', remplir_param) ;
       remplir_var = num2str( (tab_out_old_var(param-1, :)) , 5);
       fprintf(fid_var, '%d:%d: ', nb_part, param) ;
       fprintf(fid_var, '%s\n', remplir_var) ;
    end %for
    fclose(fid_moy) ;
    fclose(fid_param) ;    
    fclose(fid_var) ;    
    
    ligne_old_moy = fgets(fid_old_moy) ; %#ok<NASGU>
    ligne_old_param = fgets(fid_old_param) ; %#ok<NASGU>
    ligne_old_var = fgets(fid_old_var) ; %#ok<NASGU>
    
end%for
 
fclose(fid_old_moy) ;
fclose(fid_old_param) ;
fclose(fid_old_var) ;
   
fid_moy = fopen(nouvelle_matrice_moy, 'r+t','native') ;
fprintf(fid_moy, '%.6d %.6d', nb_part_max, nb_image_video) ;
fclose(fid_moy) ;

fid_param = fopen(nouvelle_matrice_param, 'r+t','native') ;
fprintf(fid_param, '%.6d %.6d', nb_part_max, nb_image_video) ;
fclose(fid_param) ;

fid_var = fopen(nouvelle_matrice_var, 'r+t','native') ;
fprintf(fid_var, '%.6d %.6d', nb_part_max, nb_image_video) ;
fclose(fid_var) ;
   
end