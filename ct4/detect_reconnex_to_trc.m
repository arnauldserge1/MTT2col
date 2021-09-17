function [trc, pk, t_last] = detect_reconnex_to_trc(data_in, print_out, dirname, t_first, t_last, n_first, n_last)
% function [trc, pk, t_last] = detect_reconnex_to_trc(data_in, print_out, dirname, t_first, t_last, n_first, n_last)
% réordonne les paramètres des pics selon chaque convention
%
% Estimation/Reconnexion (transposé)
%               1    2  3      4       5         6          7      8
% tab_param' = [num, t, i,     j,     alpha,     rayon,     m0,    blink]
% tab_var' =   [num, t, sig_i, sig_j, sig_alpha, sig_rayon, sig_b, blink] = STD!!! sauf sig_b!!!
%
%       1  2     3     4  5         6           7   8   9   10  11   12 13 14    15
% pk = [t, x(j), y(i), w, I(alpha), offset(m0), dx, dy, dw, dI, do , 3tests(=>0) view]
%
%        1       2 3 4 5
% trc = [num_trc t x y pipo]

include_global
pars = MTTparams_def; %% MTT_param %if isempty(sig_free) || isempty(SAVE_VAR_MOY), , end % ; ???

if nargin<2, print_out = 1 ; end% affiche % effectué
if nargin<3, dirname = '' ; end
if isempty(dirname), dirname = pars{4} ; end %% disp(['Hey! we are using ' dirname])
if nargin<5,  t_first = 1; t_last = inf; end
if nargin<7,  n_first = 1; n_last = inf; end

if ~isempty(strfind(cd,dirname)), cd .., end % nb: on se retrouve svt déjà ds output si crash ou interupt

% if ischar(data_in) && isempty(dir(dirname))
%     dirname = 'output22' ; disp('no output23?? back to version 22!!!')
% end

%% valeurs algo Nico
filename = '';
if isempty(data_in)
    pk = []; trc = []; return % disp('He ben??')
    
elseif isstruct(data_in) || isnumeric(data_in)
    if isstruct(data_in)
        tab_param = data_in.tab_param;
    else
        tab_param = data_in;
    end
    if (mod(size(tab_param, 2), N_PARAM)==1) && isequal(tab_param(:, 1)',1:size(tab_param, 1)) % tab_param(1,end)==size(tab_param,2)% test ajouté le 26/3/9
        tab_param = tab_param';
        tab_param(1,:) = []; % entete
    end
    if t_last==inf
        t_last2 = size(tab_param, 1)/N_PARAM;
    else
        t_last2 = t_last;
    end % pb pour inf
    if n_last==inf
        n_last2 = size(tab_param, 2);
    else
        n_last2 = n_last;
    end
    ind = (N_PARAM*(t_first-1)+1:N_PARAM*t_last2);
    tab_param = tab_param(ind, n_first:n_last2);
    
    if (nargout>1) && (SAVE_VAR_MOY==1)
        if isstruct(data_in)
            tab_var = data_in.tab_var;
            if (mod(size(tab_var, 2), N_PARAM)==1) && isequal(tab_var(:, 1)', 1:size(tab_var, 1))
                tab_var = tab_var';
                tab_var(1,:) = [];
            end
            tab_var = tab_var(ind, n_first:n_last2);
        end
    end
    
elseif ischar(data_in)
    filename = data_in;
    t_last = check_mem(filename, dirname, t_last, n_last);
    %% read
    if isdir(dirname)
        cd (dirname)
        tab_param = fread_params_timewindow([filename '_tab_param.dat'], print_out, t_first, t_last, n_first, n_last) ;
        if (nargout>1) && (SAVE_VAR_MOY==1)
            tab_var = fread_params_timewindow([filename '_tab_var.dat'], print_out, t_first, t_last, n_first, n_last) ;
        end
        cd ..
    end
    
end

if isempty(tab_param) && (t_last==inf) && (n_last==inf)
    disp(['ya dleau dans lgaz: ' dirname filesep filename '_tab_param.dat est vide ou inexistant...'])
    pk = []; trc = [];
    return
end


nb_t = size(tab_param, 1)/N_PARAM ; % nombre d'image (= valeur temps max. 300 en gl)
nb_part_max = size(tab_param, 2) ;

tab_t = tab_param(PARAM_T-1:N_PARAM:end,:) ;
tab_i = tab_param(PARAM_I-1:N_PARAM:end,:) ;
tab_j = tab_param(PARAM_J-1:N_PARAM:end,:) ;

if nargout>1
    tab_alpha = tab_param(PARAM_ALPHA-1:N_PARAM:end,:) ; % amplitude du signal gaussien : alpha/(sqrt(pi)*rayon)) exp()
    tab_r = tab_param(PARAM_RADIUS_I-1:N_PARAM:end,:) ; % rayon gaussien (std) : exp(-r^2/2r0^2)
    if max(tab_r(:)) == 0, tab_r = tab_r + 1.1 ; end % not estimated in 3D!!
    tab_m = tab_param(PARAM_OFFSET-1:N_PARAM:end,:) ; % mean (=offset) %% =0, non eval...
    tab_sig2 = tab_param(PARAM_SIG2-1:N_PARAM:end,:) ; % mean (=offset) %% =0, non eval...
    
    if (SAVE_VAR_MOY==1)
        tab_sigij = tab_var(PARAM_I-1:N_PARAM:end,:) ;
        tab_sigalpha = tab_var(PARAM_ALPHA-1:N_PARAM:end,:) ; % sig = standard dev
        tab_sigr = tab_var(PARAM_RADIUS_I-1:N_PARAM:end,:) ;
%         tab_sigm = tab_var(PARAM_OFFSET-1:N_PARAM:end,:) ; %sqrt??????????
    end
end

tab_blink = tab_param(PARAM_BLINK-1:N_PARAM:end,:) ;

clear tab_param tab_var % AS 28/5/7

%% correspondance trc
t = tab_t(:) ;
% if isstruct(data_in)
%     x = tab_i(:) ;
%     y = tab_j(:) ;
% else
x = tab_j(:) ; % sic
y = tab_i(:) ;
% end

if nargout>1
    w = 2*sqrt(2*log(2)) * tab_r(:) ; % FWHM = 2.3548 r...
    I = floor(2*sqrt(pi) * tab_r(:).* tab_alpha(:)) ; % I = 3.5449 r alpha (!r, alpha = vecteurs!) I~1e5 et a~3e4
    I = I + tab_alpha(:)/1e6; % signal (I) codé ds partie entiere et bruit (dalpha) codé ds partie fractionnaire !!
    % intensité/puissance intégrée sur le pic (=volume du pic, et non max!)
    % sum(G_Leid(:))=I, et sum(G_Nico(:)^2)=alpha^2 !!! def. sur amplitude / puissance (??)
    o = tab_m(:); % offset, mean...
    sig2 = tab_sig2(:);
    
    zr = zeros(nb_t*nb_part_max, 1) ; % ou boutonnière de soutane à toto...
    uns = ones(nb_t*nb_part_max, 1) ; % ...
    
    %% Erreurs resp.
    if (SAVE_VAR_MOY==1)
        dx = tab_sigij(:) ;
        dy = dx ;
        dw = 2*sqrt(2*log(2))*tab_sigr(:) ;
        %dI = 2*sqrt(pi) * sqrt(alpha.^2.*dr.^2+dalpha.^2.*r.^2);
        dI = 2*sqrt(pi) * tab_r(:) .* tab_sigalpha(:) ;
        % (r.*dalpha... formule de propagation des erreurs, approxim, négl possible crossvariance)
%         do = tab_sigm(:) ; % var_offset === bruit;
    else
        dx = sig_free*uns;
        dy = sig_free*uns;
        dw = w/2; % cf init_tab
        dI = I/2;
        %%% do = o/2; % !!! bruit bidon !!!
    end
end

bad = tab_blink(:)<=0 ;

clear tab_t tab_i tab_j tab_blink
if (nargout>1) && (SAVE_VAR_MOY==1)
    clear tab_sigij tab_alpha tab_sigalpha tab_r tab_sigr tab_m %tab_sigm
end

n = ones(nb_t,1) * (1:nb_part_max) + (n_first-1);
n = n(:) ;

%% matrice trc
%size(n),size(t),size(x)
trc =[n t x y zeros(size(n))] ;

trc(bad,:) = [] ; %% supprimes val. nulles (blink, apparitions)

%% matrice pk
if nargout>1
    pk = [t x y w I o dx dy dw dI sig2 zr zr uns zr] ; % pk(14)=>ok ou non...
    pk(bad,:) = [] ;
    clear w I o dx dy dw dI sig2 zr uns
end

clear n t x y bad
%%%