function Lbin = MTT_conf_index_rel(column_param, conf_params)

global N_PARAM PARAM_BLINK


if nargin<2
    s = which('MTT_conf_index_rel');
    param_dir = fileparts(s);
    conf_params = load([param_dir '\proba_conf_params_rel.dat']);
end

seuil = conf_params(1);
Nb = conf_params(2);
Na = conf_params(3);

r2 = calcul_r2(column_param);
r2(isnan(r2)) = 0; %%%% disp('remplace Nan par zeros; foireux???')
blink = column_param(PARAM_BLINK:N_PARAM:end)<=0;
blkr2 = (blink(1:end-1) | blink(2:end));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test = detect_conf_lent(r2,Na,Nb,seuil,blkr2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(test)
    test(1:Na) = 0;
    test(end-Nb:end) = 0; % on squizze les bords
end

Lbin = [test'; 0]; % rem r2 non def pour dernier point, pas de detection possible
Lbin = conv(Lbin, ones(1,Nb)); % 1 point allumé => Nb points 'confs' (3 => Nb+2, etc)
Lbin = Lbin(1:end-Nb+1)~=0;