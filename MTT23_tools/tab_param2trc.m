function [trc, pk] = tab_param2trc(tab_param)%, tab_var)
% cf. detect_reconnex_to_trc, idem, but simpler

global N_PARAM PARAM_ALPHA

if isempty(N_PARAM), MTTparams_def; end

%     x = Trc(3:N_PARAM:end,:); y = Trc(4:N_PARAM:end,:); I = Trc(5:N_PARAM:end,:);
%     Trc = [ones(sum(I>0),1),(1:sum(I>0))',x(I>0),y(I>0)];

    Tmax = size(tab_param,1)/N_PARAM;
    Ntrc = size(tab_param,2);
    
    nn = ones(Tmax,1) * (1:Ntrc);
    pk = [nn(:) reshape(tab_param,N_PARAM,[])']; % n t i j alpha radius offset blk
    alpha_pos = pk(:,PARAM_ALPHA)>0; % PARAM_ALPHA = 5; compensated by column shift with nn (instead of PARAM_ALPHA-1!!)
    
    trc = pk(alpha_pos,[1,2,4,3]); % n t j i