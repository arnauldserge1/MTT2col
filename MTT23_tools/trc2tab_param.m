function tab_param = trc2tab_param(trc) % pk

global N_PARAM
if isempty(N_PARAM), MTTparams_def; end

Tmax = max(trc(:,2));
Ntrc = trc(end,1);
N_param_trc = 4; % n t i j %%%size(trc,2); % <=N_PARAM!!
% if N_param_trc>N_PARAM, disp('No more than N_PARAM params allowed!!!'); end
    
tab_param = zeros(N_PARAM*Tmax,Ntrc); 

for i=1:Ntrc
    indi = (trc(:,1)==i);
    trci = trc(indi,2:N_param_trc); % t i j
    tt = trci(:,1);
    tab_trc = ones(Tmax,N_PARAM);
    tab_trc(tt,:) = [trci ones(length(tt),N_PARAM-N_param_trc+1)]; % alpha = r0 = offset = blink = 1!!
    tab_trc = tab_trc';
    tab_param(:,i) = tab_trc(:);
end