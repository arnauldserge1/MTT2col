function [PDr PFAr PDa PFAa] = test_MTT_probaconf(n, param, conf_params_abs, conf_params_rel)
% see loop_simul_conf

global N_PARAM
if isempty(N_PARAM), MTTparams_def; end

[dirname simul_name Nfiles Nimg ImSz Nppi gwidth I sig_I SNR_dB offset D ...
    do_bleach do_blink Ton Toff do_conf Tfree Tconf conf_factor do_graph] ...
    = read_simul_param(param);

% cd(dirname)
trc = load(['trcsimul\' simul_name '_file' num2str(n) '.tif.trc']);
tab_param = trc2tab_param(trc);% tab_param = fread_all_params([simul_name '_file' num2str(n) '.tif']);

if ~do_graph
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Lr = MTT_probaconf(tab_param, 'rel', 36, 0, '', conf_params_rel);
    La = MTT_probaconf(tab_param, 'abs', 36, 0, '', conf_params_abs);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Lr = Lr(:,1:Nppi); La = La(:,1:Nppi);
else
    trcsimul = load(['trcsimul\' simul_name '_file' num2str(n) '.tif.trc']);
    tab_param_simul = trc2tab_param(trcsimul);
    Ls = [tab_param_simul(11:N_PARAM:end,:); zeros(1,Nppi)]; % N_PARAM!!!
    
    for t=1:Nppi
        column_param = tab_param(:,t);
        r2 = calcul_r2(column_param);
        
        Lr = MTT_probaconf(column_param, 'rel');
        La = MTT_probaconf(column_param, 'abs');
        
        clf, hold on
        plot(r2)
        plot(Ls(:,t)-1.1)
        plot(Lr-2.2)
        plot(La-3.3)
        ylabel('abs, rel, sim, r^2'), xlabel('frames')
        pause
    end
end

% PFAr = sum(Lr(:)==Ls(:)+1)/sum(Ls(:)==0); % Note: si Ls=1 (vrai alarme), pas de stat: Lr ne peut en tous cas pas etre à 2 (0 ou 1)
val_test_PFA = (Tfree+2*Tconf:2*Tfree-Tconf); %(Tconf:Tfree-Tconf); % safe, free part of simul (3:9)

PFAr = sum(sum(Lr(val_test_PFA,:)==1))/numel(Lr(val_test_PFA,:));
PDr = sum(Lr(Tfree+1,:))/size(Lr,2);

PFAa = sum(sum(La(val_test_PFA,:)==1))/numel(La(val_test_PFA,:));
PDa = sum(La(Tfree+1,:))/size(La,2);

% if do_graph
%     subplot(1,2,1)
%     imagesc(Lr-Ls)
%     title('rel') % axis([0 10 0 100]+0.5)
%     fprintf('rel: %g%% PFA, %g%% PD\r',100*PFAr,100*PDr)
%     
%     subplot(1,2,2)
%     imagesc(La-Ls)
%     title('abs') % axis([0 10 0 100]+0.5)
%     fprintf('abs: %g%% PFA, %g%% PD\r',100*PFAa,100*PDa)
%     pause
% end