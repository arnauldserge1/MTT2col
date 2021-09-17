function PFAinPFAout(PFA1, PFA2, step)%, param_opt)

% AS 2014

clear global Nb_STK tab_num

if strfind(cd,'Claire'), param_opt = 3;
else param_opt = 0; end

p = MTTparams_def(param_opt);
filename_def = p{1};
files = dir(filename_def);

if nargin<1 % full test
    PFA1 = 16; PFA2 = 40; step = 1;
    Nb_STK = 10;
    Nfiles = length(files);
elseif (nargin==1) && (PFA1==1) % quick test
    PFA1 = 20; PFA2 = 32; step = 4;
    Nb_STK = 3;
    Nfiles = 1;
end

output_dir = 'output23_DinDout'; % ..

PFAin = PFA1:step:PFA2; % valeurs d'entrée, en pxl/frame
% % PFAoptim = zeros(Nfiles,1);

for nf = 1:Nfiles
    filename = files(nf).name;
    
    N_out = zeros(size(PFAin));
%     dN = N_out;
    
    figure('WindowStyle','docked')
    
    for i = 1:length(PFAin)
        p = MTTparams_def(param_opt);
        PFAdef = eval(p{5});
        
        p{1} = filename;
        p{2} = '1'; % refit
        p{4} = output_dir;
        
        p{5} = num2str(PFAin(i));
        p{6} = num2str(PFAin(i)+4); % final = pre-threshold + 4 (?)
        
        p{10} = num2str(Nb_STK);
        p{20} = ''; % conf method
        
        MTT23i(p)
        tab_param = fread_all_params([output_dir filesep filename '_tab_param.dat']);
        
        if ~isempty(tab_param)
            N_tab = tab_param(4:8:end,:); % alpha = mean(N_tab(N_tab>0));
            N_out(i) = sum(N_tab(:)>0); % /Nb_STK; % dN(i) = std(N_tab(N_tab>0));
        end
        dN = diff(N_out);
        
        clf
        subplot(1,2,1)
        title({cd filename 'PFA_in N_out'},'interpreter','none')
        plot(PFAin,N_out,'r.-')%,PFAin,N_out+dN,'r^','MarkerSize',8)
        hold on
        a = axis;
        plot([PFAdef PFAdef], a(3:4), 'r:')
        xlabel('PFA_{in}'); ylabel('N_{out}')
        
        subplot(1,2,2)
        plot(PFAin(1:end-1), dN, 'r.-')
        xlabel('PFA_{in}'); ylabel('slope: delta N_{out}');
        hold on
        if i>1
            axis tight
            a = axis;
            axis([PFAin(1) PFAin(end) a(3:4)])
            plot([PFAdef PFAdef], a(3:4), 'r:')
        end
        pause(.1)
    end
    
% %     PFAoptim(nf) = find(dN==max(dN),1); % slope neg for N!
% %     
% %     subplot(1,2,2), hold on
% %     plot([PFAoptim(nf) PFAoptim(nf)], a(3:4), 'g:')
% %     subplot(1,2,1), hold on
% %     str = sprintf('Optimum: PFAin=%3.2g',PFAoptim(nf));
% %     text(.1,.9,str,'Units','normalized')
% %     a = axis;
% %     plot([PFAoptim(nf) PFAoptim(nf)], a(3:4), 'g:')
    
    savefile = ['PFA_in_N_out_' filename(1:end-4)];
    cd(output_dir)
    saveas(gcf, savefile, 'fig')
    saveas(gcf, savefile, 'png')
    cd ..
end % files

% % m = mean(PFAoptim);
% % sd = std(PFAoptim);
% % fprintf('Average optimum: PFAin = %3.2g+/-%3.2g',m,sd);
%%%