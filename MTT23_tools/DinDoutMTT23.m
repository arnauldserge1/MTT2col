function DinDoutMTT23(D1, D2, nsteps, param_opt, Nb_STK)

% plot Dout versus Din

% AS 19/12/5

clear global Nb_STK tab_num

% if nargin<3
%      if strfind(cd,'Claire'), param_opt = 3; else param_opt = 0; end
% end
p = MTTparams_def(param_opt);
filename_def = p{1};

% if strfind(cd,'Claire-Lise')
%     filename_def = '*_filt1.tif';
% end

files = dir(filename_def);
Nfiles = length(files);

if nargin<1 % full test  
    D1 = 1e-2; D2 = 1; nsteps = 28;
%     D1 = 1e-3; D2 = 10; nsteps = 28;
    Nb_STK = 100;%     D1 = 1e-4, D2 = 1, nsteps = 9, param_opt = 2 % STORM!!
elseif (nargin==1) && (D1==1) % quick test
    D1 = 1e-3; D2 = 10; nsteps = 9;
    Nb_STK = 30;
    Nfiles = 1;
end

% if strfind(cd,'Avais'), Nb_STK = 11; end

output_dir = 'output23_DinDout'; % output is erased at each D value

step = (log(D2)-log(D1))/(nsteps-1);
Din = exp(log(D1):step:log(D2)); % valeurs d'entrée, en pxl/frame
Din2 = (Din(1:end-1)+Din(2:end))/2;
Doptim = zeros(Nfiles,1);

for nf = 1:Nfiles
    filename = files(nf).name;
    filename_filt = [filename(1:end-4) '_filt1.tif'];
    if ~isempty(dir(filename_filt)), filename = filename_filt; end
    Dout = zeros(size(Din));
    r2 = Dout; dr2 = Dout;
    
    figure('WindowStyle','docked')
    
    for i = 1:length(Din)
        p = MTTparams_def(param_opt);
        Ddef = eval(p{8});
        
        p{1} = filename;
        p{2} = '1'; % refit
        p{4} = output_dir;
        
        p{8} = num2str(Din(i));
        
        p{10} = num2str(Nb_STK);
        
%         if strfind(cd,'Claire-Lise')
%             p{13} = num2str(0); %Toff
%             p{16} = num2str(4000); %seuil_alpha
%         end        
        
        p{20} = ''; % conf method
        
        MTT23i(p)
        tab_param = fread_all_params([output_dir filesep filename '_tab_param.dat']);
% % %         tab_param = importdata([output_dir filesep filename '_tab_param.mat']);

        if ~isempty(tab_param)
            r2_tab = calcul_r2(tab_param);
            r2(i) = exp(mean(log(r2_tab(r2_tab>0))));
            dr2(i) = exp(std(log(r2_tab(r2_tab>0))));
        end
        
        clf
        subplot(1,2,1)
        title({cd filename 'D_in D_out'},'interpreter','none')
        Dout = r2/4;
        dDout = Dout.*dr2;
        loglog(Din,Dout,'r.-',Din,dDout,'r^','MarkerSize',8)           
        hold on
        a = axis;
        loglog([Ddef Ddef], a(3:4), 'r:')
        xlabel('D_{in} (pxl2/frame)'); ylabel('D_{out} (pxl2/frame)')
        
        subplot(1,2,2)
        deltaD = diff(log10(Dout)); % deltaD = 10.^diff(log10(Dout));
        semilogx(Din2, deltaD, 'r.-')
        xlabel('D_{in} (pxl2/frame)'); ylabel('slope: delta D_{out} (pxl2/frame)');
        hold on
        if i>1
%             axis tight
            a = axis;
            axis([Din(1) Din(end) a(3:4)])
            loglog([Ddef Ddef], a(3:4), 'r:')
        end
        pause(.1)
    end % D
    
    lDin2 = log10(Din2);
    w = max(log10(deltaD))-log10(deltaD); % weight for ponderated mean
    Doptim(nf) = 10^(sum(lDin2.*w)/sum(w)); % optim = find(deltaD==min(deltaD),1);

    subplot(1,2,2), hold on
    semilogx([Doptim(nf) Doptim(nf)], a(3:4), 'g:')
    subplot(1,2,1), hold on
    str = sprintf('Optimum: Din=%3.2g pxl2/frame',Doptim(nf));
    text(.1,.9,str,'Units','normalized')
    a = axis;
    loglog([Doptim(nf) Doptim(nf)], a(3:4), 'g:')
    
    savefile = ['DinDout_' filename(1:end-4)];
    cd(output_dir)
    saveas(gcf, savefile, 'fig')
    saveas(gcf, savefile, 'png')
    cd ..
end % files

m = exp(mean(log(Doptim)));
sd = exp(std(log(Doptim)));
fprintf('Average optimum: Din = %3.2g+/-%3.2g pxl2/frame',m,sd);
%%%