function TinToutMTT23

% TinTout - testconnect: test valeurs maxblink vs trajlength, cf DinDout (?)

clear global Nb_STK tab_num

if strfind(cd,'Claire'), param_opt = 3;
else param_opt = 0; end

T1 = 10; T2 = 100; nsteps = 15;
step = (log(T2)-log(T1))/(nsteps-1);
Tin = - [1:9, round(exp(log(T1):step:log(T2)))];% Tin = 3:7;

output_dir = 'output23_DinDout'; %..
figure('WindowStyle','docked')

p = MTTparams_def(param_opt);
files = dir(p{1});
p{2} = '1'; % refit
p{4} = output_dir;
p{10} = '300'; % Nb_STK
p{20} = ''; % conf method
Tdef = eval(p{13});

for f=1:length(files)
    Tout = zeros(size(Tin));
    dTout = Tout;
    Dout = Tout;
    dDout = Tout;
    Ntrc = Tout;
    
    filename = files(f).name;
    p{1} = filename;
    
    for i=1:length(Tin)
        p{13} = num2str(Tin(i));
        
        MTT23i(p)
        trc = detect_reconnex_to_trc(filename,1,output_dir);
        
        if size(trc)>0
            msddata = msd(trc,0);
            Dinst = calculDinst(msddata); % = D1-5
            Dvect = Dinst(Dinst(:,2)>0,2);
            Dout(i) = exp(mean(log(Dvect))); % moyenne sur les trajs
            dDout(i) = exp(std(log(Dvect))); %/sqrt(length(Diff)) %dD = mean(dDvect)
            
            TT = trc(:,2);
            dN = diff(trc(:,1));
            Tstart = [TT(1); TT(find(dN>0)+1)];
            Tend = [TT(dN>0); TT(end)];
            TrcLen = Tend-Tstart+1; % >1 ??
            Tout(i) = exp(mean(log(TrcLen))); % moyenne sur les trajs
            dTout(i)= exp(std(log(TrcLen))); %/sqrt(length())
            Ntrc(i) = trc(end,1);
        end
        
        subplot(121)
        loglog(abs(Tin),Ntrc,'g.-',abs(Tin),Tout,'b.-',abs(Tin),Tout.*dTout,'b^',...
            abs(Tin),Dout,'k.-',abs(Tin),Dout.*dDout,'k^','MarkerSize',8)
        hold on
        a = axis;
        loglog(abs([Tdef Tdef]), a(3:4),'r:')
        hold off
        xlabel('T_{blink} (frame)'); ylabel('trc length, <D_{1-5}>, N_{trc}');
        title({cd filename 'T_in T_out'},'interpreter','none')
        %     legend({'N_{trc}','trc length','sd','<D_{1-5}>','sd'})
        
        subplot(122)
        x = abs(Tin(1:i-1));
        semilogx(x,diff(Ntrc(1:i))/max(-diff(Ntrc)),'g.-',...
            x,diff(Tout(1:i))/max(diff(Tout)),'b.-',...
            x,diff(Dout(1:i))/max(diff(Dout)),'k.-')
        hold on
        a = axis;
        semilogx(abs([Tdef Tdef]), a(3:4),'r:')
        hold off
        legend({'N_{trc}','trc length','<D_{1-5}>'})
        xlabel('T_{blink} (frame)'); ylabel('delta');
        if i<length(Tin), title(['Tin(',num2str(i),') = ',num2str(Tin(i))]), end
       
        pause(.1)
    end % i
    savefile = ['TinTout_' filename(1:end-4)];
    cd(output_dir)
    saveas(gcf, savefile, 'fig')
    saveas(gcf, savefile, 'png')
    cd ..
end % f
%%%