function compil_fig

dir0 = 'Z:\Vincent Rouger\20100906 QDot';
Npuits = 10;
varname = 'z';

figure('windowstyle', 'docked', 'PaperOrientation', 'landscape', ...
    'PaperUnits', 'normalized', 'PaperPosition', [0.02 0.02 0.98 0.98])

for i = 1:Npuits % puit
    cd(sprintf('%s\\Puit %i\\Sans EGF\\carto',dir0, i))
    
    f1 = compil_dir(varname);
    
    cd(sprintf('Puit %i\\avec EGF\\carto', i))
    
    f2 = compil_dir(varname);
    
    fc = cell2mat([f1; f2]);
    imshow(fc)
    cd(dir0)
    saveas(gcf, sprintf('carto_compil%d_%s.pdf', i, varname))
    pause(.1) % figure(gcf)
end

%%%

function f = compil_dir(varname) %(Nstk, f, n)

Nstk = 7;
f = cell(1,Nstk);
disp(cd)

for j = 1:Nstk
    file = sprintf('Stack%i.stk_%s_3D.png', j, varname);
    if ~isempty(dir(file))
        f{j} = imread(sprintf('Stack%i.stk_%s_3D.png', j, varname));
        f{j} = f{j}(200:750,80:1120,:);
        fprintf('%i ',j)
    else
        f{j} = uint8(zeros(750-200+1,1120-80+1,3))+255;
    end
end
fprintf('\r')

%%%

% % % f = cell(2,6);
% % % figure('windowstyle','docked')
% % % 
% % % for i=1:6 % puit
% % %     cd('\\amenophis\matlab$\Vincent Rouger\20100713 stim EGF 3D')
% % %     cd(sprintf('Puit_%i\\Controle\\carto', i))
% % %     disp(cd)
% % %     for j = 1:6 % stk
% % %         f{1,j} = imread(sprintf('Stack%i.stk_r2_3D.png', j));
% % %         f{1,j} = f{1,j}(200:750,80:1120,:);
% % %     end
% % %     
% % %     cd('\\amenophis\matlab$\Vincent Rouger\20100713 stim EGF 3D')
% % %     cd(sprintf('Puit_%i\\Experience\\carto', i))
% % %     disp(cd)
% % %     for j = 1:6 % stk
% % %         f{2,j} = imread(sprintf('Stack%i.stk_r2_3D.png', j));
% % %         f{2,j} = f{2,j}(200:750,80:1120,:);
% % %     end
% % %     
% % %     cd('\\amenophis\matlab$\Vincent Rouger\20100713 stim EGF 3D')
% % %     fc = cell2mat(f');
% % %     imshow(fc)
% % %     saveas(gcf, sprintf('carto_compil%d_r2_.pdf', i))
% % %     figure(gcf)%, pause
% % % end