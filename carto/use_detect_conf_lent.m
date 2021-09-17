function     [Lconf lenConf lenFree R2Conf R2Free trc_conf trc_free beginConf endConf Dlocal] = ...
    use_detect_conf_lent(trc,graph)

seuil = 4.5 ; 
% Nb = 3 ;
Na = 10 ;

Lconf = cell(1,Ntraj);
lenConf = cell(1,Ntraj);
lenFree = cell(1,Ntraj);
R2Conf = cell(1,Ntraj);
R2Free = cell(1,Ntraj);
beginConf = cell(1,Ntraj);
endConf = cell(1,Ntraj);
trc_conf = [];
trc_free = [];
Dlocal = cell(1,Ntraj);

Ntrc = trc(end,1);
if graph, figure('WindowStyle','docked'), end

for i=1:Ntrc
    trci = trc(trc(:,1)==i,:); 
    if size(trci,1)>Na 
%         t = (trci(1:,2))'; % on transpose!
        dx = diff(trci(:,3))';
        dy = diff(trci(:,4))';
        r2 = dx.^2+dy.^2;
        
        [test, glrt] = detect_conf_lent(r2);%,Na,Nb,seuil,t);
        
        Lconf{i} = glrt;
        
        if graph
            plot([r2; test-2; glrt-8]')
            disp('next?')
            pause%(.1)
        end
    end
end

global USECELL
if isempty(USECELL), USECELL = 0; end

if ~USECELL
    lenConf = cell2mat(lenConf);
    lenFree = cell2mat(lenFree);
    R2Conf = cell2mat(R2Conf);
    R2Free = cell2mat(R2Free);
    beginConf = cell2mat(beginConf');
    endConf = cell2mat(endConf');
end