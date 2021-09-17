function visites_map = revisites(trc,lg,ht)
% place sur la carte les zones détectées plusieurs fois comme confinée,
% en utilisant le critère 1 pixel(160nm) ~ lambda/2

if isempty(trc), visites_map =[]; return, end
if nargin<2
    if max(trc(:,3))>128, lg = 512; 
    elseif  max(trc(:,3))<=100
    else lg = 128; 
    end
    ht = lg;
end

visites_map = zeros(lg,ht,3);
%time_visite = cell(taille);

[Lconf lenConf lenFree R2Conf R2Free trc_conf trc_free] = probaconf(trc, 0);
clear Lconf lenConf lenFree R2Conf R2Free

for k=1:3
    if k==1, dat = trc_conf;
    elseif k==2, dat = trc_free;
    elseif k==3, dat = trc;
    end
    Ntrc = dat(end,1);
    for i=1:Ntrc
        trci = dat(dat(:,1)==i,:);
        ic = round(trci(:,4));
        jc = round(trci(:,3));
        visites_map(ic,jc,k) = visites_map(ic,jc,k)+1;
        %time_visite{ic,jc} = trci(:,2);
    end
end
%%%