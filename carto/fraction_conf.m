function [res f] = fraction_conf

dd = {'M:\Arnauld Serge\2005-09-29 SPT specif',...
    'M:\Arnauld Serge\2005-06-03 LatB',...
    'M:\Arnauld Serge\2007-02-27 gamme temps acquis'};
ff = {'cell*.stk', 'ctl*.stk', 'cell*36ms 300img.stk'};

Tc = zeros(length(dd),8);
Tf = Tc; ac = Tc; 
f = cell(length(dd),8);

for i=1:length(dd)
    cd (dd{i})
    files = dir(ff{i});
    for j=1:length(files) % 8 8 5
        f{i,j} = files(j).name;
        trcj = detect_reconnex_to_trc(files(j).name);
        [Lconf lenConf lenFree ] = probaconf(trcj,0);
        Tc(i,j) = exp(mean(log(lenConf)));
        Tf(i,j) = exp(mean(log(lenFree)));
        ac(i,j) = Tc(i,j)/(Tc(i,j)+Tf(i,j));
    end
end
Tc = Tc'; Tf = Tf'; ac = ac'; f=f';
res = [Tc(:) Tf(:) ac(:)];