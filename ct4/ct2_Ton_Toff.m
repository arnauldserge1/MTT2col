function [mTon, sdTon, mToff, sdToff, fraction, Ton, Toff] = ct2_Ton_Toff(ntrcdata)

% function [mTon, sdTon, mToff, sdToff, fraction, Ton, Toff] = ct2_Ton_Toff(ntrcdata)
% 
% ct2_Ton_Toff
% -
% calcule temps allume et temps eteind moyen pour des
% molecules qui blinkent, a partir de trajs reconnectees

% -
% Version 1.0 AS 27/5/4

ntraces = ntrcdata(end, 1); % nombre de trajs
Ton = cell(1, ntraces); Toff = cell(1, ntraces); fraction = zeros(1, ntraces);

for i = 1:ntraces
    trci = ntrcdata(ntrcdata(:,1)==i, 2)'; % on ne garde que les # d'image, la 2e colonne, de la trace i
    if ~isempty(trci)
        dtrci = diff(trci); % increment: 1 = normal, >1 = blink
        blink = find(dtrci>1);
        if isempty(blink) % traj continue, sans blink ou non reconnectee
            Ton{i} = length(trci);
            fraction(i) = 1;
        else
            Ton{i} = blink-[1 blink(1:end-1)];
            Toff{i} = dtrci(dtrci>1);
            fraction(i) = sum(blink-[1 blink(1:end-1)])/length(trci);
        end
    end
end
Ton = cell2mat(Ton);
Toff = cell2mat(Toff);
fraction = fraction(fraction>0);

mTon = mean(Ton);
sdTon = std(Ton);
if isempty(Toff)
    mToff = NaN;
    sdToff = NaN;
else
    mToff = mean(Toff);
    sdToff = std(Toff);
end

% that's all folks