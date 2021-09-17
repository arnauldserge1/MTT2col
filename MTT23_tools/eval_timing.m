function timing = eval_timing(filename)
% cherche une chaine de caractere indiquant le timing utilisé pour un
% fichier donné, sur 1 ou plusieurs chiffres suivis de 'ms' pour millisecondes
% ex: eval_timing('toto 20ms.stk')=20
% ! Pour cell_track, filename = 'toto_filt1.tif' et timing = 7000 ms = 7s !

default_timing = 100;%36;
% % % if ~isempty(strfind(cd,['IMAGERIE' filesep 'Arnauld'])), default_timing = 100; end
if ~isempty(strfind(cd,['Geneve' filesep 'Alexandre'])), default_timing = 20; end

if ~ischar(filename), timing = default_timing; fprintf('using default: %i ms', default_timing), return, end

k = strfind(filename,'ms');
if ~isempty(k)
    k = k(1); % si plusieurs 'ms'...
%     digits = blanks(k-1);
    for i = 1:k-1
        digit = filename(k-i);
        if isnan(str2double(digit)), timing = str2double(filename(k-i+1:k-1)); break, end
        if (i==k-1), timing = str2double(filename(k-i:k-1)); end
    end
%     if k-3>0
%         timing = str2double(filename(k-3:k-1)); %'100 ms...
%         if isnan(timing), timing = str2double(filename(k-1)); end
%     if k-2>0
%         timing = str2double(filename(k-2:k-1)) %'_5' 10 20 ms...
%         if isnan(timing), timing = str2double(filename(k-1)); end
%     else % 1 seul chiffre
%         timing = str2double(filename(k-1)); %5 3 ms...
%     end
else
    timing = default_timing;
    fprintf('using default: %i ms for %s\r', default_timing, filename)
end
timing = abs(timing); %...

if ~isempty(strfind(cd,'HTHLAB')) || ~isempty(strfind(filename,'_filt')), timing = 7000; end % 7s pour manips calcium

% % % r = bfGetReader(filename); % get metadata assez long..
