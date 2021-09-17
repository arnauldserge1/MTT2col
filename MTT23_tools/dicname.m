function DIC = dicname(filename, verbose)

% function DIC = dicname(filename, verbose)
%%% convention usuelle : blabla cell1a.stk => DIC\blabla cell1.tif
% donc 5(par ex.) caratères à changer à la fin (1:end-5) ex:  a.stk => .tif
%%% OU blabla.stk => DIC\blabla.tif =>(1:end-4)

if nargin < 2, verbose = 1; end

if ~ischar(filename), DIC = ''; return, end

letter_removed = 0;
% if useDIC
DIC = ['dic', filesep, filename(1:end-4) ,'.tif'];
if isempty(dir(DIC))
    DIC = ['trans', filesep, filename(1:end-4) ,'.tif'];
end
if str2double(DIC(end-5))>0 && ~isempty(strfind('abcde',DIC(end-4))) && isempty(dir(DIC))
    DIC(end-4) = []; % xxx[chiffre lettre].tif, ex: toto1a.tif, pipo 3c.tif...
    letter_removed = 1;
end

% *** bas ou milieu... ***
m = strfind(DIC,'mil');
if m>0, DIC(m:m+2) = []; end
b = strfind(DIC,'bas');
if b>0, DIC(b:b+2) = []; end

m = strfind(DIC,'m.tif');
if m>0, DIC(m) = []; end
b = strfind(DIC,'b.tif');
if b>0, DIC(b) = []; end
t = strfind(DIC,'t.tif');
if t>0, DIC(t) = []; end

d = strfind(DIC,'_dessous');
if d>0, DIC(d-1:d+7) = []; end
d = strfind(DIC,'_dessus');
if d>0, DIC(d-1:d+6) = []; end
f = strfind(DIC,'_focus');
if f>0, DIC(f-1:f+5) = []; end

if ~isempty(strfind(cd,'HIV\2007-07-24')) || ~isempty(strfind(cd,'HIV\2007-07-17')) || ~isempty(strfind(cd,'HIV\2007-07-11'))
    if m>0, DIC = [DIC(1:end-4) 'm.tif']; end
    if b>0, DIC = [DIC(1:end-4) 'b.tif']; end
end

% *** '1a' pas à la fin... ***
if strcmp(cd, '\\Pcdm100906\datae\Axel\20060329 Biotine sans r_phénol')
    DIC(10) = [];
elseif strcmp(cd, '\\Pcdm100906\datae\Axel\20060602 Drogues')
    DIC = ['dic', filesep, filename(1:end-4) ,'_av.tif']; DIC(10)=[];
end

if (~isempty(strfind(cd,'2006-07-19 CO+latB')) || ...
        ~isempty(strfind(cd,'2005-06-03 LatB'))) ...
        && ~letter_removed %useDIC
    if strcmp(filename(1:2),'CO'), DIC(8) = [];
    elseif strcmp(filename(1:3),'ctl'), DIC(9) = [];
    elseif strcmp(filename(1:4),'latB'), DIC(10) = [];
    end
end

if isempty(dir(DIC))
    if strfind(DIC,'cell'),
        if strcmp(DIC(end-8),' '), DIC = [DIC(1:end-9) '.tif']; end % cellNN NNms.stk
        if strcmp(DIC(end-7),' '), DIC = [DIC(1:end-8) '.tif']; end% cellNN Nms.stk
    end
end


if isempty(dir(DIC))
    DIC = '';
    if ~isempty(dir('dic')) && verbose
        disp('DIC image not found..')
    end
else
    if verbose, disp(['DIC image : ' DIC]), end
end

%%%