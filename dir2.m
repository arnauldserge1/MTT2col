function listing = dir2(filename, dirs_only)
% function files = dir2(filename, dirs_only)
% same as dir, but get rid of files starting by '.'
% and can keep only dirs, not files (if dirs_only = 1, def 0)
% AS 28/2/2017

if (nargin < 1), filename = '*'; end
if (nargin < 2), dirs_only = 0; end

listing = dir(filename);
Nf = length(listing);
ok = true(Nf, 1);

for nf = 1:Nf
    item = listing(nf).name;
    if strcmp(item(1), '.')
        ok(nf) = false;
    end
    if dirs_only && ~isdir(item) % 3/4/2018
        ok(nf) = false;
    end
end

listing = listing(ok); % Nf = length(files);

%%%