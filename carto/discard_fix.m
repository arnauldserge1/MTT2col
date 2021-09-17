function [trc_mobiles trc_immobiles] = discard_fix(trc, filename, Dmin)
% selectionne et sauve traces mobiles, ie D>Dmin
% AS 2/9/8

% way_out = ['carto_v3_output22\' filename '_mobiles'];
% if ~isempty(dir(way_out))
%     disp('sélection déjà faite')
%     return
% end

if isempty(trc), trc = detect_reconnex_to_trc(filename); end %, 1, dirname);

timelag = eval_timing(filename);
if nargin<3, Dmin = 3.5e-3*timelag/36; end % pxl/lag % disp('Using Dmin=3.5e-3pixel/lag')

trc_mobiles = trc; % alloc
if nargout>1, trc_immobiles = trc; end% alloc
Ntrc = trc(end,1);

for i=1:Ntrc
    trci = trc(trc(:,1)==i,:);
    if size(trci,1)>5 % size msd>=5
        Di = calculDinst(msd(trci,1,0));
    else Di = [0+0 inf]; % N'EXCLUE PAS TRAJS COURTES (??)
    end
    
    if Di(2)<Dmin %%% CHECK DIFF MIN
        trc_mobiles(trc_mobiles(:,1)==i,:) = [];
%          fprintf('trc %i/%i discarded, D = %g\r', i, Ntrc, Di(2))
    else
        if nargout>1, trc_immobiles(trc_immobiles(:,1)==i,:) = []; end
    end
end

oldNtrc = sum(diff(trc(:,1))~=0); % >< Ntrc si il manque des n°...
newNtrc = sum(diff(trc_mobiles(:,1))~=0);
fprintf(' %i/%i trc discarded, D < %g\r', oldNtrc-newNtrc, oldNtrc, Dmin)

% save way_out trc_mobiles -ascii