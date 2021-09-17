function  [Msdout,FullMsdMatrix,MSDcompil] = ...
    msd (Trc, NewFullOutput, print_out, first_step_only, do_compil, MTT_format, use_sem)

% function  [Msdout,FullMsdMatrix,MSDcompil] = ...
%     msd (Trc, NewFullOutput, print_out, first_step_only, do_compil, MTT_format, use_sem)
%----------------------------------------------------------
%
% calculate the mean-square displacement from the given traces
%
%
% input:   Trc - matrix of particle traces as output from mktrace or MTT
%
%          NewFullOutput - if nonzero , FullMsdMatrix is a 
%               list instead of a matrix
%          if print_out = 1, display calcul progress (% done) (def.)
%          if first_step_only = 1, calcul only FullMsdMatrix{1}, ie 1st step (not def!)
%          AND MSDcompil
% def: msd (Trc, 1, 1, 0, 0, 0, 1)
%
% output:  Msd - mean-sqare displacement and it's standard
%                deviation for each particle
%          FullMsdMatrix - msd for all steps, nth column 
%                stands for n steps
%                If NewFullOutput is nonzero, the matrix is
%                replaced with a list, n-th entry containing
%                all steps of length n
%          MSDcompil = [n msd_n std(msd_n)/sqrt(n)] if 
%           with n=length(msd_n), number of msd elements for each step #n
%           caution: FullMsdMatrix={}, to save computer memory!
%
%
% date:    8.7.1994
% author:  ts
% version: <01.30> from <941014.0000>
% version: <02.00> from <021101.0000> by GAB reducing runtime
% version: <02.01> from <021126.0000> by GAB including full output
% version: <02.02> from <061123.0000> by AS including MSD compilation
%------------------------------------------------------------

global MSD_FRACTION
if isempty(MSD_FRACTION), MSD_FRACTION = 0.5; end % calcul only first 50%

if nargin<1 || isempty(Trc)
    Msdout = []; FullMsdMatrix = []; MSDcompil = []; %     help msd
    return 
end
if nargin<2, NewFullOutput = 1; end
if nargin<3, print_out = 1; end
if nargin<4, first_step_only = 0; end
if nargin<5, do_compil = 0; end
if nargin<6, MTT_format = 0; end
if nargin<7, use_sem = 1; end

if nargout>=3, do_compil = 1; end
% if nargin<5, if nargout>=3, do_compil = 1; else do_compil = 0; end, end

if MTT_format==1 %%%% format tab_param, modulo 7 => use msd(tab_param,1,1,0,0,1)
    Trc = detect_reconnex_to_trc(Trc);
end
if size(Trc,2)==2 %%%% just xy
    Trc = [ones(size(Trc,1),1),(1:size(Trc,1))',Trc];
    use_minimal_columns = 1;
else
    use_minimal_columns = 0;
end

if nargout>1    % only calculate full matrix if necessary
    DoFullMSD=1;
    if NewFullOutput    
        FullMsdMatrix = {};%cell(1,Nmax);???
    else
        FullMsdMatrix = [];
    end
    Nstep = 0; FMMlen = 0;
else
    DoFullMSD = 0;
end

%%% supprime trtaces orphelines (1 seul point) % AS 4/4/7
der1 = [1; diff(Trc(:,1))]; %derivee 1e colonne (ntrc)=1 si nvle trc
der2 = [diff(Trc(:,1)); 1]; %permut 
Trc(der1&der2,:) = []; %der1 et der20 => trc d'un seul point (der=[...1 1...])

if isempty(Trc) 
    Msdout = [0 0 0 0]; FullMsdMatrix = {}; MSDcompil = [0 0 0];
    help msd, return 
end

MaxPart = max(Trc(:,1));

Msdout = zeros(size(Trc,1)-Trc(end,1),4); 
currlgn = 1;
 
%------------------------------------------
%loop through all particles
if print_out
    disp('computing MSD...           ')
end
for Ipart=1:MaxPart % for each trace
    if mod(Ipart,50)==0 && print_out
        fprintf([repmat('\b',1,9) '%3.0f%% done'], 100*Ipart/MaxPart)
    end

    iTrc = Trc((Trc(:,1)==Ipart),2:4); % current trace: [image #, x, y]
    
    if ~isempty(iTrc) % if there are points in this trace
        Nlag = size(iTrc,1); % # of points in this trace
        Mlag = iTrc(Nlag,1)-iTrc(1,1); % length of trace in images
        %    H2  = []; %    H2 = zeros(Mlag*(Mlag-1)/2,Mlag);
        % 90 percent zeros => out of memory!
        H2 = zeros(Mlag,Mlag);
        indeX = ones(1,Mlag);
        
        %%% MSD from 1 to Nlag-2
        if Nlag>2
            for lag=1:Nlag-2
                L = iTrc(1+lag:Nlag,1)-iTrc(1:Nlag-lag,1); % delta T
                H = (iTrc(1:Nlag-lag,2:3)-iTrc(1+lag:Nlag,2:3)).^2; % delta x:y^2
                H = sum(H,2);
                for il=1:length(L)
                    H2(L(il), indeX(L(il))) = H(il); %if ~isnan(indeX(L(il))), end % AS 11/12/8
                    indeX(L(il)) = indeX(L(il))+1;
                end
            end
        end

        H  = (iTrc(1,2:3)-iTrc(Nlag,2:3)).^2;
        H2(Mlag,indeX(Mlag)) = sum(H);
        indeX(Mlag) = 2;
        
        if DoFullMSD || do_compil % someone really wants ALL the data (ou au moins compil..)
            if first_step_only && ~do_compil % AS 3/5/2006 -> 24/11/6 add do_compil (qui nécessite full msd, transitoirement)
                Nmax = 1;
            else
                Nmax = round(Mlag*MSD_FRACTION); % trace length * 50%
            end
            
            if NewFullOutput
                if Nstep<Nmax
                    [FullMsdMatrix{(Nstep+1):Nmax}] = deal([]); % fill them up with empties (zeros(1,Nmax)??)
                    Nstep=Nmax;
                end
                for i=1:Nmax
                    FullMsdMatrix{i} = [FullMsdMatrix{i},H2(i,1:(indeX(i)-1))]; 
                end
            else
                if Nstep<Nmax   % bugger, matrix too small
                    FullMsdMatrix = [FullMsdMatrix, zeros(FMMlen,Mlag-Nmax)];
                    Nstep = Nmax;
                end
                for i=1:Nmax
                    FullMsdMatrix((FMMlen)+(1:(indeX(i)-1)),i) = H2(i,1:(indeX(i)-1))';
                    FMMlen = FMMlen+indeX(i)-1;
                end
            end
        end        
        Lag=zeros(1,Mlag); Msd=zeros(1,Mlag); dMsd=zeros(1,Mlag); % AS 17/10/5, prealoc for speed
        ii = 1;
        for il=1:Mlag
            noZ = indeX(il)-1;
            H1 = H2(il, 1:noZ);
            if noZ>0

                Lag(ii)  = il;
                Msd(ii)  = mean(H1);

                if noZ>1
                    if use_sem
                        dMsd(ii) = std(H1)/sqrt(noZ);
                    else
                        dMsd(ii) = std(H1);
                    end
                else
                    dMsd(ii) = Msd(length(Msd));
                end
                ii=ii+1;
            end
        end
        if ii<=Mlag, Lag(ii:Mlag)=[]; Msd(ii:Mlag)=[]; dMsd(ii:Mlag)=[]; end % NB: full if ii=Mlag+1
        %MSD of that trace
        %Msdout = [Msdout; Ipart*ones(length(Lag),1),Lag',Msd',dMsd'];
        Msdout(currlgn:currlgn+length(Lag)-1,:) = [Ipart*ones(length(Lag),1),Lag',Msd',dMsd'];
        currlgn = currlgn+length(Lag);
    end   % if ~isempty(iTrc)
end % for each trace

% ***** compil *****
if do_compil && ~isempty(FullMsdMatrix)
    MSDcompil = zeros(length(FullMsdMatrix),3);
    for n=1:length(FullMsdMatrix) % =ntrc
        val = FullMsdMatrix{n};
        if use_sem
            MSDcompil(n,:) = [n, mean(val), std(val)/sqrt(length(val))];
        else
            MSDcompil(n,:) = [n, mean(val), std(val)];
        end
    end
    MSDcompil(isnan(MSDcompil(:,2)),:) = []; % AS 11/12/8
    if first_step_only 
        FullMsdMatrix = FullMsdMatrix{1}; % éco. mémoire vive
    end
end

if print_out, fprintf([repmat('\b',1,9) '100%% done\r']), end

if use_minimal_columns, Msdout = Msdout(:,3); end % just MSD, not n, dt, dMSD