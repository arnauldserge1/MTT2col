function B = SCI_Code_no_norm(trc, n)
% function B = SCI_Code_no_norm(trc, n)
% Correl_Speed_Index

if nargin<2, n = 30; end
% if nargin<3, verbose = 0; end

A = trc(:, [3 4 2]);

% Ntrc = trc(end, 1);
% B = cell(1, Ntrc);
% if verbose, fprintf('fitting Correl_Speed_Index           '), end

% for nt = 1:Ntrc
%     A = trc2(trc(:, 1) =  = nt, :);
    
    Sl = size(A, 1);
    %     Sl = S(1);    % Sc = S(2);    % Temps(1, :) = pas*A(:, 3)';
    Tab = 999*ones(n, Sl-n-1);
    
    
    for start = 0:n-1
        
        %------------------Smoothing---------------------------
        
        Inter = [];
        %     Deltax = [];    %     Deltay = [];    %     DL = [];
        
        k = start+1;
        i = 1;
        
        while k <= Sl
            Inter(i, :) = A(k, :); %#ok
            i = i+1;
            k = k+n;
        end
        
        if isempty(Inter), continue, end
        
        %----------------------Correlations------------------------
        
        vx = diff(Inter(:, 1)) ./ diff(Inter(:, 3));
        vy = diff(Inter(:, 2)) ./ diff(Inter(:, 3));
        
        Sx = size(vx);
        SxL = Sx(1);
        
        C = vx(1:SxL-1).*vx(2:SxL) + vy(1:SxL-1).*vy(2:SxL);
% % % % % %         N1 = sqrt(vx(1:SxL-1).*vx(1:SxL-1) + vy(1:SxL-1).*vy(1:SxL-1));
% % % % % %         N2 = sqrt(vx(2:SxL).*vx(2:SxL) + vy(2:SxL).*vy(2:SxL));
% % % % % %         N = N1.*N2;
% % % % % %         C = Ct./N;
        
        SDLL = size(C, 1);
        if isempty(C), SDLL = 0; end
        
        %------------------Different reading frames --------------
        
        i = 1+start;
        while i <= Sl
            ir = fix((i-1-start)/n)+1;
            if ir <= SDLL
                Tab(start+1, i) = C(ir);
                i = i+1;
            else
                i = Sl+1;
            end
        end
    end % start
    
    %-------------------Averaging--------------------
    
    B = zeros(1, Sl); %{nt} zeros(1, Sl-n-1);
    for m = 1:Sl-n-1
        Temp = 0;
        compt = 0;
        
        for j = 1:n
            tp = Tab(j, m);
            
            if tp == 999
                %            Temp = Temp;                %            compt = compt;
            else
                Temp = Temp+tp;
                compt = compt+1;
            end
            B(1, m+floor(n/2)+1) = Temp/compt; %{nt}% if n = 2m, m+1 0 before and m after, if n = 2m+1, m+1 0 before and m+1 after
            
        end
    end % m
    
    %-----------------------------------------------------------------------------
    
    %     Res(1, :) = ;    Res(2, :) = Moy;
%     B{nt} = [A(1:Sl-n-1, 3)'; Moy]; % Moy>0.6??
    %     subplot(121), plot(A(:, 1), A(:, 2)), axis equal, title(nt)
    %     subplot(122), plot(A(1:Sl-n-1, 3)', Moy), set(gca, 'ylim', [-1 1]), pause
%     if verbose, fprintf([repmat('\b', 1, 11) '%5i/%5i'], nt, Ntrc), end
end % trc