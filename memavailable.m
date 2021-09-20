function [num,str] = memavailable
%MEMAVAILABLE Returns the amount of memory available 
%   [NUM,STR] = MEMAVAILABLE returns the numerical number of bytes in NUM, 
%   and a string representation in STR. 

% Get the results of DUMPMEMMEX 

if ismac
    
    num = inf; %......
    
elseif ~strcmp(computer,'PCWIN64')
    
    smem = evalc('dumpmemmex');
    
    % Extract tokens
    m = regexp(smem,'(?<bytes>\d*)\sbytes[^\(]*\((?<mb>[^\)]*)\).*$','names');
    num = str2double(m.bytes);
    str = m.mb;
    
else % Windows 64 bits
    
    user = memory;
    num =  user.MaxPossibleArrayBytes; % MemAvailableAllArrays ;% 16*2^30; % inf; => 16 Go

end