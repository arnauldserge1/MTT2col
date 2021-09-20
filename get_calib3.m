function [pixel, timing] = get_calib3
% function [pixel, timing] = get_calib3
% read pixel (um) and timing (s) values saved in output23/MTTparams.txt

fid = fopen('output23/MTTparams.txt', 'r');

if fid == -1
    par_def = MTTparams_def;
    pixel = eval(par_def{24});
    timing = eval(par_def{25});
    fprintf('Caution, no MTTparams.txt found! Using default values: pixel = %s um, timing = %s s\r', par_def{24}, par_def{25})
    return
end

tline = fgets(fid);
n = 1;
while ischar(tline)
    if contains(tline, 'pixel size ')
        pixel = str2double(tline(length('pixel size (um) : '):end-1));
        tline2 = fgets(fid);
        timing = str2double(tline2(length('time lag (s) : '):end-1));
        break
    end
    tline = fgets(fid);
    n = n+1;
end

fclose(fid);