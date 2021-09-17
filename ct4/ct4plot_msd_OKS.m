function D = ct4plot_msd_OKS(msd_compil, timelag, clr, fit_lin, fit_ano, text_offset)
% function D = ct4plot_msd(msd_compil, timelag, clr, fit_lin, fit_ano, text_offset)
% def: D = ct4plot_msd(msd_compil, 1, 'b', 1, 1, 0)
% see also ct4

if nargin<6, text_offset = 0; end
if nargin<5, fit_lin = 1; end
if nargin<4, fit_ano = 1; end
if nargin<3, clr = 'b'; end
if nargin<2 % if timelag == -1...
    use_units = 0;% % %     unit_strings = {'(img' 'pxl²'};
    timelag = 1;
else
    use_units = 1;
end

if ~isempty(msd_compil)
    r2 = msd_compil(msd_compil(:, 2)>0, :); %in  pxl2 or um2
    tmax = size(r2, 1);
    
    nmax = round(tmax/2); % MSD_FRACTION, rem: already in msd.m => 25%!
    tt = (1:nmax)*timelag; % if use_units
    hold on
    plot(tt, r2(1:nmax, 2), 'color', clr)
%     plot(tt, r2(1:nmax, 2)+r2(1:nmax, 3), 'color', clr, 'linestyle', ':')
%     plot(tt, r2(1:nmax, 2)-r2(1:nmax, 3), 'color', clr, 'linestyle', ':')
    
    if use_units, xlabel ('Time(s)'), ylabel('MSD (\mum^2)')
    else, xlabel ('t (img)'), ylabel('MSD (\mum^2)') % um???????
    end
    hold on
    
    % calcul D et RD, AS 8/3/5
    D = calculDinst(r2); % = [n D dD o]
    ttf = (1:round(tmax/4))*timelag;
    if fit_lin
        offset = D(4);
        if use_units
            D(2:3) = D(2:3)/timelag; % if use_units, msd already in um2 => D in um2/s
            s1 = sprintf('D_{inst} = %4.2g +/- %4.2gµm²/s', D(2:3));
        else
            s1 = sprintf('D_{inst} = %4.2g +/- %4.2gpxl²/img', D(2:3));
        end
        if D(2) > 0, plot([ttf(1) ttf(end)], 4*D(2)*[ttf(1) ttf(end)] + offset, 'color', clr, 'linestyle', '--'), end % 'c'??
        
        RD = calculRD(r2(:, 2));
        if isempty(RD), s2 = 'trc too short for RD';
        else, s2 = sprintf('RD = %4.2g', RD); end
    else
        s1 = ''; s2 = '';
    end
    
    % fit anomal
    if fit_ano
        [Da, gamma, offset_a] = fit_anomal2(r2);%, offset, 0);
        
        if use_units
            Da = Da/timelag;
            s3 = sprintf('D_{anom}=%4.2gum�/s', Da);
        else
            s3 = sprintf('D_{anom}=%4.2gpxl�/img', Da);
        end
        s4 = sprintf('\\gamma=%4.2g', gamma);
        
        if Da > 0, plot(ttf, 4*Da*ttf.^gamma + offset_a, 'color', clr, 'linestyle', '--'), end % fit anomal 'm'??
    else, s3 = ''; s4 = '';
    end
    if (text_offset == 0), txtclr = 'k'; else, txtclr = clr; end
    if (text_offset >= 0), text(0.05, 0.6+text_offset, {s1 s2 s3 s4}, 'Units', 'normalized', 'color', txtclr), end %  text(0.1, 0.7, {[s1 ' | ' s2] [s3 ' | ' s4]}, 'Units', 'normalized')
    D = D(2);
else
    D = -1;
    text(0.25, 0.9, 'No Traces', 'Units', 'normalized');
end