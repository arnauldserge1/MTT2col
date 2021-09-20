function [v, f, d, t, N, corr_lin, sci_max] = compute_speed_from_SCI(input, do_plot, SCI_threshold, SCI_window, v_min)

% function [v, f, d, t, N, corr_lin, sci_max] = compute_speed_from_SCI(files, do_plot, SCI_threshold, SCI_window, v_min)
%
% velocity for each correlated episode (SCI > SCI_threshold, 0.6),
% computed in a sliding window (neighbor points or the trace) SCI_window = 30
% and discarded for v < v_min, 0.01um/s (too slow, irrelevant)
% f: fraction of time spent above threshold for each trace,
% d: duration and N: number of event/trace
% def: compute_speed_from_SCI('*.tif', 1, 0.6, 30, 0.01)
% see also SCI_Code


if nargin<1, input = '*.tif'; end
if isnumeric(input)
    tab_param = input;
    Nfiles = 1;
else files = dir(input);
    if isempty(files), files = dir('*.stk'); end
    Nfiles = length(files);
end
[pixel, timing] = get_calib3;

if nargin<2, do_plot = 1; end
if nargin<3, SCI_threshold = 0.6;  end % SCI_threshold2 = 0.6/300;SCI_threshold2 for test w/o norm, hence smaller values (about /300??)
if nargin<4, SCI_window = 30; end % times points
if nargin<5, v_min = 0.01; end % um/s

params_def = MTTparams_def; dirname = params_def{4};
min_trc_length = 5;
MSD_FRACTION = 0.5;

fit_directed = @(p, x, weight) ((p(1) + 4*p(2)*x + p(3)*x.^2).*weight); % = (2r02 + 4Dt + v2t2).*weight => yfit == fit_directed./weight
optfit = optimset('Display', 'off');
warning off Matlab:nearlySingularMatrix, warning off Matlab:SingularMatrix


v = cell(Nfiles, 1);
f = cell(Nfiles, 1);
d = cell(Nfiles, 1);
t = cell(Nfiles, 1);
N = cell(Nfiles, 1);
corr_lin = cell(Nfiles, 1);
sci_max = cell(Nfiles, 1);

if do_plot, figure('WindowStyle', 'docked'), end

for nf = 1:Nfiles
    if isnumeric(input)
        filename = '';
    else
        filename = files(nf).name;
        filename_full = [dirname filesep filename '_tab_param.mat'] ;
        tab_param = importdata(filename_full);
    end
    
    Ntrc = size(tab_param, 2);%    Ntrc = 10
    v{nf} = cell(Ntrc, 1);
    f{nf} = zeros(Ntrc, 1);
    d{nf} = cell(Ntrc, 1); % d{nf} = zeros(Ntrc, 1);
    t{nf} = cell(Ntrc, 1);
    N{nf} = zeros(Ntrc, 1);
    corr_lin{nf} = cell(Ntrc, 1);
    sci_max{nf} = cell(Ntrc, 1);
    
    fprintf('%s, trc           ', filename) %if do_plot, , end
    
    for nt = 1:Ntrc
        trc = detect_reconnex_to_trc(tab_param(:, nt));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        corr_lin{nf}{nt} = SCI_Code(trc, SCI_window);
% % %         corr_lin2{nf}{nt} = SCI_Code_no_norm(trc, SCI_window);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ind_lin = (corr_lin{nf}{nt} > SCI_threshold); % binary: linear (1) or not (0)
        trc_length = size(trc, 1);
        
        if do_plot &&  (trc_length >= SCI_window)
            clf
            subplot(221), plot(trc(:, 3)*pixel, trc(:, 4)*pixel), hold on % traj
            trc_begin = 1:floor(SCI_window/2+1);
            plot(trc(trc_begin, 3)*pixel, trc(trc_begin, 4)*pixel, 'k')
            trc_end = trc_length-floor(SCI_window/2)+1:trc_length;
            plot(trc(trc_end, 3)*pixel, trc(trc_end, 4)*pixel, 'k')
            xlabel('x (um)'), ylabel('y (um)'), axis equal, title(nt)
            
            subplot(223), plot((1:length(corr_lin{nf}{nt}))*timing, corr_lin{nf}{nt})
            xlabel('time (s)'), ylabel('Speed Correlation Index') % index SCI
            hold on, plot([trc_begin(end) trc_end(1)]*timing, [1 1]*SCI_threshold, ':')% SCI threshold
            plot(trc_begin*timing, trc_begin*0, 'k')
            plot(trc_end*timing, trc_end*0, 'k')
            axis tight
            set(gca, 'ylim', [-1 1])
            
% % %             subplot(325), plot((1:length(corr_lin2{nf}{nt}))*timing, corr_lin2{nf}{nt})
% % %             xlabel('time (s)'), ylabel('Speed Correlation Index') % index SCI
% % %             hold on, plot([trc_begin(end) trc_end(1)]*timing, [1 1]*SCI_threshold2, ':')% SCI threshold, no norm
% % %             plot(trc_begin*timing, trc_begin*0, 'k')
% % %             plot(trc_end*timing, trc_end*0, 'k')
% % %             axis tight 
        end
        
        start_lin = find(diff([0 ind_lin]) == 1);
        end_lin = find(diff([ind_lin 0]) == -1);
        Nlin = length(start_lin);
        if Nlin>0, fprintf(repmat('\b', 1, 11)), fprintf('%5i/%5i', nt, Ntrc), end % && do_plot
        vv = zeros(Nlin, 1); % speeds of linear events in current traj
        dd = zeros(Nlin, 1); % straight distance of linear events in current traj
        tt = zeros(Nlin, 1); % durations of linear events in current traj
        cc = zeros(Nlin, 1); % corr max of linear events in current traj
        
        %% loop over linear events
        for nl = 1:Nlin
            trc2 = trc(start_lin(nl):end_lin(nl), :);
            n = round((size(trc2, 1)-1)*MSD_FRACTION);
            
            if n > min_trc_length*MSD_FRACTION % 5*0.5 = 2.5
                msddata = msd(trc2, 1, 0); % = [n t r2 sd_r2] ntraj, nimage, msd, error
                tau = msddata(1:n, 2)*timing;
                r2 = msddata(1:n, 3)*pixel^2;
                weight = 1./msddata(1:n, 4); % 1/error
                if sum(weight) == 0
                    weight = ones(size(r2));
                end
                tau(weight == inf) = []; r2(weight == inf) = []; weight(weight == inf) = []; % si err = 0..
                
                n = length(r2);
                if n<=2, continue, end
                
                r0ini = r2(1) - (r2(2)-r2(1));
                p0 = [r0ini 0.1 0.1]; % p0 = [r0ini Dini vini]; dr, Din, dr/tau???
                %% fit MSD
                p = lsqcurvefit(fit_directed, p0, tau, r2.*weight, [0 0 0], [inf inf inf], optfit, weight); %         r0(nt) = sqrt(p(1))/2; % [i, sqrt(p(1))/2, err(1)];        %         D(nt) = p(2)*pxl_size^2/time_lag;        %         chi2(i) = sum((r2 - fit_directed(p, tt, weight)./weight).^2);
                vv(nl) = sqrt(p(3)); % in um/s, with r in um and t in s
                tau2 = start_lin(nl):end_lin(nl);
                tt(nl) = length(tau2)*timing; % duration, in s
                dd(nl) = sqrt((trc2(end, 3)-trc2(1, 3))^2 + (trc2(end, 4)-trc2(1, 4))^2)*pixel; % straight length d = sqrt(dx2 + dy2), in um
                cc(nl) = max(corr_lin{nf}{nt}(tau2));
                
                if do_plot
                    if vv(nl) < v_min, color = 'r'; % v slow
                    else color = 'g'; % v fast, relevant
                    end
                    subplot(221), plot(trc2(:, 3)*pixel, trc2(:, 4)*pixel, color) % linear event
                    subplot(122), plot(tau, r2, color), hold on % MSD
                    plot(tau, fit_directed(p, tau, weight)./weight, [color ':']);
                    xlabel('time (s)'), ylabel('MSD (um^2)')
                    title(['v = ' num2str(vv(nl)) ' um/s'])
                    
                    subplot(223), plot((tau2)*timing, corr_lin{nf}{nt}(tau2), color)
                    pause(0.1)
                    if vv(nl)>1, pause(1), end
                end
                
                if vv(nl) < v_min % remove events below v min: put back ind_lin at 0
                    ind_lin(tau2) = 0;
%                     vv(nl) = 0; % to be removed
                end
            end % if n > 2
        end % lin
        %         if (Nlin > 1) && (max(v{nf}{nt}) > v_min) && (min(v{nf}{nt}(v{nf}{nt} > 0)) < v_min), pause, end
        %         if (Nlin > 0), pause, end
        
%         v{nf}{nt} = v{nf}{nt}(v{nf}{nt} > 0); %keep only vals above v min
        v{nf}{nt} = vv;
        f{nf}(nt) = mean(ind_lin);
        d{nf}{nt} = dd; % d{nf}(nt) = sum(ind_lin);
        t{nf}{nt} = tt;
        start_lin2 = find(diff([0 ind_lin]) == 1); % ind_lin binaire, dc diff = 1 qd passe de 0 à 1 (=start) (& -1 pour end)
        N{nf}(nt) = numel(start_lin2);
        sci_max{nf}{nt} = cc;
         
%         pause
    end % trc
    
    fprintf('\r') % if do_plot, , end
    v{nf} = cell2mat(v{nf});
    d{nf} = cell2mat(d{nf});
    t{nf} = cell2mat(t{nf});
    sci_max{nf} = cell2mat(sci_max{nf});
    
end % files

v = cell2mat(v);
f = cell2mat(f);
d = cell2mat(d);
t = cell2mat(t);
N = cell2mat(N);
sci_max = cell2mat(sci_max);

%%%

% cd('/Volumes/IMAGERIE/Acquaviva/Data time lapse/2014/Spinning/2014-06-26/MATLAB-indiv-cell')
% [vFOP fFOP] = compute_speed_from_SCI('H4WF-FOP*.tif');
% [vGAPDH fGAPDH] = compute_speed_from_SCI('H4WF-GAPDH*.tif');
% scatter_hist({vFOP(vFOP>.01) vGAPDH(vGAPDH>.01)}, {'FOP', 'GAPDH'}, 'speed (um/s)', 1, 1, 1)

% name = cell(4,1); v = cell(4,1); f = cell(4,1); d = cell(4,1); N = cell(4,1);
% dirs=dir('2014*'); for nd = 1:4, name{nd} = dirs(nd).name; name{nd} = name{nd}(15:end); end
% for nd = 1:4, cd(dirs(nd).name), [v{nd}, f{nd}, d{nd}, N{nd}] = compute_speed_from_SCI; cd .., end
% save data v f d N
% figure, scatter_hist(v([1 4 2 3]), name([1 4 2 3]), 'speed', 1, 1, 1)
% figure, scatter_hist(f([1 4 2 3]), name([1 4 2 3]), 'f', 1, 1, 1)
% figure, scatter_hist(d([1 4 2 3]), name([1 4 2 3]), 'd', 1, 1, 1)
% figure, scatter_hist(N([1 4 2 3]), name([1 4 2 3]), 'N', 1, 1, 1)