function lt_batchsong_plotFF(DATSTRUCT, MotifsToExtract, TrainON, SwitchTimes, subtractMean, ...
    dozscore)

%% LT 11/20/17 - plots

%% INPUTS

% TrainON = '05Nov2017-1125'; WN onset

% MotifsToExtract = {'a(b)', 'j(b)', 'ab(h)', 'jb(h)',  'jbh(h)', '(g)'};


% SwitchTimes = {'05Nov2017-1235', '05Nov2017-1355', '05Nov2017-1548', ...
%     '05Nov2017-1811'}; % will places lines in plot at these times

% subtractMean = 0; (baseline mean)


if ~exist('dozscore', 'var')
    dozscore = 0; % if 1, then zscore rel to baseline std.
end

%%

TrainON_dnum = datenum(TrainON, 'ddmmmyyyy-HHMM');

tval_min = [];
for i=1:length(DATSTRUCT.motif)
    tval_min = min([tval_min min([DATSTRUCT.motif(i).rendnum.datenum_song_SecRes])]);
end

firstday = datestr(tval_min, 'ddmmmyyyy');

%%
figcount=1;
subplotrows=4;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for i = 1:length(MotifsToExtract)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(MotifsToExtract{i});
    
    % ####################################### UNDIR
    inds = [DATSTRUCT.motif(i).rendnum.isDIR]==0;
    plotcol = 'k';
    
    % ----------------- RUN
    ffvals = [DATSTRUCT.motif(i).rendnum(inds).ff];
    tvals = [DATSTRUCT.motif(i).rendnum(inds).datenum_song_SecRes];
    
    % ------- subtract baseline FF
    baseinds = tvals < TrainON_dnum;
    if subtractMean==1
        ffvals = ffvals - nanmean(ffvals(baseinds));
    end
    if dozscore ==1
        basemean = nanmean(ffvals(baseinds));
        basestd = nanstd(ffvals(baseinds));
        ffvals = (ffvals - basemean)./basestd;
    end
    
    % -- convert tvals to days from start
    tvals = lt_convert_EventTimes_to_RelTimes(firstday, tvals);
    tvals = tvals.FinalValue;
    
    % --- HITS/MISSES
    if isfield(DATSTRUCT.motif(i).rendnum(inds), 'isWNhit')
        ishit = [DATSTRUCT.motif(i).rendnum(inds).isWNhit];
    else
        ishit = [];
    end
    
    % ----- plot
    plot(tvals, ffvals, 'o', 'Color', plotcol);
    if ~isempty(ishit)
        plot(tvals(ishit), ffvals(ishit), 'or');
    end
    
    
    % ----- plot day means
    numdays = floor(max(tvals));
    for j=1:numdays
        
        indstmp = floor(tvals)==j;
        tt = tvals(indstmp);
        ff = ffvals(indstmp);
        
        if isempty(ff)
            continue
        end
        
        lt_plot(max(tt)+0.1, mean(ff), {'Errors', lt_sem(ff), 'Color', plotcol});
    end
    
    % ######################### lines
    % --- train onset
    tmp = lt_convert_EventTimes_to_RelTimes(firstday, TrainON_dnum);
    line([tmp.FinalValue tmp.FinalValue], ylim, 'Color', 'r');
    
    for j=1:length(SwitchTimes)
        tmp = lt_convert_EventTimes_to_RelTimes(firstday, datenum(SwitchTimes{j}, 'ddmmmyyyy-HHMM'));
        line([tmp.FinalValue tmp.FinalValue], ylim, 'Color', 'm');
    end
    lt_plot_zeroline;
    
    
    % ####################################### DIR
    inds = [DATSTRUCT.motif(i).rendnum.isDIR]==1;
    if any(inds)
        plotcol = 'b';
        
        % ----------------- RUN
        ffvals = [DATSTRUCT.motif(i).rendnum(inds).ff];
        tvals = [DATSTRUCT.motif(i).rendnum(inds).datenum_song_SecRes];
        
        % ------- subtract baseline FF
        if subtractMean==1
            baseinds = tvals < TrainON_dnum;
            if ~any(baseinds)
                disp('TO PLOT DIR, DO NOT NORM TO BASE!!! [lackign base data]');
            end
            ffvals = ffvals - mean(ffvals(baseinds));
        end
        if dozscore ==1
            % ------- use base mean and std for UNDIR (to maintain ability to
            % compare)
            ffvals = (ffvals - basemean)./basestd;
        end
        
        % -- convert tvals to days from start
        tvals = lt_convert_EventTimes_to_RelTimes(firstday, tvals);
        tvals = tvals.FinalValue;
        
        
        % ====== HITS/MISSES
        
        
        % ----- plot
        lt_plot(tvals, ffvals, {'Color', plotcol});
        %      plot(tvals, ffvals, 'o', 'Color', plotcol);
        
        % ----- plot day means
        numdays = floor(max(tvals));
        for j=1:numdays
            
            indstmp = floor(tvals)==j;
            tt = tvals(indstmp);
            ff = ffvals(indstmp);
            
            if isempty(ff)
                continue
            end
            
            lt_plot(max(tt)+0.15, mean(ff), {'Errors', lt_sem(ff), 'Color', plotcol});
        end
    end
    
end
linkaxes(hsplots, 'xy');

