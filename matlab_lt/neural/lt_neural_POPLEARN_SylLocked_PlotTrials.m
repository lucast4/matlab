function [DATSTRUCT_SPK, DATSTRUCT_LFP, PARAMS] = ...
    lt_neural_POPLEARN_SylLocked_PlotTrials(OUTSTRUCT, OUTSTRUCT_XCOV, SwitchStruct, ...
    SwitchCohStruct, MOTIFSTATS_pop, SummaryStruct, ...
    PARAMS, onlygoodexpt, xtoplot, plotraw, zscore_lfp, fpass, ...
    sta_wind, kernelSD, removeBadLFP)
%% lt 2/8/19 - plots all syl locked LFP and smoothed MU

epochtoplot = 'base_epoch';
% fpass = [10 200]; % for bandpass filtering LFP.

Ntrialsplot =2; %^ for each experiments/motif

% ============== SPIKE TRIGGERED LFP
% --- how much flank time to get?
% sta_wind = [-0.05 0.05]; % relative to spike, in sec % will only get spikes that are within data...
disp('note that on average spike will be in real life about 0.33ms after the time that is plotted... [made code simpler.]');
%%
if onlygoodexpt==1
    % ===== filter outstruct
    if removeBadLFP==0
        expttype = 'xcov_spikes';
    elseif removeBadLFP==1
        expttype = 'lfp_good_and_spikes';
    end
    [OUTSTRUCT] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT, SwitchStruct, expttype);
    [OUTSTRUCT_XCOV] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT_XCOV, SwitchStruct, expttype);
end


%% ============ go thru each switch [PLOT]
if plotraw==1
    [indsgrp, indsgrpU] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, OUTSTRUCT_XCOV.enum, OUTSTRUCT_XCOV.switch, ...
        OUTSTRUCT_XCOV.motifnum});
    
    % ==== go thru each motif from each switch...
    figcount=1;
    subplotrows=6;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    
    % ================ COLLECT (AND PLOT)
    for i=1:length(indsgrpU)
        disp(i);
        
        indlist= find(indsgrp == indsgrpU(i));
        
        bnum = unique(OUTSTRUCT_XCOV.bnum(indlist));
        enum = unique(OUTSTRUCT_XCOV.enum(indlist));
        ss = unique(OUTSTRUCT_XCOV.switch(indlist));
        mm = unique(OUTSTRUCT_XCOV.motifnum(indlist));
        neurlist = unique(OUTSTRUCT_XCOV.neurpair(indlist, :));
        inds_base = OUTSTRUCT_XCOV.inds_base_epoch{indlist(1)};
        inds_WN = OUTSTRUCT_XCOV.inds_WN_epoch{indlist(1)};
        
        bregionlist = {SummaryStruct.birds(bnum).neurons(neurlist).NOTE_Location};
        chanlist = [SummaryStruct.birds(bnum).neurons(neurlist).channel];
        
        bname = SummaryStruct.birds(bnum).birdname;
        ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
        motifname = unique(OUTSTRUCT_XCOV.motifname(indlist));
        
        motifID = lt_neural_QUICK_MotifID(bname, motifname);
        % ==== which trials to plot?
        if strcmp(epochtoplot, 'base_epoch')
            indstoplot = inds_base;
        elseif strcmp(epochtoplot, 'base_all')
            indstoplot = OUTSTRUCT_XCOV.inds_base_allgood{indlist(1)};
        else
            disp('WHICH TRIALS TO PLOT?')
            pause;
        end
        
        
        % ################### COLLECT DATASETS
        % ============ SPIKING
        nset = SwitchCohStruct.bird(bnum).exptnum(enum).switchlist(ss).motifnum(mm).neursetused;
        DAT = MOTIFSTATS_pop.birds(bnum).exptnum(enum).DAT.setnum(nset).motif(mm).SegExtr_neurfakeID;
        segglobal = DAT(1).SegmentsExtract;
        
        % ============= LFP
        DATLFP = SwitchCohStruct.bird(bnum).exptnum(enum).switchlist(ss).motifnum(mm);
        assert(size(DATLFP.lfpall,2)==length(DATLFP.lfpall_chans), 'just to make sure I rememebered format correctly');
        assert(size(DATLFP.lfpall,1)==length(DATLFP.tvals), 'just to make sure I rememebered format correctly');
        assert(length(segglobal) == size(DATLFP.lfpall,1), 'making sure that lfp and seg match up');
        assert(all([segglobal.song_datenum] == DATLFP.tvals), 'making sure that lfp and seg match up');
        
        % --------- GET UNIQUE CHANELS
        [~, indtmp] = unique(chanlist);
        chanlist_lfp = chanlist(indtmp);
        bregion_lfp = bregionlist(indtmp);
        
        
        
        %% ========= ONE PLOT FOR EACH TRIAL
        trialstoplot = indstoplot(randperm(length(indstoplot), Ntrialsplot));
        
        for tt = 1:length(trialstoplot)
            trialthis = trialstoplot(tt);
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title({[bname '-' ename '-sw' num2str(ss)], [motifname{1} ',trial' num2str(trialthis)]});
            
            
            % ###################### [PLOT SPIKING]
            for j=1:length(DAT)
                spktimes = DAT(j).SegmentsExtract(trialthis).spk_Times;
                
                % -- only keep times within boundaries
                spktimes = spktimes - PARAMS.motif_predur;
                indstmp = spktimes>=xtoplot(1) & spktimes<=xtoplot(2);
                spktimes = spktimes(indstmp);
                
                spkclust = {DAT(j).SegmentsExtract(trialthis).spk_Clust};
                assert(all(cellfun(@(x)(all(x==1)), spkclust)), 'not all in same spk cluster?');
                
                assert(neurlist(j)==DAT(j).neurID_orig);
                if strcmp(bregionlist{j}, 'LMAN')
                    pcol = [0.1 0.5 0.1];
                elseif strcmp(bregionlist{j}, 'RA')
                    pcol = 'r';
                end
                lt_neural_PLOT_rasterline(spktimes, j, pcol);
            end
            
            % ###################### [PLOT LFP]
            dozscore =1;
            [lfpcollect, lfpx] = lt_neural_POPLEARN_SylLocked_Over_sub1(chanlist_lfp, ...
                bregion_lfp, DATLFP, trialthis, xtoplot, fpass, 0, 0, dozscore);
            
            for cc=1:length(chanlist_lfp)
                if strcmp(bregion_lfp{cc}, 'LMAN')
                    pcol = [0.1 0.5 0.1];
                elseif strcmp(bregion_lfp{cc}, 'RA')
                    pcol = 'r';
                end
                % --- scale and shift lfp for plotting purposes
                lfpthis = lfpcollect(cc,:);
                lfpthis = (length(DAT)+cc)+lfpthis./2;
                plot(lfpx, lfpthis, 'Color', pcol);
            end
            
            
            % ####################### FORMAT
            ylabel('MU(lower), LFP(zscore/2)');
            axis tight;
            xlim(xtoplot);
            lt_plot_zeroline_vert;
            
            
            % ###################### [OVERLAY SYL ONSET/OFFSET]
            ons = segglobal(trialthis).motifsylOnsets - PARAMS.motif_predur;
            offs = segglobal(trialthis).motifsylOffsets - PARAMS.motif_predur;
            
            YLIM = ylim;
            YLIM(1) = YLIM(2)-0.5;
            lt_neural_QUICK_PlotSylPatches(ons, offs, YLIM);
            
            
            
        end
        
    end
end



%% ########################### COLLECT DATA [ALL INDIVUDAL TRIALS]

[indsgrp, indsgrpU] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, OUTSTRUCT_XCOV.enum, OUTSTRUCT_XCOV.switch, ...
    OUTSTRUCT_XCOV.motifnum});


% ----- expected inds for lfp
tmax_lfp = floor(-1 + (xtoplot(2)-xtoplot(1))*1500);

% =================== COLLECT THINGS ACROSS TRIALS
% === SPIKE
DATSTRUCT_SPK.bnum = [];
DATSTRUCT_SPK.enum = [];
DATSTRUCT_SPK.switch = [];
DATSTRUCT_SPK.motifID = [];
DATSTRUCT_SPK.motifnum = [];

DATSTRUCT_SPK.spike_dat = {};
DATSTRUCT_SPK.spike_neurID = [];
DATSTRUCT_SPK.spike_bregions = {};
DATSTRUCT_SPK.spike_chan = [];

DATSTRUCT_SPK.frmat = {};

% ====== LFP
DATSTRUCT_LFP.bnum = [];
DATSTRUCT_LFP.enum = [];
DATSTRUCT_LFP.switch = [];
DATSTRUCT_LFP.motifID = [];
DATSTRUCT_LFP.motifnum = [];

DATSTRUCT_LFP.LFP_dat = {};
DATSTRUCT_LFP.LFP_chan = [];
DATSTRUCT_LFP.LFP_bregions = {};

% === GENERAL THINGS
DATSTRUCT_SPK.sylonsets = {};
DATSTRUCT_SPK.syloffsets = {};
DATSTRUCT_SPK.inds_base_epoch = {};
DATSTRUCT_SPK.inds_base = {};
DATSTRUCT_SPK.inds_WN_epoch = {};




% ================ [COLLECT EACH TRIAL]
for i=1:length(indsgrpU)
    disp([num2str(i) '/' num2str(max(indsgrpU))]);
    indlist= find(indsgrp == indsgrpU(i));
    
    bnum = unique(OUTSTRUCT_XCOV.bnum(indlist));
    enum = unique(OUTSTRUCT_XCOV.enum(indlist));
    ss = unique(OUTSTRUCT_XCOV.switch(indlist));
    mm = unique(OUTSTRUCT_XCOV.motifnum(indlist));
    neurlist = unique(OUTSTRUCT_XCOV.neurpair(indlist, :));
    inds_base_epoch = single(OUTSTRUCT_XCOV.inds_base_epoch{indlist(1)});
    inds_base = single(OUTSTRUCT_XCOV.inds_base_allgood{indlist(1)});
    inds_WN_epoch = single(OUTSTRUCT_XCOV.inds_WN_epoch{indlist(1)});
    
    bregionlist = {SummaryStruct.birds(bnum).neurons(neurlist).NOTE_Location};
    chanlist = [SummaryStruct.birds(bnum).neurons(neurlist).channel];
    
    bname = SummaryStruct.birds(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    motifname = unique(OUTSTRUCT_XCOV.motifname(indlist));
    motifID = lt_neural_QUICK_MotifID(bname, motifname);
    
    
    % ################### COLLECT DATASETS
    % ============ SPIKING
    nset = SwitchCohStruct.bird(bnum).exptnum(enum).switchlist(ss).motifnum(mm).neursetused;
    DAT = MOTIFSTATS_pop.birds(bnum).exptnum(enum).DAT.setnum(nset).motif(mm).SegExtr_neurfakeID;
    segglobal = DAT(1).SegmentsExtract;
    
    % ============= LFP
    DATLFP = SwitchCohStruct.bird(bnum).exptnum(enum).switchlist(ss).motifnum(mm);
    assert(size(DATLFP.lfpall,2)==length(DATLFP.lfpall_chans), 'just to make sure I rememebered format correctly');
    assert(size(DATLFP.lfpall,1)==length(DATLFP.tvals), 'just to make sure I rememebered format correctly');
    assert(length(segglobal) == size(DATLFP.lfpall,1), 'making sure that lfp and seg match up');
    assert(all([segglobal.song_datenum] == DATLFP.tvals), 'making sure that lfp and seg match up');
    % --------- GET UNIQUE CHANELS
    [~, indtmp] = unique(chanlist);
    chanlist_lfp = chanlist(indtmp);
    bregion_lfp = bregionlist(indtmp);
    
    
    
    % ############################# REMOVE IF BAD TRIALS
    tvals = SwitchCohStruct.bird(bnum).exptnum(enum).switchlist(ss).motifnum(mm).tvals;
    assert(length(tvals)==length(DAT(1).SegmentsExtract));
    
    if removeBadLFP==1
        % --- check trials
        badtrials = find(lt_neural_QUICK_RemoveTrials(bname, ename, ss, tvals, 'lfp'));
        
        inds_base_epoch(ismember(inds_base_epoch, badtrials)) = [];
        inds_base(ismember(inds_base, badtrials)) = [];
        inds_WN_epoch(ismember(inds_WN_epoch, badtrials)) = [];
    end
    
    
    
    %% ========= COLLECT ALL TRIALS
    % ###################### [PLOT SPIKING]
    for j=1:length(neurlist)
        assert(neurlist(j)==DAT(j).neurID_orig);
        
        % ============= COLLECT SPIKE TIMES AND GET IN WINDOW OF INTEREST
        spktimes = {DAT(j).SegmentsExtract.spk_Times};
        
        spktimes = cellfun(@(x)(x-PARAMS.motif_predur), spktimes, 'UniformOutput', 0); % get rel to syl onset
        
        functmp = @(x)(x(x>=xtoplot(1) & x<=xtoplot(2))); % get in window of interest
        spktimes = cellfun(functmp, spktimes, 'UniformOutput', 0);
        spktimes = cellfun(@(x)single(x), spktimes, 'UniformOutput', 0);
        
        % ================ COLLECT SMOOTHED FR
        DAT(j).SegmentsExtract = ...
            lt_neural_SmoothFR(DAT(j).SegmentsExtract, [], kernelSD, [], 0, segglobal);
        frmat = [DAT(j).SegmentsExtract.FRsmooth_rate_CommonTrialDur];
        frx = DAT(j).SegmentsExtract(1).FRsmooth_xbin_CommonTrialDur - PARAMS.motif_predur;
        
        frmat = frmat(frx>xtoplot(1) & frx<xtoplot(2), :);
        frx = frx(frx>xtoplot(1) & frx<xtoplot(2));
        
        
        % ========= [OVERLAY SYL ONSET/OFFSET]
        ons = cellfun(@(x)(x-PARAMS.motif_predur), {segglobal.motifsylOnsets}, 'UniformOutput', 0);
        ons = cellfun(@(x)single(x), ons, 'UniformOutput', 0);
        offs = cellfun(@(x)(x-PARAMS.motif_predur), {segglobal.motifsylOffsets}, 'UniformOutput', 0);
        offs = cellfun(@(x)single(x), offs, 'UniformOutput', 0);
        %         ons = segglobal.motifsylOnsets - PARAMS.motif_predur;
        %         offs = segglobal.motifsylOffsets - PARAMS.motif_predur;
        
        % ============= OUTPUT
        DATSTRUCT_SPK.bnum = [DATSTRUCT_SPK.bnum; bnum];
        DATSTRUCT_SPK.enum = [DATSTRUCT_SPK.enum; enum];
        DATSTRUCT_SPK.switch = [DATSTRUCT_SPK.switch; ss];
        DATSTRUCT_SPK.motifID = [DATSTRUCT_SPK.motifID; motifID];
        DATSTRUCT_SPK.motifnum = [DATSTRUCT_SPK.motifnum; mm];
        
        DATSTRUCT_SPK.frmat = [DATSTRUCT_SPK.frmat; frmat];
        
        DATSTRUCT_SPK.spike_dat = [DATSTRUCT_SPK.spike_dat; {spktimes}];
        DATSTRUCT_SPK.spike_neurID = [DATSTRUCT_SPK.spike_neurID; neurlist(j)];
        DATSTRUCT_SPK.spike_bregions = [DATSTRUCT_SPK.spike_bregions; bregionlist{j}];
        DATSTRUCT_SPK.spike_chan = [DATSTRUCT_SPK.spike_chan; chanlist(j)];
        
        DATSTRUCT_SPK.sylonsets = [DATSTRUCT_SPK.sylonsets; {ons}];
        DATSTRUCT_SPK.syloffsets = [DATSTRUCT_SPK.syloffsets; {offs}];
        DATSTRUCT_SPK.inds_base_epoch = [DATSTRUCT_SPK.inds_base_epoch; inds_base_epoch];
        DATSTRUCT_SPK.inds_base = [DATSTRUCT_SPK.inds_base; inds_base];
        DATSTRUCT_SPK.inds_WN_epoch = [DATSTRUCT_SPK.inds_WN_epoch; inds_WN_epoch];
    end
    
    % ###################### [LFP]
    % ==== COLLECT LFP DATA
    lfpx = DATLFP.t_lfp;
    indx = find(lfpx>=xtoplot(1) & lfpx<=xtoplot(2));
    indx = indx(1:tmax_lfp);
    lfpx = lfpx(indx);
    PARAMS.THIS.lfpx = lfpx;
    for j=1:length(chanlist_lfp)
        chanthis = chanlist_lfp(j);
        bregionthis = bregion_lfp{j};
        
        datlfp = DATLFP.lfpall(:, DATLFP.lfpall_chans==chanthis);
        datlfp = squeeze(lt_neural_Coher_Cell2Mat(datlfp));
        
        % ------- filter
        datlfp = lt_neural_filter(datlfp, 1500, 0, fpass(1), fpass(2));
        
        % -- get window of interest
        datlfp = datlfp(indx,:);
        
        % --- get zscore if desired
        if zscore_lfp==1
            lfpmean = mean(datlfp(:));
            lfpstd = std(datlfp(:));
            datlfp = (datlfp - lfpmean)./lfpstd;
        end
        datlfp = single(datlfp);
        
        % ==================== SAVE OUTPUT
        DATSTRUCT_LFP.bnum = [DATSTRUCT_LFP.bnum; bnum];
        DATSTRUCT_LFP.enum = [DATSTRUCT_LFP.enum; enum];
        DATSTRUCT_LFP.switch = [DATSTRUCT_LFP.switch; ss];
        DATSTRUCT_LFP.motifID = [DATSTRUCT_LFP.motifID; motifID];
        DATSTRUCT_LFP.motifnum = [DATSTRUCT_LFP.motifnum; mm];
        
        DATSTRUCT_LFP.LFP_dat = [DATSTRUCT_LFP.LFP_dat; datlfp];
        DATSTRUCT_LFP.LFP_chan = [DATSTRUCT_LFP.LFP_chan; chanthis];
        DATSTRUCT_LFP.LFP_bregions = [DATSTRUCT_LFP.LFP_bregions; bregionthis] ;
    end
    
    
end


PARAMS.THIS.frmat_x = frx;


%% ============ [PREPROCESS]
% for each spike datapoint, find and save the ind for its channel in LFP
tmp = [];
for i=1:length(DATSTRUCT_SPK.bnum)
    
    bnum = DATSTRUCT_SPK.bnum(i);
    enum = DATSTRUCT_SPK.enum(i);
    ss = DATSTRUCT_SPK.switch(i);
    mm = DATSTRUCT_SPK.motifnum(i);
    chan = DATSTRUCT_SPK.spike_chan(i);
    
    indthis = find(DATSTRUCT_LFP.bnum==bnum & DATSTRUCT_LFP.enum==enum & ...
        DATSTRUCT_LFP.switch==ss & DATSTRUCT_LFP.motifnum==mm ...
        & DATSTRUCT_LFP.LFP_chan==chan);
    assert(length(indthis)==1); % each spike should only have one chan.
    
    tmp = [tmp; indthis];
end
DATSTRUCT_SPK.LFP_indthis = tmp;

%% ============ [PROCESS] SPIKE TRIGGERED LFP

t_lfp = PARAMS.THIS.lfpx;
STA_all = cell(length(DATSTRUCT_SPK.bnum),1);
STA_all_xregions = cell(length(DATSTRUCT_SPK.bnum),1);
for i=1:length(DATSTRUCT_SPK.bnum)
    disp(i);
    
    bnum = DATSTRUCT_SPK.bnum(i);
    enum = DATSTRUCT_SPK.enum(i);
    ss = DATSTRUCT_SPK.switch(i);
    mm = DATSTRUCT_SPK.motifnum(i);
    chan = DATSTRUCT_SPK.spike_chan(i);
    
    % ################################################ [WITHIN REGION]
    spkdat = DATSTRUCT_SPK.spike_dat{i};
    lfpdat = DATSTRUCT_LFP.LFP_dat{DATSTRUCT_SPK.LFP_indthis(i)};
    assert(length(spkdat)==size(lfpdat,2));
    
    trialstoget = DATSTRUCT_SPK.inds_base_epoch{i};
    
    % --- go thru all trials and get STA
    sta_all = lt_neural_POPLEARN_STA(spkdat, lfpdat, trialstoget, sta_wind, ...
        t_lfp);
    
    STA_all{i} = sta_all;
    
    
    % ################################################ [ACROSS REGION]
    bregion_spk = DATSTRUCT_SPK.spike_bregions{i};
    spkdat = DATSTRUCT_SPK.spike_dat{i};
    % --- find all LFP datasets using oppsoite brain region
    if strcmp(bregion_spk, 'LMAN')
        bregiontofind = 'RA';
    elseif strcmp(bregion_spk, 'RA')
        bregiontofind = 'LMAN';
    else
        disp('PROBLEM!!! NOT SURE WHAT BRAIN REGION TO GET...');
        bregiontofind = '';
        pause;
    end
    
    indsLFP = find(DATSTRUCT_LFP.bnum==bnum & DATSTRUCT_LFP.enum==enum & ...
        DATSTRUCT_LFP.switch==ss & DATSTRUCT_LFP.motifnum==mm ...
        & strcmp(DATSTRUCT_LFP.LFP_bregions, bregiontofind));
    assert(length(indsLFP)>0);
    
    sta_all_xregion = cell(length(indsLFP), 1);
    for j=1:length(indsLFP)
        lfpdat = DATSTRUCT_LFP.LFP_dat{indsLFP(j)};
        assert(length(spkdat)==size(lfpdat,2));
        
        % --- go thru all trials and get STA
        sta_all = lt_neural_POPLEARN_STA(spkdat, lfpdat, trialstoget, sta_wind, ...
            t_lfp);
        
        sta_all_xregion{j} = sta_all;
    end
    sta_all_xregion = lt_neural_Coher_Cell2Mat(sta_all_xregion);
    
    STA_all_xregions{i} = sta_all_xregion;
end

DATSTRUCT_SPK.STA_all = STA_all;
DATSTRUCT_SPK.STA_all_xregions = STA_all_xregions;


%% =========== time bins (lags) for sta
PARAMS.THIS.sta_lags = linspace(sta_wind(1), sta_wind(end), size(sta_all,2));


if all(isempty(STA_all))
    disp('NOT PLOTTING STA, not enough data in xtoplot to get STA...');
else
    %% ================= [PLOT] STA
    
    figcount=1;
    subplotrows=4;
    subplotcols=3;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    
    % ===========
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    xlabel('lag(rel spike)');
    ylabel('motif x cell (random order)');
    title('LMAN');
    indsthis = strcmp(DATSTRUCT_SPK.spike_bregions, 'LMAN');
    y = DATSTRUCT_SPK.STA_all(indsthis);
    y = cellfun(@(x)mean(x,1), y, 'UniformOutput', 0);
    y = cell2mat(y)';
    % randomize y axis order
    y = y(:, randperm(size(y,2)));
    f = 1:size(y,1);
    clim = [-0.75*max(abs(y(:))) 0.75*max(abs(y(:)))];
    lt_neural_Coher_Plot(y, PARAMS.THIS.sta_lags, f, 1, '', clim);
    
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    xlabel('lag(rel spike)');
    ylabel('motif x cell (sort by max)');
    title('LMAN');
    indsthis = strcmp(DATSTRUCT_SPK.spike_bregions, 'LMAN');
    y = DATSTRUCT_SPK.STA_all(indsthis);
    y = cellfun(@(x)mean(x,1), y, 'UniformOutput', 0);
    y = cell2mat(y)';
    [~, indtmp] = max(y);
    [~, indtmp] = sort(indtmp);
    y = y(:, indtmp);
    f = 1:size(y,1);
    clim = [-0.75*max(abs(y(:))) 0.75*max(abs(y(:)))];
    lt_neural_Coher_Plot(y, PARAMS.THIS.sta_lags, f, 1, '', clim);
    
    
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    xlabel('lag(rel spike)');
    ylabel('motif x cell (sort by min)');
    title('LMAN');
    indsthis = strcmp(DATSTRUCT_SPK.spike_bregions, 'LMAN');
    y = DATSTRUCT_SPK.STA_all(indsthis);
    y = cellfun(@(x)mean(x,1), y, 'UniformOutput', 0);
    y = cell2mat(y)';
    [~, indtmp] = min(y);
    [~, indtmp] = sort(indtmp);
    y = y(:, indtmp);
    f = 1:size(y,1);
    clim = [-0.75*max(abs(y(:))) 0.75*max(abs(y(:)))];
    lt_neural_Coher_Plot(y, PARAMS.THIS.sta_lags, f, 1, '', clim);
    
    
    % ===========
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    xlabel('lag(rel spike)');
    ylabel('lfp(V)');
    title('LMAN');
    indsthis = strcmp(DATSTRUCT_SPK.spike_bregions, 'LMAN');
    y = DATSTRUCT_SPK.STA_all(indsthis);
    y = cellfun(@(x)mean(x,1), y, 'UniformOutput', 0);
    y = cell2mat(y);
    x = PARAMS.THIS.sta_lags;
    % plot(x, y', 'Color', [0.7 0.7 0.7]);
    plot(x, y');
    
    
    
    % ===========
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    xlabel('lag(rel spike)');
    ylabel('lfp(V)');
    title('LMAN');
    indsthis = strcmp(DATSTRUCT_SPK.spike_bregions, 'LMAN');
    y = DATSTRUCT_SPK.STA_all(indsthis);
    y = cellfun(@(x)mean(x,1), y, 'UniformOutput', 0);
    y = cell2mat(y);
    x = PARAMS.THIS.sta_lags;
    ymean = mean(y,1);
    ysem = lt_sem(y);
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
    
    % ===========
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    xlabel('lag(rel spike)');
    ylabel('lfp(V)');
    title('LMAN (each bird)');
    pcols = lt_make_plot_colors(max(DATSTRUCT_SPK.bnum), 0,0);
    for i=1:max(DATSTRUCT_SPK.bnum)
        indsthis = strcmp(DATSTRUCT_SPK.spike_bregions, 'LMAN') & ...
            DATSTRUCT_SPK.bnum==i;
        if ~any(indsthis)
            continue
        end
        y = DATSTRUCT_SPK.STA_all(indsthis);
        y = cellfun(@(x)mean(x,1), y, 'UniformOutput', 0);
        y = cell2mat(y);
        x = PARAMS.THIS.sta_lags;
        ymean = mean(y,1);
        ysem = lt_sem(y);
        shadedErrorBar(x, ymean, ysem, {'Color', pcols{i}}, 1);
    end
    lt_plot_zeroline_vert;
    
    
    % ===========
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    xlabel('lag(rel spike)');
    ylabel('motif x cell (random order)');
    title('RA');
    indsthis = strcmp(DATSTRUCT_SPK.spike_bregions, 'RA');
    y = DATSTRUCT_SPK.STA_all(indsthis);
    y = cellfun(@(x)mean(x,1), y, 'UniformOutput', 0);
    y = cell2mat(y)';
    % randomize y axis order
    y = y(:, randperm(size(y,2)));
    f = 1:size(y,1);
    clim = [-0.75*max(abs(y(:))) 0.75*max(abs(y(:)))];
    lt_neural_Coher_Plot(y, PARAMS.THIS.sta_lags, f, 1, '', clim);
    
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    xlabel('lag(rel spike)');
    ylabel('motif x cell (sort by max)');
    title('RA');
    indsthis = strcmp(DATSTRUCT_SPK.spike_bregions, 'RA');
    y = DATSTRUCT_SPK.STA_all(indsthis);
    y = cellfun(@(x)mean(x,1), y, 'UniformOutput', 0);
    y = cell2mat(y)';
    [~, indtmp] = max(y);
    [~, indtmp] = sort(indtmp);
    y = y(:, indtmp);
    f = 1:size(y,1);
    clim = [-0.75*max(abs(y(:))) 0.75*max(abs(y(:)))];
    lt_neural_Coher_Plot(y, PARAMS.THIS.sta_lags, f, 1, '', clim);
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    xlabel('lag(rel spike)');
    ylabel('motif x cell (sort by min)');
    title('RA');
    indsthis = strcmp(DATSTRUCT_SPK.spike_bregions, 'RA');
    y = DATSTRUCT_SPK.STA_all(indsthis);
    y = cellfun(@(x)mean(x,1), y, 'UniformOutput', 0);
    y = cell2mat(y)';
    [~, indtmp] = min(y);
    [~, indtmp] = sort(indtmp);
    y = y(:, indtmp);
    f = 1:size(y,1);
    clim = [-0.75*max(abs(y(:))) 0.75*max(abs(y(:)))];
    lt_neural_Coher_Plot(y, PARAMS.THIS.sta_lags, f, 1, '', clim);
    
    
    % ===========
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    xlabel('lag(rel spike)');
    ylabel('lfp(V)');
    title('RA');
    indsthis = strcmp(DATSTRUCT_SPK.spike_bregions, 'RA');
    y = DATSTRUCT_SPK.STA_all(indsthis);
    y = cellfun(@(x)mean(x,1), y, 'UniformOutput', 0);
    y = cell2mat(y);
    x = PARAMS.THIS.sta_lags;
    % plot(x, y', 'Color', [0.7 0.7 0.7]);
    plot(x, y');
    
    
    % ===========
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    xlabel('lag(rel spike)');
    ylabel('lfp(V)');
    title('RA');
    indsthis = strcmp(DATSTRUCT_SPK.spike_bregions, 'RA');
    y = DATSTRUCT_SPK.STA_all(indsthis);
    y = cellfun(@(x)mean(x,1), y, 'UniformOutput', 0);
    y = cell2mat(y);
    x = PARAMS.THIS.sta_lags;
    ymean = mean(y,1);
    ysem = lt_sem(y);
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
    
    % ===========
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    xlabel('lag(rel spike)');
    ylabel('lfp(V)');
    title('RA (each bird)');
    pcols = lt_make_plot_colors(max(DATSTRUCT_SPK.bnum), 0,0);
    for i=1:max(DATSTRUCT_SPK.bnum)
        indsthis = strcmp(DATSTRUCT_SPK.spike_bregions, 'RA') & ...
            DATSTRUCT_SPK.bnum==i;
        if ~any(indsthis)
            continue
        end
        y = DATSTRUCT_SPK.STA_all(indsthis);
        y = cellfun(@(x)mean(x,1), y, 'UniformOutput', 0);
        y = cell2mat(y);
        x = PARAMS.THIS.sta_lags;
        ymean = mean(y,1);
        ysem = lt_sem(y);
        shadedErrorBar(x, ymean, ysem, {'Color', pcols{i}}, 1);
    end
    lt_plot_zeroline_vert;
    
    
    %% ================= [PLOT] STA [across regions]
    % -- get one trace for each spike dataset
    
    figcount=1;
    subplotrows=4;
    subplotcols=3;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    
    % ===========
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    xlabel('lag(rel spike)');
    ylabel('motif x cell (random order)');
    title('LMAN(spk)vs RA(LFP)[one dat per spk dataset]');
    indsthis = strcmp(DATSTRUCT_SPK.spike_bregions, 'LMAN');
    y = DATSTRUCT_SPK.STA_all_xregions(indsthis);
    y = cellfun(@(x)mean(x,3), y, 'UniformOutput', 0); % -- get one trace for each spike dataset
    y = cellfun(@(x)mean(x,1), y, 'UniformOutput', 0);
    y = cell2mat(y)';
    % randomize y axis order
    y = y(:, randperm(size(y,2)));
    f = 1:size(y,1);
    clim = [-0.75*max(abs(y(:))) 0.75*max(abs(y(:)))];
    lt_neural_Coher_Plot(y, PARAMS.THIS.sta_lags, f, 1, '', clim);
    
    
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    xlabel('lag(rel spike)');
    ylabel('motif x cell (sort by max)');
    title('LMAN(spk)vs RA(LFP)[one dat per spk dataset]');
    
    indsthis = strcmp(DATSTRUCT_SPK.spike_bregions, 'LMAN');
    y = DATSTRUCT_SPK.STA_all_xregions(indsthis);
    y = cellfun(@(x)mean(x,3), y, 'UniformOutput', 0); % -- get one trace for each spike dataset
    
    y = cellfun(@(x)mean(x,1), y, 'UniformOutput', 0);
    y = cell2mat(y)';
    [~, indtmp] = max(y);
    [~, indtmp] = sort(indtmp);
    y = y(:, indtmp);
    f = 1:size(y,1);
    clim = [-0.75*max(abs(y(:))) 0.75*max(abs(y(:)))];
    lt_neural_Coher_Plot(y, PARAMS.THIS.sta_lags, f, 1, '', clim);
    
    
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    xlabel('lag(rel spike)');
    ylabel('motif x cell (sort by min)');
    title('LMAN(spk)vs RA(LFP)[one dat per spk dataset]');
    
    indsthis = strcmp(DATSTRUCT_SPK.spike_bregions, 'LMAN');
    y = DATSTRUCT_SPK.STA_all_xregions(indsthis);
    y = cellfun(@(x)mean(x,3), y, 'UniformOutput', 0); % -- get one trace for each spike dataset
    y = cellfun(@(x)mean(x,1), y, 'UniformOutput', 0);
    y = cell2mat(y)';
    [~, indtmp] = min(y);
    [~, indtmp] = sort(indtmp);
    y = y(:, indtmp);
    f = 1:size(y,1);
    clim = [-0.75*max(abs(y(:))) 0.75*max(abs(y(:)))];
    lt_neural_Coher_Plot(y, PARAMS.THIS.sta_lags, f, 1, '', clim);
    
    
    % ===========
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    xlabel('lag(rel spike)');
    ylabel('lfp(V)');
    title('LMAN(spk)vs RA(LFP)[one dat per spk dataset]');
    
    indsthis = strcmp(DATSTRUCT_SPK.spike_bregions, 'LMAN');
    y = DATSTRUCT_SPK.STA_all_xregions(indsthis);
    y = cellfun(@(x)mean(x,3), y, 'UniformOutput', 0); % -- get one trace for each spike dataset
    y = cellfun(@(x)mean(x,1), y, 'UniformOutput', 0);
    y = cell2mat(y);
    x = PARAMS.THIS.sta_lags;
    % plot(x, y', 'Color', [0.7 0.7 0.7]);
    plot(x, y');
    
    
    
    % ===========
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    xlabel('lag(rel spike)');
    ylabel('lfp(V)');
    title('LMAN(spk)vs RA(LFP)[one dat per spk dataset]');
    
    indsthis = strcmp(DATSTRUCT_SPK.spike_bregions, 'LMAN');
    y = DATSTRUCT_SPK.STA_all_xregions(indsthis);
    y = cellfun(@(x)mean(x,3), y, 'UniformOutput', 0); % -- get one trace for each spike dataset
    y = cellfun(@(x)mean(x,1), y, 'UniformOutput', 0);
    y = cell2mat(y);
    x = PARAMS.THIS.sta_lags;
    ymean = mean(y,1);
    ysem = lt_sem(y);
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
    lt_plot_zeroline_vert;
    
    % ===========
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    xlabel('lag(rel spike)');
    ylabel('lfp(V)');
    title('LMAN(spk)vs RA(LFP)[one dat per spk dataset]');
    pcols = lt_make_plot_colors(max(DATSTRUCT_SPK.bnum), 0,0);
    for i=1:max(DATSTRUCT_SPK.bnum)
        indsthis = strcmp(DATSTRUCT_SPK.spike_bregions, 'LMAN') & ...
            DATSTRUCT_SPK.bnum==i;
        if ~any(indsthis)
            continue
        end
        y = DATSTRUCT_SPK.STA_all_xregions(indsthis);
        y = cellfun(@(x)mean(x,3), y, 'UniformOutput', 0); % -- get one trace for each spike dataset
        y = cellfun(@(x)mean(x,1), y, 'UniformOutput', 0);
        y = cell2mat(y);
        x = PARAMS.THIS.sta_lags;
        ymean = mean(y,1);
        ysem = lt_sem(y);
        shadedErrorBar(x, ymean, ysem, {'Color', pcols{i}}, 1);
    end
    lt_plot_zeroline_vert;
    
    
    
    
    
    
    
    % ===========
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    xlabel('lag(rel spike)');
    ylabel('motif x cell (random order)');
    title('RA(spk)vs LMAN(LFP)[one dat per spk dataset]');
    indsthis = strcmp(DATSTRUCT_SPK.spike_bregions, 'RA');
    y = DATSTRUCT_SPK.STA_all_xregions(indsthis);
    y = cellfun(@(x)mean(x,3), y, 'UniformOutput', 0); % -- get one trace for each spike dataset
    y = cellfun(@(x)mean(x,1), y, 'UniformOutput', 0);
    y = cell2mat(y)';
    % randomize y axis order
    y = y(:, randperm(size(y,2)));
    f = 1:size(y,1);
    clim = [-0.75*max(abs(y(:))) 0.75*max(abs(y(:)))];
    lt_neural_Coher_Plot(y, PARAMS.THIS.sta_lags, f, 1, '', clim);
    
    
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    xlabel('lag(rel spike)');
    ylabel('motif x cell (sort by max)');
    title('RA(spk)vs LMAN(LFP)[one dat per spk dataset]');
    
    indsthis = strcmp(DATSTRUCT_SPK.spike_bregions, 'RA');
    y = DATSTRUCT_SPK.STA_all_xregions(indsthis);
    y = cellfun(@(x)mean(x,3), y, 'UniformOutput', 0); % -- get one trace for each spike dataset
    
    y = cellfun(@(x)mean(x,1), y, 'UniformOutput', 0);
    y = cell2mat(y)';
    [~, indtmp] = max(y);
    [~, indtmp] = sort(indtmp);
    y = y(:, indtmp);
    f = 1:size(y,1);
    clim = [-0.75*max(abs(y(:))) 0.75*max(abs(y(:)))];
    lt_neural_Coher_Plot(y, PARAMS.THIS.sta_lags, f, 1, '', clim);
    
    
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    xlabel('lag(rel spike)');
    ylabel('motif x cell (sort by min)');
    title('RA(spk)vs LMAN(LFP)[one dat per spk dataset]');
    
    indsthis = strcmp(DATSTRUCT_SPK.spike_bregions, 'RA');
    y = DATSTRUCT_SPK.STA_all_xregions(indsthis);
    y = cellfun(@(x)mean(x,3), y, 'UniformOutput', 0); % -- get one trace for each spike dataset
    y = cellfun(@(x)mean(x,1), y, 'UniformOutput', 0);
    y = cell2mat(y)';
    [~, indtmp] = min(y);
    [~, indtmp] = sort(indtmp);
    y = y(:, indtmp);
    f = 1:size(y,1);
    clim = [-0.75*max(abs(y(:))) 0.75*max(abs(y(:)))];
    lt_neural_Coher_Plot(y, PARAMS.THIS.sta_lags, f, 1, '', clim);
    
    
    % ===========
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    xlabel('lag(rel spike)');
    ylabel('lfp(V)');
    title('RA(spk)vs LMAN(LFP)[one dat per spk dataset]');
    
    indsthis = strcmp(DATSTRUCT_SPK.spike_bregions, 'RA');
    y = DATSTRUCT_SPK.STA_all_xregions(indsthis);
    y = cellfun(@(x)mean(x,3), y, 'UniformOutput', 0); % -- get one trace for each spike dataset
    y = cellfun(@(x)mean(x,1), y, 'UniformOutput', 0);
    y = cell2mat(y);
    x = PARAMS.THIS.sta_lags;
    % plot(x, y', 'Color', [0.7 0.7 0.7]);
    plot(x, y');
    
    
    
    % ===========
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    xlabel('lag(rel spike)');
    ylabel('lfp(V)');
    title('RA(spk)vs LMAN(LFP)[one dat per spk dataset]');
    
    indsthis = strcmp(DATSTRUCT_SPK.spike_bregions, 'RA');
    y = DATSTRUCT_SPK.STA_all_xregions(indsthis);
    y = cellfun(@(x)mean(x,3), y, 'UniformOutput', 0); % -- get one trace for each spike dataset
    y = cellfun(@(x)mean(x,1), y, 'UniformOutput', 0);
    y = cell2mat(y);
    x = PARAMS.THIS.sta_lags;
    ymean = mean(y,1);
    ysem = lt_sem(y);
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
    lt_plot_zeroline_vert;
    
    % ===========
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    xlabel('lag(rel spike)');
    ylabel('lfp(V)');
    title('RA(spk)vs LMAN(LFP)[one dat per spk dataset]');
    pcols = lt_make_plot_colors(max(DATSTRUCT_SPK.bnum), 0,0);
    for i=1:max(DATSTRUCT_SPK.bnum)
        indsthis = strcmp(DATSTRUCT_SPK.spike_bregions, 'RA') & ...
            DATSTRUCT_SPK.bnum==i;
        if ~any(indsthis)
            continue
        end
        y = DATSTRUCT_SPK.STA_all_xregions(indsthis);
        y = cellfun(@(x)mean(x,3), y, 'UniformOutput', 0); % -- get one trace for each spike dataset
        y = cellfun(@(x)mean(x,1), y, 'UniformOutput', 0);
        y = cell2mat(y);
        x = PARAMS.THIS.sta_lags;
        ymean = mean(y,1);
        ysem = lt_sem(y);
        shadedErrorBar(x, ymean, ysem, {'Color', pcols{i}}, 1);
    end
    lt_plot_zeroline_vert;
end