function OUTSTRUCT_XCOV = lt_neural_POPLEARN_Xcov_ExtrLearn(OUTSTRUCT, OUTSTRUCT_XCOV, ... 
    SwitchStruct, PARAMS, SwitchCohStruct, MOTIFSTATS_pop, windowprem)
%%  lt 3/4/19 - extracts learning into time bins... EPOCHS
% motifpredur = PARAMS.motif_predur;

%% ======= first downsample OUTSRRUCT to get only the experiments in xcov
% [OUTSTRUCT, OUTSTRUCT_XCOV, indsXcov] = lt_neural_POPLEARN_MatchOutstructs(OUTSTRUCT, OUTSTRUCT_XCOV, ...
%     SwitchStruct, 0);

% %% ====== for each case compute z-scored laerning
%
% LearnTargDir_Z = cell(length(OUTSTRUCT.bnum),1);
% % LearnTargDir_Slope = nan(size(OUTSTRUCT.bnum));
% % TrainDuration = nan(size(OUTSTRUCT.bnum)); % time from mid of base bin to mid of training
% LearnTargDir_Z_XCOV = cell(length(OUTSTRUCT_XCOV.bnum),1);
% for i=1:length(OUTSTRUCT.bnum)
%
%     ffvals = OUTSTRUCT.ffvals{i};
%     tvals = OUTSTRUCT.tvals{i};
%
%     if all(isnan(ffvals))
%         continue
%     end
%
%     % ---- USE THE INDS FROM XCOV
%     %     indsbase = OUTSTRUCT_XCOV.inds_base_epoch{indsXcov{i}(1)};
%     indsbase = OUTSTRUCT_XCOV.inds_base_epoch{indsXcov{i}(1)};
%     indsWN = OUTSTRUCT_XCOV.inds_WN_allgood{indsXcov{i}(1)};
%
%     % --- sanity check....
%     assert(all(OUTSTRUCT.indsXcov_all{i} == indsXcov{i}));
%
%     % =========== GET EDGES OF TRIAL BINS
%     N = size(OUTSTRUCT_XCOV.XcovgramWN_epochs{1},3); % num epochs to get
%     % --- divide up WN into different bins
%     binsize = length(indsWN)/N;
%
%     trialedges = round(linspace(1, length(indsWN), N+1));
%     trialedges = indsWN(trialedges);
%     trialedges(end)=trialedges(end)+1;
%
%
%     % =========== GET ZSCORE LEARNING IN EACH BIN
%     basemean = mean(ffvals(indsbase));
%     basestd = std(ffvals(indsbase));
%
%     ffzscore = [];
%     for j=1:N
%         trialsthis = trialedges(j):(trialedges(j+1)-1);
%
%         ffz = (mean(ffvals(trialsthis))-basemean)./basestd;
%         ffzscore = [ffzscore; ffz];
%     end
%
%
%     % ========= flip if training is down
%     traindir = OUTSTRUCT_XCOV.learndirTarg(indsXcov{i}(1));
%     ffzscore = ffzscore*traindir;
%
%
%     % ============ SAVE TO OUTPUT
%     LearnTargDir_Z{i} = ffzscore;
%
%     for j=indsXcov{i}'
%         if ~isempty(LearnTargDir_Z_XCOV{j})
%             assert(all(LearnTargDir_Z_XCOV{j}==ffzscore));
%         end
%         LearnTargDir_Z_XCOV{j} = ffzscore;
%     end
%
%     %     %  ################### TRAINING DURATION
%     %     traindur = (mean(tvals(indsWN)) - mean(tvals(indsbase)))*24;
%     %     TrainDuration(i) = traindur;
%     %
%     %     % ############################### LEARNING RATE (SLOPE OF FF, DURING
%     %     % TRAINING)
%     %     indsWN_all = OUTSTRUCT.indsWN{i};
%     %
%     %     ffthis = ffvals(indsWN_all);
%     %     tthis = tvals(indsWN_all);
%     %     tthis = (tthis-tthis(1))*24; % conver to hours
%     %
%     %     [~,~,~,~,~, stats] = lt_regress(ffthis, tthis, 0);
%     %     learnslope = stats.slope*traindir; % flip slope so up is in direction of training.
%     %
%     %     LearnTargDir_Slope(i) = learnslope;
% end
%
% % OUTSTRUCT.LearnTargDir_Z = LearnTargDir_Z;
% % OUTSTRUCT.TrainDuration = TrainDuration;
% % OUTSTRUCT.LearnTargDir_Slope = LearnTargDir_Slope;
%
% OUTSTRUCT_XCOV.LearnTargDir_Z_Epochs = LearnTargDir_Z_XCOV;


%% ==========
LearnTargDir_Z_XCOV = cell(length(OUTSTRUCT_XCOV.bnum), 1);
LearnTargDir_Z_XCOV_Split_Adaptive = cell(length(OUTSTRUCT_XCOV.bnum), 1);
LearnTargDir_Z_XCOV_Split_NonAdap = cell(length(OUTSTRUCT_XCOV.bnum), 1);
TimeBins_Base_Epochs = cell(length(OUTSTRUCT_XCOV.bnum), 1);
Nall = nan(length(OUTSTRUCT_XCOV.bnum),1);

epochSplitStatsAll = cell(length(OUTSTRUCT_XCOV.bnum), 1);

bregionpairAll = cell(length(OUTSTRUCT_XCOV.bnum), 1);

for i=1:length(OUTSTRUCT_XCOV.bnum)
    
    bnum = OUTSTRUCT_XCOV.bnum(i);
    enum = OUTSTRUCT_XCOV.enum(i);
    sw = OUTSTRUCT_XCOV.switch(i);
    mm = OUTSTRUCT_XCOV.motifnum(i);
    trialedges = OUTSTRUCT_XCOV.trialedges_epoch{i};
    neurpair = OUTSTRUCT_XCOV.neurpair(i,:);
    
    inds_lfp = find(OUTSTRUCT.bnum==bnum & OUTSTRUCT.enum==enum & ...
        OUTSTRUCT.switch==sw & OUTSTRUCT.motifnum==mm);
    
    bregionpair = OUTSTRUCT.bregionpair{inds_lfp(1)};
    ffvals = OUTSTRUCT.ffvals{inds_lfp(1)};
    tvals = OUTSTRUCT.tvals{inds_lfp(1)};
    tvals_hours = (tvals-tvals(1))*24;
    
    
    % =============== CONFIRM TVALS AND FFVALS ARE CORRECT
    assert(all(tvals== ...
        SwitchCohStruct.bird(bnum).exptnum(enum).switchlist(sw).motifnum(mm).tvals), 'why no match?');
    if ~any(isnan(ffvals))
    assert(all(ffvals== ...
        SwitchCohStruct.bird(bnum).exptnum(enum).switchlist(sw).motifnum(mm).ffvals), 'why no match?');
    end
    
    % ---- USE THE INDS FROM XCOV
    %     indsbase = OUTSTRUCT_XCOV.inds_base_epoch{indsXcov{i}(1)};
    indsbase = OUTSTRUCT_XCOV.inds_base_epoch{i};
    indsWN = OUTSTRUCT_XCOV.inds_WN_allgood{i};
    
    % --- sanity check....
    %     assert(all(OUTSTRUCT.indsXcov_all{i} == indsXcov{i}));
    
    % =========== GET EDGES OF TRIAL BINS
    N = size(OUTSTRUCT_XCOV.XcovgramWN_epochs{1},3); % num epochs to get
    % --- divide up WN into different bins
    binsize = length(indsWN)/N;
    
    %     trialedges = round(linspace(1, length(indsWN), N+1));
    %     Nthis = median(diff(trialedges));
    %     trialedges = indsWN(trialedges);
    %     trialedges(end)=trialedges(end)+1;
    
    Nthis = median(diff(round(linspace(1, length(indsWN), N+1))));
    
    
    % ===== samle size
    Nall(i) = Nthis;
    
    % =========== GET ZSCORE LEARNING IN EACH BIN
    basemean = mean(ffvals(indsbase));
    basestd = std(ffvals(indsbase));
    timebase = median(tvals_hours(indsbase));
    
    ffzscore = [];
    ffhi_z_all = [];
    fflo_z_all = [];
    timemedian = [];
    for j=1:N
        trialsthis = trialedges(j):(trialedges(j+1)-1);
        if any(isnan(trialsthis))
            
            ffzscore = [ffzscore; nan];
            ffhi_z_all = [ffhi_z_all; nan];
            fflo_z_all = [fflo_z_all; nan];
            timemedian = [timemedian; nan];
            continue
        end
        
        
        ffz = (mean(ffvals(trialsthis))-basemean)./basestd;
        ffzscore = [ffzscore; ffz];
        
        % ======= get mean time within bins
        timemedian = [timemedian; nanmedian(tvals_hours(trialsthis))];
        
        % ======= get high and low pitch within bin
        ffmid = median(ffvals(trialsthis));
        ffhi = ffvals(trialsthis(ffvals(trialsthis)>ffmid));
        fflo = ffvals(trialsthis(ffvals(trialsthis)<ffmid));
        
        ffhi_z = (mean(ffhi)-basemean)./basestd;
        fflo_z = (mean(fflo)-basemean)./basestd;
        
        ffhi_z_all = [ffhi_z_all; ffhi_z];
        fflo_z_all = [fflo_z_all; fflo_z];
    end
    
    
    % ========= flip if training is down
    traindir = OUTSTRUCT_XCOV.learndirTarg(i);
    ffzscore = ffzscore*traindir;
    ffhi_z_all =ffhi_z_all*traindir;
    fflo_z_all =fflo_z_all*traindir;
    
    % -- flip hi and lo z if traindir neg
    if traindir==-1
        tmp = ffhi_z_all;
        tmp2 = fflo_z_all;
        
        ffhi_z_all = tmp2;
        fflo_z_all = tmp;
    end
    
    % ============ SAVE TO OUTPUT
    %     LearnTargDir_Z{i} = ffzscore;
    LearnTargDir_Z_XCOV{i} = ffzscore;
    
    LearnTargDir_Z_XCOV_Split_Adaptive{i} = ffhi_z_all; % i.e. ff (z) in half of trials in which ff is furthest from baseline (i.e. in adaptive direction fo trainig_)
    LearnTargDir_Z_XCOV_Split_NonAdap{i} = fflo_z_all;
    
    TimeBins_Base_Epochs{i} = [timebase timemedian'];
    
    bregionpairAll{i} = bregionpair;

    %% =================== ASK HOW FF HIGH AND LOW ARE DISTRIBUTED
    % 1) across epoch
    % 2) across song bout
    % 3) mean frate difference
    
    % ================ COMBINE BASELINE INTO TRIAL EDGES
    if isnan(trialedges(1))
        % tack on baseline. have to replace first in
        assert(isnan(trialedges(2)), 'otherwise will get a first traiing epoch, even though should not (since flanking inds are both nan)');
        trialedges = [min(indsbase) trialedges];
        trialedges(2) = max(indsbase)+1;
    else
    assert(max(indsbase)+1 == trialedges(1), 'this lets me combine base with trial edges easily. if not true then come upw ith different metod to include baseline');
    trialedges = [min(indsbase) trialedges] ;
    end
    
    % ============= EXTRACT FRATE OVER ALL TRIALS
    neurset = SwitchCohStruct.bird(bnum).exptnum(enum).switchlist(sw).motifnum(mm).neursetused;
    segdat = MOTIFSTATS_pop.birds(bnum).exptnum(enum).DAT.setnum(neurset).motif(mm).SegExtr_neurfakeID;
%     SwitchCohStruct.bird(bnum).exptnum(enum).switchlist(sw).motifnum(mm).bregionpair_originalorder
    
    % ========== STATS - time/rendition in song
    posinbout = [segdat(1).SegmentsExtract.BOUT_RendInBout]; % bout defined as sounds separated by 1s. rend in bout is 1, 2, ...
    
    
    % ============= extract mean spike rate on each trial for each neuron
    % in this pair
    NspksByNeuron = cell(1,2);
    for j=1:length(neurpair)
        nn = neurpair(j);
        indtmp = find([segdat.neurID_orig]==nn);
        assert(length(indtmp)==1, 'need to find original neruon id');
        
        segthis = segdat(indtmp).SegmentsExtract;
        
        % ===== count spikes in some premotor window
        Nspks = lt_neural_QUICK_SpkCnts(segthis, PARAMS.motif_predur, windowprem);        
        NspksByNeuron{j} = single(Nspks);
    end
    
    epochSplitStats = struct;
    for j=1:length(trialedges)-1
        trialsthis = trialedges(j):(trialedges(j+1)-1);
        
        if any(isnan(trialsthis))
            % ========= OUTUPT
            epochSplitStats(j).inds_hi = nan;
            epochSplitStats(j).inds_lo = nan;
            epochSplitStats(j).ffthis = nan;
            epochSplitStats(j).songboutID = nan;
            epochSplitStats(j).nspksByNeuron = nan;
            epochSplitStats(j).rendnumInBout_1secIBI= nan;
           continue
        end
        
        % ============= HI AND LOW TRIALS
        ffthis = ffvals(trialsthis);
        ffmed = median(ffthis);
        
        inds_hi = ffthis>ffmed;
        inds_lo = ffthis<ffmed;
        
        % =========== GROUP BY SONG (i.e. file time)
        tthis = tvals(trialsthis);
        songboutID = grp2idx(tthis);
        
        % ========== Nspks
        nspksthis = cellfun(@(x)x(trialsthis), NspksByNeuron, 'UniformOutput', 0);
        
        
        % ==== pos in bout
        pos = posinbout(trialsthis);
        
        % ========= OUTUPT
        epochSplitStats(j).inds_hi = inds_hi;
        epochSplitStats(j).inds_lo = inds_lo;
        epochSplitStats(j).ffthis = ffthis;
        epochSplitStats(j).songboutID = single(songboutID);
        epochSplitStats(j).nspksByNeuron = nspksthis;
        epochSplitStats(j).rendnumInBout_1secIBI= pos;
    end
    
        
    epochSplitStatsAll{i} = epochSplitStats;
end

OUTSTRUCT_XCOV.LearnTargDir_Z_Epochs = LearnTargDir_Z_XCOV;
OUTSTRUCT_XCOV.LearnTargDir_Z_XCOV_Split_Adaptive = LearnTargDir_Z_XCOV_Split_Adaptive;
OUTSTRUCT_XCOV.LearnTargDir_Z_XCOV_Split_NonAdap = LearnTargDir_Z_XCOV_Split_NonAdap;
OUTSTRUCT_XCOV.NperBin_Epochs = Nall;
OUTSTRUCT_XCOV.TimeBins_Base_Epochs = TimeBins_Base_Epochs;
OUTSTRUCT_XCOV.epochSplitStatsAll = epochSplitStatsAll;
OUTSTRUCT_XCOV.bregionpair = bregionpairAll;
