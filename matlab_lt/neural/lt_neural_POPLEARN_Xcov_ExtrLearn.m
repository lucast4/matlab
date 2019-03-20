function OUTSTRUCT_XCOV = lt_neural_POPLEARN_Xcov_ExtrLearn(OUTSTRUCT, OUTSTRUCT_XCOV, SwitchStruct, PARAMS, ...
    onlygoodexpt)
%%  lt 3/4/19 - extracts learning into time bins... EPOCHS

%% ======= first downsample OUTSRRUCT to get only the experiments in xcov
[OUTSTRUCT, OUTSTRUCT_XCOV, indsXcov] = lt_neural_POPLEARN_MatchOutstructs(OUTSTRUCT, OUTSTRUCT_XCOV, ...
    SwitchStruct, 0);

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
for i=1:length(OUTSTRUCT_XCOV.bnum)
    
    bnum = OUTSTRUCT_XCOV.bnum(i);
    enum = OUTSTRUCT_XCOV.enum(i);
    sw = OUTSTRUCT_XCOV.switch(i);
    mm = OUTSTRUCT_XCOV.motifnum(i);
    
    inds_lfp = find(OUTSTRUCT.bnum==bnum & OUTSTRUCT.enum==enum & ...
        OUTSTRUCT.switch==sw & OUTSTRUCT.motifnum==mm);
    
    
    ffvals = OUTSTRUCT.ffvals{inds_lfp(1)};
    tvals = OUTSTRUCT.tvals{inds_lfp(1)};
    tvals_hours = (tvals-tvals(1))*24;
    
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
    
    trialedges = round(linspace(1, length(indsWN), N+1));
    Nthis = median(diff(trialedges));
    trialedges = indsWN(trialedges);
    trialedges(end)=trialedges(end)+1;
    
    
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
end

OUTSTRUCT_XCOV.LearnTargDir_Z_Epochs = LearnTargDir_Z_XCOV;
OUTSTRUCT_XCOV.LearnTargDir_Z_XCOV_Split_Adaptive = LearnTargDir_Z_XCOV_Split_Adaptive;
OUTSTRUCT_XCOV.LearnTargDir_Z_XCOV_Split_NonAdap = LearnTargDir_Z_XCOV_Split_NonAdap;
OUTSTRUCT_XCOV.NperBin_Epochs = Nall;
OUTSTRUCT_XCOV.TimeBins_Base_Epochs = TimeBins_Base_Epochs;
