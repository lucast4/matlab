function [OUTSTRUCT, PARAMS] = lt_neural_Coher_RecalcMat(SwitchStruct, OUTSTRUCT, OUTSTRUCT_CohMatOnly, ...
    SwitchCohStruct, PARAMS, twind, fwind, wntouse, useWNtiming, ...
    prewind_relWN, COHSTRUCT, RemoveIfTooFewTrials, removebadtrialtype, extractLFP, lfpUseMedian, ...
    specScaleType, cohUseMedian, ...
    removebadsyl, normtoshuff, normtype, OUTSTRUCT_CohMatOnly_shift)
%% lt 2/21/19 - also extract mean LFP power...
disp('NOTE: normalization doesnt apply for LFP (i.e. normn to shuffle.)');


% disp('NOTE: !!! check line 319 - using median');
% pause;

%% lt 2/15/19 - uses new set of inds to recalcalate coherence matrices
% NOTE: saves old version as backup
Nmin = 10; % n trials... otherwise will ignore... [TO DO]


% == for getting scal rel wn
% - linearly interpolate coherence?
interpol =0; % across time  [NOT YET DONE].
onlyusegoodtargsyls = 1; % default = 1, i.e. gets timing only from that that will actually analyze.


%% ==== make sure datasets match up
if length(OUTSTRUCT.bnum) ~= size(OUTSTRUCT_CohMatOnly,1)
    disp('NOTE!!!! coh stuff will be incorrect - datasets are not matched...');
    pause;
end

%% === if extract LFP, have to first load lfp data
if extractLFP==1
    disp('LOADING LFP DATA');
    sdir = ['/bluejay0/bluejay2/lucas/analyses/neural/COHERENCE/Learn_Extr/' PARAMS.savemarker];
    load([sdir '/OUTSTRUCT_S1MatOnly.mat']);
    load([sdir '/OUTSTRUCT_S2MatOnly.mat']);
    
    assert(length(OUTSTRUCT.bnum) == size(OUTSTRUCT_S1MatOnly,1));
    assert(length(OUTSTRUCT.bnum) == size(OUTSTRUCT_S2MatOnly,1));
    
    % =========== VARIOUS SANITY CHECKS
    
    try
        assert(length(PARAMS.tbins) == size(OUTSTRUCT_S2MatOnly{1},1));
        assert(length(PARAMS.ffbins) == size(OUTSTRUCT_S2MatOnly{1},2));
        assert(length(PARAMS.tbins) == size(OUTSTRUCT_S1MatOnly{1},1));
        assert(length(PARAMS.ffbins) == size(OUTSTRUCT_S1MatOnly{1},2));
    catch err
        try
            % ================
            settmp = 1;
            if isempty(COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif)
                settmp=2;
                assert(~isempty(COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif));
            end
            PARAMS.tbins = COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif(1).t_relons;
            
            assert(length(PARAMS.tbins) == size(OUTSTRUCT_S2MatOnly{1},1));
            assert(length(PARAMS.ffbins) == size(OUTSTRUCT_S2MatOnly{1},2));
            assert(length(PARAMS.tbins) == size(OUTSTRUCT_S1MatOnly{1},1));
            assert(length(PARAMS.ffbins) == size(OUTSTRUCT_S1MatOnly{1},2));
        catch err
            PARAMS.tbins = PARAMS.tbins_old;
            PARAMS.ffbins = PARAMS.ffbins_old;
            
            assert(length(PARAMS.tbins) == size(OUTSTRUCT_S2MatOnly{1},1));
            assert(length(PARAMS.ffbins) == size(OUTSTRUCT_S2MatOnly{1},2));
            assert(length(PARAMS.tbins) == size(OUTSTRUCT_S1MatOnly{1},1));
            assert(length(PARAMS.ffbins) == size(OUTSTRUCT_S1MatOnly{1},2));
        end
    end
end


%% ==== make sure datasets match up
if length(OUTSTRUCT.bnum) == size(OUTSTRUCT_CohMatOnly,1)
    try
    assert(length(PARAMS.tbins)==size(OUTSTRUCT_CohMatOnly{1},1), 'tbins not corrent unsure why..');
    assert(length(PARAMS.ffbins)==size(OUTSTRUCT_CohMatOnly{1},2), 'fbins not corrent unsure why..');
    catch err
                    settmp = 1;
            if isempty(COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif)
                settmp=2;
                assert(~isempty(COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif));
            end
            
            PARAMS.tbins = COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif(1).t_relons;
            PARAMS.ffbins = COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif(1).ffbins;
    
            
            assert(length(PARAMS.tbins)==size(OUTSTRUCT_CohMatOnly{1},1), 'tbins not corrent unsure why..');
    assert(length(PARAMS.ffbins)==size(OUTSTRUCT_CohMatOnly{1},2), 'fbins not corrent unsure why..');
    end
end

%%

Cohmean_base = {};
Cohmean_wn = {};

cohscal_all = {};
cohscal_diff_all = [];

specscal_chan1_all = {};
specdiff_chan1_all = [];

specscal_chan2_all = {};
specdiff_chan2_all = [];

SpecMean_Base_Chan1 = {};
SpecMean_WN_Chan1 = {};

SpecMean_Base_Chan2 = {};
SpecMean_WN_Chan2 = {};

% ========= params for scalars
indT = PARAMS.tbins>=twind(1) & PARAMS.tbins<=twind(2);
indF = PARAMS.ffbins>=fwind(1) & PARAMS.ffbins<=fwind(2);

% if length(PARAMS.tbins)~=size(OUTSTRUCT_CohMatOnly{1},1)
%     disp('REEXTRACTING tbine: been overwritten with wrong one (wn realign...)');
%     % reextract correct tbins
%     % ================
%     settmp = 1;
%     if isempty(COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif)
%         settmp=2;
%         assert(~isempty(COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif));
%     end
%     PARAMS.tbins = COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif(1).t_relons;
% end
%

% ===============
for i=1:length(OUTSTRUCT.bnum)
    disp(i);
    
    indsbase = find(OUTSTRUCT.indsbase{i});
    indswn = find(OUTSTRUCT.indsWN{i});
    ntrial = length(OUTSTRUCT.ffvals);
    
    
    % ====================== GET COHERENCE DATA
    bnum = OUTSTRUCT.bnum(i);
    enum = OUTSTRUCT.enum(i);
    sw = OUTSTRUCT.switch(i);
    mm = OUTSTRUCT.motifnum(i);
    mID = OUTSTRUCT.motifID_unique(i);
    mname = OUTSTRUCT.motifname{i};
    bname = SwitchStruct.bird(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;

    
    % =============== REMOVE BAD TRIALS IF DESIRED.
    tvals = OUTSTRUCT.tvals{i};
    if ~isempty(removebadtrialtype)
        badtrials = lt_neural_QUICK_RemoveTrials(bname, ename, sw, tvals, ...
            removebadtrialtype);
        
        badtrials=find(badtrials);
        
        % -- remove bad trials
        indsbase(ismember(indsbase, badtrials)) = [];
        indswn(ismember(indswn, badtrials)) = [];
    end
    
   
    
    % ======================
    indsbase_epoch = indsbase(end-round(length(indsbase)/2):end);
    
    if strcmp(wntouse, 'half')
        indswn_epoch = indswn(end-round(length(indswn)/2):end);
    elseif strcmp(wntouse, 'quarter')
        indswn_epoch = indswn(end-round(length(indswn)/4):end);
    elseif strcmp(wntouse, 'third')
        indswn_epoch = indswn(end-round(length(indswn)/3):end);
    elseif strcmp(wntouse, 'firsthalf')
        indswn_epoch = indswn(1:round(length(indswn)/2));
    end
    
    
    
    
    
   if removebadsyl==1
        sylbad = lt_neural_QUICK_LearnRemoveBadSyl(bname, ename, sw, mname, ...
            {'wn'});
        if sylbad==1
            Cohmean_base = [Cohmean_base; nan(length(indT), length(indF))];
            Cohmean_wn = [Cohmean_wn; nan(length(indT), length(indF))];
            
            SpecMean_Base_Chan1 = [SpecMean_Base_Chan1; nan(length(indT), length(indF))];
            SpecMean_WN_Chan1 = [SpecMean_WN_Chan1; nan(length(indT), length(indF))];
            SpecMean_Base_Chan2 = [SpecMean_Base_Chan2; nan(length(indT), length(indF))];
            SpecMean_WN_Chan2 = [SpecMean_WN_Chan2; nan(length(indT), length(indF))];
            
            cohscal_all = [cohscal_all; nan(1,ntrial)];
            cohscal_diff_all = [cohscal_diff_all; nan];
            specscal_chan1_all = [specscal_chan1_all; nan];
            specdiff_chan1_all = [specdiff_chan1_all; nan];
            specscal_chan2_all = [specscal_chan2_all; nan];
            specdiff_chan2_all = [specdiff_chan2_all; nan];
            continue
        end
    end
    
    
   
    
    if RemoveIfTooFewTrials==1
        if length(indsbase_epoch)<Nmin | length(indswn_epoch)<Nmin
            disp(['SKIP SYL (N too low): '  num2str(length(indsbase_epoch)) '-' num2str(length(indswn_epoch))]);
            Cohmean_base = [Cohmean_base; nan(length(indT), length(indF))];
            Cohmean_wn = [Cohmean_wn; nan(length(indT), length(indF))];
            SpecMean_Base_Chan1 = [SpecMean_Base_Chan1; nan(length(indT), length(indF))];
            SpecMean_WN_Chan1 = [SpecMean_WN_Chan1; nan(length(indT), length(indF))];
            SpecMean_Base_Chan2 = [SpecMean_Base_Chan2; nan(length(indT), length(indF))];
            SpecMean_WN_Chan2 = [SpecMean_WN_Chan2; nan(length(indT), length(indF))];
            
            cohscal_all = [cohscal_all; nan(1,ntrial)];
            cohscal_diff_all = [cohscal_diff_all; nan];
            specscal_chan1_all = [specscal_chan1_all; nan];
            specdiff_chan1_all = [specdiff_chan1_all; nan];
            specscal_chan2_all = [specscal_chan2_all; nan];
            specdiff_chan2_all = [specdiff_chan2_all; nan];
            continue
        end
    end
    
    
    %% ###################### LFP STUFF
    if extractLFP==1
        
        % =============== SITE 1
        SmatThis = OUTSTRUCT_S1MatOnly;
        
        % --- run
        specmat = SmatThis{i};
        % --- 1) RESCALE SPECMAT
        [~, specmat] = lt_neural_LFP_ScaleSpecgram(specmat(:,:, indsbase_epoch), ...
            specmat, specScaleType);
        
        % -- 1) extract trial vector of power scalar
        specscal = squeeze(mean(mean(specmat(indT, indF, :), 1),2));
        % -- 2) scalar difference (train minus base)
        if lfpUseMedian==1
            specdiff = median(specscal(indswn_epoch)) - median(specscal(indsbase_epoch));
        else
            specdiff = mean(specscal(indswn_epoch)) - mean(specscal(indsbase_epoch));
        end
        specscal_chan1_all = [specscal_chan1_all; specscal];
        specdiff_chan1_all = [specdiff_chan1_all; specdiff];
        
        % --------------- EXTRACT BASELINE AND WN MEAN SPECTROGRAM
        if lfpUseMedian==1
            sgram_base = median(specmat(:,:, indsbase_epoch),3);
            sgram_WN = median(specmat(:,:, indswn_epoch),3);
        else
            sgram_base = mean(specmat(:,:, indsbase_epoch),3);
            sgram_WN = mean(specmat(:,:, indswn_epoch),3);
        end
        % -- save
        SpecMean_Base_Chan1 = [SpecMean_Base_Chan1; sgram_base];
        SpecMean_WN_Chan1 = [SpecMean_WN_Chan1; sgram_WN];
        
        
        % =============== SITE 2
        SmatThis = OUTSTRUCT_S2MatOnly;
        
        % --- run
        specmat = SmatThis{i};
        % --- 1) RESCALE SPECMAT
        [~, specmat] = lt_neural_LFP_ScaleSpecgram(specmat(:,:, indsbase_epoch), ...
            specmat, specScaleType);
        % -- 1) extract trial vector of power scalar
        specscal = squeeze(mean(mean(specmat(indT, indF, :), 1),2));
        % -- 2) scalar difference (train minus base)
        if lfpUseMedian==1
            specdiff = median(specscal(indswn_epoch)) - median(specscal(indsbase_epoch));
        else
            specdiff = mean(specscal(indswn_epoch)) - mean(specscal(indsbase_epoch));
        end
        specscal_chan2_all = [specscal_chan2_all; specscal];
        specdiff_chan2_all = [specdiff_chan2_all; specdiff];
        
        % --------------- EXTRACT BASELINE AND WN MEAN SPECTROGRAM
        if lfpUseMedian==1
            sgram_base = median(specmat(:,:, indsbase_epoch),3);
            sgram_WN = median(specmat(:,:, indswn_epoch),3);
        else
            sgram_base = mean(specmat(:,:, indsbase_epoch),3);
            sgram_WN = mean(specmat(:,:, indswn_epoch),3);
        end
        
        % -- save
        SpecMean_Base_Chan2 = [SpecMean_Base_Chan2; sgram_base];
        SpecMean_WN_Chan2 = [SpecMean_WN_Chan2; sgram_WN];
    end
    
    
    
    %% =============== recalc coherence mat
    if isempty(OUTSTRUCT_CohMatOnly)
        Cohmean_base = [Cohmean_base; nan(length(indT), length(indF))];
        Cohmean_wn = [Cohmean_wn; nan(length(indT), length(indF))];
        cohscal_all = [cohscal_all; nan(1,ntrial)];
        cohscal_diff_all = [cohscal_diff_all; nan];
        continue
    end
    
    % ================== GET RAW COHMAT
    if (1)
        cohmat = OUTSTRUCT_CohMatOnly{i};
    else
        % === olde version, go back and load - tajkes a while since does
        % once for each case.
        chpair = OUTSTRUCT.chanpair(i,:);
        DAT = SwitchCohStruct.bird(bnum).exptnum(enum).switchlist(sw).motifnum(mm);
        fpre = DAT.fileprefix;
        fpost = DAT.filesuffix;
        
        tmp = ismember(DAT.chanpair, chpair, 'rows');
        chind = DAT.chanpairstokeep(tmp);
        cohmat = lt_neural_LFP_loadProcessDat([fpre '/Coh' fpost], chind);
    end
    
    % ==== skip if cohmat is empyt.. [previuosly skipped during recalc...]
    if isempty(cohmat)
        Cohmean_base = [Cohmean_base; nan(length(indT), length(indF))];
        Cohmean_wn = [Cohmean_wn; nan(length(indT), length(indF))];
        cohscal_all = [cohscal_all; nan(1,ntrial)];
        cohscal_diff_all = [cohscal_diff_all; nan];
        %         specscal_chan1_all = [specscal_chan1_all; nan];
        %         specdiff_chan1_all = [specdiff_chan1_all; nan];
        %         specscal_chan2_all = [specscal_chan2_all; nan];
        %         specdiff_chan2_all = [specdiff_chan2_all; nan];
        continue
    end
    
    
    
    cohbase = [];
    cohwn = [];
    if normtoshuff==1
        %% norm to shuffle?
        cohmat_shuff = OUTSTRUCT_CohMatOnly_shift{i};
        
        % ---- 1) GET SHUFFLE TRIALS THAT MATCH DATA TRIALS
        ntrials = size(cohmat,3);
        tmp = sort([1 2:ntrials-1 2:ntrials-1 ntrials]); % tmp(2) tells you the ind (in real dat) that shuff trial 2 corresponds to
        assert(length(tmp)==size(cohmat_shuff,3));
        shuffinds = tmp;
        
        inds_base_shuff = find(ismember(tmp, indsbase_epoch));
        inds_WN_shuff = find(ismember(tmp, indswn_epoch));
        
        % ---- 2) GET Z=SCORE FOR EACH T/F BIN, OVER ALL TRIALS
        shuffmean = mean(cohmat_shuff(:,:, [inds_base_shuff inds_WN_shuff]),3);
        
        % BASE
        cohtmp = cohmat_shuff(:,:, inds_base_shuff);
        ymean = mean(cohtmp,3);
        ystd = std(cohtmp, [], 3);
        if strcmp(normtype, 'zscore')
            cohbase = (mean(cohmat(:,:, indsbase_epoch), 3) - ymean)./ystd;
        elseif strcmp(normtype, 'minus')
            cohbase = mean(cohmat(:,:, indsbase_epoch), 3) - ymean;
        elseif strcmp(normtype, 'minusglob')
            cohbase = mean(cohmat(:,:, indsbase_epoch), 3) - shuffmean;
        end
        
        % WN
        cohtmp = cohmat_shuff(:,:, inds_WN_shuff);
        ymean = mean(cohtmp,3);
        ystd = std(cohtmp, [], 3);
        
        if strcmp(normtype, 'zscore')
            cohwn = (mean(cohmat(:,:, indswn_epoch), 3) - ymean)./ystd;
        elseif strcmp(normtype, 'minus')
            cohwn = mean(cohmat(:,:, indswn_epoch), 3) - ymean;
        elseif strcmp(normtype, 'minusglob')
            cohwn = mean(cohmat(:,:, indswn_epoch), 3) - shuffmean;
        end
    else
        % =================== GET PRE AND POST
        cohbase = mean(cohmat(:,:, indsbase_epoch), 3);
        cohwn = mean(cohmat(:,:, indswn_epoch), 3);
    end
    
    % =========== SAVE IN OUTSTRUCT
    Cohmean_base = [Cohmean_base; cohbase];
    Cohmean_wn = [Cohmean_wn; cohwn];
    
    
    
    %% ========== recalc coherence scalar
    if useWNtiming==1
        % -- recalculate indT
        wnontime = OUTSTRUCT.WNonset(i);
        
        twindtmp = wnontime + prewind_relWN;
        indT = PARAMS.tbins>=twindtmp(1) & PARAMS.tbins<=twindtmp(2);
    end
    
    if normtoshuff==1 & strcmp(normtype, 'zscore')
        cohscal = {[]}; % ignore, since must figure out how to normalize locally..
        % --- then must get mean scalar directly from the normalized base
        % and WN mean matrices
        cohscal_base = squeeze(mean(mean(cohbase(indT, indF, :), 1),2));
        cohscal_wn = squeeze(mean(mean(cohwn(indT, indF, :), 1),2));
    elseif normtoshuff==1 & strcmp(normtype, 'minus')
        % --- then get mean scalar from previous means, buit still can get
        % scalar timecourse...
        cohscal_base = squeeze(mean(mean(cohbase(indT, indF, :), 1),2));
        cohscal_wn = squeeze(mean(mean(cohwn(indT, indF, :), 1),2));
        
        cohscal_prenorm = squeeze(mean(mean(cohmat(indT, indF, :), 1),2));
        cohscal_shuff= squeeze(mean(mean(cohmat_shuff(indT, indF, :), 1),2));
        cohscal_shuff = grpstats(cohscal_shuff, shuffinds); % take means to get same number of trials a data.
        cohscal = cohscal_prenorm - cohscal_shuff;
        
    elseif normtoshuff==1 & strcmp(normtype, 'minusglob')
        % --- then get mean scalar from previous means, buit still can get
        % scalar timecourse...
        cohscal_base = squeeze(mean(mean(cohbase(indT, indF, :), 1),2));
        cohscal_wn = squeeze(mean(mean(cohwn(indT, indF, :), 1),2));
        
        cohscal_shuff = squeeze(mean(mean(shuffmean(indT, indF), 1),2));
        cohscal_prenorm = squeeze(mean(mean(cohmat(indT, indF, :), 1),2));
        cohscal = cohscal_prenorm - cohscal_shuff;
        
        
    elseif normtoshuff==0 & strcmp(normtype, 'zscorerelbase')
        % -- each trial get as zscore relative to baseline
        cohscal = squeeze(mean(mean(cohmat(indT, indF, :), 1),2));
        basemean = mean(cohscal(indsbase_epoch));
        basestd = std(cohscal(indsbase_epoch));
        cohscal = (cohscal-basemean)./basestd;
        
        cohscal_base = mean(cohscal(indsbase_epoch));
        cohscal_wn = mean(cohscal(indswn_epoch));
    elseif normtoshuff==0
        % if no norming, or if norming but just subtraacting shiffled, then
        % this should be fine.
        cohscal = squeeze(mean(mean(cohmat(indT, indF, :), 1),2));
        % --- then can get scalars from timecourse of cohscalars
        cohscal_base = mean(cohscal(indsbase_epoch));
        cohscal_wn = mean(cohscal(indswn_epoch));
        %         cohscal_base = median(cohscal(indsbase_epoch));
        %         cohscal_wn = median(cohscal(indswn_epoch));
    end
    cohscal_diff = cohscal_wn - cohscal_base;
    
    % ====================== OUTPUT
    cohscal_all = [cohscal_all; cohscal'];
    cohscal_diff_all = [cohscal_diff_all; cohscal_diff];
    
    
    
end

% === save old versions
OUTSTRUCT.CohMean_Base_orig = OUTSTRUCT.CohMean_Base;
OUTSTRUCT.CohMean_WN_orig = OUTSTRUCT.CohMean_WN;

OUTSTRUCT.cohscal_orig = OUTSTRUCT.cohscal;
OUTSTRUCT.cohscal_diff_orig= OUTSTRUCT.cohscal_diff;

% ===== save new versions
OUTSTRUCT.CohMean_Base = Cohmean_base;
OUTSTRUCT.CohMean_WN = Cohmean_wn;

OUTSTRUCT.cohscal = cohscal_all;
OUTSTRUCT.cohscal_diff = cohscal_diff_all;


if extractLFP==1
    OUTSTRUCT.Spec1Mean_Base = SpecMean_Base_Chan1;
    OUTSTRUCT.Spec1Mean_WN = SpecMean_WN_Chan1;
    
    OUTSTRUCT.Spec2Mean_Base = SpecMean_Base_Chan2;
    OUTSTRUCT.Spec2Mean_WN = SpecMean_WN_Chan2;
    
    % =======
    OUTSTRUCT.specscal_chan1_all = specscal_chan1_all;
    OUTSTRUCT.specscal_chan2_all = specscal_chan2_all;
    
    OUTSTRUCT.specdiff_chan1_all = specdiff_chan1_all;
    OUTSTRUCT.specdiff_chan2_all = specdiff_chan2_all;
end


