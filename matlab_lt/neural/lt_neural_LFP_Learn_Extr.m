function [OUTSTRUCT, OUTSTRUCT_CohMatOnly, OUTSTRUCT_LFPOnly] = lt_neural_LFP_Learn_Extr(SwitchStruct, SwitchCohStruct, ...
    plotON, averagechanpairs, PARAMS, onlyfirstswitch, removeBadSyls, ...
    collectAllProcess, zscoreLFP, collectDiffMats, removeBadChans, typesToRemove, ...
    saveSpec)
%%

if ~exist('saveSpec', 'var')
    saveSpec = 0;
end
%% lt 11/1/18 - added option to collect also phi and psds

if collectAllProcess==1
    % then cannot get all trials
    collectonlyMeans=1;
else
    collectonlyMeans=0;
end
% NOTE: regardless, will put all cohmat into a structure OUTSTRUCT_CohMatOnly


%% lt 10/9/18 - % THIS 1) plots and 2) collects to do stats

% removeBadSyls. leave at 1. this removes syls that follow target (i.e.
% acute effects. and also things that are "same-type" but actually a
% previous target..

%%
% plotON=0;
tbins = PARAMS.tbins;
ffbins = PARAMS.ffbins;
numbirds = length(SwitchCohStruct.bird);

%%
% =============== INITIATE
OUTSTRUCT.S1Mat = {};
OUTSTRUCT.S2Mat = {};
OUTSTRUCT.CohMat = {};
OUTSTRUCT.PhiMat = {};
OUTSTRUCT.bnum = [];
OUTSTRUCT.enum = [];
OUTSTRUCT.switch = [];
OUTSTRUCT.motifnum = [];
OUTSTRUCT.issame = [];
OUTSTRUCT.istarg = [];
OUTSTRUCT.tvals = {};
OUTSTRUCT.indsbase = {};
OUTSTRUCT.indsWN = {};

OUTSTRUCT.chanpair = [];
OUTSTRUCT.bregionpair = {};
OUTSTRUCT.bregionpairs_originalorder = {};
OUTSTRUCT.motifname = {};

OUTSTRUCT.CohMean_WNminusBase = {};
OUTSTRUCT.CohMean_WN = {};
OUTSTRUCT.CohMean_Base = {};

OUTSTRUCT.Spec1Mean_WNminusBase = {};
OUTSTRUCT.Spec1Mean_WN= {};
OUTSTRUCT.Spec1Mean_Base = {};

OUTSTRUCT.Spec2Mean_WNminusBase = {};
OUTSTRUCT.Spec2Mean_WN= {};
OUTSTRUCT.Spec2Mean_Base = {};

OUTSTRUCT.indsbase_epoch = {};
OUTSTRUCT.indsWN_epoch = {};

OUTSTRUCT.cohscal = {};
OUTSTRUCT.ffvals = {};

scalarcounter = 1;

% =============== RUN
for i=1:numbirds
    bthis = SwitchStruct.bird(i).birdname;
    numexpts = length(SwitchCohStruct.bird(i).exptnum);
    for ii=1:numexpts
        ename = SwitchStruct.bird(i).exptnum(ii).exptname;
        numswitch = length(SwitchCohStruct.bird(i).exptnum(ii).switchlist);
        foundatarg =0;
        for ss=1:numswitch
            if onlyfirstswitch==1
                if ss>1
                    continue
                end
            end
            
            nummotifs = length(SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum);
            datthis = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss);
            if isempty(datthis.motifnum)
                continue
            end
            % =========== for this switch, get edge times and switch time
            tstart = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum_previous;
            tswitch = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum;
            tend = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum_next;
            targsyls = {SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningDirs{1:2:end}};
            targdirs = [SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningDirs{2:2:end}];
            motiflist = {datthis.motifnum.motifname};
            motiflist = motiflist(~cellfun(@isempty, motiflist));
            sylssame = lt_neural_QUICK_ExtractSameType(motiflist, targsyls);
            
            % =========== one figure for each switch
            figcount=1;
            subplotrows=3;
            subplotcols=6;
            fignums_alreadyused=[];
            hfigs=[];
            hsplots = [];
            if isempty(targsyls)
                keyboard
            end
            
            %             disp(['TARG SYLS: ' num2str(targsyls)]);
            for mm=1:nummotifs
                motifthis = datthis.motifnum(mm).motifname;
                
                if removeBadSyls==1
                    sylbad = lt_neural_QUICK_LearnRemoveBadSyl(bthis, ename, ...
                        ss, motifthis, typesToRemove);
                    if sylbad==1
                        continue
                    end
                end
                
                if isempty(datthis.motifnum(mm).tvals)
                    continue
                end
                
                %% ========= [COHERENCE] COLLECT MEAN COHERENCE PRE AND POST SWITCH
                %% ================= EXTRACT COHERENCE DATA
                if isfield(datthis.motifnum(mm), 'cohmat')
                    cohmat = datthis.motifnum(mm).cohmat;
                else
                    % --- then go and load
                    filename = [datthis.motifnum(mm).fileprefix '/Coh' datthis.motifnum(mm).filesuffix];
                    pairstoget = datthis.motifnum(mm).chanpairstokeep;
                    cohmat = lt_neural_LFP_loadProcessDat(filename, pairstoget);
                end
                
                if isempty(cohmat)
                    continue
                end
                
                % shoudl save one series for each channel pair
                chanpairs = datthis.motifnum(mm).chanpair;
                bregionpairs = datthis.motifnum(mm).bregionpair;
                bregionpairs_originalorder = datthis.motifnum(mm).bregionpair_originalorder;
                
                
                
                % ###################### [OTHER PROCESSED DATA]
                if collectAllProcess==1
                    pairstoget = datthis.motifnum(mm).chanpairstokeep;
                    
                    % ======== 1) PHI
                    filename = [datthis.motifnum(mm).fileprefix '/phi' datthis.motifnum(mm).filesuffix];
                    phimat = lt_neural_LFP_loadProcessDat(filename, pairstoget);
                    
                    % ====== 2) S1
                    filename = [datthis.motifnum(mm).fileprefix '/S1' datthis.motifnum(mm).filesuffix];
                    S1mat = lt_neural_LFP_loadProcessDat(filename, pairstoget);
                    
                    % ====== 2) S2
                    filename = [datthis.motifnum(mm).fileprefix '/S2' datthis.motifnum(mm).filesuffix];
                    S2mat = lt_neural_LFP_loadProcessDat(filename, pairstoget);
                end
                
                
                
                % ######################## REMOVE ANY PAIRS THAT INCLUDE A
                % --- go thru each channel pair and ask if it is good.
                numpairs = size(chanpairs,1);
                if removeBadChans==1
                    chanbad_all = [];
                    for k=1:numpairs
                        chanbad = lt_neural_QUICK_RemoveBadChans(bthis, ename, ss, chanpairs(k,:));
                        chanbad_all = [chanbad_all; chanbad];
                    end
                    % ---- DO REMOVAL
                    chanpairs = chanpairs(~chanbad_all, :);
                    bregionpairs = bregionpairs(~chanbad_all);
                    bregionpairs_originalorder = bregionpairs_originalorder(~chanbad_all);
                    cohmat = cohmat(:,:,:,~chanbad_all);
                    phimat = phimat(:,:,:,~chanbad_all);
                    S1mat = S1mat(:,:,:,~chanbad_all);
                    S2mat = S2mat(:,:,:,~chanbad_all);
                end
                
                %% ############### take mean across all channels paiors
                if averagechanpairs ==1
                    cohmat = nanmean(cohmat, 4);
                    S1mat = nanmean(S1mat, 4);
                    S2mat = nanmean(S2mat, 4);
                end
                
                
                %% =============== CALCULATE DEVIATION (WN MINUS BASE)
                tvals = datthis.motifnum(mm).tvals;
                indsbase = tvals>tstart & tvals<tswitch;
                indsWN = tvals>tswitch & tvals<tend;

                % ======================== REMOVE BAD TRIALS?

                
                % ---------------- EPOCHS
                indsbase_epoch = datthis.motifnum(mm).indsbase_epoch;
                indsWN_epoch = datthis.motifnum(mm).indsWN_epoch;
                
                
                %                 indsbase_epoch = find(indsbase);
                %                 indsbase_epoch = indsbase_epoch(round(length(indsbase_epoch)/2):end);
                
                %                 indsWN_epoch = find(indsWN);
                %                 indsWN_epoch = indsWN_epoch(round(length(indsWN_epoch)/2):end);
                npairs = size(cohmat,4);
                
                % ------- COHERENCE
                cohmean_diff = nan(size(cohmat,1), size(cohmat,2), npairs);
                cohmean_base = nan(size(cohmat,1), size(cohmat,2), npairs);
                cohmean_wn = nan(size(cohmat,1), size(cohmat,2), npairs);
                
                for j=1:npairs
                    cohmean_base_this = nanmean(cohmat(:,:, indsbase_epoch, j),3);
                    cohmean_WN_this = nanmean(cohmat(:,:, indsWN_epoch, j),3);
                    cohmean_diff(:,:,j) = cohmean_WN_this - cohmean_base_this;
                    
                    cohmean_base(:,:,j) = cohmean_base_this;
                    cohmean_wn(:,:,j) = cohmean_WN_this;
                end
                
                
                % ------ PHI
                % don't do anything here yet (since want to keep entire
                % distribution
                
                
                % ------ S1
                S1mean_diff = nan(size(cohmat,1), size(cohmat,2), npairs);
                S1mean_base = nan(size(cohmat,1), size(cohmat,2), npairs);
                S1mean_WN = nan(size(cohmat,1), size(cohmat,2), npairs);
                for j=1:npairs
                    
                    S_base = S1mat(:,:, indsbase_epoch, j);
                    S_WN = S1mat(:,:, indsWN_epoch, j);
                    
                    if zscoreLFP ==1
                        % --- do zscore for each t,ff bin (relative to baseline
                        % distribution)
                        [S_base, S_WN] = fn_zscorespec(S_base, S_WN);
                    elseif zscoreLFP ==2
                        % --- for each slice and trial normalize as power
                        % proportion (i.e. divide by sum across
                        % frequencies)
                        S_base = lt_neural_LFP_NormSpec(S_base);
                        S_WN = lt_neural_LFP_NormSpec(S_WN);
                    elseif zscoreLFP ==3
                        % each ff window z-scored
                        [S_base, S_WN] = fn_zscorespec_withinFF(S_base, S_WN);
                    elseif zscoreLFP==4
                        [S_base, S_WN] = fn_zscorespec_withinFF(S_base, S_WN);
                        S_base = lt_neural_LFP_NormSpec(S_base);
                        S_WN = lt_neural_LFP_NormSpec(S_WN);
                    end
                    
                    % get means
                    Specmean_base = nanmean(S_base,3);
                    Specmean_WN = nanmean(S_WN,3);
                    
                    % save
                    S1mean_diff(:,:,j) = Specmean_WN - Specmean_base;
                    S1mean_base(:,:,j) = Specmean_base;
                    S1mean_WN(:,:,j) = Specmean_WN;
                end
                
                
                % ------ S2
                S2mean_diff = nan(size(cohmat,1), size(cohmat,2), npairs);
                S2mean_base = nan(size(cohmat,1), size(cohmat,2), npairs);
                S2mean_WN = nan(size(cohmat,1), size(cohmat,2), npairs);
                for j=1:npairs
                    
                    S_base = S2mat(:,:, indsbase_epoch, j);
                    S_WN = S2mat(:,:, indsWN_epoch, j);
                    
                    if zscoreLFP ==1
                        % --- do zscore for each t,ff bin (relative to baseline
                        % distribution)
                        [S_base, S_WN] = fn_zscorespec(S_base, S_WN);
                    elseif zscoreLFP ==2
                        % --- for each slice and trial normalize as power
                        % proportion (i.e. divide by sum across
                        % frequencies)
                        S_base = lt_neural_LFP_NormSpec(S_base);
                        S_WN = lt_neural_LFP_NormSpec(S_WN);
                    elseif zscoreLFP ==3
                        % each ff window z-scored
                        [S_base, S_WN] = fn_zscorespec_withinFF(S_base, S_WN);
                    elseif zscoreLFP==4
                        [S_base, S_WN] = fn_zscorespec_withinFF(S_base, S_WN);
                        S_base = lt_neural_LFP_NormSpec(S_base);
                        S_WN = lt_neural_LFP_NormSpec(S_WN);
                    end
                    
                    % get means
                    Specmean_base = nanmean(S_base,3);
                    Specmean_WN = nanmean(S_WN,3);
                    
                    S2mean_diff(:,:,j) = Specmean_WN - Specmean_base;
                    S2mean_base(:,:,j) = Specmean_base;
                    S2mean_WN(:,:,j) = Specmean_WN;
                end
                
                %% ============== GET SCALAR VALUES OF COHERENCE
                if isfield(datthis.motifnum(mm), 'cohscalar')
                    cohscalar = datthis.motifnum(mm).cohscalar';
                else
                    if scalarcounter==1
                        scalarcounter=2;
                        disp('NOTE: skipping coh scalar as not extracted previously!!');
                        pause
                    end
                end
                
                %% ============= get ff
                ffvals = datthis.motifnum(mm).ffvals;
                
                %% ============== GET SPECTRUM POWER BY Z-SCORING
                % z-score separtely for each (t, ff) bin
                
                
                
                %% OTHER PARAMS FOR THIS SWITCH
                
                
                % --- TARG? what type of syl is this?
                if any(strcmp(targsyls, motifthis))
                    istarg =1;
                    learndir = targdirs(find(strcmp(targsyls, motifthis)));
                    titlecol = 'r';
                else
                    istarg = 0;
                    learndir = [];
                    titlecol = 'k';
                end
                
                % --- SAME?
                issame = any(strcmp(sylssame, motifthis));
                
                if istarg==1
                    foundatarg=1;
                end
                
                %% sanity check plots
                if plotON==1
                    % ========= PLOT COHEROGRAMS
                    % -- baseline
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    title(['[BASE]' motifthis], 'Color', titlecol);
                    ylabel([bthis '-' ename '-sw' num2str(ss)]);
                    cohthis = cohmat(:,:,indsbase);
                    lt_neural_Coher_Plot(cohthis, tbins, ffbins, 1);
                    % -- WN
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    title(['[WN]'], 'Color', titlecol);
                    cohthis = cohmat(:,:,indsWN);
                    lt_neural_Coher_Plot(cohthis, tbins, ffbins, 1);
                    % =========== OVERLAY BASELINE AND WN FREQUENCY BANDS
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    title(['[base(dash), WN(solid)]'], 'Color',titlecol);
                    
                    cohthis = cohmat(:,:,indsbase);
                    lt_neural_Coher_Plot(cohthis, tbins, ffbins, 2, ':');
                    
                    cohthis = cohmat(:,:,indsWN);
                    lt_neural_Coher_Plot(cohthis, tbins, ffbins, 2, '-');
                    
                    % ========= NOTE DOWN SYLTYPE
                end
                
                %% ######################## EXTRACT
                if averagechanpairs==1
                    % then there is one series of cohmat, averaged over
                    % chanel p[airs
                    OUTSTRUCT.CohMat = [OUTSTRUCT.CohMat; cohmat];
                    OUTSTRUCT.S1Mat = [OUTSTRUCT.S1Mat; S1mat];
                    OUTSTRUCT.S2Mat = [OUTSTRUCT.S2Mat; S2mat];
                    OUTSTRUCT.bnum = [OUTSTRUCT.bnum; i];
                    OUTSTRUCT.enum = [OUTSTRUCT.enum; ii];
                    OUTSTRUCT.switch = [OUTSTRUCT.switch; ss];
                    OUTSTRUCT.motifnum = [OUTSTRUCT.motifnum; mm];
                    OUTSTRUCT.motifname = [OUTSTRUCT.motifname; motifthis];
                    OUTSTRUCT.issame = [OUTSTRUCT.issame; issame];
                    OUTSTRUCT.istarg = [OUTSTRUCT.istarg; istarg];
                    OUTSTRUCT.tvals = [OUTSTRUCT.tvals; tvals];
                    OUTSTRUCT.indsbase = [OUTSTRUCT.indsbase; indsbase];
                    OUTSTRUCT.indsWN = [OUTSTRUCT.indsWN; indsWN];
                else
                    
                    numpairs = size(chanpairs,1);
                    for k=1:numpairs
                        % -- save for this pair
                        OUTSTRUCT.CohMat = [OUTSTRUCT.CohMat; single(cohmat(:,:,:,k))];

                        OUTSTRUCT.S1Mat = [OUTSTRUCT.S1Mat; single(S1mat(:,:,:,k))];
                        OUTSTRUCT.S2Mat = [OUTSTRUCT.S2Mat; single(S2mat(:,:,:,k))];
                        
                        OUTSTRUCT.PhiMat = [OUTSTRUCT.PhiMat; single(phimat(:,:,:,k))];
                        
                        OUTSTRUCT.bnum = [OUTSTRUCT.bnum; i];
                        OUTSTRUCT.enum = [OUTSTRUCT.enum; ii];
                        OUTSTRUCT.switch = [OUTSTRUCT.switch; ss];
                        OUTSTRUCT.motifnum = [OUTSTRUCT.motifnum; mm];
                        OUTSTRUCT.motifname = [OUTSTRUCT.motifname; motifthis];
                        OUTSTRUCT.issame = [OUTSTRUCT.issame; issame];
                        OUTSTRUCT.istarg = [OUTSTRUCT.istarg; istarg];
                        OUTSTRUCT.tvals = [OUTSTRUCT.tvals; tvals];
                        OUTSTRUCT.ffvals = [OUTSTRUCT.ffvals; single(ffvals)];
                        
                        OUTSTRUCT.indsbase = [OUTSTRUCT.indsbase; indsbase];
                        OUTSTRUCT.indsWN = [OUTSTRUCT.indsWN; indsWN];
                        
                        OUTSTRUCT.chanpair = [OUTSTRUCT.chanpair; chanpairs(k,:)];
                        OUTSTRUCT.bregionpair = [OUTSTRUCT.bregionpair; bregionpairs{k}];
                        
                        OUTSTRUCT.bregionpairs_originalorder = ...
                            [OUTSTRUCT.bregionpairs_originalorder; bregionpairs_originalorder{k}];
                        
                        
                        OUTSTRUCT.CohMean_WN = [OUTSTRUCT.CohMean_WN; single(cohmean_wn(:,:,k))];
                        OUTSTRUCT.CohMean_Base= [OUTSTRUCT.CohMean_Base; single(cohmean_base(:,:,k))];
                        
                        OUTSTRUCT.Spec1Mean_Base = [OUTSTRUCT.Spec1Mean_Base; single(S1mean_base(:,:,k))];
                        OUTSTRUCT.Spec1Mean_WN = [OUTSTRUCT.Spec1Mean_WN; single(S1mean_WN(:,:,k))];
                        
                        OUTSTRUCT.Spec2Mean_Base = [OUTSTRUCT.Spec2Mean_Base; single(S2mean_base(:,:,k))];
                        OUTSTRUCT.Spec2Mean_WN = [OUTSTRUCT.Spec2Mean_WN; single(S2mean_WN(:,:,k))];
                        
                        OUTSTRUCT.indsbase_epoch = [OUTSTRUCT.indsbase_epoch; single(indsbase_epoch)];
                        OUTSTRUCT.indsWN_epoch = [OUTSTRUCT.indsWN_epoch; single(indsWN_epoch)];
                        
                        if scalarcounter==1
                            OUTSTRUCT.cohscal = [OUTSTRUCT.cohscal; cohscalar(k,:)];
                        end
                        
                        if collectDiffMats ==1
                            OUTSTRUCT.CohMean_WNminusBase = [OUTSTRUCT.CohMean_WNminusBase; cohmean_diff(:,:,k)];
                            OUTSTRUCT.Spec1Mean_WNminusBase = [OUTSTRUCT.Spec1Mean_WNminusBase; S1mean_diff(:,:,k)];
                            OUTSTRUCT.Spec2Mean_WNminusBase = [OUTSTRUCT.Spec2Mean_WNminusBase; S2mean_diff(:,:,k)];
                        end
                        
                    end
                end
                
            end
        end
        if foundatarg==0 & numswitch>0
            % PROBLEM: this expt is missing motifs... see previous
            % extractio function, is likely eb cuase is throwing out
            % anything with nans.
            keyboard
        end
        if plotON==1
            pause; close all;
        end
    end
end


%% ======== 2) CALCULATE TRAIN MINUS BASE COHERENCE
if (0)
    cohdiff_all = {};
    indsbase_epoch = {};
    indsWN_epoch = {};
    % --- takes 2nd half of training and base
    for j=1:length(OUTSTRUCT.bnum)
        
        cohmat = OUTSTRUCT.CohMat{j};
        indsbase = find(OUTSTRUCT.indsbase{j});
        indsWN = find(OUTSTRUCT.indsWN{j});
        
        % -------------- take second half of WN inds
        indsbase = indsbase(round(length(indsbase)/2):end);
        %     try
        %     indsbase = indsbase(end-25:end);
        %     catch err
        %     end
        indsWN = indsWN(round(length(indsWN)/2):end);
        
        % -----
        cohmean_base = nanmean(cohmat(:,:, indsbase),3);
        cohmean_WN = nanmean(cohmat(:,:, indsWN),3);
        cohmean_diff = cohmean_WN - cohmean_base;
        
        if any(isnan(cohmean_diff(:)))
            keyboard
        end
        
        cohdiff_all = [cohdiff_all; cohmean_diff];
        indsbase_epoch = [indsbase_epoch; indsbase];
        indsWN_epoch = [indsWN_epoch; indsWN];
    end
    OUTSTRUCT.CohMean_WNminusBase = cohdiff_all;
    OUTSTRUCT.indsbase_epoch = indsbase_epoch;
    OUTSTRUCT.indsWN_epoch = indsWN_epoch;
end
%% ======== remove individual trials of cohernece>?
OUTSTRUCT_CohMatOnly = OUTSTRUCT.CohMat;
if collectonlyMeans==1
    
    OUTSTRUCT = rmfield(OUTSTRUCT, 'CohMat');
end

%% ============ SAVE SPECTROGRAMS
OUTSTRUCT_S1MatOnly = OUTSTRUCT.S1Mat;
OUTSTRUCT_S2MatOnly = OUTSTRUCT.S2Mat;
OUTSTRUCT = rmfield(OUTSTRUCT, 'S1Mat');
OUTSTRUCT = rmfield(OUTSTRUCT, 'S2Mat');

if saveSpec==1
    sdir = '/bluejay0/bluejay2/lucas/analyses/neural/COHERENCE/Learn_Extr';
    if ~exist([sdir '/' PARAMS.savemarker], 'dir')
        mkdir([sdir '/' PARAMS.savemarker]);
    end
    save([sdir '/' PARAMS.savemarker '/OUTSTRUCT_S1MatOnly.mat'], 'OUTSTRUCT_S1MatOnly');
    save([sdir '/' PARAMS.savemarker '/OUTSTRUCT_S2MatOnly.mat'], 'OUTSTRUCT_S2MatOnly');
end

%% for each case get wn minus base, mean coh scalar
cohscal_diff = [];
for i=1:length(OUTSTRUCT.bnum)
    
    indsbase = OUTSTRUCT.indsbase_epoch{i};
    indswn = OUTSTRUCT.indsWN_epoch{i};
    
    cohscal = OUTSTRUCT.cohscal{i};
    
    cohdiff = mean(cohscal(indswn)) - mean(cohscal(indsbase));
    
    cohscal_diff = [cohscal_diff; cohdiff];
end

OUTSTRUCT.cohscal_diff = cohscal_diff;


%% ====== COLLECT LEARN DIR AT TARGET

OUTSTRUCT = lt_neural_LFP_GetLearnDir(OUTSTRUCT, SwitchStruct);

%% ====== GET UNIQUE MOTIF ID FOR EACH CASE
motifID_unique = nan(length(OUTSTRUCT.bnum), 1);

for j=1:length(OUTSTRUCT.bnum)
    
    bnum = unique(OUTSTRUCT.bnum(j));
    bname = SwitchStruct.bird(bnum).birdname;
    
    enum = unique(OUTSTRUCT.enum(j));
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    
    motifname = OUTSTRUCT.motifname{j};
    
    try
        motifID = lt_neural_QUICK_MotifID(bname, {motifname});
        assert(~isempty(motifID));
    catch err
        motifID = nan;
    end
    motifID_unique(j) = motifID;
    
end

OUTSTRUCT.motifID_unique = motifID_unique;


end

function [S_base_z, S_WN_z] = fn_zscorespec(S_base, S_WN)
% will zscore each t,ff bin based on baseline
% input sould be power (i.e not log). code will take log for you.

% -- take log [makes distrubtions more gaussina]
S_base = 10*log10(S_base);
S_WN = 10*log10(S_WN);

% get baseline mean and std
S_mean = mean(S_base,3);
S_std = std(S_base, [], 3);

% perform zscoring
S_base = S_base - repmat(S_mean, 1, 1, size(S_base,3));
S_base = S_base./repmat(S_std, 1, 1, size(S_base,3));
S_base_z = S_base;

S_WN = S_WN - repmat(S_mean, 1, 1, size(S_WN,3));
S_WN = S_WN./repmat(S_std, 1, 1, size(S_WN,3));
S_WN_z = S_WN;
end

function [S_base_z, S_WN_z] = fn_zscorespec_withinFF(S_base, S_WN)
% will zscore each t,ff bin based on baseline
% input sould be power (i.e not log). code will take log for you.

% -- take log [makes distrubtions more gaussina]
S_base = 10*log10(S_base);
S_WN = 10*log10(S_WN);

% get baseline mean and std
nffbin = size(S_base,2);
Smeans = nan(1, nffbin);
Sstds = nan(1, nffbin);
for n=1:nffbin
    tmp = S_base(:,n,:);
    s_mean = mean(tmp(:));
    s_std = std(tmp(:));
    
    Smeans(n) = s_mean;
    Sstds(n) = s_std;
end

S_base = S_base - repmat(Smeans, size(S_base,1), 1, size(S_base,3));
S_base_z = S_base./repmat(Sstds, size(S_base,1), 1, size(S_base,3));

S_WN = S_WN - repmat(Smeans, size(S_WN,1), 1, size(S_WN,3));
S_WN_z = S_WN./repmat(Sstds, size(S_WN,1), 1, size(S_WN,3));
end






