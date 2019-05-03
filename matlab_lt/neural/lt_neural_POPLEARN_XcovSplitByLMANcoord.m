function [OUTSTRUCT_MeanRhoSplit, PARAMS, NanCountAll] = lt_neural_POPLEARN_XcovSplitByLMANcoord(...
    SwitchXCovStruct, SwitchStruct, PARAMS, SwitchCohStruct, dosmooth, dosmooth_sigma, ...
    xcovver, getxgram, removebadchans, getHiLoFFSplit, hilosplit_shuffver, OUTSTRUCT_XCOV, ...
    DATSTRUCT_POP, bregionToSplitBy)
%% 4/5/19 - like FF split, but split trials by LMAN-LMAN  ensemble coordination
% like from Hamish.
% Must first run lt_neural_POPLEARN_SylLocked_PlotTrials to extract

PARAMS.Xcov_ccLags = PARAMS.Xcov_ccLags_beforesmooth;

%% what do use as shuffle trials if do ffsplit

% hilosplit_shuffver=2;
% then use the exact same shuffle trials for both low and high (i.e. take data without splitting)

%% NOTE: if any trials have nan, then just ignroes in calcualting dat minus shuff
% i.e. uses nanmean. I do not expect there to be any nan for xcorr
% ("unbiased");

%% lt 2/2019 - extracts OUTSTRUCT_MeanRhoSplit for xcov analysis during learning.
% MinTotRends = 10; % skips if etiher base or Wn (all, not epoch) has fewer trials than this

plotraw=0;

%% EXTRACT PARAMS
getxgram_epochbins = size(OUTSTRUCT_XCOV.Xcovslice_epochs{1},1);

%% % ===================================
OUTSTRUCT_MeanRhoSplit = struct;

OUTSTRUCT_MeanRhoSplit.XcovBase = [];
OUTSTRUCT_MeanRhoSplit.XcovWN = [];

OUTSTRUCT_MeanRhoSplit.XcovgramBase = {};
OUTSTRUCT_MeanRhoSplit.XcovgramWN = {};

OUTSTRUCT_MeanRhoSplit.XcovgramWN_epochs = {};

OUTSTRUCT_MeanRhoSplit.neurpairnum = [];
OUTSTRUCT_MeanRhoSplit.neurpair = [];
OUTSTRUCT_MeanRhoSplit.bnum = [];
OUTSTRUCT_MeanRhoSplit.enum =[];
OUTSTRUCT_MeanRhoSplit.swnum = [];
OUTSTRUCT_MeanRhoSplit.motifnum = [];
OUTSTRUCT_MeanRhoSplit.issame = [];
OUTSTRUCT_MeanRhoSplit.istarg = [];
OUTSTRUCT_MeanRhoSplit.learndirTarg = [];

OUTSTRUCT_MeanRhoSplit.XcovBase_NoMinShuff = ...
    [];
OUTSTRUCT_MeanRhoSplit.XcovWN_NoMinShuff = ...
    [];

OUTSTRUCT_MeanRhoSplit.inds_base_epoch = {}; % used for analysis
OUTSTRUCT_MeanRhoSplit.inds_WN_epoch = {}; % used for analysis
OUTSTRUCT_MeanRhoSplit.inds_base_allgood = {}; % all, after remove bad trials
OUTSTRUCT_MeanRhoSplit.inds_WN_allgood = {}; % all, after remove bad trials
OUTSTRUCT_MeanRhoSplit.trialedges_epoch = {};


OUTSTRUCT_MeanRhoSplit.XcovgramWN_FFsplits_Base = {};
OUTSTRUCT_MeanRhoSplit.XcovgramWN_FFsplits_Epochs = {};

OUTSTRUCT_MeanRhoSplit.Xcov_DotProd_trials = {};

NanCountAll = []; % if Nan then means: (1) if "unbiased" normalikzation, then should not ever be nan.
% (2) if using "coeff" then means that one of the trials had variance of 0.
for i=1:length(SwitchXCovStruct.bird)
    bname = SwitchStruct.bird(i).birdname;
    for ii=1:length(SwitchXCovStruct.bird(i).exptnum)
        ename = SwitchStruct.bird(i).exptnum(ii).exptname;
        disp([i ii]);
        for iii=1:length(SwitchXCovStruct.bird(i).exptnum(ii).switchlist)
            nmotifs = length(SwitchXCovStruct.bird(i).exptnum(ii).switchlist(iii).motif);
            
            
            for mm=1:nmotifs
                motifthis = SwitchCohStruct.bird(i).exptnum(ii).switchlist(iii).motifnum(mm).motifname;
                if isempty(motifthis)
                    continue
                end
                
                datthis = SwitchXCovStruct.bird(i).exptnum(ii).switchlist(iii).motif(mm);
                
                %                 if i==1 & ii==2 & iii==1
                %                     keyboard
                %                 end
                
                %% =========== SKIP IF DOESNT HAVE WITHIN-REGION ENSEMBLE DATA
                indstmp = find(DATSTRUCT_POP.bnum==i & DATSTRUCT_POP.enum==ii & ...
                    DATSTRUCT_POP.switch==iii & strcmp(DATSTRUCT_POP.bregion, bregionToSplitBy) ...
                    & DATSTRUCT_POP.motifnum == mm);
                if isempty(indstmp)
                    disp('SKIP (no pop data)');
                    continue
                end
                
                
                %% ============ get data from OUTSTRUCT_XCOV
                indsOUT = find(OUTSTRUCT_XCOV.bnum==i & OUTSTRUCT_XCOV.enum==ii & OUTSTRUCT_XCOV.switch==iii ...
                    & OUTSTRUCT_XCOV.motifnum==mm);
                indsOUT = indsOUT(1);
                
                
                inds_base = OUTSTRUCT_XCOV.inds_base_epoch{indsOUT};
                inds_WN = OUTSTRUCT_XCOV.inds_WN_epoch{indsOUT};
                inds_base_all = OUTSTRUCT_XCOV.inds_base_allgood{indsOUT};
                inds_WN_all = OUTSTRUCT_XCOV.inds_WN_allgood{indsOUT};
                trialedges = OUTSTRUCT_XCOV.trialedges_epoch{indsOUT};
                
                
                istarg = OUTSTRUCT_XCOV.istarg(indsOUT);
                issame = OUTSTRUCT_XCOV.issame(indsOUT);
                learndir = OUTSTRUCT_XCOV.learndirTarg(indsOUT);
                
                %% ================= COLLEC THINGS ACROSS TRILAS (E.G. FF)
                FFvals = SwitchCohStruct.bird(i).exptnum(ii).switchlist(iii).motifnum(mm).ffvals;
                tvals = SwitchCohStruct.bird(i).exptnum(ii).switchlist(iii).motifnum(mm).tvals;
                
                
                %% ================== things about this syl
                
                
                % =============== GET BASELINE AND WN FOR ALL PAIRS
                npairs = length(datthis.ccRealAllPair);
                for np=1:npairs
                    datmat_real = datthis.ccRealAllPair{np};
                    datmat_shuff = datthis.ccShiftAllPair{np};
                    neurpair = datthis.neurPair(np, :);
                    
                    if removebadchans==1
                        chanbad = lt_neural_QUICK_RemoveBadChans(bname, ename, iii, datthis.chanPair(np,:));
                        if chanbad==1
                            continue
                        end
                    end
                    
                    
                    % ===== get autocovariance too
                    if strcmp(xcovver, 'coherency')
                        datcell_auto_real = datthis.ccRealAuto(np,:);
                        datcell_auto_shift = datthis.ccShiftAuto(np,:);
                    else
                        datcell_auto_real = {};
                        datcell_auto_shift = {};
                    end
                    
                    
                    % ==== process (e.g. smooth) and get xcov output
                    [datbase, datWN, datWN_notminshuff, datbase_notminshuff, Xq, NanCount, ...
                        dattrials] = ...
                        lt_neural_POPLEARN_XCov_sub1(datmat_real, datmat_shuff, dosmooth, ...
                        dosmooth_sigma, inds_base, inds_WN, PARAMS.Xcov_ccLags, plotraw, ...
                        xcovver, datcell_auto_real, datcell_auto_shift);
                    
                    
                    %% ===== get xcov-gram... [BASE AND WN]
                    if getxgram==1
                        % ====== BASE AND WN
                        [xcovgram_base, xcovgram_wn] = lt_neural_POPLEARN_XcovExtr_sub1(datthis, ...
                            dosmooth, dosmooth_sigma, inds_base, inds_WN, PARAMS, xcovver, datbase, ...
                            np);
                        
                    else
                        nwinds = size(datthis.ccRealAllPair_allwind,2);
                        xcovgram_base = nan(nwinds, length(datbase)); % win x lags
                        xcovgram_wn= nan(nwinds, length(datbase));
                    end
                    
                    %% ==== SPLIT BY ENSEMBLE COORD
                    % ------- get trial by trial values for ensemble coord
                    indsPOP = find(DATSTRUCT_POP.bnum==i & DATSTRUCT_POP.enum==ii & ...
                        DATSTRUCT_POP.switch==iii & strcmp(DATSTRUCT_POP.bregion, bregionToSplitBy) ...
                        & DATSTRUCT_POP.motifnum == mm);
                    assert(length(indsPOP)==1);
                    
                    rho_allpair = DATSTRUCT_POP.fr_RhoPairwise{indsPOP};
                    assert(length(rho_allpair) == length(FFvals));
                    
                    xcovgram_base_ffsplits = cell(1,2);
                    
                    if getHiLoFFSplit==1
                        rhothis = rho_allpair(inds_base);
                        rhomid = median(rhothis);
                        
                        %                        % ========= THEN get, for baseline, high and low FF.
                        %                        ffthis = FFvals(inds_base);
                        %                        ffmid = median(ffthis);
                        %
                        % -------------- HI CORR BASELINE TRIALS
                        inds_base_hi = inds_base(rhothis>rhomid);
                        [~, xcovgram_wn] = lt_neural_POPLEARN_XcovExtr_sub1(datthis, ...
                            dosmooth, dosmooth_sigma, inds_base, inds_base_hi, PARAMS, xcovver, datbase, ...
                            np, hilosplit_shuffver);
                        
                        xcovgram_base_ffsplits{2} = xcovgram_wn;
                        
                        % -------- LO PITCH
                        inds_base_lo = inds_base(rhothis<rhomid);
                        [~, xcovgram_wn] = lt_neural_POPLEARN_XcovExtr_sub1(datthis, ...
                            dosmooth, dosmooth_sigma, inds_base, inds_base_lo, PARAMS, xcovver, datbase, ...
                            np, hilosplit_shuffver);
                        
                        xcovgram_base_ffsplits{1} = xcovgram_wn;
                        
                    end
                    
                    
                    %                     xcenters = mean(windlist,2);
                    % ============= COUNT How many trials had nan - i
                    %                     NanCountAll = [NanCountAll; NanCount];
                    
                    
                    %% ======= get xcov in different time bins during learning
                    
                    N = getxgram_epochbins;
                    XcovgramWN_epochs = nan([size(xcovgram_base) N]);
                    
                    XcovgramWN_epochs_hiFF = nan([size(xcovgram_base) N]);
                    XcovgramWN_epochs_loFF = nan([size(xcovgram_base) N]);
                    
                    % ============== divide up WN into different bins
                    %                     binsize = length(inds_WN_all)/N;
                    
                    trialedges2 = inds_WN_all(round(linspace(1, length(inds_WN_all), N+1)));
                    %                         trialedges = inds_WN_all(trialedges);
                    trialedges2(end)=trialedges2(end)+1;
                    assert(all(trialedges==trialedges2));
                    
                    % ====== go thru all bins
                    for bb=1:N
                        
                        trialsthis = trialedges(bb):(trialedges(bb+1)-1);
                        
                        [~, xcovgram_wn] = lt_neural_POPLEARN_XcovExtr_sub1(datthis, ...
                            dosmooth, dosmooth_sigma, inds_base, trialsthis, PARAMS, xcovver, datbase, ...
                            np);
                        
                        XcovgramWN_epochs(:,:, bb) = xcovgram_wn;
                        
                        
                        % ########################## GET HI LO SPLIT?
                        rhothis = rho_allpair(trialsthis);
                        rhomid = median(rhothis);
                        
                        % ============ HGIH PITCH TRIALS
                        trialsthis_hi = trialsthis(rhothis>rhomid);
                        [~, xcovgram_hi] = lt_neural_POPLEARN_XcovExtr_sub1(datthis, ...
                            dosmooth, dosmooth_sigma, inds_base, trialsthis_hi, PARAMS, xcovver, datbase, ...
                            np, hilosplit_shuffver);
                        
                        XcovgramWN_epochs_hiFF(:,:, bb) = xcovgram_hi;
                        
                        % ============ LOW PITCH TRIALS
                        trialsthis_lo = trialsthis(rhothis<rhomid);
                        [~, xcovgram_lo] = lt_neural_POPLEARN_XcovExtr_sub1(datthis, ...
                            dosmooth, dosmooth_sigma, inds_base, trialsthis_lo, PARAMS, xcovver, datbase, ...
                            np, hilosplit_shuffver);
                        
                        XcovgramWN_epochs_loFF(:,:, bb) = xcovgram_lo;
                        
                    end
                    
                    % ======== save ff split into cell
                    XcovgramWN_epochs_FFsplits = cell(1,2);
                    XcovgramWN_epochs_FFsplits{1} = XcovgramWN_epochs_loFF;
                    XcovgramWN_epochs_FFsplits{2} = XcovgramWN_epochs_hiFF;
                    
                    %%
                    
                    OUTSTRUCT_MeanRhoSplit.XcovgramBase = [OUTSTRUCT_MeanRhoSplit.XcovgramBase; single(xcovgram_base)];
                    OUTSTRUCT_MeanRhoSplit.XcovgramWN = [OUTSTRUCT_MeanRhoSplit.XcovgramWN; single(xcovgram_wn)];
                    
                    OUTSTRUCT_MeanRhoSplit.XcovBase = [OUTSTRUCT_MeanRhoSplit.XcovBase; datbase];
                    OUTSTRUCT_MeanRhoSplit.XcovWN = [OUTSTRUCT_MeanRhoSplit.XcovWN; datWN];
                    
                    OUTSTRUCT_MeanRhoSplit.Xcov_DotProd_trials = [OUTSTRUCT_MeanRhoSplit.Xcov_DotProd_trials; dattrials];
                    
                    OUTSTRUCT_MeanRhoSplit.XcovBase_NoMinShuff = ...
                        [OUTSTRUCT_MeanRhoSplit.XcovBase_NoMinShuff; datbase_notminshuff];
                    OUTSTRUCT_MeanRhoSplit.XcovWN_NoMinShuff = ...
                        [OUTSTRUCT_MeanRhoSplit.XcovWN_NoMinShuff; datWN_notminshuff];
                    
                    OUTSTRUCT_MeanRhoSplit.XcovgramWN_epochs = [OUTSTRUCT_MeanRhoSplit.XcovgramWN_epochs; single(XcovgramWN_epochs)];
                    
                    OUTSTRUCT_MeanRhoSplit.neurpairnum = [OUTSTRUCT_MeanRhoSplit.neurpairnum; np];
                    OUTSTRUCT_MeanRhoSplit.neurpair = [OUTSTRUCT_MeanRhoSplit.neurpair; neurpair];
                    OUTSTRUCT_MeanRhoSplit.bnum = [OUTSTRUCT_MeanRhoSplit.bnum; i];
                    OUTSTRUCT_MeanRhoSplit.enum = [OUTSTRUCT_MeanRhoSplit.enum; ii];
                    OUTSTRUCT_MeanRhoSplit.swnum = [OUTSTRUCT_MeanRhoSplit.swnum; iii];
                    OUTSTRUCT_MeanRhoSplit.motifnum = [OUTSTRUCT_MeanRhoSplit.motifnum; mm];
                    OUTSTRUCT_MeanRhoSplit.issame = [OUTSTRUCT_MeanRhoSplit.issame; issame];
                    OUTSTRUCT_MeanRhoSplit.istarg = [OUTSTRUCT_MeanRhoSplit.istarg; istarg];
                    OUTSTRUCT_MeanRhoSplit.learndirTarg = [OUTSTRUCT_MeanRhoSplit.learndirTarg; learndir];
                    
                    OUTSTRUCT_MeanRhoSplit.inds_base_epoch = [OUTSTRUCT_MeanRhoSplit.inds_base_epoch; inds_base]; % used for analysis
                    OUTSTRUCT_MeanRhoSplit.inds_WN_epoch = [OUTSTRUCT_MeanRhoSplit.inds_WN_epoch; inds_WN]; % used for analysis
                    OUTSTRUCT_MeanRhoSplit.inds_base_allgood = [OUTSTRUCT_MeanRhoSplit.inds_base_allgood; inds_base_all]; % all, after remove bad trials
                    OUTSTRUCT_MeanRhoSplit.inds_WN_allgood = [OUTSTRUCT_MeanRhoSplit.inds_WN_allgood; inds_WN_all]; % all, after remove bad trials
                    OUTSTRUCT_MeanRhoSplit.trialedges_epoch = [OUTSTRUCT_MeanRhoSplit.trialedges_epoch; trialedges]; % all, after remove bad trials
                    
                    
                    
                    % ============== SAVE HI LO
                    if getHiLoFFSplit==1
                        OUTSTRUCT_MeanRhoSplit.XcovgramWN_FFsplits_Base = ...
                            [OUTSTRUCT_MeanRhoSplit.XcovgramWN_FFsplits_Base; xcovgram_base_ffsplits];
                        
                        OUTSTRUCT_MeanRhoSplit.XcovgramWN_FFsplits_Epochs = ...
                            [OUTSTRUCT_MeanRhoSplit.XcovgramWN_FFsplits_Epochs; XcovgramWN_epochs_FFsplits];
                    end
                end
            end
        end
    end
end

if dosmooth==1
    PARAMS.Xcov_ccLags_beforesmooth = PARAMS.Xcov_ccLags;
    PARAMS.Xcov_ccLags = Xq;
else
    PARAMS.Xcov_ccLags = datthis.ccLags;
end
% PARAMS.xcenters_gram = xcenters;