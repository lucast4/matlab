function [OUTSTRUCT_XCOV, PARAMS, NanCountAll] = lt_neural_POPLEARN_XcovExtr(SwitchXCovStruct, ...
    SwitchStruct, PARAMS, SwitchCohStruct, OUTSTRUCT, usealltrials, ...
    useallbase, dosmooth, dosmooth_sigma, removebadsyl, ...
    windlist, plotrawtrials, xcovver, wntouse, removebadtrials, getxgram, removebadchans, ...
    getxgram_epochbins, getHiLoFFSplit)
%% NOTE: if any trials have nan, then just ignroes in calcualting dat minus shuff
% i.e. uses nanmean. I do not expect there to be any nan for xcorr
% ("unbiased");
%% lt 2/2019 - extracts OUTSTRUCT_XCOV for xcov analysis during learning.

%% % ===================================
OUTSTRUCT_XCOV = struct;

OUTSTRUCT_XCOV.XcovBase = [];
OUTSTRUCT_XCOV.XcovWN = [];

OUTSTRUCT_XCOV.XcovgramBase = {};
OUTSTRUCT_XCOV.XcovgramWN = {};

OUTSTRUCT_XCOV.XcovgramWN_epochs = {};

OUTSTRUCT_XCOV.neurpairnum = [];
OUTSTRUCT_XCOV.neurpair = [];
OUTSTRUCT_XCOV.bnum = [];
OUTSTRUCT_XCOV.enum =[];
OUTSTRUCT_XCOV.swnum = [];
OUTSTRUCT_XCOV.motifnum = [];
OUTSTRUCT_XCOV.issame = [];
OUTSTRUCT_XCOV.istarg = [];
OUTSTRUCT_XCOV.learndirTarg = [];

OUTSTRUCT_XCOV.XcovBase_NoMinShuff = ...
    [];
OUTSTRUCT_XCOV.XcovWN_NoMinShuff = ...
    [];

OUTSTRUCT_XCOV.inds_base_epoch = {}; % used for analysis
OUTSTRUCT_XCOV.inds_WN_epoch = {}; % used for analysis
OUTSTRUCT_XCOV.inds_base_allgood = {}; % all, after remove bad trials
OUTSTRUCT_XCOV.inds_WN_allgood = {}; % all, after remove bad trials
OUTSTRUCT_XCOV.trialedges_epoch = {};


OUTSTRUCT_XCOV.XcovgramWN_FFsplits_Base = {};
OUTSTRUCT_XCOV.XcovgramWN_FFsplits_Epochs = {};

OUTSTRUCT_XCOV.Xcov_DotProd_trials = {};

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
                
                % ========== CHECK IF IS BAD SYL?
                if removebadsyl==1
                    sylbad = lt_neural_QUICK_LearnRemoveBadSyl(bname, ename, iii, motifthis);
                    if sylbad==1
                        continue
                    end
                end
                datthis = SwitchXCovStruct.bird(i).exptnum(ii).switchlist(iii).motif(mm);
                %                 if usealltrials==1
                % %                     inds_base = SwitchCohStruct.bird(i).exptnum(ii).switchlist(iii).motifnum(mm).indsbase;
                % %                     inds_WN = SwitchCohStruct.bird(i).exptnum(ii).switchlist(iii).motifnum(mm).indsWN;
                %                 elseif usealltrials==0
                %                     inds_base = SwitchCohStruct.bird(i).exptnum(ii).switchlist(iii).motifnum(mm).indsbase_epoch;
                %                     inds_WN = SwitchCohStruct.bird(i).exptnum(ii).switchlist(iii).motifnum(mm).indsWN_epoch;
                %                 end
                %
                %                 if useallbase==1
                % %                     inds_base = SwitchCohStruct.bird(i).exptnum(ii).switchlist(iii).motifnum(mm).indsbase;
                %                 end
                
                %                 if useallwn==1
                %                     inds_WN = SwitchCohStruct.bird(i).exptnum(ii).switchlist(iii).motifnum(mm).indsWN;
                %                 end
                
                %                 if usehalfwn==1
                %                     inds_WN = SwitchCohStruct.bird(i).exptnum(ii).switchlist(iii).motifnum(mm).indsWN_epoch;
                %                     inds_WN = inds_WN(round(length(inds_WN)/2):end);
                %                 end
                
                
                %% SAVE 2 SETS OF INDS:
                % 1) all base and WN inds (not removed)
                inds_base_all = SwitchCohStruct.bird(i).exptnum(ii).switchlist(iii).motifnum(mm).indsbase;
                inds_WN_all = SwitchCohStruct.bird(i).exptnum(ii).switchlist(iii).motifnum(mm).indsWN;
                tvals = SwitchCohStruct.bird(i).exptnum(ii).switchlist(iii).motifnum(mm).tvals;
                
                % ---- filter out the good trials
                if removebadtrials==1
                    badtrials = lt_neural_QUICK_RemoveTrials(bname, ename, iii, tvals(inds_base_all));
                    inds_base_all(badtrials) = [];
                    
                    badtrials = lt_neural_QUICK_RemoveTrials(bname, ename, iii, tvals(inds_WN_all));
                    inds_WN_all(badtrials) = [];
                end
                
                if isempty(inds_base_all) | isempty(inds_WN_all)
                    continue
                end
                
                % 2) epoch base and WN inds
                inds_base = [];
                inds_WN = [];
                if usealltrials==1
                    inds_base = inds_base_all;
                    inds_WN = inds_WN_all;
                elseif usealltrials ==0
                    % ---- BASE
                    if useallbase==1
                        inds_base = inds_base_all;
                    elseif useallbase==0
                        inds_base = inds_base_all(end-round(length(inds_base_all)/2):end);
                    end
                    
                    % ----- WN
                    if strcmp(wntouse, 'all')
                        inds_WN = inds_WN_all;
                    elseif strcmp(wntouse, 'half')
                        inds_WN = inds_WN_all(end-round(length(inds_WN_all)/2):end);
                    elseif strcmp(wntouse, 'third')
                        inds_WN = inds_WN_all(end-round(length(inds_WN_all)/3):end);
                    elseif strcmp(wntouse, 'quarter')
                        inds_WN = inds_WN_all(end-round(length(inds_WN_all)/4):end);
                    end
                end
                if isempty(inds_base) | isempty(inds_WN)
                    disp('WHY EMPTYY')
                    keyboard
                end
                
                % ======== AD HOC, one experiment only has epoch of good
                % trials (at endpoint of leajrning) so always use all those
                % trials...
                if strcmp(bname, 'pu69wh78') & strcmp(ename, 'RALMANlearn2') ...
                        & iii==1
                    inds_WN = inds_WN_all;
                end
                
                
                %% ================= COLLEC THINGS ACROSS TRILAS (E.G. FF)
                FFvals = SwitchCohStruct.bird(i).exptnum(ii).switchlist(iii).motifnum(mm).ffvals;
                
                
                %% ================== things about this syl
                %                 motifname = SwitchCohStruct.bird(i).exptnum(ii).switchlist(iii).motifnum(mm);
                indsOUT = OUTSTRUCT.bnum==i & OUTSTRUCT.enum==ii & OUTSTRUCT.switch==iii & ...
                    OUTSTRUCT.motifnum==mm;
                if ~any(indsOUT)
                    % skip, since not part of oringial outsrtuct. (e.g., is
                    % bad syl
                    continue
                end
                istarg = unique(OUTSTRUCT.istarg(indsOUT));
                issame = unique(OUTSTRUCT.issame(indsOUT));
                learndir = unique(OUTSTRUCT.learndirTarg(indsOUT));
                
                
                
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
                    
                    plotraw=0;
                    if istarg==1 & iii==1 & i>3 & np<3 & plotrawtrials==1
                        disp([bname '-' ename '-sw' num2str(iii)]);
                        plotraw = 1;
                        keyboard
                    else
                        plotraw=0;
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
                    
                    xcovgram_base_ffsplits = cell(1,2);
                    if getHiLoFFSplit==1
                       % ========= THEN get, for baseline, high and low FF.
                       ffthis = FFvals(inds_base);
                       ffmid = median(ffthis);
                       
                       % -------------- HI PITCH BASELINE TRIALS
                       inds_base_hi = inds_base(ffthis>ffmid);
                       [~, xcovgram_wn] = lt_neural_POPLEARN_XcovExtr_sub1(datthis, ...
                           dosmooth, dosmooth_sigma, inds_base, inds_base_hi, PARAMS, xcovver, datbase, ...
                           np);
                       
                       xcovgram_base_ffsplits{2} = xcovgram_wn;
                       
                       % -------- LO PITCH
                       inds_base_lo = inds_base(ffthis<ffmid);
                       [~, xcovgram_wn] = lt_neural_POPLEARN_XcovExtr_sub1(datthis, ...
                           dosmooth, dosmooth_sigma, inds_base, inds_base_lo, PARAMS, xcovver, datbase, ...
                           np);
                       
                       xcovgram_base_ffsplits{1} = xcovgram_wn;
                        
                    end
                    
                    
                    xcenters = mean(windlist,2);
                    % ============= COUNT How many trials had nan - i
                    NanCountAll = [NanCountAll; NanCount];
                    
                    
                    %% ======= get xcov in different time bins during learning
                    if ~isempty(getxgram_epochbins)
                        N = getxgram_epochbins;
                        XcovgramWN_epochs = nan([size(xcovgram_base) N]);
                        
                        XcovgramWN_epochs_hiFF = nan([size(xcovgram_base) N]);
                        XcovgramWN_epochs_loFF = nan([size(xcovgram_base) N]);
                        
                        % ============== divide up WN into different bins
                        binsize = length(inds_WN_all)/N;
                        
                        trialedges = inds_WN_all(round(linspace(1, length(inds_WN_all), N+1)));
                        %                         trialedges = inds_WN_all(trialedges);
                        trialedges(end)=trialedges(end)+1;
                        
                        % ====== go thru all bins
                        for bb=1:N
                            
                            trialsthis = trialedges(bb):(trialedges(bb+1)-1);
                            
                            [~, xcovgram_wn] = lt_neural_POPLEARN_XcovExtr_sub1(datthis, ...
                                dosmooth, dosmooth_sigma, inds_base, trialsthis, PARAMS, xcovver, datbase, ...
                                np);
                            
                            XcovgramWN_epochs(:,:, bb) = xcovgram_wn;
                            
                            
                            % ########################## GET HI LO SPLIT?
                            if getHiLoFFSplit ==1
                                
                                ffthis = FFvals(trialsthis);
                                midff = median(ffthis);
                                
                                % ============ HGIH PITCH TRIALS
                                trialsthis_hi = trialsthis(ffthis>midff);
                                [~, xcovgram_hi] = lt_neural_POPLEARN_XcovExtr_sub1(datthis, ...
                                    dosmooth, dosmooth_sigma, inds_base, trialsthis_hi, PARAMS, xcovver, datbase, ...
                                    np);
                                
                                XcovgramWN_epochs_hiFF(:,:, bb) = xcovgram_hi;
                            
                                % ============ LOW PITCH TRIALS
                                trialsthis_lo = trialsthis(ffthis<midff);
                                 [~, xcovgram_lo] = lt_neural_POPLEARN_XcovExtr_sub1(datthis, ...
                                    dosmooth, dosmooth_sigma, inds_base, trialsthis_lo, PARAMS, xcovver, datbase, ...
                                    np);
                                
                                XcovgramWN_epochs_loFF(:,:, bb) = xcovgram_lo;
                                
                           end
                        end
                        
                        % ======== save ff split into cell
                        XcovgramWN_epochs_FFsplits = cell(1,2);
                        XcovgramWN_epochs_FFsplits{1} = XcovgramWN_epochs_loFF;
                        XcovgramWN_epochs_FFsplits{2} = XcovgramWN_epochs_hiFF;
                        
                    end
                    
                    %%
                    
                    OUTSTRUCT_XCOV.XcovgramBase = [OUTSTRUCT_XCOV.XcovgramBase; single(xcovgram_base)];
                    OUTSTRUCT_XCOV.XcovgramWN = [OUTSTRUCT_XCOV.XcovgramWN; single(xcovgram_wn)];
                    
                    OUTSTRUCT_XCOV.XcovBase = [OUTSTRUCT_XCOV.XcovBase; datbase];
                    OUTSTRUCT_XCOV.XcovWN = [OUTSTRUCT_XCOV.XcovWN; datWN];
                    
                    OUTSTRUCT_XCOV.Xcov_DotProd_trials = [OUTSTRUCT_XCOV.Xcov_DotProd_trials; dattrials];
                    
                    OUTSTRUCT_XCOV.XcovBase_NoMinShuff = ...
                        [OUTSTRUCT_XCOV.XcovBase_NoMinShuff; datbase_notminshuff];
                    OUTSTRUCT_XCOV.XcovWN_NoMinShuff = ...
                        [OUTSTRUCT_XCOV.XcovWN_NoMinShuff; datWN_notminshuff];
                    
                    OUTSTRUCT_XCOV.XcovgramWN_epochs = [OUTSTRUCT_XCOV.XcovgramWN_epochs; single(XcovgramWN_epochs)];
                    
                    OUTSTRUCT_XCOV.neurpairnum = [OUTSTRUCT_XCOV.neurpairnum; np];
                    OUTSTRUCT_XCOV.neurpair = [OUTSTRUCT_XCOV.neurpair; neurpair];
                    OUTSTRUCT_XCOV.bnum = [OUTSTRUCT_XCOV.bnum; i];
                    OUTSTRUCT_XCOV.enum = [OUTSTRUCT_XCOV.enum; ii];
                    OUTSTRUCT_XCOV.swnum = [OUTSTRUCT_XCOV.swnum; iii];
                    OUTSTRUCT_XCOV.motifnum = [OUTSTRUCT_XCOV.motifnum; mm];
                    OUTSTRUCT_XCOV.issame = [OUTSTRUCT_XCOV.issame; issame];
                    OUTSTRUCT_XCOV.istarg = [OUTSTRUCT_XCOV.istarg; istarg];
                    OUTSTRUCT_XCOV.learndirTarg = [OUTSTRUCT_XCOV.learndirTarg; learndir];
                    
                    OUTSTRUCT_XCOV.inds_base_epoch = [OUTSTRUCT_XCOV.inds_base_epoch; inds_base]; % used for analysis
                    OUTSTRUCT_XCOV.inds_WN_epoch = [OUTSTRUCT_XCOV.inds_WN_epoch; inds_WN]; % used for analysis
                    OUTSTRUCT_XCOV.inds_base_allgood = [OUTSTRUCT_XCOV.inds_base_allgood; inds_base_all]; % all, after remove bad trials
                    OUTSTRUCT_XCOV.inds_WN_allgood = [OUTSTRUCT_XCOV.inds_WN_allgood; inds_WN_all]; % all, after remove bad trials
                    OUTSTRUCT_XCOV.trialedges_epoch = [OUTSTRUCT_XCOV.trialedges_epoch; trialedges]; % all, after remove bad trials
                    
                    
                    
                    % ============== SAVE HI LO
                    if getHiLoFFSplit==1
                    OUTSTRUCT_XCOV.XcovgramWN_FFsplits_Base = ...
                        [OUTSTRUCT_XCOV.XcovgramWN_FFsplits_Base; xcovgram_base_ffsplits];
                    
                    OUTSTRUCT_XCOV.XcovgramWN_FFsplits_Epochs = ...
                        [OUTSTRUCT_XCOV.XcovgramWN_FFsplits_Epochs; XcovgramWN_epochs_FFsplits];
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
PARAMS.xcenters_gram = xcenters;