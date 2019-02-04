function OutStruct = lt_neural_v2_ANALY_FRsmooth_PredLearn(OUTDAT, MOTIFSTATS_Compiled, ...
    SwitchStruct, corrwindow, syltypesneeded, epochtoplot, analytoplot, ...
    syltoplot, docorrvsdiff, onlyPlotIfAllTargSameDir, onlyIfLearnCorrectDir, ...
    doregression, plotRaw, plotSummary)

%% lt - based on baseline relationship between smoothed FR and FF, predict change in neural?
% 9/12/18


% prctile_divs = [33 66 100]; % percentiles to divide up data by
% epochtoplot = 3; % i.e. out of the epochs decided by prctile_divs
% prctile_divs = [50 100]; % percentiles to divide up data by
% epochtoplot = 2; % i.e. out of the epochs decided by prctile_divs

% usepercent = 0; % if 1, then gets fr percent deviation from baseline. otherwise uses hz diff
% nbasetime = 60; % 60 minutes bnefore first WN trial is min time for baseline
% nbasetime = []; % 60 minutes bnefore first WN trial is min time for baseline

% analytype = 'AllMinusAll_FRsmooth';
% % AllMinusAll_FRsmooth
% % AllOnlyMinusDiff_FRsmooth
% % AllMinusAllMinusDiff_FRsmooth
% % AllOnlyMinusBase_FRsmooth
%
%
% syltypesneeded = [1 1 1] means needs minimum 1 targ, 1 same, 1 diff.
% [1 0 1] means doesnt care if has same type

% docorrvsdiff=1, means first get diff in fr, then corr vs that diff
% if 0, then do corr separately vs FR, then subtract those corrs.


%% ###################### LIMIT TO NEURONS THAT CONTAIN ALL SYL TYPES?

% =-==== go thru all switches. if bad then throw ou
maxneur = max(OUTDAT.All_neurnum);
numbirds = length(SwitchStruct.bird);

% ======================== SHUFFLE SYL TYPE? % within each neuron
indstokeep = [];
for i=1:numbirds
    numexpts = length(SwitchStruct.bird(i).exptnum);
    for ii=1:numexpts
        numswitch = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        for ss = 1:numswitch
            % ################### SECOND, PLOT SMOOTH FR FOR ALL SYLS, SUBTRACT BASE
            for nn=1:maxneur
                
                indsthis = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss ...
                    & OUTDAT.All_neurnum==nn);
                
                if isempty(indsthis)
                    continue
                end
                
                % ==== count how many of each syl type there exists
                numtarg = sum(OUTDAT.All_istarg(indsthis)==1);
                numsame = sum(OUTDAT.All_istarg(indsthis)==0 & OUTDAT.All_issame(indsthis)==1);
                numdiff = sum(OUTDAT.All_istarg(indsthis)==0 & OUTDAT.All_issame(indsthis)==0);
                
                if all([numtarg numsame numdiff] >= syltypesneeded)
                    % then keep
                    disp([numtarg numsame numdiff]);
                    disp(indsthis)
                    indstokeep = [indstokeep; indsthis];
                end
                
            end
        end
    end
end

disp(['Keeping ' num2str(length(indstokeep)) '/' num2str(length(OUTDAT.All_birdnum)) ' datapoints, passes syltypes required criterion']);

OUTDAT = lt_structure_subsample_all_fields(OUTDAT, indstokeep, 1);

if plotRaw==1
    %% ============= FOR EACH SWITCH, PLOT ALL AS FUNCTION OF LEARN MAG
    %% ############## go thru all switches, for each neuron COLLECT THINGS ABOUT PREDICT LEARNING
    
    % analytoplot = 'AllDevDiff_NotAbs';
    % syltoplot = 'targ';
    
    for i=1:numbirds
        numexpts = length(SwitchStruct.bird(i).exptnum);
        
        figcount=1;
        subplotrows=5;
        subplotcols=4;
        fignums_alreadyused=[];
        hfigs=[];
        hsplots = [];
        
        
        for ii=1:numexpts
            numswitch = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
            for ss = 1:numswitch
                
                if strcmp(syltoplot, 'targ')
                    indsthis = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & ...
                        OUTDAT.All_swnum==ss & OUTDAT.All_istarg==1);
                elseif strcmp(syltoplot, 'same')
                    indsthis = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & ...
                        OUTDAT.All_swnum==ss & OUTDAT.All_istarg==0 & OUTDAT.All_issame==1);
                elseif strcmp(syltoplot, 'diff')
                    indsthis = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & ...
                        OUTDAT.All_swnum==ss & OUTDAT.All_istarg==0 & OUTDAT.All_issame==0);
                end
                
                if isempty(indsthis)
                    continue
                end
                
                
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
                bname = SwitchStruct.bird(i).birdname;
                ename = SwitchStruct.bird(i).exptnum(ii).exptname;
                title([bname '-' ename '-sw' num2str(ss)]);
                ylabel('corr (frdiff) vs (basehi - baselo) [r = learn up]');
                xlabel('pitch from base (z)');
                
                
                % ########################## CALCUALTE THINGS
                for j=indsthis'
                    
                    % ============== skip if bad syl
                    
%                 if removeBadSyls==1
%                     sylbad = lt_neural_QUICK_LearnRemoveBadSyl(bthis, ename, ss, motifthis);
%                     if sylbad==1
%                         continue
%                     end
%                 end
                    
                    % ----- get time window
                    t = OUTDAT.AllMinusBase_tbinAll{j};
                    indstokeep_t = t>=corrwindow(1) & t<corrwindow(2);
                    t = t(indstokeep_t);
                    
                    % ---- SMOOTHED FR, baseline
                    ybase_lo = OUTDAT.AllBase_FRsmooth_Lo{j}(indstokeep_t);
                    ybase_hi = OUTDAT.AllBase_FRsmooth_Hi{j}(indstokeep_t);
                    ybase_himinuslo = ybase_hi - ybase_lo;
                    
                    if all(isnan(ybase_himinuslo))
                        assert(all(isnan(OUTDAT.All_FF{j,1})));
                        continue
                    end
                    
                    % ---- SMOOTHED FR, ACOUNT FOR DRIFT
                    if strcmp(analytoplot, 'AllMinusBase_FRmeanAll')
                        ydev = OUTDAT.AllMinusBase_FRmeanAll{j}{epochtoplot}(indstokeep_t);
                    else
                        ydev = OUTDAT.(analytoplot){j}(indstokeep_t);
                    end
                    
                    % ------ PITCH, DEVIATION FROM BASELINE
                    pitchdiffZ = OUTDAT.AllMinusBase_PitchZ(j, epochtoplot);
                    
                    
                    % ================== CALCULATE THINGS
                    frcorr_hi = corr(ybase_hi, ydev);
                    frcorr_lo = corr(ybase_lo, ydev);
                    
                    if docorrvsdiff ==1 % 2 methods for getting difference in correlations
                        frcorr_himinuslo = corr(ybase_himinuslo, ydev);
                    elseif docorrvsdiff==0
                        % then get corr vs. hi and lo, and then subtract
                        frcorr_himinuslo = frcorr_hi - frcorr_lo;
                    end
                    
                    
                    learndir = unique([SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningDirs{2:2:end}]);
                    if length(learndir)>1
                        if onlyPlotIfAllTargSameDir==1
                            disp(learndir)
                            continue
                        end
                        learndir = nan;
                    end
                    istarg = OUTDAT.All_istarg(j);
                    issame = OUTDAT.All_issame(j);
                    
                    %                 if i==1 & ii==1 & ss==7 & istarg==1
                    %                     keyboard
                    %                 end
                    %
                    
                    % ########################### PLOT
                    if learndir==1
                        pcol = 'r';
                    elseif learndir==-1
                        pcol = 'b';
                    end
                    
                    % ----
                    plot(pitchdiffZ, frcorr_himinuslo, 'o', 'Color', pcol);
                    %                 plot(ss, frcorr_himinuslo, 'o', 'Color', pcol);
                    
                    % --- note down if this is a switch that turns off WN
                    functmp = @(x)x(2);
                    postConts = cellfun(functmp, ...
                        SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningContingencies(2:2:end));
                    turnsoff = any(postConts==0);
                    tmp = [SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningContingencies{2:2:end}];
                    lt_plot_text(0, 0.9, num2str(tmp), 'k')
                    
                end
                axis tight;
                lt_plot_zeroline;
                lt_plot_zeroline_vert;
                ylim([-1 1]);
                
            end
            linkaxes(hsplots, 'xy');
        end
    end
    
    
    %% ############## go thru all switches, for each neuron COLLECT THINGS ABOUT PREDICT LEARNING
    
    % analytoplot = 'AllDevDiff_NotAbs';
    % syltoplot = 'targ';
    
    LearnDirAll = [];
    CorrAll = [];
    BirdnumAll = [];
    EnumAll = [];
    EcounterAll = [];
    PitchDiffZ = [];
    ecount =1 ;
    for i=1:numbirds
        numexpts = length(SwitchStruct.bird(i).exptnum);
        
        figcount=1;
        subplotrows=2;
        subplotcols=3;
        fignums_alreadyused=[];
        hfigs=[];
        hsplots = [];
        
        
        for ii=1:numexpts
            numswitch = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            bname = SwitchStruct.bird(i).birdname;
            ename = SwitchStruct.bird(i).exptnum(ii).exptname;
            title([bname '-' ename]);
            ylabel('corr (frdiff) vs (basehi - baselo) [r = learn up]');
            xlabel('swnum');
            
            for ss = 1:numswitch
                
                if strcmp(syltoplot, 'targ')
                    indsthis = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & ...
                        OUTDAT.All_swnum==ss & OUTDAT.All_istarg==1);
                elseif strcmp(syltoplot, 'same')
                    indsthis = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & ...
                        OUTDAT.All_swnum==ss & OUTDAT.All_istarg==0 & OUTDAT.All_issame==1);
                elseif strcmp(syltoplot, 'diff')
                    indsthis = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & ...
                        OUTDAT.All_swnum==ss & OUTDAT.All_istarg==0 & OUTDAT.All_issame==0);
                end
                
                if isempty(indsthis)
                    continue
                end
                
                % ########################## CALCUALTE THINGS
                for j=indsthis'
                    
                    % ----- get time window
                    t = OUTDAT.AllMinusBase_tbinAll{j};
                    indstokeep_t = t>=corrwindow(1) & t<corrwindow(2);
                    t = t(indstokeep_t);
                    
                    % ---- SMOOTHED FR, baseline
                    ybase_lo = OUTDAT.AllBase_FRsmooth_Lo{j}(indstokeep_t);
                    ybase_hi = OUTDAT.AllBase_FRsmooth_Hi{j}(indstokeep_t);
                    ybase_himinuslo = ybase_hi - ybase_lo;
                    if all(isnan(ybase_himinuslo))
                        assert(all(isnan(OUTDAT.All_FF{j})));
                        continue
                    end
                    
                    % ---- SMOOTHED FR, ACOUNT FOR DRIFT
                    if strcmp(analytoplot, 'AllMinusBase_FRmeanAll')
                        ydev = OUTDAT.AllMinusBase_FRmeanAll{j}{epochtoplot}(indstokeep_t);
                    else
                        ydev = OUTDAT.(analytoplot){j}(indstokeep_t);
                    end
                    
                    % ------ PITCH, DEVIATION FROM BASELINE
                    pitchdiffZ = OUTDAT.AllMinusBase_PitchZ(j, epochtoplot);
                    
                    % ================== CALCULATE THINGS
                    frcorr_hi = corr(ybase_hi, ydev);
                    frcorr_lo = corr(ybase_lo, ydev);
                    
                    if docorrvsdiff ==1
                        frcorr_himinuslo = corr(ybase_himinuslo, ydev);
                    elseif docorrvsdiff==0
                        % then get corr vs. hi and lo, and then subtract
                        frcorr_himinuslo = frcorr_hi - frcorr_lo;
                    end
                    
                    
                    learndir = unique([SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningDirs{2:2:end}]);
                    if length(learndir)>1
                        if onlyPlotIfAllTargSameDir==1
                            disp(learndir)
                            continue
                        end
                        learndir = nan;
                    end
                    istarg = OUTDAT.All_istarg(j);
                    issame = OUTDAT.All_issame(j);
                    
                    %                 if i==1 & ii==1 & ss==7 & istarg==1
                    %                     keyboard
                    %                 end
                    %
                    
                    % ########################### PLOT
                    if learndir==1
                        pcol = 'r';
                    elseif learndir==-1
                        pcol = 'b';
                    end
                    
                    % ----
                    plot(ss, frcorr_himinuslo, 'o', 'Color', pcol);
                    
                    % --- note down if this is a switch that turns off WN
                    functmp = @(x)x(2);
                    postConts = cellfun(functmp, ...
                        SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningContingencies(2:2:end));
                    turnsoff = any(postConts==0);
                    tmp = [SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningContingencies{2:2:end}];
                    lt_plot_text(ss-0.1, 0.9, num2str(tmp), 'k')
                    
                    % ================= do not collect if not showing signs of
                    % learning
                    if onlyIfLearnCorrectDir==1
                        if sign(pitchdiffZ)~=learndir
                            disp('skipped, since no sign of learning');
                            continue
                        end
                    end
                    
                    % === collect
                    PitchDiffZ = [PitchDiffZ; pitchdiffZ];
                    LearnDirAll = [LearnDirAll; learndir];
                    CorrAll = [CorrAll; frcorr_himinuslo];
                    EcounterAll = [EcounterAll; ecount];
                    
                end
            end
            lt_plot_zeroline;
            ylim([-1 1]);
            xlim([0 numswitch+1]);
            
            ecount = ecount+1;
        end
    end
end


%%  ################################ COLLECT SUMMARY DATA
%% ############## go thru all switches, for each neuron COLLECT THINGS ABOUT PREDICT LEARNING

% analytoplot = 'AllDevDiff_NotAbs';
% syltoplot = 'targ';

LearnDirAll = [];
CorrAll = [];
BirdnumAll = [];
EnumAll = [];
EcounterAll = [];
PitchDiffZ = [];
ecount =1 ;
for i=1:numbirds
    numexpts = length(SwitchStruct.bird(i).exptnum);
    
    for ii=1:numexpts
        numswitch = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        
        for ss = 1:numswitch
            
            % ======= which syl class
            if strcmp(syltoplot, 'targ')
                indsthis = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & ...
                    OUTDAT.All_swnum==ss & OUTDAT.All_istarg==1);
            elseif strcmp(syltoplot, 'same')
                indsthis = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & ...
                    OUTDAT.All_swnum==ss & OUTDAT.All_istarg==0 & OUTDAT.All_issame==1);
            elseif strcmp(syltoplot, 'diff')
                indsthis = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & ...
                    OUTDAT.All_swnum==ss & OUTDAT.All_istarg==0 & OUTDAT.All_issame==0);
            end
            
            if isempty(indsthis)
                continue
            end
            
            % ########################## CALCUALTE THINGS
            for j=indsthis'
                
                % ----- get time window
                t = OUTDAT.AllMinusBase_tbinAll{j};
                indstokeep_t = t>=corrwindow(1) & t<corrwindow(2);
                %                 t = t(indstokeep_t);
                
                % ---- SMOOTHED FR, baseline
                ybase_lo = OUTDAT.AllBase_FRsmooth_Lo{j}(indstokeep_t);
                ybase_hi = OUTDAT.AllBase_FRsmooth_Hi{j}(indstokeep_t);
                ybase_himinuslo = ybase_hi - ybase_lo;
                if all(isnan(ybase_himinuslo))
                    assert(all(isnan(OUTDAT.All_FF{j, 1})));
                    continue
                end
                
                % ---- SMOOTHED FR, ACOUNT FOR DRIFT
                if strcmp(analytoplot, 'AllMinusBase_FRmeanAll')
                    ydev = OUTDAT.(analytoplot){j}{epochtoplot}(indstokeep_t);
                else
                    ydev = OUTDAT.(analytoplot){j}(indstokeep_t);
                end
                
                % ------ PITCH, DEVIATION FROM BASELINE
                pitchdiffZ = OUTDAT.AllMinusBase_PitchZ(j, epochtoplot);
                
                % ================== CALCULATE THINGS
                frcorr_hi = corr(ybase_hi, ydev);
                frcorr_lo = corr(ybase_lo, ydev);
                
                if docorrvsdiff ==1
                    frcorr_himinuslo = corr(ybase_himinuslo, ydev);
                elseif docorrvsdiff==0
                    % then get corr vs. hi and lo, and then subtract
                    frcorr_himinuslo = frcorr_hi - frcorr_lo;
                end
                
                % ------------------ LEARN DIR
                learndir = unique([SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningDirs{2:2:end}]);
                if length(learndir)>1
                    if onlyPlotIfAllTargSameDir==1
                        disp(learndir)
                        continue
                    end
                    learndir = nan;
                end
                %                 istarg = OUTDAT.All_istarg(j);
                %                 issame = OUTDAT.All_issame(j);
                
                % ================= do not collect if not showing signs of
                % learning
                if onlyIfLearnCorrectDir==1
                    if sign(pitchdiffZ)~=learndir
                        disp('skipped, since no sign of learning');
                        continue
                    end
                end
                
                % === collect
                PitchDiffZ = [PitchDiffZ; pitchdiffZ];
                LearnDirAll = [LearnDirAll; learndir];
                CorrAll = [CorrAll; frcorr_himinuslo];
                EcounterAll = [EcounterAll; ecount];
                
            end
        end
        ecount = ecount+1;
    end
end


%% ======== STORE OUTPUTS
clear OutStruct;
OutStruct.PitchDiffZ = PitchDiffZ;
OutStruct.LearnDirAll = [LearnDirAll];
OutStruct.CorrAll = [CorrAll];
OutStruct.EcounterAll = [EcounterAll];

%% =========
lt_figure; hold on;

% plot(CorrAll(LearnDirAll==1), 'ok');
plot(LearnDirAll, CorrAll, 'ok');
xlim([-2 2])



%% ##################################### PLOT, BASED ON DIR OF LEARNING
if plotSummary==1
    figcount=1;
    subplotrows=2;
    subplotcols=3;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    % (INSTRUCTION)
    % ====== SUMMARY, ONE DATAPOINT FOR EACH EXPERIMENT (MEAN OVER UNITS,
    % SWITCHES, MOTIFS)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('each expt one line');
    maxecount = max(EcounterAll);
    Yall = single(nan(max(EcounterAll), 2));
    Ysemall = nan(max(EcounterAll), 2);
    for ee=1:maxecount
        
        indsthis = EcounterAll==ee;
        pcol = 0.2 + 0.8*[rand rand rand];
        [ymean, ysem] = grpstats(CorrAll(indsthis), LearnDirAll(indsthis), {'mean', 'sem'});
        x = unique(LearnDirAll(indsthis));
        lt_plot(x+rand, ymean, {'Errors', ysem, 'Color', pcol, 'LineStyle', '-'});
        
        if any(x==-1)
            Yall(ee, 1) = ymean(find(x==-1));
        end
        if any(x==1)
            Yall(ee, 2) = ymean(find(x==1));
        end
        
    end
    lt_plot_zeroline;
    
    try
        indstmp = ~any(isnan(Yall)');
        p = signrank(Yall(indstmp, 1), Yall(indstmp,2));
        lt_plot_pvalue(p, 'signrank (only those paired)');
        p = ranksum(Yall(:,1), Yall(:,2));
        lt_plot_text(0, 0, ['ranksum:' num2str(p)], 'b');
    catch err
    end
    
    
    % ================ SUMMARY ACROSS ALL NEURONS
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([syltoplot ', all neurons']);
    
    [ymean, ysem] = grpstats(CorrAll, LearnDirAll, {'mean', 'sem'});
    x = unique(LearnDirAll);
    
    plot(LearnDirAll, CorrAll, 'ok');
    lt_plot(x, ymean, {'Errors', ysem, 'Color', 'r'});
    
    % ---- significance
    try
        p = ranksum(CorrAll(LearnDirAll==1), CorrAll(LearnDirAll==-1));
        lt_plot_pvalue(p, 'ranksum');
    catch err
    end
    
    % --- significance (each one)
    try
        p = signrank(CorrAll(LearnDirAll==1));
        if p<0.2
            lt_plot_text(1, 1, [num2str(p)], 'r');
        end
        p = signrank(CorrAll(LearnDirAll==-1));
        if p<0.2
            lt_plot_text(-1, 1, [num2str(p)], 'r');
        end
    catch err
    end
    
    xlim([-2 2]);
    lt_plot_zeroline;
    
end


%% ================ REGRESSION MODEL WITH RANDOM EFFECTS
if doregression==1
    CorrAll = double(CorrAll);
    tbl = table(CorrAll, LearnDirAll, EcounterAll, PitchDiffZ);
    mdl = 'CorrAll ~ LearnDirAll + (LearnDirAll|EcounterAll)';
    mdl = 'CorrAll ~ LearnDirAll';
    lme = fitlme(tbl, mdl);
    
    
    % ################################## PLOTS, BASED ON MAGNITUDE OF LEARNING
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('all datapoints');
    xlabel('pitch from base (z)');
    ylabel('corr of fr diff (up is high FF)');
    
    x = PitchDiffZ;
    y = CorrAll;
    
    lt_regress(y, x, 1, 0, 1, 1, 'k');
    lt_plot_makesquare_plot45line(gca, 'b');
    % plot(x, y, 'ok');
    
    
    % ===================== REGRESSION MODEL
    mdl = 'CorrAll ~ PitchDiffZ + (PitchDiffZ|EcounterAll)';
    % mdl = 'CorrAll ~ PitchDiffZ';
    lme = fitlme(tbl, mdl)
    
    % ===================== REGRESSION MODEL
    mdl = 'CorrAll ~ PitchDiffZ + LearnDirAll + (PitchDiffZ|EcounterAll) + (LearnDirAll|EcounterAll)';
    % mdl = 'CorrAll ~ PitchDiffZ';
    lme = fitlme(tbl, mdl)
end

try
    % ========== ALWAYS DO THIS REGRESSION FOR OUTPUT
    CorrAll = double(CorrAll);
    tbl = table(CorrAll, LearnDirAll, EcounterAll, PitchDiffZ);
    mdl = 'CorrAll ~ LearnDirAll + (LearnDirAll|EcounterAll)';
    lme = fitlme(tbl, mdl);
    assert(strcmp(lme.Coefficients.Name{2}, 'LearnDirAll'));
    beta = lme.Coefficients.Estimate(2);
    beta_CI = [lme.Coefficients.Lower(2) lme.Coefficients.Upper(2)];
    beta_p = lme.Coefficients.pValue(2);
    
    OutStruct.FITLME.beta = beta;
    OutStruct.FITLME.beta_CI = beta_CI;
    OutStruct.FITLME.beta_p = beta_p;
catch err
end

