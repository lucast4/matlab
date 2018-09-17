function lt_neural_v2_ANALY_FRsmooth_PredLMult(OUTDAT, MOTIFSTATS_Compiled, ...
    SwitchStruct, syltoplot, timeshift, windsize, syltypesneeded, ...
    epochtoplot,analytoplot, onlyIfLearnCorrectDir)
%% inputs
% timeshift = 0.01;
% windsize = 0.1;
% syltoplot = 'targ';
%
%     syltypesneeded = [1 0 1];

% default params:
%     onlyIfLearnCorrectDir = 1; % onoy applies to summary plots
%     analytoplot = 'AllDevDiff_NotAbs';

%% =========== figure out time windows to iterate over.

predur = MOTIFSTATS_Compiled.birds(1).exptnum(1).MOTIFSTATS.params.motif_predur;
postdur = MOTIFSTATS_Compiled.birds(1).exptnum(1).MOTIFSTATS.params.motif_postdur;

tmp1 = [(-predur+0.02):timeshift:(postdur-windsize-0.02)];
tmp2 = [(-predur+windsize+0.02):timeshift:(postdur-0.02)];
TimeWindList = [tmp1' tmp2'];


%% =========== iterate over time windows and collect summary stats
PREDSTRUCT = [];
for j=1:size(TimeWindList, 1)
    % corrwindow = [-0.025 0.04];
    corrwindow = TimeWindList(j,:);
    OutStruct = lt_neural_v2_ANALY_FRsmooth_PredLearn(OUTDAT, MOTIFSTATS_Compiled, ...
        SwitchStruct, corrwindow, syltypesneeded, epochtoplot, analytoplot, ...
        syltoplot, 1, 1, onlyIfLearnCorrectDir, ...
        0, 0, 0);
    PREDSTRUCT(j).OutStruct = OutStruct;
end


%% ----------- PLOT ALL TIME WINDOWS
figcount=1;
subplotrows=4;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];
DiffMeanAll = [];
LME_beta = [];
LME_CI = [];

for j=1:length(PREDSTRUCT)
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([num2str(TimeWindList(j,:))]);
    
    EcounterAll = PREDSTRUCT(j).OutStruct.EcounterAll;
    CorrAll = PREDSTRUCT(j).OutStruct.CorrAll;
    LearnDirAll = PREDSTRUCT(j).OutStruct.LearnDirAll;
    
    maxecount = max(EcounterAll);
    Yall = single(nan(max(EcounterAll), 2));
    Ysemall = nan(max(EcounterAll), 2);
    for ee=1:maxecount
        
        indsthis = EcounterAll==ee;
        pcol = 0.2 + 0.8*[rand rand rand];
        [ymean, ysem] = grpstats(CorrAll(indsthis), LearnDirAll(indsthis), {'mean', 'sem'});
        x = unique(LearnDirAll(indsthis));
        lt_plot(x-0.5+rand, ymean, {'Errors', ysem, 'Color', pcol, 'LineStyle', '-'});
        
        if any(x==-1)
            Yall(ee, 1) = ymean(find(x==-1));
        end
        if any(x==1)
            Yall(ee, 2) = ymean(find(x==1));
        end
        
    end
    lt_plot_zeroline;
    
    indstmp = ~any(isnan(Yall)');
    p = signrank(Yall(indstmp, 1), Yall(indstmp,2));
    lt_plot_pvalue(p, 'signrank (only those paired)');
    p = ranksum(Yall(:,1), Yall(:,2));
    lt_plot_text(0, 0, ['ranksum:' num2str(p)], 'b');
    ylim([-1 1]);
    
    % ======= collect measure of effect size (first just take diff of
    % means)
    lt_plot_bar(unique(LearnDirAll), nanmean(Yall,1), {'Errors', lt_sem(Yall), ...
        'Color', 'k', 'FaceAlpha', 0.8});
    
    diffmean = nanmean(Yall(:,2)) - nanmean(Yall(:,1));
    DiffMeanAll = [DiffMeanAll; diffmean];
    
    LME_beta = [LME_beta; PREDSTRUCT(j).OutStruct.FITLME.beta];
    LME_CI = [LME_CI; PREDSTRUCT(j).OutStruct.FITLME.beta_CI];
    
end
lt_figure; hold on;

lt_subplot(2,2,1); hold on
title('mean of each expt');
ylabel('mean');
xlabel('cetner of time wind');
plot(mean(TimeWindList,2), DiffMeanAll, '-ok');
lt_plot_zeroline;

lt_subplot(2,2,2); hold on
ylabel('lme, est and CI');
lt_plot(mean(TimeWindList,2), LME_beta, {'LineStyle', '-', 'LineWidth', 2});
plot(mean(TimeWindList,2), LME_beta, '-k');
plot(mean(TimeWindList,2), LME_CI(:,1), '-k');
plot(mean(TimeWindList,2), LME_CI(:,2), '-k');
lt_plot_zeroline;
