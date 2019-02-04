close all;
% ---- what to plot?
birdtoplot = 'pu69wh78';
expttoplot = 'RALMANOvernightLearn1';
% expttoplot = 'RALMANlearn1';
swnum = 1;


% ##################################### 1) ======= PLOT RAW (EACH MOTIF AND CHAN PAIR)
% ----
averagechanpairs= 0; % for each motif, average over all chan pairs
removeBadSyls = 1; % LEAVE AT 1.
lt_neural_Coher_PlotChanMotif(SwitchStruct, SwitchCohStruct, OUTSTRUCT, ...
    birdtoplot, expttoplot, swnum, averagechanpairs, PARAMS, removeBadSyls)


% ##################################### 1) ======= PLOT RAW (EACH MOTIF, AVERAGE CHAN PAIR)
% ----
ffbinsedges = [20 35 80 120];
averagechanpairs= 1; % for each motif, average over all chan pairs
removeBadSyls = 1; % LEAVE AT 1.
lt_neural_Coher_PlotChanMotif(SwitchStruct, SwitchCohStruct, OUTSTRUCT, ...
    birdtoplot, expttoplot, swnum, averagechanpairs, PARAMS, removeBadSyls, ...
    ffbinsedges)



% 3) ################### ONE PLOT FOR EACH CHANPAIR (AVERAGES OVER MOTIFS)
plotrawcohdiff = 1; % if 0, then takes average over motifs. if 1, then plots each raw
swtoplot = swnum;

% assert(averagechanpairs==0, 'following needs each pair broken out');
tbins = PARAMS.tbins;
ffbins = PARAMS.ffbins;
clim = [-0.15 0.15];

% ======= for each channel pair, get mean for targ, nontarg
% -- gets motifs for each channel pair
indsgrp = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch, OUTSTRUCT.chanpair});
indsgrp_unique = unique(indsgrp);

figcount=1;
subplotrows=3;
subplotcols=10;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for j=indsgrp_unique'
    
    % ----- info...
    bname = SwitchStruct.bird(unique(OUTSTRUCT.bnum(indsgrp==j))).birdname;
    ename = SwitchStruct.bird(unique(OUTSTRUCT.bnum(indsgrp==j))).exptnum(unique(OUTSTRUCT.enum(indsgrp==j))).exptname;
    swnum = unique(OUTSTRUCT.switch(indsgrp==j));
    chpair = OUTSTRUCT.chanpair(indsgrp==j, :);
    chpair = chpair(1,:);
    
    if ~(strcmp(birdtoplot, bname) & strcmp(ename, expttoplot) & swnum==swtoplot)
        continue
    end
    
    % ---------- target
    indsthis = indsgrp==j & OUTSTRUCT.istarg==1;
    if ~any(indsthis)
        continue
    end
    cohmat = lt_neural_Coher_Cell2Mat(OUTSTRUCT.CohMean_WNminusBase(indsthis));
    
    if isempty(cohmat)
        continue
    end
    
    % 1. cohgram of diff
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 1, '', clim);
    title('TARG');
    ylabel({[bname '-' ename '-sw' num2str(ss)], ['ch:' num2str(chpair)]});
    % 2. ffband differences
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1, plotrawcohdiff);
    lt_plot_zeroline;
    lt_plot_text(0, clim(2)-0.05, ['n=' num2str(size(cohmat,3))], 'b');
    
    % ---------- same-type
    indsthis = indsgrp==j & OUTSTRUCT.istarg==0 & OUTSTRUCT.issame==1;
    if ~any(indsthis)
        continue
    end
    cohmat = lt_neural_Coher_Cell2Mat(OUTSTRUCT.CohMean_WNminusBase(indsthis));
    % 1. cohgram of diff
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('SAME');
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 1, '', clim);
    % 2. ffband differences
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1, plotrawcohdiff);
    lt_plot_zeroline;
    lt_plot_text(0, clim(2)-0.05, ['n=' num2str(size(cohmat,3))], 'b');
    
    % -------- diff type
    indsthis = indsgrp==j & OUTSTRUCT.istarg==0 & OUTSTRUCT.issame==0;
    if ~any(indsthis)
        continue
    end
    cohmat = lt_neural_Coher_Cell2Mat(OUTSTRUCT.CohMean_WNminusBase(indsthis));
    % 1. cohgram of diff
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('DIFF');
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 1, '', clim);
    % 2. ffband differences
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1, plotrawcohdiff);
    lt_plot_zeroline;
    lt_plot_text(0, clim(2)-0.05, ['n=' num2str(size(cohmat,3))], 'b');
    
    % ------------ TARGET MINUS SAMETYPE
    indsthis = indsgrp==j & OUTSTRUCT.istarg==1;
    cohmat1 = lt_neural_Coher_Cell2Mat(OUTSTRUCT.CohMean_WNminusBase(indsthis));
    indsthis = indsgrp==j & OUTSTRUCT.istarg==0 & OUTSTRUCT.issame==1;
    cohmat2 = lt_neural_Coher_Cell2Mat(OUTSTRUCT.CohMean_WNminusBase(indsthis));
    cohmat = nanmean(cohmat1,3)-nanmean(cohmat2,3);
    
    % 1. cohgram of diff
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('TARG-SAME');
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 1, '', clim);
    % 2. ffband differences
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1);
    lt_plot_zeroline;
    lt_plot_text(0, clim(2)-0.05, ['n=' num2str(size(cohmat,3))], 'b');
    
    
    % ------------ TARGET MINUS DIFFTYPE
    indsthis = indsgrp==j & OUTSTRUCT.istarg==1;
    cohmat1 = lt_neural_Coher_Cell2Mat(OUTSTRUCT.CohMean_WNminusBase(indsthis));
    indsthis = indsgrp==j & OUTSTRUCT.istarg==0 & OUTSTRUCT.issame==0;
    cohmat2 = lt_neural_Coher_Cell2Mat(OUTSTRUCT.CohMean_WNminusBase(indsthis));
    cohmat = nanmean(cohmat1,3)-nanmean(cohmat2,3);
    
    % 1. cohgram of diff
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('TARG-DIFF');
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 1, '', clim);
    % 2. ffband differences
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1);
    lt_plot_zeroline;
    lt_plot_text(0, clim(2)-0.05, ['n=' num2str(size(cohmat,3))], 'b');
end



%% ##################################### 2) PLOT TIMECOURSE OF COHERENCE SCALAR CHANGES
swnum = swtoplot;

% =================== RUN
bnumthis = find(strcmp(birdtoplot, {MOTIFSTATS_pop.birds.birdname}));
enumthis = find(strcmp(expttoplot, {MOTIFSTATS_pop.birds(bnumthis).exptnum.exptname}));

indsthis = OUTSTRUCT.bnum==bnumthis & OUTSTRUCT.enum==enumthis ...
    & OUTSTRUCT.switch==swnum;

% PLOT EACH MOTIF SEPARATELY
motifnum_unique = unique(OUTSTRUCT.motifnum(indsthis));
figcount=1;
subplotrows=6;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for mm=motifnum_unique'
    indsthis = OUTSTRUCT.bnum==bnumthis & OUTSTRUCT.enum==enumthis ...
        & OUTSTRUCT.switch==swnum & OUTSTRUCT.motifnum==mm;
    
    motifname = SwitchCohStruct.bird(bnumthis).exptnum(enumthis).switchlist(swnum).motifnum(mm).motifname;
    istarg = unique(OUTSTRUCT.istarg(indsthis));
    issame = unique(OUTSTRUCT.issame(indsthis));
    if istarg==1
        pcoltit = 'r';
    elseif issame ==1
        pcoltit = 'b';
    else
        pcoltit = 'k';
    end
    
    t = OUTSTRUCT.tvals(indsthis); assert(length(unique(cellfun(@mean, t)))==1, 'should all be same /./');
    t = t{1};
    ff = SwitchCohStruct.bird(bnumthis).exptnum(enumthis).switchlist(swnum).motifnum(mm).ffvals;
    cohscalall = OUTSTRUCT.CohMat_scalar(indsthis);
    cohscalall = cell2mat(cellfun(@transpose, cohscalall, 'UniformOutput', 0));
    
    % ------ only get inds during base and WN
    indsbase = OUTSTRUCT.indsbase(indsthis); assert(length(unique(cellfun(@sum, indsbase)))==1, 'should all be same /./');
    indsWN = OUTSTRUCT.indsWN(indsthis); assert(length(unique(cellfun(@sum, indsWN)))==1, 'should all be same /./');
    indsbase = indsbase{1};
    indsWN = indsWN{1};
    
    t = t(indsbase | indsWN);
    ff = ff(indsbase | indsWN);
    cohscalall = cohscalall(:, indsbase | indsWN);
    
    % ================ PLOT T
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([birdtoplot '-' expttoplot '-sw' num2str(swnum)]);
    ylabel([motifname], 'Color', pcoltit);
    xlabel('t');
    plot(t, ff, 'xk');
    axis tight
    % -- line for WN on
    tswitch = mean(t(sum(indsbase):sum(indsbase)+1));
    line([tswitch tswitch], ylim, 'Color', 'r');
    
    % ================ PLOT COHERENCE, EACH CHANNEL
    for j=1:size(cohscalall,1)
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['chanpair: ']);
        ylabel('coherence (in t/f wind)');
        xlabel('t');
        cohscal = cohscalall(j,:);
        plot(t, cohscal, 'xk');
        axis tight;
        ylim([0.1 0.9]);
        % -- line for WN on
        tswitch = mean(t(sum(indsbase):sum(indsbase)+1));
        line([tswitch tswitch], ylim, 'Color', 'r');
    end
    
    
    
end
