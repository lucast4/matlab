function lt_neural_POPLEARN_Xcov_Epochs_Distr(OUTSTRUCT_XCOV, PARAMS, ...
    onlygoodexpt, SwitchStruct, plotRaw, epochWN)
%% lt 3/26/19 - distributions of trials within each epoch

%% ONE PLOT FOR EACH SWITCH (I.E. EACH TARGET SYL)

%%
assert(all(strcmp(OUTSTRUCT_XCOV.bregionpair, 'LMAN-RA')));

%%
if onlygoodexpt==1
    
    % ===== filter outstruct
    [OUTSTRUCT_XCOV, indstokeep] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT_XCOV, ...
        SwitchStruct, 'xcov_spikes');
end



%% go thru each expriemnts


[indsgrp, indsgrpU] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, ...
    OUTSTRUCT_XCOV.enum, OUTSTRUCT_XCOV.switch});

if plotRaw==1
    for i=1:length(indsgrpU)
        
        figcount=1;
        subplotrows=6;
        subplotcols=4;
        fignums_alreadyused=[];
        hfigs=[];
        hsplots = [];
        
        indthis = find(indsgrp==indsgrpU(i) & OUTSTRUCT_XCOV.istarg==1);
        %     assert(sum(indthis)==1, 'multipel targs?');
        
        bnum = unique(OUTSTRUCT_XCOV.bnum(indthis));
        enum = unique(OUTSTRUCT_XCOV.enum(indthis));
        sw = unique(OUTSTRUCT_XCOV.switch(indthis));
        motifname = unique(OUTSTRUCT_XCOV.motifname(indthis)); assert(length(motifname)==1);
        %     bregionpair = unique(OUTSTRUCT_XCOV.b (indthis)); assert(length(motifname)==1);
        %                     hsplots = [];
        
        for j=1:length(indthis)
            
            epochsplit = OUTSTRUCT_XCOV.epochSplitStatsAll{indthis(j)};
            ldir = OUTSTRUCT_XCOV.learndirTarg(indthis(j));
            
            neurpair = OUTSTRUCT_XCOV.neurpair(indthis(j), :);
            
            if j==1
                % ====== PLOT DISTRIBUTION OF FF OVER ALL TRIALS
                for k=1:length(epochsplit)
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    hsplots = [hsplots hsplot];
                    title([num2str(bnum) '-' num2str(enum) '-' num2str(sw) '-' motifname{1}]);
                    ylabel(['laerndir: ' num2str(ldir)]);
                    xlabel('trial num')
                    ff = epochsplit(k).ffthis;
                    if isnan(ff)
                        lt_plot_annotation(1, 'no data, so skip');
                        continue
                    end
                    
                    indshi = epochsplit(k).inds_hi;
                    indslo = epochsplit(k).inds_lo;
                    x = 1:length(ff);
                    
                    plot(x(indshi), ff(indshi), 'or');
                    plot(x(indslo), ff(indslo), 'ok');
                    
                end
                
                % ====== PLOT DISTRIBUTION OF FF WITHIN SONG BOUTS
                for k=1:length(epochsplit)
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    hsplots = [hsplots hsplot];
                    title(['epoch num ' num2str(k)]);
                    xlabel('rendition in bout (10syl min, 1sec IBI min)');
                    ff = epochsplit(k).ffthis;
                    if isnan(ff)
                        lt_plot_annotation(1, 'no data, so skip');
                        continue
                    end
                    indshi = epochsplit(k).inds_hi;
                    indslo = epochsplit(k).inds_lo;
                    songID = epochsplit(k).songboutID;
                    
                    rendinbout = epochsplit(k).rendnumInBout_1secIBI; % if nan, then it is not part of bout with at least 10 syls.
                    
                    plot(rendinbout(indshi), ff(indshi), 'xr');
                    plot(rendinbout(indslo), ff(indslo), 'xk');
                    XLIM = xlim;
                    xlim([XLIM(1)-1 XLIM(2)+1]);
                end
                
                linkaxes(hsplots, 'y');
            end
            
            
            % ====== PLOT DISTRIBUTION OF SPIKE COUNTS
            hsplots = [];
            for k=1:length(epochsplit)
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
                title(['LMAN, ' num2str(neurpair)]);
                ylabel('spikecount');
                xlabel('LO-HI');
                
                %                 ff = epochsplit(k).ffthis;
                if isnan(epochsplit(k).ffthis)
                    lt_plot_annotation(1, 'no data, so skip');
                    continue
                end
                indshi = epochsplit(k).inds_hi;
                indslo = epochsplit(k).inds_lo;
                nspksByNeuron = epochsplit(k).nspksByNeuron;
                
                % --- neuron 1
                nn=1;
                x = [1 2];
                
                Y = {};
                Y{1} = nspksByNeuron{nn}(indslo);
                Y{2} = nspksByNeuron{nn}(indshi);
                lt_plot_MultDist(Y, x, 1, 'k');
                
            end
            linkaxes(hsplots, 'xy');
            
            hsplots = [];
            for k=1:length(epochsplit)
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
                title(['RA, ' num2str(neurpair)]);
                ylabel('spikecount');
                xlabel('LO-HI');
                
                %                 ff = epochsplit(k).ffthis;
                if isnan(epochsplit(k).ffthis)
                    lt_plot_annotation(1, 'no data, so skip');
                    continue
                end
                indshi = epochsplit(k).inds_hi;
                indslo = epochsplit(k).inds_lo;
                nspksByNeuron = epochsplit(k).nspksByNeuron;
                
                % --- neuron 1
                nn=2;
                x = [1 2];
                
                Y = {};
                Y{1} = nspksByNeuron{nn}(indslo);
                Y{2} = nspksByNeuron{nn}(indshi);
                lt_plot_MultDist(Y, x, 1, 'k');
                
            end
            linkaxes(hsplots, 'xy');
        end
    end
end

%% ==================== PLOT SUMMARY IS THERE OVERALL INCREASE IN FIRING RATE?

% ================ BASELINE
epochtoplot = 0;
[spkmeanNonadAdapt_LMAN_BASE, spkmeanNonadAdapt_RA_BASE] = lt_neural_POPLEARN_Xcov_Epochs_DistrSub(...
    OUTSTRUCT_XCOV, epochtoplot);

% ================= WN EPOCH, LAST
for j=1:length(epochtoplot)
    %     epochtoplot = epochWN;
    [spkmeanNonadAdapt_LMAN, spkmeanNonadAdapt_RA] = lt_neural_POPLEARN_Xcov_Epochs_DistrSub(...
        OUTSTRUCT_XCOV, epochWN(j));
end
% ============= PLOT
lt_figure; hold on;
indsthis = OUTSTRUCT_XCOV.istarg==1;

lt_subplot(3,4,1); hold on;
title('LMAN');
xlabel('NONADAPT (BASE - WN) ===== ADAPT (BASE - WN)');

Y = [spkmeanNonadAdapt_LMAN_BASE(indsthis,1) spkmeanNonadAdapt_LMAN(indsthis,1)];
x = [1 2];
plot(x, Y', '-k');
lt_plot(x+0.2, nanmean(Y,1), {'Errors', lt_sem(Y), 'Color', 'r'});

Y = [spkmeanNonadAdapt_LMAN_BASE(indsthis,2) spkmeanNonadAdapt_LMAN(indsthis,2)];
x = [4 5];
plot(x, Y', '-k');
lt_plot(x+0.2, nanmean(Y,1), {'Errors', lt_sem(Y), 'Color', 'r'});

xlim([0 6]);

% ----------------
lt_subplot(3,4,2); hold on;
title('LMAN');
xlabel('BASE(NONAD-AD) ===== WN(NONAD-AD)');

Y = [spkmeanNonadAdapt_LMAN_BASE(indsthis,1) spkmeanNonadAdapt_LMAN_BASE(indsthis,2)];
x = [1 2];
plot(x, Y', '-k');
lt_plot(x+0.2, nanmean(Y,1), {'Errors', lt_sem(Y), 'Color', 'r'});

Y = [spkmeanNonadAdapt_LMAN(indsthis,1) spkmeanNonadAdapt_LMAN(indsthis,2)];
x = [4 5];
plot(x, Y', '-k');
lt_plot(x+0.2, nanmean(Y,1), {'Errors', lt_sem(Y), 'Color', 'r'});

xlim([0 6]);

% -------
lt_subplot(3,4,3); hold on;
title('LMAN');
xlabel('NONADAPT (WN minus base) ===== ADAPT (WN minus base)');

Y =  [spkmeanNonadAdapt_LMAN(indsthis,1) - spkmeanNonadAdapt_LMAN_BASE(indsthis,1), ...
    spkmeanNonadAdapt_LMAN(indsthis,2) - spkmeanNonadAdapt_LMAN_BASE(indsthis,2)];
x = [1 2];

plot(x, Y', '-k');
lt_plot(x+0.2, nanmean(Y,1), {'Errors', lt_sem(Y), 'Color', 'r'});

p = signrank(Y(:,1), Y(:,2));
xlim([0 3]);


% --
lt_subplot(3,4,4); hold on;
title('LMAN (one dat per site)');
xlabel('NONADAPT (WN minus base) ===== ADAPT (WN minus base)');

Y =  [spkmeanNonadAdapt_LMAN(indsthis,1) - spkmeanNonadAdapt_LMAN_BASE(indsthis,1), ...
    spkmeanNonadAdapt_LMAN(indsthis,2) - spkmeanNonadAdapt_LMAN_BASE(indsthis,2)];
[indsgrp] = lt_tools_grp2idx({Y(:,1), Y(:,2)});
Y = grpstats(Y, indsgrp);
x = [1 2];

plot(x, Y', '-k');
lt_plot(x+0.2, nanmean(Y,1), {'Errors', lt_sem(Y), 'Color', 'r'});

p = signrank(Y(:,1), Y(:,2));
% p = lt_plot_zeroline;
lt_plot_pvalue(p, 'srank', 1)
xlim([0 3]);
N = size(Y,1);
lt_plot_annotation(1, ['N=' num2str(N)], 'm');

% #######################

lt_subplot(3,4,5); hold on;
title('RA');
xlabel('NONADAPT (BASE - WN) ===== ADAPT (BASE - WN)');

Y = [spkmeanNonadAdapt_RA_BASE(indsthis,1) spkmeanNonadAdapt_RA(indsthis,1)];
x = [1 2];
plot(x, Y', '-k');
lt_plot(x+0.2, nanmean(Y,1), {'Errors', lt_sem(Y), 'Color', 'r'});

Y = [spkmeanNonadAdapt_RA_BASE(indsthis,2) spkmeanNonadAdapt_RA(indsthis,2)];
x = [4 5];
plot(x, Y', '-k');
lt_plot(x+0.2, nanmean(Y,1), {'Errors', lt_sem(Y), 'Color', 'r'});

xlim([0 6]);

% ----------------
lt_subplot(3,4,6); hold on;
title('RA');
xlabel('BASE(NONAD-AD) ===== WN(NONAD-AD)');

Y = [spkmeanNonadAdapt_RA_BASE(indsthis,1) spkmeanNonadAdapt_RA_BASE(indsthis,2)];
x = [1 2];
plot(x, Y', '-k');
lt_plot(x+0.2, nanmean(Y,1), {'Errors', lt_sem(Y), 'Color', 'r'});

Y = [spkmeanNonadAdapt_RA(indsthis,1) spkmeanNonadAdapt_RA(indsthis,2)];
x = [4 5];
plot(x, Y', '-k');
lt_plot(x+0.2, nanmean(Y,1), {'Errors', lt_sem(Y), 'Color', 'r'});

xlim([0 6]);

% -------
lt_subplot(3,4,7); hold on;
title('RA');
xlabel('NONADAPT (WN minus base) ===== ADAPT (WN minus base)');

Y =  [spkmeanNonadAdapt_RA(indsthis,1) - spkmeanNonadAdapt_RA_BASE(indsthis,1), ...
    spkmeanNonadAdapt_RA(indsthis,2) - spkmeanNonadAdapt_RA_BASE(indsthis,2)];
x = [1 2];

plot(x, Y', '-k');
lt_plot(x+0.2, nanmean(Y,1), {'Errors', lt_sem(Y), 'Color', 'r'});

p = signrank(Y(:,1), Y(:,2));
xlim([0 3]);


% --
lt_subplot(3,4,8); hold on;
title('RA (one dat per site)');
xlabel('NONADAPT (WN minus base) ===== ADAPT (WN minus base)');

Y =  [spkmeanNonadAdapt_RA(indsthis,1) - spkmeanNonadAdapt_RA_BASE(indsthis,1), ...
    spkmeanNonadAdapt_RA(indsthis,2) - spkmeanNonadAdapt_RA_BASE(indsthis,2)];
[indsgrp] = lt_tools_grp2idx({Y(:,1), Y(:,2)});
Y = grpstats(Y, indsgrp);
x = [1 2];

plot(x, Y', '-k');
lt_plot(x+0.2, nanmean(Y,1), {'Errors', lt_sem(Y), 'Color', 'r'});

p = signrank(Y(:,1), Y(:,2));
% p = lt_plot_zeroline;
lt_plot_pvalue(p, 'srank', 1)
xlim([0 3]);
N = size(Y,1);
lt_plot_annotation(1, ['N=' num2str(N)], 'm');


% #######################
spkmeanNonadAdapt_BOTH_BASE = spkmeanNonadAdapt_RA_BASE + spkmeanNonadAdapt_LMAN_BASE;
spkmeanNonadAdapt_BOTH = spkmeanNonadAdapt_RA + spkmeanNonadAdapt_LMAN;

lt_subplot(3,4,9); hold on;
title('RA+LMAN');
xlabel('NONADAPT (BASE - WN) ===== ADAPT (BASE - WN)');

Y = [spkmeanNonadAdapt_BOTH_BASE(indsthis,1) spkmeanNonadAdapt_BOTH(indsthis,1)];
x = [1 2];
plot(x, Y', '-k');
lt_plot(x+0.2, nanmean(Y,1), {'Errors', lt_sem(Y), 'Color', 'r'});

Y = [spkmeanNonadAdapt_BOTH_BASE(indsthis,2) spkmeanNonadAdapt_BOTH(indsthis,2)];
x = [4 5];
plot(x, Y', '-k');
lt_plot(x+0.2, nanmean(Y,1), {'Errors', lt_sem(Y), 'Color', 'r'});

xlim([0 6]);

% ----------------
lt_subplot(3,4,10); hold on;
title('RA+LMAN');
xlabel('BASE(NONAD-AD) ===== WN(NONAD-AD)');

Y = [spkmeanNonadAdapt_BOTH_BASE(indsthis,1) spkmeanNonadAdapt_BOTH_BASE(indsthis,2)];
x = [1 2];
plot(x, Y', '-k');
lt_plot(x+0.2, nanmean(Y,1), {'Errors', lt_sem(Y), 'Color', 'r'});

Y = [spkmeanNonadAdapt_BOTH(indsthis,1) spkmeanNonadAdapt_BOTH(indsthis,2)];
x = [4 5];
plot(x, Y', '-k');
lt_plot(x+0.2, nanmean(Y,1), {'Errors', lt_sem(Y), 'Color', 'r'});

xlim([0 6]);

% -------
lt_subplot(3,4,11); hold on;
title('RA+LMAN');
xlabel('NONADAPT (WN minus base) ===== ADAPT (WN minus base)');

Y =  [spkmeanNonadAdapt_BOTH(indsthis,1) - spkmeanNonadAdapt_BOTH_BASE(indsthis,1), ...
    spkmeanNonadAdapt_BOTH(indsthis,2) - spkmeanNonadAdapt_BOTH_BASE(indsthis,2)];
x = [1 2];

plot(x, Y', '-k');
lt_plot(x+0.2, nanmean(Y,1), {'Errors', lt_sem(Y), 'Color', 'r'});

p = signrank(Y(:,1), Y(:,2));
xlim([0 3]);


% --
lt_subplot(3,4,12); hold on;
title('RA+LMAN (one dat per site)');
xlabel('NONADAPT (WN minus base) ===== ADAPT (WN minus base)');

Y =  [spkmeanNonadAdapt_BOTH(indsthis,1) - spkmeanNonadAdapt_BOTH_BASE(indsthis,1), ...
    spkmeanNonadAdapt_BOTH(indsthis,2) - spkmeanNonadAdapt_BOTH_BASE(indsthis,2)];
[indsgrp] = lt_tools_grp2idx({Y(:,1), Y(:,2)});
Y = grpstats(Y, indsgrp);
x = [1 2];

plot(x, Y', '-k');
lt_plot(x+0.2, nanmean(Y,1), {'Errors', lt_sem(Y), 'Color', 'r'});

p = signrank(Y(:,1), Y(:,2));
% p = lt_plot_zeroline;
lt_plot_pvalue(p, 'srank', 1)
xlim([0 3]);

N = size(Y,1);
lt_plot_annotation(1, ['N=' num2str(N)], 'm');


