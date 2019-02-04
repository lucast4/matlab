function lt_neural_POP_PlotSummary2(MOTIFSTATS_pop, SummaryStruct, DATSTRUCT, ...
    plotRaw)
%% lt 1/22/18 - plots output of lt_neural_POP_PlotSummary

consolidateSamePairs=1; % if 1, then for any pair of neurons and a given motif,
% takes average of pair
bregionWanted = 'LMAN-RA';

%% PARAMS

smoothBeforeSubtr = 0; % smooth before subtract psth or shuffle corrector - see Kass

%% ============= PLOT DISTIRBUTIONS ACROSS MOTIFS/BIRDS

numbirds = max(DATSTRUCT.Pairs.birdID);
CCallbirds = [];
for i=1:numbirds
    
    nummotifs = max(DATSTRUCT.Pairs.motifnum(DATSTRUCT.Pairs.birdID==i));
    
    CCrealCell = cell(1,nummotifs);
    CCshiftCell =cell(1,nummotifs);
    CCCell = cell(1,nummotifs);
    MotifNames = cell(1,nummotifs);
    FFcvMat = nan(1,nummotifs);
    
    birdname = MOTIFSTATS_pop.birds(i).birdname;
    for mm =1:nummotifs
        
        inds = find(DATSTRUCT.Pairs.birdID==i & DATSTRUCT.Pairs.motifnum==mm ...
            & strcmp(DATSTRUCT.Pairs.bregion_string	, bregionWanted));
        
        if isempty(inds)
            continue
        end
        
        
        % ===== SANITY CHECKS
        motifname = unique(DATSTRUCT.Pairs.motifregexp(inds));
        assert(length(motifname) ==1, 'diff motifs, problem ...');
        
        % ##################### OVERLAY ALL XCOV ACROSS PAIRS
        functmp = @(x) mean(x,1);
        ccreal_all = cell2mat(cellfun(functmp, DATSTRUCT.Pairs.ccRealAll(inds), 'UniformOutput', 0));
        ccshift_all = cell2mat(cellfun(functmp, DATSTRUCT.Pairs.ccShiftAll(inds), 'UniformOutput', 0));
        xlags = unique(DATSTRUCT.Pairs.xlags(inds, :));
        
        % ===================== CONSOLIDATE SAME PAIRS (TAKE MEAN)
        if consolidateSamePairs==1
            neurID_all = DATSTRUCT.Pairs.neurIDreal(inds,:);
            neurID_str = [num2str(neurID_all(:,1)) num2str(neurID_all(:,2))];
            neurID_str = mat2cell(neurID_str, ones(size(neurID_str,1),1));
            
            strclasses = tabulate(neurID_str);
            if any(cell2mat(strclasses(:,2))>1)
                % then there are repeats ...
                badclasses = find(cell2mat(strclasses(:,2))>1);

                % ============ AVERAGE OVER THE REPEATS
                ccreal_all_orig = ccreal_all;
                ccshift_all_orig = ccshift_all;
                
                indstoremove_all = [];
                for j=badclasses
                    strid = strclasses{j,1};
                    indstoaverage = find(strcmp(neurID_str, strid));
                    indstoremove_all = [indstoremove_all; indstoaverage];
                    
                    % ==== average over these inds and append to end
                    ccreal_all = [ccreal_all; ...
                        mean(ccreal_all_orig(indstoaverage, :),1)];
                    
                    ccshift_all = [ccshift_all; ...
                        mean(ccshift_all_orig(indstoaverage, :),1)];
                end
                % ====== REMOVE original repeats
                ccreal_all(indstoremove_all, :) = [];
                ccshift_all(indstoremove_all, :) = [];
            end
        end
        cc_all = ccreal_all - ccshift_all;
        
        
        % ========== cv of FF
        if all(isnan(DATSTRUCT.Pairs.FF_cv(inds)))
            ffcv = nan;
        else
            ffcv = mean(DATSTRUCT.Pairs.FF_cv(inds));
        end
        
        
        % =========== save across motifs
        CCrealCell{mm} = ccreal_all;
        CCshiftCell{mm} = ccshift_all;
        CCCell{mm} = cc_all;
        MotifNames{mm} = motifname{1};
        FFcvMat(mm) = ffcv;
        
    end
    
    if isempty(CCrealCell)
        continue
    end
    
    % ################### PLOT FOR ALL MOTIFS [this bird]
    figcount=1;
    subplotrows=3;
    subplotcols=nummotifs;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    % ==== FIRST, raw CC
    hsplots = [];
    for mm=1:nummotifs
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(MotifNames{mm});
        if mm==1
            ylabel([birdname]);
        end
        %         title(['-' motifstr '-nn' num2str(neurons)]);
        
        plot(xlags, CCrealCell{mm}', '-', 'Color', [0.7 0.7 0.7]);
        plot(xlags, mean(CCrealCell{mm}), '-k', 'LineWidth', 3);
        lt_plot_zeroline_vert;
    end
    
    % ==== SECOND, shifted CC
    for mm=1:nummotifs
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        %         title(['-' motifstr '-nn' num2str(neurons)]);
        
        plot(xlags, CCshiftCell{mm}', '-', 'Color', [0.7 0.7 0.7]);
        plot(xlags, mean(CCshiftCell{mm}), '-r', 'LineWidth', 3);
        lt_plot_zeroline_vert;
    end
    linkaxes(hsplots, 'xy');
    
    % ==== THIRD, difference
    hsplots = [];
    for mm=1:nummotifs
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        %         title(['-' motifstr '-nn' num2str(neurons)]);
        
        plot(xlags, CCCell{mm}', '-', 'Color', [0.7 0.7 0.7]);
        plot(xlags, mean(CCCell{mm}), '-b', 'LineWidth', 3);
        lt_plot_zeroline_vert;
        lt_plot_zeroline;
    end
    linkaxes(hsplots, 'xy');
    
    
    % #################### PLOT JUST AVERAGE XCOV (DIFF), FOR ALL MOTIFS
    figcount=1;
    subplotrows=1;
    subplotcols=nummotifs;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    for mm=1:nummotifs
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(MotifNames{mm});
        if mm==1
            ylabel([birdname]);
        end
        
        plot(xlags, mean(CCCell{mm}), '-b', 'LineWidth', 3);
        shadedErrorBar(xlags, mean(CCCell{mm}), lt_sem(CCCell{mm}), {'Color','b'},1);
        lt_plot_zeroline_vert;
        lt_plot_zeroline;
    end
    linkaxes(hsplots, 'xy');
    
    % ########### PLOT CV OVER MOTIFS
    lt_figure; hold on;
    
    % --- 1) cv vs. position along motif
    lt_subplot(2,1,1); hold on;
    ylabel('ff cv');
    plot(1:length(FFcvMat), FFcvMat, '-ok');
    lt_plot_zeroline;
    
    % --- 2) correlate cv with xcov)
    lt_subplot(2,1,2); hold on;
    
    
    % #################### GET AVERAGE ACROSS ALL MOTIFS (THIS BIRD)
    lt_figure; hold on;
    title([birdname ', grand mean, pairs and syls']);
    CCall = [];
    for mm=1:nummotifs
        ccmat = CCCell{mm};
        CCall = [CCall; ccmat];
    end
    
    plot(xlags, CCall', '-', 'Color', [0.7 0.7 0.7]);
    shadedErrorBar(xlags, mean(CCall,1), lt_sem(CCall), {'Color', 'b'}, 1);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    % ============= COLLECT TO PLOT ACROSS BIRDS
    CCallbirds = [CCallbirds; CCall];
    
end

% ==== plot across all birds
lt_figure; hold on;
title(['grand mean, all birds pairs and syls']);
plot(xlags, CCallbirds', '-', 'Color', [0.7 0.7 0.7]);
shadedErrorBar(xlags, mean(CCallbirds,1), lt_sem(CCallbirds), {'Color', 'b'}, 1);
lt_plot_zeroline;
lt_plot_zeroline_vert;



%% ====== PLOT RAW,. EACH MOTIF
if plotRaw==1

    numbirds = max(DATSTRUCT.Pairs.birdID);
    numexpts = max(DATSTRUCT.Pairs.exptID);
    nummotifs = max(DATSTRUCT.Pairs.motifnum);
    
    
    
    Ccallall = [];
    for i=1:numbirds
        bname = MOTIFSTATS_pop.birds(i).birdname;
        for ii=1:numexpts
            exptname = MOTIFSTATS_pop.birds(i).exptnum(ii).exptname;
            for mm =1:nummotifs
                
                motifstr = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(1).motif(mm).regexpstr;
                
                inds = find(DATSTRUCT.Pairs.birdID==i & DATSTRUCT.Pairs.exptID==ii ...
                    & DATSTRUCT.Pairs.motifnum==mm ...
                    & strcmp(DATSTRUCT.Pairs.bregion_string	, bregionWanted));
                
                if isempty(inds)
                    continue
                end
                
                % one figure for each bird/expt/motif
                figcount=1;
                subplotrows=2;
                subplotcols=3;
                fignums_alreadyused=[];
                hfigs=[];
                % #################### PLOT EACH PAIR SEPARATELY
                CCfinalAll = [];
                for jj = inds' % each pair of neuron
                    neurons = DATSTRUCT.Pairs.neurIDreal(jj,:);
                    
                    % ------------------ 0) PLOT CROSS-COV
                    ccRealAll = DATSTRUCT.Pairs.ccRealAll{jj};
                    ccShiftAll = DATSTRUCT.Pairs.ccShiftAll{jj};
                    x = DATSTRUCT.Pairs.xlags(jj,:);
                    
                    % --- get means
                    ccreal = mean(ccRealAll,1);
                    ccreal_sem = lt_sem(ccRealAll);
                    ccneg = mean(ccShiftAll,1);
                    ccneg_sem = lt_sem(ccShiftAll);
                    
                    ccreal_minus =  ccreal - ccneg;
                    %
                    %                 ccreal_minus_all = ccRealAll - repmat(ccneg, size(ccRealAll,1),1);
                    %
                    %% ==========
                    if smoothBeforeSubtr==1
                        % -------- get gaussian window
                        filtdur = 0.015;
                        binsize = x(2)-x(1);
                        N = ceil(filtdur/binsize);
                        if mod(N,2)==0
                            N = N+1;
                        end
                        wind = gausswin(N);
                        wind = wind./sum(wind);
                        
                        % ------- smooth
                        ccreal_sm = conv(ccreal, wind);
                        edge = (length(ccreal_sm)-length(ccreal))/2;
                        ccreal_sm = ccreal_sm(edge+1:end-edge);
                        
                        ccneg_sm = conv(ccneg, wind);
                        edge = (length(ccneg_sm)-length(ccneg))/2;
                        ccneg_sm = ccneg_sm(edge+1:end-edge);
                        
                        % ------- get new final
                        ccreal_minus = ccreal_sm - ccneg_sm;
                        
                        if (0)
                            lt_figure; hold on ;
                            plot(x, ccreal, 'k');
                            plot(x, ccreal_sm, 'r');
                        end
                    end
                    
                    %% ======== auto covariance
                    if (0)
                        ccAuto1minus = DATSTRUCT.Pairs.ccAuto1_minus(jj, :);
                        ccAuto2minus = DATSTRUCT.Pairs.ccAuto2_minus(jj, :);
                        
                        ccAuto1 = DATSTRUCT.Pairs.ccAuto1_mean(jj,:);
                        ccAuto2 = DATSTRUCT.Pairs.ccAuto2_mean(jj,:);
                        
                        
                        ccreal_minus = ccreal./sqrt(ccAuto1.*ccAuto2);
                    end
                    %%
                    if plotRaw==1
                        % ---- plot mean raw vs. shuffle
                        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                        title([bname '-' exptname '-' motifstr '-nn' num2str(neurons)]);
                        shadedErrorBar(x, ccneg, ccneg_sem, {'Color', 'r', 'LineStyle', '--'}, 1)
                        shadedErrorBar(x, ccreal, ccreal_sem, {'Color', 'k'}, 1)
                        axis tight;
                        lt_plot_zeroline;
                        lt_plot_zeroline_vert;
                        
                        % ---- plot autocorr
                        if (0)
                            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                            title('auto, 1(k) and 2(b)');
                            plot(x, ccAuto1, 'k');
                            plot(x, ccAuto2, 'b');
                        end
                        
                        % ---- plot shuffle subtracted
                        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                        %                 title(['b- ' ]);
                        plot(x, ccreal_minus, 'b');
                        axis tight;
                        lt_plot_zeroline;
                        lt_plot_zeroline_vert;
                    end
                    
                    
                    % ------------------- 1) PLOT SMOOTHED FR
                    
                    % ------------------- 2) OVERLAY SYL SEGMENTS
                    
                    % ========================= COLLECT
                    CCfinalAll = [CCfinalAll; ccreal_minus];
                    Ccallall = [Ccallall; ccreal_minus];
                    
                end
                
                % ================== PLOT ALL PAIRS FOR THIS MOTIF
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title('all neuron pairs');
                
                plot(x, CCfinalAll, 'b');
                axis tight;
                lt_plot_zeroline;
                lt_plot_zeroline_vert;
                
                
            end
        end
    end
    
    
    %% =========== plot all
    lt_figure; hold on;
    title('all');
    
    plot(x, Ccallall, 'Color', [0.7 0.7 0.7]);
    ymean = mean(Ccallall);
    ysem = lt_sem(Ccallall);
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'},1);
    
    axis tight;
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    
    
    
end

