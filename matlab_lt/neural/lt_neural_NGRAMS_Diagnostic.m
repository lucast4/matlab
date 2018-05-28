%% =========== CHECK OF RESAMPLING PROCEDURE
% =========== [CHECK] Compare resampled to old data, separated by pairtype
% === one plot per neuron, grouping by pairtype, comparing oold vs. new
% === only plots fr diff (for both data and shuffle).
close all;
maxbirds = max(OUTSTRUCT.All_birdnum);
maxneur = max(OUTSTRUCT.All_neurnum);

figcount=1;
subplotrows=6;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for i=1:maxbirds
    birdname = SummaryStruct.birds(i).birdname;
    for ii=1:maxneur
        
        inds = OUTSTRUCT.All_birdnum==i & OUTSTRUCT.All_neurnum==ii;
        if ~any(inds)
            continue
        end
        
        % =====================
        PairTypes = OUTSTRUCT.All_diffsyl_PairType(inds);
        
        Y_old = OUTSTRUCT.All_AbsFRdiff_ORIG(inds);
        Yneg_old = OUTSTRUCT.All_AbsFRdiff_NEG_ORIG(inds);
        
        Y_new = OUTSTRUCT.All_AbsFRdiff(inds);
        Yneg_new = OUTSTRUCT.All_AbsFRdiff_NEG(inds);
        
        % ==================================== PLOTS [DAT]
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title([birdname ',' num2str(i) '-' num2str(ii)]);
        xlabel('old(gr) - new(bk)');
        ylabel('fr dist');
        
        % ==== plot old
        plot(PairTypes-0.2, Y_old, 'x', 'Color', [0.7 0.7 0.7]);
        [ymean, ystd] = grpstats(Y_old, PairTypes, {'mean', 'std'}); % -- mean, std
        x = unique(PairTypes);
        lt_plot(x-0.1, ymean, {'Errors', ystd, 'Color', [0.7 0.7 0.7]});
        
        % ==== plot new
        plot(PairTypes+0.2, Y_new, 'xk');
        [ymean, ystd] = grpstats(Y_new, PairTypes, {'mean', 'std'}); % -- mean, std
        x = unique(PairTypes);
        lt_plot(x+0.3, ymean, {'Errors', ystd, 'Color', 'k'});
        
        % ----
        set(gca, 'XTick', 1:max(PairTypes));
        
        % ==================================== PLOTS [NEG SHUFF]
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['[SHUFF]' birdname ',' num2str(i) '-' num2str(ii)]);
        xlabel('old(gr) - new(bk)');
        ylabel('fr dist');
        
        % ==== plot old
        plot(PairTypes-0.2, Yneg_old, 'x', 'Color', [0.7 0.7 0.7]);
        [ymean, ystd] = grpstats(Yneg_old, PairTypes, {'mean', 'std'}); % -- mean, std
        x = unique(PairTypes);
        lt_plot(x-0.1, ymean, {'Errors', ystd, 'Color', [0.7 0.7 0.7]});
        
        % ==== plot new
        plot(PairTypes+0.2, Yneg_new, 'xk');
        [ymean, ystd] = grpstats(Yneg_new, PairTypes, {'mean', 'std'}); % -- mean, std
        x = unique(PairTypes);
        lt_plot(x+0.3, ymean, {'Errors', ystd, 'Color', 'k'});
        
        % ----
        axis tight;
        set(gca, 'XTick', 1:max(PairTypes));
        lt_plot_zeroline;
    end
    
    if mod(i,5)==0
        disp('PRESS ANYTHING TO GO TO NEXT SET OF 5 BIRDS')
        pause
        close all;
    end
end


%% ==== [CHECK] compare sample size and negative control decoding
% ===== for two pairtypes, does scatterplot, with each datapoint one
% neuron.
% -== one plot per animal.
% === can compare things like sample size.

PairtypesToplot = {...
    '1  1  1', ... % xaxis
    '1  0  0'}; % yaxis
Indpaircomp = find(ismember(OUTSTRUCT.PairTypesInOrder, PairtypesToplot));

plottype = 'N_new';
% ---- options:
% frdiff_shuff_new
% frdiff_shuff_old
% N_old
% N_new

% ============
maxbirds = max(OUTSTRUCT.All_birdnum);
maxneur = max(OUTSTRUCT.All_neurnum);

figcount=1;
subplotrows=6;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for i=1:maxbirds
    birdname = SummaryStruct.birds(i).birdname;
    
    % ================== one subplot per bird
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['[' plottype ']' birdname ',' num2str(i)]);
    xlabel([OUTSTRUCT.PairTypesInOrder{Indpaircomp(1)}]);
    ylabel([OUTSTRUCT.PairTypesInOrder{Indpaircomp(2)}]);
    
    YYall = [];
    YYstdall = [];
    for ii=1:maxneur
        
        inds = OUTSTRUCT.All_birdnum==i & OUTSTRUCT.All_neurnum==ii;
        if ~any(inds)
            continue
        end
        
        % =====================
        PairTypes = OUTSTRUCT.All_diffsyl_PairType(inds);
        
        if strcmp(plottype, 'frdiff_shuff_new')
            Ythis = OUTSTRUCT.All_AbsFRdiff_NEG(inds);
        elseif strcmp(plottype, 'frdiff_shuff_old')
            Ythis = OUTSTRUCT.All_AbsFRdiff_NEG_ORIG(inds);
        elseif strcmp(plottype, 'N_old')
            Ythis = OUTSTRUCT.All_N_ORIG(inds, :);
            Ythis = mean(Ythis,2); % combine across ngrams
        elseif strcmp(plottype, 'N_new')
            Ythis = OUTSTRUCT.All_N(inds, :);
            Ythis = mean(Ythis,2); % combine across ngrams
        end
        %
        %         Y_old = OUTSTRUCT.All_AbsFRdiff_ORIG(inds);
        %         Yneg_old = OUTSTRUCT.All_AbsFRdiff_NEG_ORIG(inds);
        %
        %         Y_new = OUTSTRUCT.All_AbsFRdiff(inds);
        %         Yneg_new = OUTSTRUCT.All_AbsFRdiff_NEG(inds);
        %
        %         N_old = OUTSTRUCT.All_N(inds,:);
        %         N_new = OUTSTRUCT.All_N_ORIG(inds,:);
        %
        
        % ========================= PLOTS [neg shuff distribition];
        
        % ========== NEW
        Y_all = {};
        Ystd_all = {};
        
        % --------------------------- first pairtype (x)
        indtmp = PairTypes == Indpaircomp(1);
        y = Ythis(indtmp);
        
        % -- collect
        Y_all{1} = mean(y);
        Ystd_all{1} = std(y);
        
        % --------------------------- second pairtype (x)
        indtmp = PairTypes == Indpaircomp(2);
        y = Ythis(indtmp);
        
        % -- collect
        Y_all{2} = mean(y);
        Ystd_all{2} = std(y);
        
        % ======================
        YYall = [YYall; [Y_all{1} Y_all{2}]];
        YYstdall = [YYstdall ; [Ystd_all{1} Ystd_all{2}]];
        
        % ==================== PLOT
        %         lt_plot(Y_all{1}, Y_all{2});
        line([Y_all{1} Y_all{1}], [Y_all{2}-Ystd_all{2} Y_all{2}+Ystd_all{2}], 'Color', 'k');
        line([Y_all{1}-Ystd_all{1} Y_all{1}+Ystd_all{1}], [Y_all{2} Y_all{2}], 'Color', 'k');
    end
    
    % =========== PLOT
    %     lt_plot_45degScatter(YYall(:,1), YYall(:,2), 'k', 1);
    lt_plot(YYall(:,1), YYall(:,2), {'Color', 'k'});
    YLIM = ylim;
    XLIM = xlim;
    line([0 max([XLIM(2) YLIM(2)])], [0 max([XLIM(2) YLIM(2)])]);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
end



%% ============ VARIOUS THINGS HERE [PLOTS]
close all;
dispNgramStrings = 0; % then for 10% of neurons will disp.
plotRawSampSize = 0; % then plots number of trials.
compareNtoFRdist = 1; % MANY PLOTS - compares samp size to fr dist. neuron by neuron.
lt_neural_NGRAMS_DIAGNOSTIC(OUTSTRUCT, SummaryStruct, Params, dispNgramStrings, ...
    plotRawSampSize, compareNtoFRdist);
