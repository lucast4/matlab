%% ============ [DOWNSAMPLING]
% QUESTION: RA vs. LMAN, what if downsample RA so that equalize positive
% control values between RA and LMAN?

% === for each neuron, regardless of brain region, get values as function
% of sample size.

measure_to_recalc = 'absfrdiff';
downfactorlist = [0.025 0.05 0.1 0.2 0.5 0.75 1]; % fractin of samples between N and Nmin.
[OUTSTRUCT_subsamp, downfactorlist] = lt_neural_NGRAMS_Downsample(OUTSTRUCT, SummaryStruct, Params, ...
    measure_to_recalc, downfactorlist);
Params.downfactorlist = downfactorlist;

%% ============= [DOWNSAMPLE] === plot each neuron, fr diff
% as function of downsample.

close all;

PairtypesToplot = {...
    '1  1  1', ... % xaxis
    '1  0  0'}; % yaxis
Indpaircomp = find(ismember(OUTSTRUCT.PairTypesInOrder, PairtypesToplot));

% ============
maxbirds = max(OUTSTRUCT.All_birdnum);
maxneur = max(OUTSTRUCT.All_neurnum);


for i=1:maxbirds
    birdname = SummaryStruct.birds(i).birdname;
    
    
    figcount=1;
    subplotrows=6;
    subplotcols=3;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    for ii=1:maxneur
        
        inds = OUTSTRUCT.All_birdnum==i & OUTSTRUCT.All_neurnum==ii;
        if ~any(inds)
            continue
        end
        
        bregion = SummaryStruct.birds(i).neurons(ii).NOTE_Location;
        
        % =====================
        PairTypes = OUTSTRUCT.All_diffsyl_PairType(inds);
        Ydat = OUTSTRUCT_subsamp.FRdiffDAT(inds,:);
        Yshuff = OUTSTRUCT_subsamp.FRdiffShuff(inds,:);
        Nmean = OUTSTRUCT_subsamp.Nboth(inds, :, :);
        Nmean = squeeze(mean(Nmean, 2));
        
        % --------------------------- first pairtype (x)
        indtype = 1;
        
        indtmp = PairTypes == Indpaircomp(indtype);
        ydat = Ydat(indtmp,:);
        yshuff = Yshuff(indtmp, :);
        %         x = Nmean(indtmp,:); % sample size, mean between 2 motifs
        x = Params.downfactorlist;
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title([birdname ',' num2str(i) '-n' num2str(ii) '-' bregion]);
        %         xlabel('subsample N (mean 2 motifs)');
        xlabel('downfactor');
        ylabel([OUTSTRUCT.PairTypesInOrder{Indpaircomp(indtype)}]);
        
        if ~isempty(ydat)
        % -- dat
        if sum(indtmp)<100
            plot(x, ydat, '-', 'Color', [0.7 0.7 0.7]);
            plot(x, yshuff, '-', 'Color', [0.7 0.2 0.2]);
        end
        % -- means
        plot(x, mean(yshuff,1), '-r', 'LineWidth', 2);
        plot(x, mean(ydat,1), '-k', 'LineWidth', 2);
        end
        
        % --------------------------- second pairtype (x)
        indtype = 2;
        
        indtmp = PairTypes == Indpaircomp(indtype);
        ydat = Ydat(indtmp,:);
        yshuff = Yshuff(indtmp, :);
        %         x = Nmean(indtmp,:); % sample size, mean between 2 motifs
        x = Params.downfactorlist;
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title([birdname ',' num2str(i) '-n' num2str(ii)]);
        %         xlabel('subsample N (mean 2 motifs)');
        xlabel('downfactor');
        ylabel([OUTSTRUCT.PairTypesInOrder{Indpaircomp(indtype)}]);
        if ~isempty(ydat)
        % -- dat
        if sum(indtmp)<100
            plot(x, ydat, '-', 'Color', [0.7 0.7 0.7]);
            plot(x, yshuff, '-', 'Color', [0.7 0.2 0.2]);
        end
        % -- means
        plot(x, mean(yshuff,1), '-r', 'LineWidth', 2);
        plot(x, mean(ydat,1), '-k', 'LineWidth', 2);
        end
        
        % ------------------------ COMBINE (EACH TYPE, SUBTRACT GLOBAL NEG)
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title([birdname ',' num2str(i) '-n' num2str(ii)]);
        %         xlabel('subsample N (mean 2 motifs)');
        xlabel('downfactor');
        ylabel(['frdiff']);
        
        Yall = {}; % dat, dat, shuff
        
        % -- type 1 (dat)
        indtype = 1;
        indtmp = PairTypes == Indpaircomp(indtype);
        ydat = mean(Ydat(indtmp,:),1);
        plot(x, ydat, '--r');
        lt_plot_text(x(1), ydat(1), OUTSTRUCT.PairTypesInOrder{Indpaircomp(indtype)});
        
        Yall{indtype} = ydat;
        
        % -- type 1 (dat)
        indtype = 2;
        indtmp = PairTypes == Indpaircomp(indtype);
        ydat = mean(Ydat(indtmp,:),1);
        plot(x, ydat, '--b');
        lt_plot_text(x(1), ydat(1), OUTSTRUCT.PairTypesInOrder{Indpaircomp(indtype)});
        
        Yall{indtype} = ydat;
        
        % -- negative
        indtmp = ismember(PairTypes, Indpaircomp);
        yshuff = mean(Yshuff(indtmp,:),1);
        
        plot(x, yshuff, '--k')
        
        % ============= OVERLAY DIFFERENCES
        plot(x, Yall{1}-yshuff, '-r');
        plot(x, Yall{2}-yshuff, '-b');
        lt_plot_zeroline;
    end
end


%% =========== [DOWNSAMPLE] - replace by certain factor of downsample
% then redo scatter plot analysis
dsampvalue_take = 0.05;
dsampbregion = {'RA'}; % if empty, downsamples none. e.g. {'LMAN', 'RA'};

OUTSTRUCT_tmp = OUTSTRUCT;

dsampind = find(Params.downfactorlist ==  dsampvalue_take);
assert(length(dsampind)==1, 'this dsamp value doesnt exist...');

% ================= REPLACE DATA WITH DOWNSAMPLED DATA
indstodo = ismember(OUTSTRUCT_tmp.All_Bregion, dsampbregion);

OUTSTRUCT_tmp.All_AbsFRdiff(indstodo) = ...
    OUTSTRUCT_subsamp.FRdiffDAT(indstodo, dsampind);
OUTSTRUCT_tmp.All_AbsFRdiff_NEG(indstodo) = ...
    OUTSTRUCT_subsamp.FRdiffShuff(indstodo, dsampind);
OUTSTRUCT_tmp.All_AbsFRdiff_Zrelshuff(indstodo) = ...
    OUTSTRUCT_subsamp.FRdiff_Z(indstodo, dsampind);
OUTSTRUCT_tmp.All_N(indstodo, :) = ...
    OUTSTRUCT_subsamp.Nboth(indstodo, :, dsampind);

% ================= RUN SCATTERPLOT
close all;

plottype = 'absfrdiff_typediff'; % oneminusrho or absfrdiff or absfrdiff_globZ or absfrdiff_typediff
usemedian = 0; % only workds for absfrdiff, absfrdiff_globZ or absfrdiff_typediff
plotON=0; % only works for absfrdiff
if strcmp(Params.strtype, 'xaa')
    PairtypesToplot = {...
        '1  1  1', ... % xaxis
        '1  0  0'}; % yaxis
elseif strcmp(Params.strtype, 'xaaa')
    PairtypesToplot = {...
        '1  1  1  1', ... % xaxis
        '1  0  0  0'}; % yaxis
end
plotRawGood = 0; % histogram for pos, negative, dat (only works for plottype = absfrdiff_globZ);

% ----------- params for one minus rho, specifically
dosubtractcontrol = 1; % then subtracts negative control before plotting [if 0, then overlays neg]
sanitycheckuseneg = 0; % uses negative control data instead of data
removeBadSyls = 1; % i.e. badly labeled...

[AllPairs_Means, AllPairs_Birdnum, AllPairs_Bregions] = ...
    lt_neural_NGRAMS_PlotScatter(OUTSTRUCT_tmp, SummaryStruct, plottype, plotON, ...
    PairtypesToplot, dosubtractcontrol, sanitycheckuseneg, plotRawGood, usemedian, ...
    removeBadSyls);

