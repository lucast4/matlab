function lt_neural_NGRAMS_PlotByPairtype(OUTSTRUCT, SummaryStruct, plottype, plotON, ...
    dosubtractcontrol)

%% ================== plot all data for all neurons/birds
% plottype = 'absfrdiff'; % oneminusrho or absfrdiff
% plotON=0; % raw plots? only works for absfrdiff
% dosubtractcontrol = 0; % then subtracts negative control before plotting

%% === pull out variables

% ========================== FOR COMPATIBILITY WITH OLD CODE, EXTRACT ALL
% FIELDS
fnamesthis = fieldnames(OUTSTRUCT);
for j=1:length(fnamesthis)
   eval([fnamesthis{j} ' = OUTSTRUCT.' fnamesthis{j} ';']); 
end

%% 

maxbirds = max(All_birdnum);
maxneur = max(All_neurnum);


%% COLLECT AND PLOT

figcount=1;
subplotrows=6;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

Yall = [];
Yall_NEG = [];
Bregionall = {};
for i=1:maxbirds
    birdname = SummaryStruct.birds(i).birdname;
    for ii=1:maxneur
        
        inds = All_birdnum==i & All_neurnum==ii;
        
        if ~any(inds)
            continue
        end
        
        % -- brainregion
        bregion = SummaryStruct.birds(i).neurons(ii).NOTE_Location;
        if strcmp(bregion, 'RA')
            pcol = 'r';
        elseif strcmp(bregion, 'LMAN')
            pcol = 'g';
        end
        
        
        pairtype = All_diffsyl_PairType(inds);
        x = unique(pairtype);
        
        if strcmp(plottype, 'oneminusrho')
            % ==================== Extract data
            oneminusrho_neg = All_OneMinusRho_NEG(inds);
            oneminusrho = All_OneMinusRho(inds);
            
            % ----- subtract negative control?
            if dosubtractcontrol==1
                oneminusrho = oneminusrho - oneminusrho_neg;
            end
            
            % ========= plot
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([birdname '-n' num2str(ii) '[' bregion ']']);
            ylabel('one minus rho');
            
            % --- overlay mean of negative control
            [ymean, ystd] = grpstats(oneminusrho_neg, pairtype, {'mean', 'std'});
            ymeanNEG = ymean;
            shadedErrorBar(x, ymean, ystd, {'Color', [0.7 0.7 0.7]}, 1);
            
            % ---- ALL datapoints
            plot(pairtype, oneminusrho, 'o', 'Color', pcol);
            
            % --- means actual data
            [ymean, ystd] = grpstats(oneminusrho, pairtype, {'mean', 'std'});
            lt_plot(x+0.2,ymean, {'Errors', ystd});
            
            ylim([-0.1 2.1]);
            lt_plot_zeroline;
            xlim([0 length(PairTypesInOrder)+1]);
            
            
            % ==== output for means
            ytmp = nan(1, length(PairTypesInOrder));
            ytmp(x) = ymean;
            Yall = [Yall; ytmp];
            Bregionall = [Bregionall; bregion];
            
            
            ytmp = nan(1, length(PairTypesInOrder));
            ytmp(x) = ymeanNEG;
            Yall_NEG = [Yall_NEG; ytmp];
            
        elseif strcmp(plottype, 'absfrdiff')
            Y = All_AbsFRdiff_Zrelshuff(inds);
            [ymean, ystd] = grpstats(Y, pairtype, {'mean', 'std'});
            if plotON==1
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title([birdname '-n' num2str(ii) '[' bregion ']']);
                ylabel('mean abs(fr diff), zscored');
                
                % ---- plot datapoionts and mean
                plot(pairtype, Y, 'o', 'Color', pcol);
                lt_plot(x+0.2,ymean, {'Errors', ystd});
                
                % --- formating
                lt_plot_zeroline;
                xlim([0 length(PairTypesInOrder)+1]);
                
            end
            
            % ------- output means
            ytmp = nan(1, length(PairTypesInOrder));
            ytmp(x) = ymean;
            Yall = [Yall; ytmp];
            Bregionall = [Bregionall; bregion];
        else
            asdfasdfasdfasddf;
        end
    end
end


%% QUICK AND DIRTY - plot mean/distributions over all data
removeNeurWithNan=1; % only works for absfrdiff

if strcmp(plottype, 'oneminusrho')
    lt_figure; hold on;
    indstokeep = ~any(isnan(Yall)');
    Yall = Yall(indstokeep,:);
    Bregionall = Bregionall(indstokeep);
    Yall_NEG = Yall_NEG(indstokeep,:);
    
    % -- lman
    lt_subplot(2,2,1); hold on;
    title('lman (only complete data)');
    inds = strcmp(Bregionall, 'LMAN');
    y = Yall(inds,:);
    x = 1:size(Yall,2);
    plot(1:size(Yall,2), y, '-k');
    ymean = nanmean(y,1);
    ysem = lt_sem(y);
    lt_plot(x+0.1, ymean, {'Errors', ysem})
    
    yneg = Yall_NEG(inds,:);
    plot(1:size(Yall,2), yneg, '-', 'Color', [0.7 0.3 0.3]);
    plot(1:size(Yall,2), mean(yneg,1), '-', 'Color', [0.7 0.3 0.3], 'LineWidth', 2);
    % -- ra
    lt_subplot(2,2,2); hold on;
    title('ra');
    inds = strcmp(Bregionall, 'RA');
    y = Yall(inds,:);
    x = 1:size(Yall,2);
    plot(1:size(Yall,2), y, '-k');
    ymean = nanmean(y,1);
    ysem = lt_sem(y);
    lt_plot(x+0.1, ymean, {'Errors', ysem})
    
    yneg = Yall_NEG(inds,:);
    plot(1:size(Yall,2), yneg, '-', 'Color', [0.7 0.3 0.3]);
    plot(1:size(Yall,2), mean(yneg,1), '-', 'Color', [0.7 0.3 0.3], 'LineWidth', 2);
    
elseif strcmp(plottype, 'absfrdiff')
    lt_figure; hold on;
    
    if removeNeurWithNan==1
        indstokeep = ~any(isnan(Yall)');
        YallPLOT = Yall(indstokeep,:);
        BregionallPLOT = Bregionall(indstokeep);
    else
        YallPLOT = Yall;
        BregionallPLOT = Bregionall;
    end
    
    % -- lman
    lt_subplot(2,2,1); hold on;
    title('lman');
    inds = strcmp(Bregionall, 'LMAN');
    y = Yall(inds,:);
    x = 1:size(Yall,2);
    plot(1:size(Yall,2), y, '-k');
    ymean = nanmean(y,1);
    ysem = lt_sem(y);
    lt_plot(x+0.1, ymean, {'Errors', ysem, 'Color', 'r'})
    
    % -- ra
    lt_subplot(2,2,2); hold on;
    title('ra');
    inds = strcmp(Bregionall, 'RA');
    y = Yall(inds,:);
    x = 1:size(Yall,2);
    plot(1:size(Yall,2), y, '-k');
    ymean = nanmean(y,1);
    ysem = lt_sem(y);
    lt_plot(x+0.1, ymean, {'Errors', ysem, 'Color', 'r'})
    
    
else
    asdfasdfasdfasddf;
end