function lt_neural_v2_ANALY_FRsmooth_BaseMinu(OUTDAT, SwitchStruct, ...
    epochtoplot, analytype, doShuff, syltypesneeded)
%% ####################### GETTING DEVIATION FROM BASELINE SMOOTHED FR

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

%% ####################### SUMMARY - collect, for each subtract global mean

% ======= ACTYAL DAT (NOT SHUFF)
shuffSylType =0;
plotOn = 1;
[OUTDAT_tmp] = lt_neural_v2_ANALY_FRsmooth_Comps(OUTDAT, SwitchStruct, shuffSylType, ...
    epochtoplot, plotOn);

% [OUTDAT_tmp] = fn_summaryplot(OUTDAT, SwitchStruct, shuffSylType, ...
%     epochtoplot, plotOn);
Ytarg_dat = mean(cell2mat(OUTDAT_tmp.(analytype)(OUTDAT.All_istarg==1)),1);
Ysame_dat = mean(cell2mat(OUTDAT_tmp.(analytype)(OUTDAT.All_istarg==0 & OUTDAT.All_issame==1)),1);
Ydiff_dat = mean(cell2mat(OUTDAT_tmp.(analytype)(OUTDAT.All_istarg==0 & OUTDAT.All_issame==0)),1);


%% ===== also compare to shuffle analysis?
if doShuff==1
    % ====== many shuffles,
    % each shuffle, permute motif type
    nshuffs = 1000;
    
    Ytarg_shuff = [];
    Ysame_shuff = [];
    Ydiff_shuff = [];
    for n=1:nshuffs
        disp(['shuff ' num2str(n)]);
        
        shuffSylType = 1;
        plotOn = 0;
        [OUTDAT_tmp] = lt_neural_v2_ANALY_FRsmooth_Comps(OUTDAT, SwitchStruct, shuffSylType, ...
            epochtoplot, plotOn);
%         [OUTDAT_tmp] = fn_summaryplot(OUTDAT, SwitchStruct, shuffSylType, ...
%             epochtoplot, plotOn);
        
        
        % ---- save means for diff syl types
        Ytarg_shuff = [Ytarg_shuff; mean(cell2mat(OUTDAT_tmp.(analytype)(OUTDAT_tmp.All_istarg==1)),1)];
        Ysame_shuff = [Ysame_shuff; mean(cell2mat(OUTDAT_tmp.(analytype)(OUTDAT_tmp.All_istarg==0 & OUTDAT_tmp.All_issame==1)),1)];
        Ydiff_shuff = [Ydiff_shuff; mean(cell2mat(OUTDAT_tmp.(analytype)(OUTDAT_tmp.All_istarg==0 & OUTDAT_tmp.All_issame==0)),1)];
        
        
        %     Ytarg_shuff = [Ytarg_shuff; mean(cell2mat(OUTDAT_tmp.AllMinusAll_FRsmooth(OUTDAT_tmp.All_istarg==1)),1)];
        %     Ysame_shuff = [Ysame_shuff; mean(cell2mat(OUTDAT_tmp.AllMinusAll_FRsmooth(OUTDAT_tmp.All_istarg==0 & OUTDAT_tmp.All_issame==1)),1)];
        %     Ydiff_shuff = [Ydiff_shuff; mean(cell2mat(OUTDAT_tmp.AllMinusAll_FRsmooth(OUTDAT_tmp.All_istarg==0 & OUTDAT_tmp.All_issame==0)),1)];
    end
    
    
    % ################################ PLOT SHUFFLE OUTCOMES
    figcount=1;
    subplotrows=4;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    tbin = OUTDAT.All_FRsmooth_t{1};
    
    % ===== 1) overlay distribution of data and shuffle
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('TARG');
    datreal = Ytarg_dat;
    datshuff = Ytarg_shuff;
    % --- run
    yCI = prctile(datshuff, [2.5 97.5], 1);
    ymean = mean(datshuff,1);
    ydat = datreal;
    
    plot(tbin, mean(ymean,1), 'k-', 'LineWidth', 2);
    plot(tbin, yCI(1,:), 'k-', 'LineWidth', 1);
    plot(tbin, yCI(2,:), 'k-', 'LineWidth', 1);
    plot(tbin, ydat, '-r', 'LineWidth', 2);
    
    % -- 2) p value at each time bin
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    ylabel('log10(p)');
    
    tmp = datshuff>datreal;
    p = log10(sum(tmp,1)./size(tmp,1)+1/nshuffs);
    plot(tbin, p, '-k');
    line(xlim, log10([0.05 0.05]));
    
    
    
    % ===== 1) overlay distribution of data and shuffle
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('SAME');
    datreal = Ysame_dat;
    datshuff = Ysame_shuff;
    % --- run
    yCI = prctile(datshuff, [2.5 97.5], 1);
    ymean = mean(datshuff,1);
    ydat = datreal;
    
    plot(tbin, mean(ymean,1), 'k-', 'LineWidth', 2);
    plot(tbin, yCI(1,:), 'k-', 'LineWidth', 1);
    plot(tbin, yCI(2,:), 'k-', 'LineWidth', 1);
    plot(tbin, ydat, '-r', 'LineWidth', 2);
    
    % -- 2) p value at each time bin
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    ylabel('log10(p)');
    
    tmp = datshuff>datreal;
    p = log10(sum(tmp,1)./size(tmp,1)+1/nshuffs);
    plot(tbin, p, '-k');
    line(xlim, log10([0.05 0.05]));
    
    
    % ===== 1) overlay distribution of data and shuffle
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('diff');
    datreal = Ydiff_dat;
    datshuff = Ydiff_shuff;
    % --- run
    yCI = prctile(datshuff, [2.5 97.5], 1);
    ymean = mean(datshuff,1);
    ydat = datreal;
    
    plot(tbin, mean(ymean,1), 'k-', 'LineWidth', 2);
    plot(tbin, yCI(1,:), 'k-', 'LineWidth', 1);
    plot(tbin, yCI(2,:), 'k-', 'LineWidth', 1);
    plot(tbin, ydat, '-r', 'LineWidth', 2);
    
    % -- 2) p value at each time bin
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    ylabel('log10(p)');
    
    tmp = datshuff>datreal;
    p = log10(sum(tmp,1)./size(tmp,1)+1/nshuffs);
    plot(tbin, p, '-k');
    line(xlim, log10([0.05 0.05]));
end
end
