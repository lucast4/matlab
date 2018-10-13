function lt_neural_v2_ANALY_FRsmooth_BaseMinu(OUTDAT, SwitchStruct, ...
    epochtoplot, analytype, doShuff, syltypesneeded, premotorwind, nshuffs)
%%

dattype = 'neuron'; % for breaking down shuffle into lower level data.
% dattype = , bird, 'switch', 'expt', 'neuron'
%     nshuffs = 500;

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

% ============ 1) DUBTRACT BASE, 2) SUBTRACT MEAN OF DIFF, 3) TAKE ABSOLUTE
% VALUE
tmp = cellfun(@transpose, OUTDAT_tmp.AllDevDiff_NotAbs, 'UniformOutput', 0);
Ytarg_v2_dat = mean(abs(cell2mat(tmp(OUTDAT.All_istarg==1))),1);

tmp = cellfun(@transpose, OUTDAT_tmp.AllDevAll_NotAbs, 'UniformOutput', 0);
Ytarg_v2_minusall_dat = mean(abs(cell2mat(tmp(OUTDAT.All_istarg==1))),1);


% ============= SAVE WITHIN TIME WINDOW, VALUE FOR ALL YTARG
[Ytimewindmetric_DAT, allgrpvals_b_e_s_n] = ...
    fn_ExtrWithinTimeWind(OUTDAT_tmp, analytype, premotorwind, dattype);

%% ======================== PLOT ACROSS EXPERIMENTS THE ABSOLUTE DEVIATION FROM BASELINE
% i.e. use scalar
% NOTE: This is absolute deviation from baseline. Then compare targ and
% nontarg. This DIFFERS from the other analyses here (shuffle below), or
% 'AllOnlyMinusDiff_FRsmooth', which is absolute deviation from the mean of
% different type (after subtract baseline). 
% NOTE: I am not using this becuase it does not capture the diff between
% targ and nontarg the same way the analyses above and below do.


indsgrp = lt_tools_grp2idx({OUTDAT.All_birdnum, OUTDAT.All_exptnum, OUTDAT.All_swnum, OUTDAT.All_neurnum});
% indsgrp = lt_tools_grp2idx({OUTDAT.All_birdnum, OUTDAT.All_exptnum, OUTDAT.All_swnum});
% ---- collect all absolute deviations from baseline
tbin = OUTDAT.All_FRsmooth_t{1};
indstime = tbin>premotorwind(1) & tbin<premotorwind(2);

tmp = cell2mat(OUTDAT.AllOnlyMinusBase_FRsmooth); % only minus base, no abs.
AllDiffBase_abs = mean(abs(tmp(:, indstime)),2);

Yallmat = [];
Xallmat = [];
% --- target
indsthis = OUTDAT.All_istarg==1;
ytmp = grpstats(AllDiffBase_abs(indsthis), indsgrp(indsthis), {'mean'});
grpids = unique(indsgrp(indsthis));

Yallmat = [Yallmat; ytmp'];
Xallmat = [Xallmat; grpids'];
% --- diff type
indsthis = OUTDAT.All_istarg==0 & OUTDAT.All_issame==0;
ytmp = grpstats(AllDiffBase_abs(indsthis), indsgrp(indsthis), {'mean'});
grpids = unique(indsgrp(indsthis));

Yallmat = [Yallmat; ytmp'];
Xallmat = [Xallmat; grpids'];

% ---- check they are aligned
tmp = diff(Xallmat);
assert(all(tmp(:)==0));

lt_figure; hold on;
lt_subplot(2,2,1); hold on;
title('data = neurons');
xlabel('targ -- diff');
ylabel('abs diff from base, prem wind (hz)');
x = [1 2];
y = Yallmat;
plot(x, y, '-k');
lt_plot(x, mean(y,2)', {'Errors', lt_sem(y'), 'color', 'r'});

% ===== overlay entire distributions over all motifs
lt_subplot(2,2,4); hold on;
title('targ(k) and diff(r)');
ytarg = AllDiffBase_abs(OUTDAT.All_istarg==1);
ydiff = AllDiffBase_abs(OUTDAT.All_istarg==0 & OUTDAT.All_issame==0);

xcenters = min(AllDiffBase_abs)-1:5:max(AllDiffBase_abs)+1;
lt_plot_histogram(ytarg, xcenters, 1, 1, [], 1, 'k');
lt_plot_histogram(ydiff, xcenters, 1, 1, [], 1, 'r');

%% ===== also compare to shuffle analysis?
if doShuff==1
        tbin = OUTDAT.All_FRsmooth_t{1};
        indstime = tbin>premotorwind(1) & tbin<premotorwind(2);
    % ====== many shuffles,
    % each shuffle, permute motif type
    
    Ytarg_shuff = [];
    Ysame_shuff = [];
    Ydiff_shuff = [];
%     Ywithinwind_all_shuff = [];
    Ytarg_v2_shuffle = [];
    Ytarg_v2_minusall_shuffle = [];
    Ytimewindmetric_SHUFF = [];
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
        
        
        tmp = cellfun(@transpose, OUTDAT_tmp.AllDevDiff_NotAbs, 'UniformOutput', 0);
        Ytarg_v2_shuffle = [Ytarg_v2_shuffle; mean(abs(cell2mat(tmp(OUTDAT_tmp.All_istarg==1))),1)];

        tmp = cellfun(@transpose, OUTDAT_tmp.AllDevAll_NotAbs, 'UniformOutput', 0);
        Ytarg_v2_minusall_shuffle = [Ytarg_v2_minusall_shuffle; mean(abs(cell2mat(tmp(OUTDAT_tmp.All_istarg==1))),1)];
        
        %     Ytarg_shuff = [Ytarg_shuff; mean(cell2mat(OUTDAT_tmp.AllMinusAll_FRsmooth(OUTDAT_tmp.All_istarg==1)),1)];
        %     Ysame_shuff = [Ysame_shuff; mean(cell2mat(OUTDAT_tmp.AllMinusAll_FRsmooth(OUTDAT_tmp.All_istarg==0 & OUTDAT_tmp.All_issame==1)),1)];
        %     Ydiff_shuff = [Ydiff_shuff; mean(cell2mat(OUTDAT_tmp.AllMinusAll_FRsmooth(OUTDAT_tmp.All_istarg==0 & OUTDAT_tmp.All_issame==0)),1)];
        
        % =================== GET WITHIN TIMEWIND VALUES FOR ALL SYLS
%         Ywithinwind_all = cell2mat(OUTDAT_tmp.(analytype));
%         Ywithinwind_all = mean(Ywithinwind_all(:,indstime),2); % get mean over time bins.
%         Ywithinwind_all_shuff = [Ywithinwind_all_shuff; Ywithinwind_all'];
        ytimewindmetric = fn_ExtrWithinTimeWind(OUTDAT_tmp, analytype, premotorwind, dattype);
        Ytimewindmetric_SHUFF = [Ytimewindmetric_SHUFF; ytimewindmetric'];
    end
    
    
    % ################################ PLOT SHUFFLE OUTCOMES
    figcount=1;
    subplotrows=5;
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
    ylabel('absolute deviation from the mean of diff type (all after subtr base)');
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
    
    % -- line for premotor winodw
    line([premotorwind(1) premotorwind(1)], ylim, 'Color', 'b')
    line([premotorwind(2) premotorwind(2)], ylim, 'Color', 'b')
    
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
    
    
    
    
    
    % ===== 1) v2 (see above)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('TARG');
    ylabel('minus base->minus diff->abs');
    datreal = Ytarg_v2_dat;
    datshuff = Ytarg_v2_shuffle;
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
    
    
    
    
    
    
    % ===== 1) v2 (see above)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('TARG');
    ylabel('minus base->minus ALL->abs');
    datreal = Ytarg_v2_minusall_dat;
    datshuff = Ytarg_v2_minusall_shuffle;
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
    
    
    
    
    
    % ############################# PLOTS USING PREMOTOR WINDOW
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('TARG');
    xlabel('minus base, minus diff, abs');
    indstime = tbin>premotorwind(1) & tbin<premotorwind(2);
    ydat = mean(Ytarg_dat(indstime));
    yshuff = mean(Ytarg_shuff(:, indstime),2);
    lt_plot_shuffresult(ydat, yshuff);
    
    
    % ######################## USING PREMOTOR WINODW, BUT PLOT SEPARATELY
    % FOR EACH SWITCH
    lt_figure; hold on;
    title('');
    xlabel('bird, expt, sw, neur');
    ylabel('scalar, within time wind');
    
    for j=1:size(Ytimewindmetric_SHUFF,2) % each ind group (e.g. switches)
    % ------ 1) plot all shuffles
    % use mult dist plot
    yshuff = Ytimewindmetric_SHUFF(:, j);
    lt_plot_MultDist({yshuff}, j, 0, [0.6 0.6 0.6], 1, 0, 0)
    
    % ------ 2 ) overlay all data
    ydat = Ytimewindmetric_DAT(j);
    line([j-0.35 j+0.35], [ydat ydat], 'LineWidth', 3)
    end
    
    set(gca, 'XTick', 1:size(Ytimewindmetric_SHUFF,2));
    set(gca, 'XTickLabel', num2str(allgrpvals_b_e_s_n));
    rotateXLabels(gca, 90)
end



%% ###################### FOR EACH EXPERIMENT, 

tmp = cell2mat(OUTDAT.AllDevDiff_NotAbs(:))

end

function [ytimewindmetric, allgrpvals_b_e_s_n] = fn_ExtrWithinTimeWind(OUTDAT_tmp, analytype, ...
    premotorwind, dattype)

% dattype = 'switch', 'expt', 'neuron'

% ytimewindmetric = scalars, one for each switch (averaged across neurona nd 
% targets. takes average of analytype
% allgrpvals_b_e_s_n = nunique groups x (bird, expt, switch, neuron)
%%
% ---- 1) convert all smoothed fr diff to one scalar in time window
Ywithinwind_all_DAT = cell2mat(OUTDAT_tmp.(analytype));
tbin = OUTDAT_tmp.All_FRsmooth_t{1};
indstime = tbin>premotorwind(1) & tbin<premotorwind(2);
Ywithinwind_all_DAT = mean(Ywithinwind_all_DAT(:,indstime),2); % get mean over time bins.

% ---- 2) extract one value for each switch
if strcmp(dattype, 'bird')
indsgrp = lt_tools_grp2idx({OUTDAT_tmp.All_birdnum});    
elseif strcmp(dattype, 'switch')
indsgrp = lt_tools_grp2idx({OUTDAT_tmp.All_birdnum, OUTDAT_tmp.All_exptnum, OUTDAT_tmp.All_swnum});
elseif strcmp(dattype, 'expt')
    indsgrp = lt_tools_grp2idx({OUTDAT_tmp.All_birdnum, OUTDAT_tmp.All_exptnum});
elseif strcmp(dattype, 'neuron')
    indsgrp = lt_tools_grp2idx({OUTDAT_tmp.All_birdnum, OUTDAT_tmp.All_exptnum, OUTDAT_tmp.All_swnum, OUTDAT_tmp.All_neurnum});
end

% === create legend, one val for each unique indsgrp
[~, ~, indstmp] = intersect(unique(indsgrp), indsgrp);
allgrpvals_b_e_s_n = [OUTDAT_tmp.All_birdnum(indstmp) OUTDAT_tmp.All_exptnum(indstmp) OUTDAT_tmp.All_swnum(indstmp) OUTDAT_tmp.All_neurnum(indstmp)];

% === FOR EACH IND IN GROUP, EXTRACT mean time window diff over all
% neurons and target motifs (for this switch)
ytimewindmetric = [];
for j=unique(indsgrp)'
    indthis = indsgrp==j & OUTDAT_tmp.All_istarg==1; % all targets for this
    ytimewindmetric = [ytimewindmetric; mean(Ywithinwind_all_DAT(indthis))];
end
end
