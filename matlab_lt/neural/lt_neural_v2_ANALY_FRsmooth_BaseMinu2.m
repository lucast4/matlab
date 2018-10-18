function lt_neural_v2_ANALY_FRsmooth_BaseMinu2(OUTDAT, SwitchStruct, ...
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

% ============= SAVE WITHIN TIME WINDOW, VALUE FOR ALL YTARG
[Ytimewindmetric_DAT, allgrpvals_b_e_s_n] = ...
    fn_ExtrWithinTimeWind(OUTDAT_tmp, analytype, premotorwind, dattype);


%% ===== also compare to shuffle analysis?
if doShuff==1

    % ====== many shuffles,
    % each shuffle, permute motif type
    
    Ytimewindmetric_SHUFF = [];
    for n=1:nshuffs
        disp(['shuff ' num2str(n)]);
        
        % ------ shuffle, within each neuron, all motifs
        syltypestoshuffle = [1 0 1]; % same and diff.
%         syltypestoshuffle = []; % same and diff.
        [OUTDAT_tmp] = lt_neural_v2_ANALY_FRsmooth_Comps(OUTDAT, SwitchStruct, 1, ...
            epochtoplot, 0, 0, syltypestoshuffle);
        ytimewindmetric = fn_ExtrWithinTimeWind(OUTDAT_tmp, analytype, premotorwind, dattype);
        Ytimewindmetric_SHUFF = [Ytimewindmetric_SHUFF; ytimewindmetric'];
    end
    
    
    %% =========== PLOTS
    % ======================= Z-SCORE DATA VS. SHUFFLE
    yshuffmean = mean(Ytimewindmetric_SHUFF,1);
    yshuffstd = std(Ytimewindmetric_SHUFF, [], 1);
    Yzscore = (Ytimewindmetric_DAT - yshuffmean')./yshuffstd';
        
    lt_figure; hold on;
    % =========== 1) DATA VS SHUFFLE
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
    
    % ==================== zscored
    lt_figure; hold on;
    lt_subplot(2,1,1); hold on;
    title('');
    xlabel('bird, expt, sw, neur');
    ylabel('scalar, within time wind');
    x = 1:length(Yzscore);
    plot(x, Yzscore, 'ok');
    lt_plot_zeroline;
    
    set(gca, 'XTick', 1:size(Ytimewindmetric_SHUFF,2));
    set(gca, 'XTickLabel', num2str(allgrpvals_b_e_s_n));
    rotateXLabels(gca, 90)
    % == z-scored
    lt_figure; hold on;
    
    % ============ 
    lt_subplot(2,3,4); hold on;
    title('histograms of z-scores')
    xlabel('z-score');
    xcenters = min(Yzscore)-0.1:0.1:max(Yzscore)+0.1;
    nbirds = max(allgrpvals_b_e_s_n(:,1));
    pcols = lt_make_plot_colors(nbirds, 0, 0);
    for j=1:nbirds
    indstmp = allgrpvals_b_e_s_n(:,1)==j;
    lt_plot_histogram(Yzscore(indstmp), xcenters, 1, 0, '', 1, pcols{j});
    end
    lt_plot_zeroline_vert;
    
    
    % ============ 
    lt_subplot(2,3,5); hold on;
    title('mean for each bird, relative own shuffle')
    xlabel('bird');
    xcenters = min(Yzscore)-0.1:0.1:max(Yzscore)+0.1;
    nbirds = max(allgrpvals_b_e_s_n(:,1));
    pcols = lt_make_plot_colors(nbirds, 0, 0);
    Yallshuff = {};
    Yall = [];
    for j=1:nbirds
    indstmp = allgrpvals_b_e_s_n(:,1)==j;
    
    ydat = Ytimewindmetric_DAT(indstmp);
    yshuff = Ytimewindmetric_SHUFF(:, indstmp);
    
    ymean = mean(ydat);
    ymean_shuff = mean(yshuff,2);
    
    lt_plot_MultDist({ymean_shuff}, j, 0, pcols{j}, 1, 0, 0);
    line([j-0.4 j+0.4], [ymean ymean], 'LineWidth', 3, 'Color', 'k');
%     Yall = [Yall; ymean];
%     Yallshuff = [Yallshuff; ymean_shuff];

    p = (sum(ymean_shuff>ymean)+1)./(length(ymean_shuff)+1);
    lt_plot_text(j, max(ymean_shuff), ['p=' num2str(p)], 'k');
    end
    
    
        
    
    
end




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
