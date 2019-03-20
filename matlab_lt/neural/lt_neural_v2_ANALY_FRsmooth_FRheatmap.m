function lt_neural_v2_ANALY_FRsmooth_FRheatmap(OUTDAT, SwitchStruct, ...
    epochtoplot, analytype, syltypesneeded, minmotifs, prewind)
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


%% =================== ONLY KEEP SWITCHES THAT HAVE A MINIMUM NUMBER OF MOTIFS

[indsgrp, indsgrpU] = lt_tools_grp2idx({OUTDAT.All_birdnum, OUTDAT.All_exptnum, OUTDAT.All_swnum, OUTDAT.All_neurnum});
indstoremove = [];
for i=indsgrpU'
    
    indsthis = indsgrp==i;
    
    nmot = length(unique(OUTDAT.All_motifnum(indsthis)));
    if nmot < minmotifs
        indstoremove = [indstoremove; find(indsthis)];
    end
    
end

indstokeep = ~ismember(1:length(OUTDAT.All_birdnum), indstoremove');
disp(['Keeping ' num2str(sum(indstokeep)) '/' num2str(length(indstokeep)) 'neurons (NOT ENOUGHT MOTIFS)']);
OUTDAT = lt_structure_subsample_all_fields(OUTDAT, indstokeep, 1);


%% =================== PLOT HEAT MAPS

[indsgrp, indsgrpU] = lt_tools_grp2idx({OUTDAT.All_birdnum, ...
    OUTDAT.All_exptnum, OUTDAT.All_swnum, OUTDAT.All_neurnum});
figcount=1;
subplotrows=5;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


% ==== t indices to get scalar
t = OUTDAT.All_FRsmooth_t{1};
indst = t>prewind(1) & t<prewind(2);

for i=indsgrpU'
    
    indsthis = indsgrp==i;
    
    Yall = OUTDAT.(analytype)(indsthis);
    t = OUTDAT.All_FRsmooth_t(indsthis); t = t{1};
    istarg = OUTDAT.All_istarg(indsthis);
    issame = OUTDAT.All_issame(indsthis);
    
    motifnum = OUTDAT.All_motifnum(indsthis);
    bnum = unique(OUTDAT.All_birdnum(indsthis));
    enum = unique(OUTDAT.All_exptnum(indsthis));
    sw = unique(OUTDAT.All_swnum(indsthis));
    neur = unique(OUTDAT.All_neurnum(indsthis));
    motifnames = {SwitchStruct.bird(bnum).exptnum(enum).switchlist(sw).STATS_motifsyl(motifnum).sylname};
    bname = SwitchStruct.bird(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    
    Yall = cell2mat(Yall);
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename]);
    ylabel(['-sw' num2str(sw) '-n' num2str(neur)]);
    if strcmp(analytype, 'AllMinusBase_FRmeanAll_epoch')
        clim = [-prctile(abs(Yall(:)), [97.5]) prctile(abs(Yall(:)), [97.5])];
        cmap = 'centered';
    else
        clim = prctile(Yall(:), [2.5 97.5]);
        cmap = 'pval';
    end
    y = 1:size(Yall,1);
    lt_neural_Coher_Plot(Yall', t, y, 1, '', clim);
    lt_plot_colormap(cmap);
    set(gca, 'YTick', 1:size(Yall,1), 'YTickLabel', motifnames);
    % --- mark target, same type
    XLIM = xlim;
    lt_plot(XLIM(1), y(istarg==1), {'Color', 'r'});
    if any(istarg==0 & issame==1)
        lt_plot(XLIM(1), y(istarg==0 & issame==1), {'Color', 'b'});
    end
    xlabel(['heat=' analytype]);
    line([prewind(1) prewind(1)], ylim, 'Color', 'm');
    line([prewind(2) prewind(2)], ylim, 'Color', 'm');
    
    % ================= PLOT DISTRIBUTION OF SCALARS FROM PREMOTOR WINDOW
    Yscal = mean(Yall(:, indst),2);
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename]);
    ylabel(['-sw' num2str(sw) '-n' num2str(neur)]);
    x = 1:length(Yscal);
     
    plot(x(istarg==1), Yscal(istarg==1), 'or');
    plot(x(istarg==0 & issame==1), Yscal(istarg==0 & issame==1), 'ob');
    plot(x(istarg==0 & issame==0), Yscal(istarg==0 & issame==0), 'ok');
    set(gca, 'XTick', x, 'XTickLabel', motifnames);
    rotateXLabels(gca, 90);
end





%% === plot all expts in same plot
[indsgrp, indsgrpU] = lt_tools_grp2idx({OUTDAT.All_birdnum, ...
    OUTDAT.All_exptnum, OUTDAT.All_swnum, OUTDAT.All_neurnum});


lt_figure; hold on;

lt_subplot(2,1,2); hold on;
% ==== t indices to get scalar
t = OUTDAT.All_FRsmooth_t{1};
indst = t>prewind(1) & t<prewind(2);

Xlaball = {};
for i=1:length(indsgrpU')
    
    indsthis = indsgrp==indsgrpU(i);
    
    Yall = OUTDAT.(analytype)(indsthis);
    t = OUTDAT.All_FRsmooth_t(indsthis); t = t{1};
    istarg = OUTDAT.All_istarg(indsthis);
    issame = OUTDAT.All_issame(indsthis);
    
    motifnum = OUTDAT.All_motifnum(indsthis);
    bnum = unique(OUTDAT.All_birdnum(indsthis));
    enum = unique(OUTDAT.All_exptnum(indsthis));
    sw = unique(OUTDAT.All_swnum(indsthis));
    neur = unique(OUTDAT.All_neurnum(indsthis));
    motifnames = {SwitchStruct.bird(bnum).exptnum(enum).switchlist(sw).STATS_motifsyl(motifnum).sylname};
    bname = SwitchStruct.bird(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    
    Yall = cell2mat(Yall);
    
    Xlaball = [Xlaball; [bname(1:4) '-' ename(end-7:end) '-sw' num2str(sw) '-n' num2str(neur)]];
    
%     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%     title([bname '-' ename]);
%     ylabel(['-sw' num2str(sw) '-n' num2str(neur)]);
%     if strcmp(analytype, 'AllMinusBase_FRmeanAll_epoch')
%         clim = [-prctile(abs(Yall(:)), [97.5]) prctile(abs(Yall(:)), [97.5])];
%         cmap = 'centered';
%     else
%         clim = prctile(Yall(:), [2.5 97.5]);
%         cmap = 'pval';
%     end

    y = 1:size(Yall,1);
%     lt_neural_Coher_Plot(Yall', t, y, 1, '', clim);
%     lt_plot_colormap(cmap);
%     set(gca, 'YTick', 1:size(Yall,1), 'YTickLabel', motifnames);
    % --- mark target, same type
%     XLIM = xlim;
%     lt_plot(XLIM(1), y(istarg==1), {'Color', 'r'});
%     if any(istarg==0 & issame==1)
%         lt_plot(XLIM(1), y(istarg==0 & issame==1), {'Color', 'b'});
%     end
%     xlabel(['heat=' analytype]);
%     line([prewind(1) prewind(1)], ylim, 'Color', 'm');
%     line([prewind(2) prewind(2)], ylim, 'Color', 'm');
    
    % ================= PLOT DISTRIBUTION OF SCALARS FROM PREMOTOR WINDOW
    Yscal = mean(Yall(:, indst),2);
    
%     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%     title([bname '-' ename]);
%     ylabel(['-sw' num2str(sw) '-n' num2str(neur)]);
%     x = 1:length(Yscal);
     
    lt_plot(i, Yscal(istarg==1), {'Color', 'r'});
    try
        lt_plot(i, Yscal(istarg==0 & issame==1), {'Color', 'b'});
% plot(i, Yscal(istarg==0 & issame==1), 'ob');
    catch err
    end
    plot(i, Yscal(istarg==0 & issame==0), 'ok');
%     set(gca, 'XTick', x, 'XTickLabel', motifnames);
%     rotateXLabels(gca, 90);
end

set(gca, 'XTick', 1:length(indsgrpU'), 'XTickLabel', Xlaball);
rotateXLabels(gca, 90);


% ################################ PLOT MAGNITUDE OF LEARNING.
lt_subplot(2,1,1); hold on;
ylabel('learning (z), targ dir');

% ==== t indices to get scalar
t = OUTDAT.All_FRsmooth_t{1};
indst = t>prewind(1) & t<prewind(2);

Xlaball = {};
for i=1:length(indsgrpU')
    
    indsthis = indsgrp==indsgrpU(i);
    
    Yall = OUTDAT.(analytype)(indsthis);
    t = OUTDAT.All_FRsmooth_t(indsthis); t = t{1};
    istarg = OUTDAT.All_istarg(indsthis);
    issame = OUTDAT.All_issame(indsthis);
    learn = OUTDAT.AllMinusBase_PitchZ_Ldir(indsthis);
    
    motifnum = OUTDAT.All_motifnum(indsthis);
    bnum = unique(OUTDAT.All_birdnum(indsthis));
    enum = unique(OUTDAT.All_exptnum(indsthis));
    sw = unique(OUTDAT.All_swnum(indsthis));
    neur = unique(OUTDAT.All_neurnum(indsthis));
    motifnames = {SwitchStruct.bird(bnum).exptnum(enum).switchlist(sw).STATS_motifsyl(motifnum).sylname};
    bname = SwitchStruct.bird(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    
    Yall = cell2mat(Yall);
    
    Xlaball = [Xlaball; [bname(1:4) '-' ename(end-7:end) '-sw' num2str(sw) '-n' num2str(neur)]];

    y = 1:size(Yall,1);
    
    % ================= PLOT DISTRIBUTION OF SCALARS FROM PREMOTOR WINDOW
    Yscal = mean(Yall(:, indst),2);
    
    lt_plot(i, learn(istarg==1), {'Color', 'r'});
    try
        lt_plot(i, learn(istarg==0 & issame==1), {'Color', 'b'});
        % plot(i, Yscal(istarg==0 & issame==1), 'ob');
    catch err
    end
    plot(i, learn(istarg==0 & issame==0), 'ok');
    %     set(gca, 'XTick', x, 'XTickLabel', motifnames);
    %     rotateXLabels(gca, 90);
end
lt_plot_zeroline;

% set(gca, 'XTick', 1:length(indsgrpU'), 'XTickLabel', Xlaball);
% rotateXLabels(gca, 90);
