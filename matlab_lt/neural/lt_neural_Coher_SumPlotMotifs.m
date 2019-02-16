function lt_neural_Coher_SumPlotMotifs(OUTSTRUCT, SwitchStruct,MOTIFSTATS_Compiled,PARAMS, indtoget_b_e_s)

%% lt 12/17/18 - plots cioherence scalar change during laerning, relative to posotion of motif
% i.e. for all motifs plot change in coh scalar.
clim = [-0.2 0.2];

% %% for each case get wn minus base, mean coh scalar
% cohscal_diff = [];
% for i=1:length(OUTSTRUCT.bnum)
%
%     indsbase = OUTSTRUCT.indsbase_epoch{i};
%     indswn = OUTSTRUCT.indsWN_epoch{i};
%
%     cohscal = OUTSTRUCT.cohscal{i};
%
%     cohdiff = mean(cohscal(indswn)) - mean(cohscal(indsbase));
%
%     cohscal_diff = [cohscal_diff; cohdiff];
% end
%
% OUTSTRUCT.cohscal_diff = cohscal_diff;

%% ############ ONE PLOT FOR EACH SWITCH, SHOW ALL MOTIFS.

[indsgrp_switch, indsgrp_switch_unique] = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch});
[indsgrp_chanpair, indsgrp_chanpair_unique] = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch, ...
    OUTSTRUCT.chanpair});

figcount=1;
subplotrows=6;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% indsgrp_switch_unique = unique(indsgrp_switch);
% indsgrp_chanpair_unique = unique(indsgrp_chanpair);

Yall = []; % targ, same, diff, nontarg (means across chans, motifs)
All_bname ={};
All_bnum = [];
All_ename = {};
All_swnum = [];
All_learndir = [];

for i=1:length(indsgrp_switch_unique)
    swgrpthis = indsgrp_switch_unique(i);
    bnum = unique(OUTSTRUCT.bnum(indsgrp_switch==swgrpthis));
    enum = unique(OUTSTRUCT.enum(indsgrp_switch==swgrpthis));
    swnum = unique(OUTSTRUCT.switch(indsgrp_switch==swgrpthis));
    bname = MOTIFSTATS_Compiled.birds(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(swnum)]);
    ylabel(['WN - base, coh']);
    
    % ====== plot each channel pair its own line
    for chanpair = indsgrp_chanpair_unique'
        indsthis = indsgrp_switch==swgrpthis & indsgrp_chanpair==chanpair;
        if ~any(indsthis)
            continue
        end
        motifID = lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)); % ---- get positions within global motif
        chnums = OUTSTRUCT.chanpair(indsthis,:);
        chnums = unique(chnums)'; assert(length(chnums)==2);
        cohscal = OUTSTRUCT.cohscal_diff(indsthis);
        %         cohscal = OUTSTRUCT.CohMean_WNminusBase_scalar(indsthis);
        
        assert(length(unique(motifID)) == length(motifID), 'each motif shoudl have at most 1 datapoint');
        % ---- sort in order of motifs
        [~, indsort] = sort(motifID);
        motifID = motifID(indsort);
        cohscal = cohscal(indsort);
        plot(motifID, cohscal, 'o-k');
        lt_plot_text(motifID(end)+0.3, cohscal(end), num2str(chnums), 'm', 9);
    end
    lt_plot_zeroline;
    
    % ====== overall
    indsthis = indsgrp_switch==swgrpthis;
    
    istarg = OUTSTRUCT.istarg(indsthis);
    issame = OUTSTRUCT.issame(indsthis);
    motifs = OUTSTRUCT.motifname(indsthis);
    motifID = lt_neural_QUICK_MotifID(bname, motifs); % ---- get positions within global motif
    %     cohscal = OUTSTRUCT.CohMean_WNminusBase_scalar(indsthis);
    cohscal = OUTSTRUCT.cohscal_diff(indsthis);
    learndir = unique(OUTSTRUCT.learndirTarg(indsthis));
    assert(length(learndir)==1);
    
    
    [ymean, ysem] = grpstats(cohscal, motifID, {'mean', 'sem'});
    [istarg] = grpstats(istarg, motifID, {'mean'});
    [issame] = grpstats(issame, motifID, {'mean'});    
    
    x = unique(motifID);
    lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
    
    % -------- NOTE DOWN POSITION OF TARGET SYSL
    indsthis = indsgrp_switch==swgrpthis & OUTSTRUCT.istarg==1;
    xtarg = unique(lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)));
    plot(xtarg, clim(1)+0.02, '^r');
    
    % ------- NOTE POSITION OF SAME_TYPES
    indsthis = indsgrp_switch==swgrpthis & OUTSTRUCT.istarg==0 & OUTSTRUCT.issame==1;
    xtarg = unique(lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)));
    if ~isempty(xtarg)
        plot(xtarg, clim(1)+0.02, '^b');
    end
    
    
    
    % ----- labels
    [~, motiflabels] = lt_neural_QUICK_MotifID(bname);
    set(gca, 'XTick', 1:length(motiflabels));
    set(gca, 'XTickLabel', motiflabels);
    rotateXLabels(gca, 90);
    ylim(clim);
    
    
    % ========= NOTE DOWN IF THIS IS AN EXPERIMENT TO GET
    if ~isempty(indtoget_b_e_s)
        if ismember([bnum enum swnum], indtoget_b_e_s, 'rows')
            lt_plot_annotation(3, 'expttoget', 'm');
            
            % ============= ALSO COLLECT DATA FOR SUMMARY PLOT
            Y = nan(1,4);
            
            % - targ
            Y(1) = nanmean(ymean(istarg==1));
            
            % - same
            Y(2) = nanmean(ymean(istarg==0 & issame==1));
            
            % - diff
            Y(3) = nanmean(ymean(istarg==0 & issame==0));
            
            % -- nontarg
            Y(4) = nanmean(ymean(istarg==0));
            
            Yall = [Yall; Y];
            
            All_bname = [All_bname; bname];
            All_ename = [All_ename; ename];
            All_swnum = [All_swnum; swnum];
            All_bnum = [All_bnum; bnum];
            All_learndir = [All_learndir; learndir];
        end
    else
        % ============= ALSO COLLECT DATA FOR SUMMARY PLOT
        Y = nan(1,4);
        
        % - targ
        Y(1) = nanmean(ymean(istarg==1));
        
        % - same
        Y(2) = nanmean(ymean(istarg==0 & issame==1));
        
        % - diff
        Y(3) = nanmean(ymean(istarg==0 & issame==0));
        
        % -- nontarg
        Y(4) = nanmean(ymean(istarg==0));
        
        Yall = [Yall; Y];
        
        All_bname = [All_bname; bname];
        All_ename = [All_ename; ename];
        All_swnum = [All_swnum; swnum];
        All_bnum = [All_bnum; bnum];
        All_learndir = [All_learndir; learndir];
    end
    
end


%% ========= SUMMARY ACROSS EXPERIMENTS

lt_figure; hold on;
pcols = lt_make_plot_colors(max(All_bnum), 0,0);

% ===== only those with all 3
lt_subplot(3,2,1); hold on;
xlabel('TARG - SAME - DIFF');
ylabel('coh (WN - base)');
colsthis = 1:3;

indsthis = all(~isnan(Yall(:,colsthis))');
ythis = Yall(indsthis,colsthis);
xthis = colsthis(1):colsthis(end);
plot(xthis, ythis', 'o-k');
xlim([0 5]);
lt_plot(colsthis+0.1, mean(ythis,1), {'Errors', lt_sem(ythis), 'Color', 'r'});
    [p] = signrank(ythis(:,1), ythis(:,2));
    lt_plot_pvalue(p, 'srank(1vs2)', 1);
bnumthis = All_bnum(indsthis);
for i=1:size(ythis,1)
   plot(xthis, ythis(i,:), '-o', 'Color', pcols{bnumthis(i)}); 
end
% -- plot each bird
[ymean, ysem] = grpstats(ythis, bnumthis, {'mean', 'sem'});
bnumtmp = unique(bnumthis);
for i=1:size(ymean,1)
    lt_plot(xthis+0.6*rand-0.3, ymean(i,:), {'Errors', ysem(i,:), 'LineStyle', '-'});
    lt_plot_text(xthis(end)+1, ymean(i, end), ['bnum ' num2str(bnumtmp(i))]);
end
xlim([0 5]);


% ===== only those with 2 (targ diff)
lt_subplot(3,2,2); hold on;
xlabel('TARG - DIFF');
ylabel('coh (WN - base)');

colsthis = [1 3];

indsthis = all(~isnan(Yall(:,colsthis))');
ythis = Yall(indsthis,colsthis);
xthis = colsthis;
plot(xthis, ythis', 'o-k');
xlim([0 5]);
lt_plot(colsthis+0.1, mean(ythis,1), {'Errors', lt_sem(ythis), 'Color', 'r'});
if length(colsthis)==2
    [p] = signrank(ythis(:,1), ythis(:,2));
    lt_plot_pvalue(p, 'srank', 1);
end
bnumthis = All_bnum(indsthis);
for i=1:size(ythis,1)
   plot(xthis, ythis(i,:), '-o', 'Color', pcols{bnumthis(i)}); 
end


% ===== 
lt_subplot(3,2,3); hold on;
xlabel('TARG - NONTARG');
ylabel('coh (WN - base)');

colsthis = [1 4];

indsthis = all(~isnan(Yall(:,colsthis))');
ythis = Yall(indsthis,colsthis);
xthis = colsthis;
plot(xthis, ythis', 'o-k');
xlim([0 5]);
lt_plot(colsthis+0.1, mean(ythis,1), {'Errors', lt_sem(ythis), 'Color', 'r'});
if length(colsthis)==2
%     [~, p] = ttest(ythis(:,1), ythis(:,2));
    [p] = signrank(ythis(:,1), ythis(:,2));
    lt_plot_pvalue(p, 'srank', 1);
end
bnumthis = All_bnum(indsthis);
for i=1:size(ythis,1)
   plot(xthis, ythis(i,:), '-o', 'Color', pcols{bnumthis(i)}); 
end


% ===== 
lt_subplot(3,2,4); hold on;
xlabel('TRAIN UP/DOWN [lines: TARG--NONTARG]');
ylabel('coh (WN - base)');

% ==== UP
indsthis = all(~isnan(Yall(:,colsthis))')' & All_learndir==1;
XSHIFT = 0;

% ---
colsthis = [1 4];

ythis = Yall(indsthis,colsthis);
xthis = colsthis + XSHIFT;
plot(xthis, ythis', 'o-k');
xlim([0 5+XSHIFT]);
if size(ythis,1)>1
lt_plot(xthis+0.1, mean(ythis,1), {'Errors', lt_sem(ythis), 'Color', 'r'});
end
if length(colsthis)==2
%     [~, p] = ttest(ythis(:,1), ythis(:,2));
    [p] = signrank(ythis(:,1), ythis(:,2));
    lt_plot_pvalue(p, 'srank', 1);
end
bnumthis = All_bnum(indsthis);
for i=1:size(ythis,1)
   plot(xthis, ythis(i,:), '-o', 'Color', pcols{bnumthis(i)}); 
end

% === DOWN
indsthis = all(~isnan(Yall(:,colsthis))')' & All_learndir==-1;
XSHIFT = 6;

% ---
colsthis = [1 4];

ythis = Yall(indsthis,colsthis);
xthis = colsthis + XSHIFT;
plot(xthis, ythis', 'o-k');
xlim([0 5+XSHIFT]);
if size(ythis,1)>1
lt_plot(xthis+0.1, mean(ythis,1), {'Errors', lt_sem(ythis), 'Color', 'r'});
end
if length(colsthis)==2
%     [~, p] = ttest(ythis(:,1), ythis(:,2));
    [p] = signrank(ythis(:,1), ythis(:,2));
    lt_plot_pvalue(p, 'srank', 1);
end
bnumthis = All_bnum(indsthis);
for i=1:size(ythis,1)
   plot(xthis, ythis(i,:), '-o', 'Color', pcols{bnumthis(i)}); 
end


% ======== EACH BIRD ONE DOT
lt_subplot(3,2,5); hold on;
title('dat = bird');
xlabel('TARG - DIFF');
ylabel('coh (WN - base)');

colsthis = [1 3];

indsthis = all(~isnan(Yall(:,colsthis))');
ythis = Yall(indsthis,colsthis);
xthis = colsthis;
bnumthis = All_bnum(indsthis);

[ymean, ysem] = grpstats(ythis, bnumthis, {'mean', 'sem'});
bnumtmp = unique(bnumthis);
for i=1:size(ymean,1)
    lt_plot(xthis+0.6*rand-0.3, ymean(i,:), {'Errors', ysem(i,:), 'LineStyle', '-'});
    lt_plot_text(xthis(end)+1, ymean(i, end), ['bnum ' num2str(i)]);
end
xlim([0 5]);



%% ========= OUTPUT, FORMATED FOR SAVING


