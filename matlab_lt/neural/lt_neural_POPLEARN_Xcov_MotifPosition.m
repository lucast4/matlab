function lt_neural_POPLEARN_XCov_PlotScal(Yscalar, OUTSTRUCT_XCOV, SwitchCohStruct, SwitchStruct, ...
    MOTIFSTATS_Compiled, PARAMS, onlygoodexpt, expttype)
%% lt 12/20/18 - plots cioherence scalar change during laerning, whether is same or diff motif.

%% 2/1/19 - Plots each scalar, each experiement...

% clim = [-0.1 0.1];
useglobalmotifname = 0; % then all expt across a bird will be aligned... motif names will not necessarily be correct for each experiement though;
minDiffN = 3; % will only noramlize to diff ytpe if diff type has N this or larger.
takeindsWNsecondhalf = 0; % then essentially takes last quartile to compute learning.
% LEAVE AT 0 - if 1, then inds will not necessariyl be exactly matched
% between xcov, coh, and laerning.


%% ======== filter experiments
if onlygoodexpt==1
    % ===== filter OUTSTRUCT_XCOV
    [OUTSTRUCT_XCOV, indstokeep] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT_XCOV, SwitchStruct, expttype);
    % ===== filter the scalar array
    Yscalar = Yscalar(indstokeep);
end

% ====== SANITY CHECK.
assert(length(Yscalar) == length(OUTSTRUCT_XCOV.bnum));


%% ============= 1) GET UNIQUE INDS FOR EACH SWITCH, AND EACH CHAN PAIR
[indsgrp_switch, indsgrp_switch_unique] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, ...
    OUTSTRUCT_XCOV.enum, OUTSTRUCT_XCOV.switch});

[indsgrp_chanpair, indsgrp_chanpair_unique] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, ...
    OUTSTRUCT_XCOV.enum, OUTSTRUCT_XCOV.switch, OUTSTRUCT_XCOV.chanpair});

%% NOTE
% if multipel targets on different motifs, then will skip that case. if on same motif
% then will take  the earliest one.
%% ############

All_istarg = [];
All_cohscal = [];
All_posreltarget = [];
All_swunique = [];
All_issame = [];
All_bnum = [];
All_enum = [];
All_channum = [];
for i=1:length(indsgrp_chanpair_unique)
    
    swgrpthis = indsgrp_chanpair_unique(i);
    
    indsthis = indsgrp_chanpair==swgrpthis;
    
    bnum = unique(OUTSTRUCT_XCOV.bnum(indsthis));
    enum = unique(OUTSTRUCT_XCOV.enum(indsthis));
    swnum = unique(OUTSTRUCT_XCOV.switch(indsthis));
    bname = MOTIFSTATS_Compiled.birds(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    %
    %     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    %     title([bname '-' ename '-sw' num2str(swnum)]);
    %     ylabel(['WN - base, coh']);
    
    %     % ====== plot each channel pair its own line
    %     for chanpair = indsgrp_chanpair_unique'
    %         indsthis = indsgrp_switch==swgrpthis & indsgrp_chanpair==chanpair;
    %         if ~any(indsthis)
    %             continue
    %         end
    %         motifID = lt_neural_QUICK_MotifID(bname, OUTSTRUCT_XCOV.motifname(indsthis)); % ---- get positions within global motif
    %
    %         cohscal = OUTSTRUCT_XCOV.cohscal_diff(indsthis);
    %         %         cohscal = OUTSTRUCT_XCOV.CohMean_WNminusBase_scalar(indsthis);
    %
    %         assert(length(unique(motifID)) == length(motifID), 'each motif shoudl have at most 1 datapoint');
    %         % ---- sort in order of motifs
    %         [~, indsort] = sort(motifID);
    %         motifID = motifID(indsort);
    %         cohscal = cohscal(indsort);
    %         plot(motifID, cohscal, 'o-k');
    %     end
    %     lt_plot_zeroline;
    
    % ====== overall
    
    istarg = OUTSTRUCT_XCOV.istarg(indsthis);
    isame = OUTSTRUCT_XCOV.issame(indsthis);
    motifs = OUTSTRUCT_XCOV.motifname(indsthis);
    [motifID, ~, motifposition] = lt_neural_QUICK_MotifID(bname, motifs); % ---- get positions within global motif
    %     cohscal = OUTSTRUCT_XCOV.CohMean_WNminusBase_scalar(indsthis);
    cohscal = Yscalar(indsthis);
    
    [ymean, ysem] = grpstats(cohscal, motifID, {'mean', 'sem'});
    
    %     x = unique(motifID);
    %     lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
    %
    %     % -------- NOTE DOWN POSITION OF TARGET SYSL
    %     indsthis = indsgrp_switch==swgrpthis & OUTSTRUCT_XCOV.istarg==1;
    %     xtarg = unique(lt_neural_QUICK_MotifID(bname, OUTSTRUCT_XCOV.motifname(indsthis)));
    %     plot(xtarg, clim(1)+0.02, '^r');
    %
    %     % ------- NOTE POSITION OF SAME_TYPES
    %     indsthis = indsgrp_switch==swgrpthis & OUTSTRUCT_XCOV.istarg==0 & OUTSTRUCT_XCOV.issame==1;
    %     xtarg = unique(lt_neural_QUICK_MotifID(bname, OUTSTRUCT_XCOV.motifname(indsthis)));
    %     if ~isempty(xtarg)
    %         plot(xtarg, clim(1)+0.02, '^b');
    %     end
    %
    %     % ----- labels
    %     [~, motiflabels] = lt_neural_QUICK_MotifID(bname);
    %     set(gca, 'XTick', 1:length(motiflabels));
    %     set(gca, 'XTickLabel', motiflabels);
    %     rotateXLabels(gca, 90);
    %     ylim(clim);
    
    % ============== OUTPUT, FOR EACH SYL GET POSOTION RELATIVE TO TARGET
    targmotif = motifposition(istarg==1,1);
    targposition = motifposition(istarg==1,2);
    
    % -- if there are multipel targets, then take the earliest tareget
    % (asserts that they are on same motif)
    if length(unique(lt_tools_grp2idx({targmotif targposition})))>1
        % then multpe targets, check that they are on same motif
        %         assert(length(unique(targmotif))==1, 'problem, targets are on diff motifs...');
        
        if length(unique(targmotif))==1
            % then good, there is only one motif, take earliset
        elseif length(unique(targmotif))>1
            % bad, multipel motifs - skip this case
            continue
        end
    end
    
    % --- earliest target
    targposition_good = [targmotif(1) min(targposition)];
    
    % ==== for each nontarget, compare to taregh
    posreltarget = nan(size(cohscal)); %
    
    indstmp = motifposition(:,1)==targposition_good(1); % cases on same motifa s targ
    posreltarget(indstmp) = motifposition(indstmp,2) - targposition_good(2);
    
    
    % ============= SAVE OUTUOT
    All_istarg = [All_istarg; istarg];
    All_issame = [All_issame; isame];
    All_cohscal = [All_cohscal; cohscal];
    All_posreltarget = [All_posreltarget; posreltarget];
    All_bnum = [All_bnum; bnum*ones(size(cohscal))];
    All_swunique = [All_swunique; i*ones(size(cohscal))];
    All_enum = [All_enum; enum*ones(size(cohscal))];
end

assert(max(All_posreltarget)<1, 'if there are some follolwing, then want to modify code to plot those. code would be fine if not, but would be intersting to look at.');

%% ============= [PUBLIATION] ADJACENT
lt_figure; hold on;


% ================ ALL
lt_subplot(2,3,1); hold on;
title('all syl (col=bird)');
xlabel('TARG -- ADJACENT (preceding) -- REST(same motif) -- diff motif');
ylabel('coh change');

pcolall = lt_make_plot_colors(max(unique(All_bnum)), 0,0);
Xall = [];
Yall =[];
Bnum =[];
Enum = [];
for i=unique(All_swunique)'
    
    indsthis = All_swunique==i;
    bnum = All_bnum(indsthis);
    bnum = unique(bnum); assert(length(bnum)==1);
    enum = All_enum(indsthis);
    enum = unique(enum); assert(length(enum)==1);    
    pos = All_posreltarget(indsthis);
    cohscal = All_cohscal(indsthis);
    issame = All_issame(indsthis);
    istarg = All_istarg(indsthis);
    
    % -- any nan, convert to large number
    pos(isnan(pos)) = min(All_posreltarget)-2;
    diffmotifID = min(All_posreltarget)-2;
    
    % --- average over chanels
    cohscal = grpstats(cohscal, pos);
    issame = grpstats(issame, pos);
    istarg = grpstats(istarg, pos);
    pos = unique(pos);
    
    % --- sort
    [~, indsort] = sort(pos);
    pos = pos(indsort);
    cohscal = cohscal(indsort);
    
    
    % === plot
    pcol = pcolall{unique(All_bnum(indsthis))};
    disp(unique(All_bnum(indsthis)));
    X = 1:4;
    Y = nan(1,4);
    
    % -- targ
    Y(1) = mean(cohscal(pos==0));
    % -- adjacnet (-preceding)
    Y(2) = mean(cohscal(pos==-1));
    % -- rest
    Y(3) = mean(cohscal(pos<-1 & pos>diffmotifID));
    % -- adjacent, folliwi
    % -- diff motif
    Y(4) = mean(cohscal(pos==diffmotifID));
    
    if any(pos>0)
        disp('NOTE: in some case ther eare syls folliowujng - ignoreing those for "adjavent');
    end
    %     assert(all(pos<1), 'assumed that all are precedeing...');
    
    plot(X+0.4*rand-0.2, Y, '-o', 'Color', pcol);
    
    Yall =[Yall; Y];
    Bnum = [Bnum; bnum]; 
    Enum = [Enum; enum];
end



% ====================== PLOTS

% lt_plot(X+0.3, nanmean(Yall,1), {'Errors', lt_sem(Yall), 'Color', 'k'});
xlim([0 X(end)+1]);
lt_plot_zeroline;
N = sum(~isnan(Yall));
lt_plot_annotation(1, num2str(N), 'm')

% -- overlay bird means
[Ymean, Ysem] = grpstats(Yall, Bnum, {"nanmean", @(x)lt_sem(x)});
for i=1:length(unique(Bnum))
    b = unique(Bnum);
    pcol = pcolall{b(i)};
    ym = Ymean(i,:);
    ys = Ysem(i,:);
    x = 1:length(ym);
    lt_plot(x+0.05*i, ym, {'Errors', ys, 'Color', pcol, 'LineStyle', '-'});
end

lt_subplot(2,3,2); hold on;
lt_plot_distributionPlot(Yall)
title("all syls")
% -- overlay bird means
[Ymean, Ysem] = grpstats(Yall, Bnum, {"nanmean", @(x)lt_sem(x)});
for i=1:length(unique(Bnum))
    b = unique(Bnum);
    pcol = pcolall{b(i)};
    ym = Ymean(i,:);
    ys = Ysem(i,:);
    x = 1:length(ym);
    lt_plot(x+0.05*i, ym, {'Errors', ys, 'Color', pcol, 'LineStyle', '-'});
end


lt_subplot(2,3,3); hold on;
title('all syl (only if has targ and adjacent)');
% xlabel('TARG -- ADJACENT (preceding) -- REST');
ylabel('coh change');

Bnum(any(isnan(Yall(:,[1 2]))')) = [];
Enum(any(isnan(Yall(:,[1 2]))')) = [];
Yall = Yall(:,[1 2]);
Yall(any(isnan(Yall(:,[1 2]))'),:) = [];
X = 1:size(Yall,2);
plot(X, Yall', '-k');
xlim([0 max(X)+1]);
% lt_plot(X+0.2, nanmean(Yall,1), {'Errors', lt_sem(Yall), 'Color', 'r'});
lt_plot_zeroline;
% -- overlay bird means
[Ymean, Ysem] = grpstats(Yall, Bnum, {"nanmean", @(x)lt_sem(x)});
for i=1:length(unique(Bnum))
    b = unique(Bnum);
    pcol = pcolall{b(i)};
    ym = Ymean(i,:);
    ys = Ysem(i,:);
    x = 1:length(ym);
    lt_plot(x+0.05*i, ym, {'Errors', ys, 'Color', pcol, 'LineStyle', '-'});
end

xlim([0 X(end)+1]);
lt_plot_zeroline;
N = sum(~isnan(Yall));
lt_plot_annotation(1, num2str(N), 'm')



lt_subplot(2,3,4); hold on;
lt_plot_distributionPlot(Yall)
title("all syls (only if have all cases)")
% -- overlay bird means
[Ymean, Ysem] = grpstats(Yall, Bnum, {"nanmean", @(x)lt_sem(x)});
for i=1:length(unique(Bnum))
    b = unique(Bnum);
    pcol = pcolall{b(i)};
    ym = Ymean(i,:);
    ys = Ysem(i,:);
    x = 1:length(ym);
    lt_plot(x+0.05*i, ym, {'Errors', ys, 'Color', pcol, 'LineStyle', '-'});
end



%% ==================== LME
% Yallcell = mat2cell(Yall, size(Yall,1), ones(size(Yall,2),1));
Yallvec = Yall(:);
Bnumvec = repmat(Bnum, 1, size(Yall,2));
Bnumvec = Bnumvec(:);
Enumvec = lt_tools_grp2idx({Bnum, Enum});
Enumvec = repmat(Enumvec, 1, size(Yall,2));
Enumvec = Enumvec(:);
Syltype = repmat(size(Yall,2)-1:-1:0, size(Yall,1), 1);
Syltype = Syltype(:);
dat = table(Yallvec, categorical(Bnumvec), categorical(Enumvec), ...
    Syltype, 'VariableNames', {'Yresponse', 'bnum', 'exptnum', 'syltype'});
% dat=dat(dat.syltype(:)==categorical(1) | dat.syltype(:)==categorical(2), :); % only care about target and adjacent
formula = 'Yresponse ~ syltype + (1|exptnum)';
formula = 'Yresponse ~ syltype';
lme = fitlme(dat, formula, 'StartMethod', 'random')
disp("LME: syltype 0 (intercept) is adjacents. syltype 1 is targ. only cases where they have both targ and adjacent.")
disp("N unique expts after only keeping cases with adjacent: ")
disp(length(unique(Enumvec)))

%% CONITNUED
pcolall = lt_make_plot_colors(max(unique(All_bnum)), 0,0);
Xall = [];
Yall =[];
Bnum =[];
Enum = [];
for i=unique(All_swunique)'
    
    indsthis = All_swunique==i;
    bnum = All_bnum(indsthis);
    bnum = unique(bnum); assert(length(bnum)==1);
    enum = All_enum(indsthis);
    enum = unique(enum); assert(length(enum)==1);    
    pos = All_posreltarget(indsthis);
    cohscal = All_cohscal(indsthis);
    issame = All_issame(indsthis);
    istarg = All_istarg(indsthis);
    
    % -- any nan, convert to large number
    pos(isnan(pos)) = min(All_posreltarget)-2;
    diffmotifID = min(All_posreltarget)-2;
    
    % --- average over chanels
    cohscal = grpstats(cohscal, pos);
    issame = grpstats(issame, pos);
    istarg = grpstats(istarg, pos);
    pos = unique(pos);
    
    % --- sort
    [~, indsort] = sort(pos);
    pos = pos(indsort);
    cohscal = cohscal(indsort);
    
    
    % === plot
    pcol = pcolall{unique(All_bnum(indsthis))};
    disp(unique(All_bnum(indsthis)));
    X = 1:4;
    Y = nan(1,4);
    
    % -- targ
    Y(1) = mean(cohscal(pos==0));
    % -- adjacnet (-preceding)
    Y(2) = mean(cohscal(pos==-1));
    % -- rest
    Y(3) = mean(cohscal(pos<-1 & pos>diffmotifID));
    % -- adjacent, folliwi
    % -- diff motif
    Y(4) = mean(cohscal(pos==diffmotifID));
    
    if any(pos>0)
        disp('NOTE: in some case ther eare syls folliowujng - ignoreing those for "adjavent');
    end
    %     assert(all(pos<1), 'assumed that all are precedeing...');
    
    plot(X+0.4*rand-0.2, Y, '-o', 'Color', pcol);
    
    Yall =[Yall; Y];
    Bnum = [Bnum; bnum]; 
    Enum = [Enum; enum];
end




lt_subplot(2,3,5); hold on;
title('all syl (only if has all cases)');
% xlabel('TARG -- ADJACENT (preceding) -- REST');
ylabel('coh change');

Bnum(any(isnan(Yall)')) = [];
Enum(any(isnan(Yall)')) = [];
Yall(any(isnan(Yall)'),:) = [];
plot(X, Yall', '-k');
xlim([0 max(X)+1]);
% lt_plot(X+0.2, nanmean(Yall,1), {'Errors', lt_sem(Yall), 'Color', 'r'});
lt_plot_zeroline;
% -- overlay bird means
[Ymean, Ysem] = grpstats(Yall, Bnum, {"nanmean", @(x)lt_sem(x)});
for i=1:length(unique(Bnum))
    b = unique(Bnum);
    pcol = pcolall{b(i)};
    ym = Ymean(i,:);
    ys = Ysem(i,:);
    x = 1:length(ym);
    lt_plot(x+0.05*i, ym, {'Errors', ys, 'Color', pcol, 'LineStyle', '-'});
end

xlim([0 X(end)+1]);
lt_plot_zeroline;
N = sum(~isnan(Yall));
lt_plot_annotation(1, num2str(N), 'm')


lt_subplot(2,3,6); hold on;
lt_plot_distributionPlot(Yall)
title("all syls (only if have all cases)")
% -- overlay bird means
[Ymean, Ysem] = grpstats(Yall, Bnum, {"nanmean", @(x)lt_sem(x)});
for i=1:length(unique(Bnum))
    b = unique(Bnum);
    pcol = pcolall{b(i)};
    ym = Ymean(i,:);
    ys = Ysem(i,:);
    x = 1:length(ym);
    lt_plot(x+0.05*i, ym, {'Errors', ys, 'Color', pcol, 'LineStyle', '-'});
end

%% LME MODEL?



%% =============== [PLOT] -- TARG, ADJACENT, NOT ADJACENT
lt_figure; hold on;


% ================ ALL
lt_subplot(2,2,1); hold on;
title('all syl (col=bird)');
xlabel('TARG -- ADJACENT (preceding) -- REST(same motif) -- diff motif');
ylabel('coh change');

pcolall = lt_make_plot_colors(max(unique(All_bnum)), 0,0);
Xall = [];
Yall =[];
for i=unique(All_swunique)'
    
    indsthis = All_swunique==i;
    
    pos = All_posreltarget(indsthis);
    cohscal = All_cohscal(indsthis);
    issame = All_issame(indsthis);
    istarg = All_istarg(indsthis);
    
    % -- any nan, convert to large number
    pos(isnan(pos)) = min(All_posreltarget)-2;
    diffmotifID = min(All_posreltarget)-2;
    
    % --- average over chanels
    cohscal = grpstats(cohscal, pos);
    issame = grpstats(issame, pos);
    istarg = grpstats(istarg, pos);
    pos = unique(pos);
    
    % --- sort
    [~, indsort] = sort(pos);
    pos = pos(indsort);
    cohscal = cohscal(indsort);
    
    
    % === plot
    pcol = pcolall{unique(All_bnum(indsthis))};
    disp(unique(All_bnum(indsthis)));
    X = 1:4;
    Y = nan(1,4);
    
    % -- targ
    Y(1) = mean(cohscal(pos==0));
    % -- adjacnet (-preceding)
    Y(2) = mean(cohscal(pos==-1));
    % -- rest
    Y(3) = mean(cohscal(pos<-1 & pos>diffmotifID));
    % -- adjacent, folliwi
    % -- diff motif
    Y(4) = mean(cohscal(pos==diffmotifID));
    
    if any(pos>0)
        disp('NOTE: in some case ther eare syls folliowujng - ignoreing those for "adjavent');
    end
    %     assert(all(pos<1), 'assumed that all are precedeing...');
    
    plot(X+0.4*rand-0.2, Y, '-o', 'Color', pcol);
    
    Yall =[Yall; Y];
end


lt_plot(X+0.3, nanmean(Yall,1), {'Errors', lt_sem(Yall), 'Color', 'k'});
xlim([0 X(end)+1]);
lt_plot_zeroline;
N = sum(~isnan(Yall));
lt_plot_annotation(1, num2str(N), 'm')


% =========== REPLOT, BUT USING DIFFERENCE (i.e paired structure included)
lt_subplot(2,2,2); hold on;
title('all syl (only if has all cases, cols1:3)');
% xlabel('TARG -- ADJACENT (preceding) -- REST');
ylabel('coh change');

Yall(any(isnan(Yall(:,1:3))'),:) = [];
plot(X(1:3), Yall(:,1:3)', '-k');
xlim([0 max(X)+1]);
lt_plot(X(1:3)+0.2, nanmean(Yall(:,1:3),1), {'Errors', lt_sem(Yall(:,1:3)), 'Color', 'r'});
lt_plot_zeroline;


% =========== REPLOT, BUT USING DIFFERENCE (i.e paired structure included)
lt_subplot(2,2,3); hold on;
title('all syl (only if has all cases)');
% xlabel('TARG -- ADJACENT (preceding) -- REST');
ylabel('coh change');

Yall(any(isnan(Yall)'),:) = [];
plot(X, Yall', '-k');
xlim([0 max(X)+1]);
lt_plot(X+0.2, nanmean(Yall,1), {'Errors', lt_sem(Yall), 'Color', 'r'});
lt_plot_zeroline;



% =========== REPLOT, BUT USING DIFFERENCE (i.e paired structure included)
lt_subplot(2,2,4); hold on;
title('all syl (only if has all cases)');
% xlabel('TARG -- ADJACENT (preceding) -- REST');
ylabel('coh change (minus "REST")');

Yall = Yall - repmat(Yall(:,3),1,length(X));
Yall(any(isnan(Yall)'),:) = [];
plot(X, Yall', '-k');
xlim([0 max(X)+1]);
lt_plot(X+0.2, nanmean(Yall,1), {'Errors', lt_sem(Yall), 'Color', 'r'});
lt_plot_zeroline;


%% ================ [PLOT]
lt_figure; hold on;


% ================ SAME
lt_subplot(2,2,1); hold on;
title('same syl');
xlabel('TARG - SAME MOTIF - DIFF MOTIF');
getsame = 1;

for i=unique(All_swunique)'
    
    Y = nan(1,3);
    
    % --- TARG
    indsthis = All_swunique == i & All_istarg==1;
    pos = All_posreltarget(indsthis); assert(all(pos==0));
    coh = All_cohscal(indsthis);
    Y(1) = mean(coh);
    
    % --- SAME MOTIF
    indsthis = All_swunique == i & All_istarg==0 & All_issame ==getsame & ...
        ~isnan(All_posreltarget);
    %     pos = All_posreltarget(indsthis); assert(all(pos==0));
    coh = All_cohscal(indsthis);
    Y(2) = mean(coh);
    
    % --- DIFF MOTIF
    indsthis = All_swunique == i & All_istarg==0 & All_issame ==getsame & ...
        isnan(All_posreltarget);
    %     pos = All_posreltarget(indsthis); assert(all(pos==0));
    coh = All_cohscal(indsthis);
    Y(3) = mean(coh);
    
    
    % ====
    X = 1:3;
    plot(X, Y, '-ok');
    
end
lt_plot_zeroline;
xlim([0 4]);



% ================ DIFF
lt_subplot(2,2,2); hold on;
title('diff syl');
xlabel('TARG - SAME MOTIF - DIFF MOTIF');
getsame = 0;

for i=unique(All_swunique)'
    
    Y = nan(1,3);
    
    % --- TARG
    indsthis = All_swunique == i & All_istarg==1;
    pos = All_posreltarget(indsthis); assert(all(pos==0));
    coh = All_cohscal(indsthis);
    Y(1) = mean(coh);
    
    % --- SAME MOTIF
    indsthis = All_swunique == i & All_istarg==0 & All_issame ==getsame & ...
        ~isnan(All_posreltarget);
    %     pos = All_posreltarget(indsthis); assert(all(pos==0));
    coh = All_cohscal(indsthis);
    Y(2) = mean(coh);
    
    % --- DIFF MOTIF
    indsthis = All_swunique == i & All_istarg==0 & All_issame ==getsame & ...
        isnan(All_posreltarget);
    %     pos = All_posreltarget(indsthis); assert(all(pos==0));
    coh = All_cohscal(indsthis);
    Y(3) = mean(coh);
    
    
    % ====
    X = 1:3;
    plot(X, Y, '-ok');
    
end
lt_plot_zeroline;
xlim([0 4]);



% ================ COMBINED
lt_subplot(2,2,3); hold on;
title('all syls');
xlabel('TARG - SAME MOTIF - DIFF MOTIF');
ylabel('change in coherence');

for i=unique(All_swunique)'
    
    Y = nan(1,3);
    
    % --- TARG
    indsthis = All_swunique == i & All_istarg==1;
    pos = All_posreltarget(indsthis); assert(all(pos==0));
    coh = All_cohscal(indsthis);
    Y(1) = mean(coh);
    
    % --- SAME MOTIF
    indsthis = All_swunique == i & All_istarg==0 & ...
        ~isnan(All_posreltarget);
    %     pos = All_posreltarget(indsthis); assert(all(pos==0));
    coh = All_cohscal(indsthis);
    Y(2) = mean(coh);
    
    % --- DIFF MOTIF
    indsthis = All_swunique == i & All_istarg==0  & ...
        isnan(All_posreltarget);
    %     pos = All_posreltarget(indsthis); assert(all(pos==0));
    coh = All_cohscal(indsthis);
    Y(3) = mean(coh);
    
    
    % ====
    X = 1:3;
    plot(X, Y, '-ok');
    
end
lt_plot_zeroline;
xlim([0 4]);


%%  === plot coherence chnage a function fo distance from target syl
lt_figure; hold on;


% ================ ALL
lt_subplot(2,2,1); hold on;
title('all syl');
xlabel('pos rel target');
ylabel('coh change');
getsame = 1;

pcolall = lt_make_plot_colors(max(unique(All_bnum)), 0,0);
Xall = [];
Yall =[];
for i=unique(All_swunique)'
    
    indsthis = All_swunique==i;
    
    pos = All_posreltarget(indsthis);
    cohscal = All_cohscal(indsthis);
    issame = All_issame(indsthis);
    istarg = All_istarg(indsthis);
    
    %     keyboard
    % -- any nan, convert to large number
    pos(isnan(pos)) = min(All_posreltarget)-2;
    
    % --- average over chanels
    cohscal = grpstats(cohscal, pos);
    issame = grpstats(issame, pos);
    istarg = grpstats(istarg, pos);
    pos = unique(pos);
    
    % --- sort
    [~, indsort] = sort(pos);
    pos = pos(indsort);
    cohscal = cohscal(indsort);
    issame = issame(indsort);
    istarg = istarg(indsort);
    
    
    % === plot
    pcol = pcolall{unique(All_bnum(indsthis))};
    plot(pos+0.4*rand-0.2, cohscal, '-o', 'Color', pcol);
    
    Xall = [Xall; pos];
    Yall =[Yall; cohscal];
end
[ymean, ysem] = grpstats(Yall, Xall, {'mean', 'sem'});
lt_plot(unique(Xall)+0.2, ymean, {'Errors', ysem});
axis tight
lt_plot_zeroline;

% xlim([0 4]);



% ================ SAME
lt_subplot(2,2,2); hold on;
title('same syl');
xlabel('pos rel target');
ylabel('coh change');
getsame = 1;

Xall = [];
Yall =[];
for i=unique(All_swunique)'
    
    indsthis = All_swunique==i & All_issame==getsame;
    
    if ~any(indsthis)
        continue
    end
    pos = All_posreltarget(indsthis);
    cohscal = All_cohscal(indsthis);
    issame = All_issame(indsthis);
    istarg = All_istarg(indsthis);
    
    %     keyboard
    % -- any nan, convert to large number
    pos(isnan(pos)) = min(All_posreltarget)-2;
    
    % --- average over chanels
    cohscal = grpstats(cohscal, pos);
    issame = grpstats(issame, pos);
    istarg = grpstats(istarg, pos);
    pos = unique(pos);
    
    % --- sort
    [~, indsort] = sort(pos);
    pos = pos(indsort);
    cohscal = cohscal(indsort);
    issame = issame(indsort);
    istarg = istarg(indsort);
    
    
    % === plot
    pcol = pcolall{unique(All_bnum(indsthis))};
    plot(pos+0.4*rand-0.2, cohscal, '-o', 'Color', pcol);
    
    Xall = [Xall; pos];
    Yall =[Yall; cohscal];
end
[ymean, ysem] = grpstats(Yall, Xall, {'mean', 'sem'});
lt_plot(unique(Xall)+0.2, ymean, {'Errors', ysem});
axis tight
lt_plot_zeroline;

% xlim([0 4]);



% ================ DIFF
lt_subplot(2,2,3); hold on;
title('diff syl');
xlabel('pos rel target');
ylabel('coh change');
getsame = 0;

Xall = [];
Yall =[];
for i=unique(All_swunique)'
    
    indsthis = All_swunique==i & All_issame==getsame;
    
    pos = All_posreltarget(indsthis);
    cohscal = All_cohscal(indsthis);
    issame = All_issame(indsthis);
    istarg = All_istarg(indsthis);
    
    %     keyboard
    % -- any nan, convert to large number
    pos(isnan(pos)) = min(All_posreltarget)-2;
    
    % --- average over chanels
    cohscal = grpstats(cohscal, pos);
    issame = grpstats(issame, pos);
    istarg = grpstats(istarg, pos);
    pos = unique(pos);
    
    % --- sort
    [~, indsort] = sort(pos);
    pos = pos(indsort);
    cohscal = cohscal(indsort);
    issame = issame(indsort);
    istarg = istarg(indsort);
    
    
    % === plot
    pcol = pcolall{unique(All_bnum(indsthis))};
    plot(pos+0.4*rand-0.2, cohscal, '-o', 'Color', pcol);
    
    Xall = [Xall; pos];
    Yall =[Yall; cohscal];
end
[ymean, ysem] = grpstats(Yall, Xall, {'mean', 'sem'});
lt_plot(unique(Xall)+0.2, ymean, {'Errors', ysem});
axis tight
lt_plot_zeroline;

% xlim([0 4]);

