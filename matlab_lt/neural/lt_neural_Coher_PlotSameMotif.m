function lt_neural_Coher_PlotSameMotif(OUTSTRUCT, SwitchStruct,MOTIFSTATS_Compiled,...
    PARAMS, indtoget_b_e_s)
%% lt 12/20/18 - plots cioherence scalar change during laerning, whether is same or diff motif.

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

%% filter? by type of switch?

OUTSTRUCT = lt_structure_RmvEmptyField(OUTSTRUCT); % so followniog code works...

% % === filter each dimension independenyl
% if ~isempty(birdstoplot)
%     indstokeep = ismember(OUTSTRUCT.bnum, birdstoplot);
%     OUTSTRUCT = lt_structure_subsample_all_fields(OUTSTRUCT, indstokeep, 1);
% end
% 
% if ~isempty(expttoplot)
%     indstokeep = ismember(OUTSTRUCT.enum, expttoplot);
%     OUTSTRUCT = lt_structure_subsample_all_fields(OUTSTRUCT, indstokeep, 1);
% end
% 
% if ~isempty(swtoplot)
%     indstokeep = ismember(OUTSTRUCT.switch, swtoplot);
%     OUTSTRUCT = lt_structure_subsample_all_fields(OUTSTRUCT, indstokeep, 1);
% end

% == filter by specific types of switches.
if ~isempty(indtoget_b_e_s)
indstokeep = ismember([OUTSTRUCT.bnum OUTSTRUCT.enum OUTSTRUCT.switch], indtoget_b_e_s, 'rows');
OUTSTRUCT = lt_structure_subsample_all_fields(OUTSTRUCT, indstokeep, 1);
end

%% NOTE
% if multipel targets on different motifs, then will skip that case. if on same motif
% then will take  the earliest one.
%% ############

[indsgrp_switch, indsgrp_switch_unique] = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch});
[indsgrp_chanpair, indsgrp_chanpair_unique] = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch, ...
    OUTSTRUCT.chanpair});

% figcount=1;
% subplotrows=6;
% subplotcols=1;
% fignums_alreadyused=[];
% hfigs=[];
% hsplots = [];
%
% indsgrp_switch_unique = unique(indsgrp_switch);
% indsgrp_chanpair_unique = unique(indsgrp_chanpair);

All_istarg = [];
All_cohscal = [];
All_posreltarget = [];
All_swunique = [];
All_issame = [];
All_bnum = [];
for i=1:length(indsgrp_switch_unique)
    swgrpthis = indsgrp_switch_unique(i);
    bnum = unique(OUTSTRUCT.bnum(indsgrp_switch==swgrpthis));
    enum = unique(OUTSTRUCT.enum(indsgrp_switch==swgrpthis));
    swnum = unique(OUTSTRUCT.switch(indsgrp_switch==swgrpthis));
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
    %         motifID = lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)); % ---- get positions within global motif
    %
    %         cohscal = OUTSTRUCT.cohscal_diff(indsthis);
    %         %         cohscal = OUTSTRUCT.CohMean_WNminusBase_scalar(indsthis);
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
    indsthis = indsgrp_switch==swgrpthis;
    
    istarg = OUTSTRUCT.istarg(indsthis);
    isame = OUTSTRUCT.issame(indsthis);
    motifs = OUTSTRUCT.motifname(indsthis);
    [motifID, ~, motifposition] = lt_neural_QUICK_MotifID(bname, motifs); % ---- get positions within global motif
    %     cohscal = OUTSTRUCT.CohMean_WNminusBase_scalar(indsthis);
    cohscal = OUTSTRUCT.cohscal_diff(indsthis);
    
    [ymean, ysem] = grpstats(cohscal, motifID, {'mean', 'sem'});
    %     x = unique(motifID);
    %     lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
    %
    %     % -------- NOTE DOWN POSITION OF TARGET SYSL
    %     indsthis = indsgrp_switch==swgrpthis & OUTSTRUCT.istarg==1;
    %     xtarg = unique(lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)));
    %     plot(xtarg, clim(1)+0.02, '^r');
    %
    %     % ------- NOTE POSITION OF SAME_TYPES
    %     indsthis = indsgrp_switch==swgrpthis & OUTSTRUCT.istarg==0 & OUTSTRUCT.issame==1;
    %     xtarg = unique(lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)));
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
end


%% =============== [PLOT] -- TARG, ADJACENT, NOT ADJACENT
lt_figure; hold on;


% ================ ALL
lt_subplot(2,2,1); hold on;
title('all syl');
xlabel('TARG -- ADJACENT (preceding) -- REST');
ylabel('coh change');

pcolall = lt_make_plot_colors(length(unique(All_bnum)), 0,0);
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
    
    X = 1:3;
    Y = nan(1,3);
    
    % -- targ
    Y(1) = mean(cohscal(pos==0));
    % -- adjacnet (-preceding)
    Y(2) = mean(cohscal(pos==-1));
    % -- rest
    Y(3) = mean(cohscal(pos<-1));
    
    if any(pos>0)
        disp('NOTE: in some case ther eare syls folliowujng - ignoreing those for "adjavent');
    end
%     assert(all(pos<1), 'assumed that all are precedeing...');
    
    plot(X+0.4*rand-0.2, Y, '-o', 'Color', pcol);
   
    Yall =[Yall; Y];
end


lt_plot([1:3]+0.3, nanmean(Yall,1), {'Errors', lt_sem(Yall), 'Color', 'k'});
xlim([0 4]);
lt_plot_zeroline;


% =========== REPLOT, BUT USING DIFFERENCE (i.e paired structure included)
lt_subplot(2,2,2); hold on;
title('all syl');
xlabel('TARG -- ADJACENT (preceding) -- REST');
ylabel('coh change (minus "REST")');

Yall = Yall - repmat(Yall(:,3),1,3);
plot(1:3, Yall', '-k');
xlim([0 4]);
lt_plot([1:3]+0.2, nanmean(Yall,1), {'Errors', lt_sem(Yall), 'Color', 'r'});
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

pcolall = lt_make_plot_colors(length(unique(All_bnum)), 0,0);
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

pcolall = lt_make_plot_colors(length(unique(All_bnum)), 0,0);
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

pcolall = lt_make_plot_colors(length(unique(All_bnum)), 0,0);
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

%%  === plot coherence chnage a function fo distance from target syl
lt_figure; hold on;


% ================ ALL
lt_subplot(2,2,1); hold on;
title('all syl');
xlabel('pos rel target');
ylabel('coh change');
getsame = 1;

pcolall = lt_make_plot_colors(length(unique(All_bnum)), 0,0);
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

pcolall = lt_make_plot_colors(length(unique(All_bnum)), 0,0);
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

pcolall = lt_make_plot_colors(length(unique(All_bnum)), 0,0);
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
