function lt_neural_POPLEARN_Xcov_PitchLearn(OUTSTRUCT, OUTSTRUCT_XCOV, SwitchStruct, PARAMS, ...
    onlygoodexpt)

%% lt 2/20/19 - summarizes learning trajectories. only experiemtns in outstructxcov

%% ======= first downsample OUTSRRUCT to get only the experiments in xcov
[OUTSTRUCT, OUTSTRUCT_XCOV, indsXcov] = lt_neural_POPLEARN_MatchOutstructs(OUTSTRUCT, OUTSTRUCT_XCOV, ...
    SwitchStruct, onlygoodexpt);

%% ====== for each case compute z-scored laerning

LearnTargDir_Z = nan(size(OUTSTRUCT.bnum));
LearnTargDir_Slope = nan(size(OUTSTRUCT.bnum));
TrainDuration = nan(size(OUTSTRUCT.bnum)); % time from mid of base bin to mid of training
for i=1:length(OUTSTRUCT.bnum)
    
    ffvals = OUTSTRUCT.ffvals{i};
    tvals = OUTSTRUCT.tvals{i};
    
    if all(isnan(ffvals))
        continue
    end
    
    % ---- USE THE INDS FROM XCOV
%     indsbase = OUTSTRUCT_XCOV.inds_base_epoch{indsXcov{i}(1)};
    indsbase = OUTSTRUCT_XCOV.inds_base_allgood{indsXcov{i}(1)};
    indsWN = OUTSTRUCT_XCOV.inds_WN_epoch{indsXcov{i}(1)};
    
    % --- sanity check....
    assert(all(OUTSTRUCT.indsXcov_all{i} == indsXcov{i}));
    
    
    % ---- get zscore
    basemean = mean(ffvals(indsbase));
    basestd = std(ffvals(indsbase));
    
    ffz = (mean(ffvals(indsWN)) - basemean)./basestd;
    
    
    % ---- flip if training is down
    traindir = OUTSTRUCT_XCOV.learndirTarg(indsXcov{i}(1));        
    ffz = ffz*traindir;
    
    
    % ============ SAVE TO OUTPUT
    LearnTargDir_Z(i) = ffz;
    
    
    %  ################### TRAINING DURATION
    traindur = (mean(tvals(indsWN)) - mean(tvals(indsbase)))*24;
    TrainDuration(i) = traindur;
    
    % ############################### LEARNING RATE (SLOPE OF FF, DURING
    % TRAINING)
    indsWN_all = OUTSTRUCT.indsWN{i};
    
    ffthis = ffvals(indsWN_all);
    tthis = tvals(indsWN_all);
    tthis = (tthis-tthis(1))*24; % conver to hours
    
    [~,~,~,~,~, stats] = lt_regress(ffthis, tthis, 0);
    learnslope = stats.slope*traindir; % flip slope so up is in direction of training.
    
    LearnTargDir_Slope(i) = learnslope;
end
    
OUTSTRUCT.LearnTargDir_Z = LearnTargDir_Z;
OUTSTRUCT.TrainDuration = TrainDuration;
OUTSTRUCT.LearnTargDir_Slope = LearnTargDir_Slope;


%% ======= SUKMARIZE LEANRING
fieldtoget = 'LearnTargDir_Z';
[~, ~, ~, ~, ...
    allswitch_bnum, allswitch_enum, allswitch_swnum, allswitch_dat, ~, allswitch_Nmotifs] = ...
    lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);

allswitch_dat = squeeze(allswitch_dat)';


fieldtoget = 'TrainDuration';
[~, ~, ~, ~, ...
    allswitch_bnum, allswitch_enum, allswitch_swnum, allswitch_dat_traindur, ~, allswitch_Nmotifs] = ...
    lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);

allswitch_dat_traindur = squeeze(allswitch_dat_traindur)';


fieldtoget = 'LearnTargDir_Slope';
[~, ~, ~, ~, ...
    allswitch_bnum, allswitch_enum, allswitch_swnum, allswitch_dat_slope, ~, allswitch_Nmotifs] = ...
    lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);

allswitch_dat_slope = squeeze(allswitch_dat_slope)';


%% ======= ONE PLOT FOR EACH EXPT, PLOT RAW DATA.
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for i=1:length(allswitch_bnum)
   bb = allswitch_bnum(i);
   ee = allswitch_enum(i);
   sw = allswitch_swnum(i);
   bname = SwitchStruct.bird(bb).birdname;
   ename = SwitchStruct.bird(bb).exptnum(ee).exptname;
   
   indsthis = find(OUTSTRUCT.bnum==bb & OUTSTRUCT.enum==ee & OUTSTRUCT.switch==sw ...
       & OUTSTRUCT.istarg==1);
   
   % == get ff over al cases. check that same num trials. then take first
   % ind
   ff = OUTSTRUCT.ffvals(indsthis);
   assert(length(unique(cellfun(@(x)length(x), ff)))==1);
   
   % ========= get data and plot
   indsthis = indsthis(1);
   ff = OUTSTRUCT.ffvals{indsthis};
   t = OUTSTRUCT.tvals{indsthis};
   t = (t-min(t))*24; % convert to hr from first trial'
   traindir = OUTSTRUCT.learndirTarg(indsthis);

   % ---- get inds from xcov
   indxcov = OUTSTRUCT.indsXcov_all{indsthis}(1);
   indsbase = OUTSTRUCT_XCOV.inds_base_allgood{indxcov};
   indswn = OUTSTRUCT_XCOV.inds_WN_allgood{indxcov};
   indsbase_epoch = OUTSTRUCT_XCOV.inds_base_epoch{indxcov};
   indswn_epoch = OUTSTRUCT_XCOV.inds_WN_epoch{indxcov};
    
   % =========== PLOT
   [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
   title([bname '-' ename '-sw' num2str(sw)]);
   plot(t, ff, 'xk');
   ylabel(['TRAINDIR: '  num2str(traindir)]);
   axis tight;
   
   % --- note down that are epoch trials
   plot(t(indsbase_epoch), ff(indsbase_epoch), 'xb');
   plot(t(indswn_epoch), ff(indswn_epoch), 'xr');
   
   % -- lines
   line([t(max(indsbase)) t(max(indsbase))], ylim, 'Color', 'r');
   ybase = mean(ff(indsbase));
   line(xlim, [ybase ybase], 'Color', [0.7 0.7 0.7]);
   
   % -- lines for epoch boundaries
   tredges = OUTSTRUCT_XCOV.trialedges_epoch{indxcov};
   for j=1:length(tredges)
       if isnan(tredges(j))
           continue
       end
       if j==length(tredges)
       line([t(tredges(j)-1) t(tredges(j)-1)]+0.25/60, ylim);    
       else
       line([t(tredges(j)) t(tredges(j))]-0.25/60, ylim);
       end
%       line([t(tredges(j+1)-1) t(tredges(j+1)-1)]+0.25/60, ylim);
   end
end


%% ======= ONE PLOT FOR EACH EXPT, PLOT RAW DATA. [SAME TYPE]
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


for i=1:length(allswitch_bnum)
    bb = allswitch_bnum(i);
    ee = allswitch_enum(i);
    sw = allswitch_swnum(i);
        bname = SwitchStruct.bird(bb).birdname;
        ename = SwitchStruct.bird(bb).exptnum(ee).exptname;
    
    indsthis = find(OUTSTRUCT.bnum==bb & OUTSTRUCT.enum==ee & OUTSTRUCT.switch==sw ...
        & OUTSTRUCT.istarg==0 & OUTSTRUCT.issame==1);
    
    if ~any(indsthis)
        continue
    end
    
    motifsdone = {};
    for ii=1:length(indsthis)
        indsthisthis = indsthis(ii);
        
        % ===== get info
        motifname = OUTSTRUCT.motifname(indsthisthis);
        if ismember(motifname, motifsdone)
            continue
        end
        motifsdone = [motifsdone motifname];
                
        % ========= get data and plot
        ff = OUTSTRUCT.ffvals{indsthisthis};
        t = OUTSTRUCT.tvals{indsthisthis};
        t = (t-min(t))*24; % convert to hr from first trial'
        traindir = OUTSTRUCT.learndirTarg(indsthisthis);
        
        % ---- get inds from xcov
        indxcov = OUTSTRUCT.indsXcov_all{indsthisthis}(1);
        indsbase = OUTSTRUCT_XCOV.inds_base_allgood{indxcov};
        indswn = OUTSTRUCT_XCOV.inds_WN_allgood{indxcov};
        indsbase_epoch = OUTSTRUCT_XCOV.inds_base_epoch{indxcov};
        indswn_epoch = OUTSTRUCT_XCOV.inds_WN_epoch{indxcov};
        
        % =========== PLOT
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title([bname '-' ename '-sw' num2str(sw)]);
        xlabel(['[SAME] ' motifname]);
        plot(t, ff, 'xk');
        ylabel(['TRAINDIR: '  num2str(traindir)]);
        axis tight;
        
        % --- note down that are epoch trials
        plot(t(indsbase_epoch), ff(indsbase_epoch), 'xb');
        plot(t(indswn_epoch), ff(indswn_epoch), 'xr');
        
        % -- lines
        line([t(max(indsbase)) t(max(indsbase))], ylim, 'Color', 'r');
        ybase = mean(ff(indsbase));
        line(xlim, [ybase ybase], 'Color', [0.7 0.7 0.7]);
        
        % -- lines for epoch boundaries
        tredges = OUTSTRUCT_XCOV.trialedges_epoch{indxcov};
        for j=1:length(tredges)
            if isnan(tredges(j))
                continue
            end
            if j==length(tredges)
                line([t(tredges(j)-1) t(tredges(j)-1)]+0.25/60, ylim);
            else
                line([t(tredges(j)) t(tredges(j))]-0.25/60, ylim);
            end
            %       line([t(tredges(j+1)-1) t(tredges(j+1)-1)]+0.25/60, ylim);
        end
        
    end
end


%% ======= PLOT (ONE LINE FOR EACH EXPERIMENT, AVERAGE OVER ALL SYLS OF A GIVBEN TYPE)

figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];



% =================== 1) 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('target learning');
xlabel('z, rel baseline');

x = allswitch_dat(:,1);
lt_plot_histogram(x);


% =================== 1) 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('target learning');
xlabel('train dur (bin means, hr');
ylabel('learning(z)');
x = [zeros(size(allswitch_dat_traindur,1),1) allswitch_dat_traindur(:,1)];
y = [zeros(size(allswitch_dat,1),1) allswitch_dat(:,1)];
plot(x', y', '-ok');
% -- note down expt name
for i=1:size(allswitch_dat,1)
   bb = allswitch_bnum(i);
   ee = allswitch_enum(i);
   sw = allswitch_swnum(i);
   bname = SwitchStruct.bird(bb).birdname;
   ename = SwitchStruct.bird(bb).exptnum(ee).exptname;
   lt_plot_text(allswitch_dat_traindur(i,1), allswitch_dat(i,1), [bname '-' ename '-sw' num2str(sw)], 'r', 9);
end
lt_plot_zeroline;

% =================== 1) 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('target learning');
ylabel('learning(z)');
y = allswitch_dat(:,1);
lt_plot_bar(1, mean(y), {'Errors', lt_sem(y)})
lt_plot_MultDist({y}, 1, 0, 'r', 1);
% --- color each bird
scatter(2*ones(size(y)), y, [], allswitch_bnum);
lt_plot_zeroline;
% [~, p] = ttest(y);
[p] = signrank(y);
lt_plot_pvalue(p, 'srank');

% =================== 1) 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('target learning');
ylabel('learning slope (hz/hr)');
y = allswitch_dat_slope(:,1);
lt_plot_bar(1, mean(y), {'Errors', lt_sem(y)})
lt_plot_MultDist({y}, 1, 0, 'r', 1);
lt_plot_zeroline;


% =================== 1) 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('target learning');
xlabel('slope, during WN (hz/hr)');

x = allswitch_dat_slope(:,1);
lt_plot_histogram(x);


%% ======== PLOT [TARGET LAERNING VS. SAME TYPE]

figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


% =================== 1) 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same type learning');
xlabel('train dur (bin means, hr');
ylabel('learning(z)');
indscol = 2;

x = [zeros(size(allswitch_dat_traindur,1),1) allswitch_dat_traindur(:,indscol)];
y = [zeros(size(allswitch_dat,1),1) allswitch_dat(:,indscol)];
indsrows = ~any(isnan(x'));
x = x(indsrows,:);
y = y(indsrows,:);
plot(x', y', '-ok');
% -- note down expt name
indsrows = find(indsrows);
for i=indsrows
   bb = allswitch_bnum(i);
   ee = allswitch_enum(i);
   sw = allswitch_swnum(i);
   bname = SwitchStruct.bird(bb).birdname;
   ename = SwitchStruct.bird(bb).exptnum(ee).exptname;
   lt_plot_text(allswitch_dat_traindur(i,indscol), allswitch_dat(i,indscol), [bname '-' ename '-sw' num2str(sw)], 'r', 9);
end
lt_plot_zeroline;


% =================== 1) 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('DIFF type learning');
xlabel('train dur (bin means, hr');
ylabel('learning(z)');
indscol = 3;

x = [zeros(size(allswitch_dat_traindur,1),1) allswitch_dat_traindur(:,indscol)];
y = [zeros(size(allswitch_dat,1),1) allswitch_dat(:,indscol)];
indsrows = ~any(isnan(x'));
x = x(indsrows,:);
y = y(indsrows,:);
plot(x', y', '-ok');
% -- note down expt name
indsrows = find(indsrows);
for i=indsrows
   bb = allswitch_bnum(i);
   ee = allswitch_enum(i);
   sw = allswitch_swnum(i);
   bname = SwitchStruct.bird(bb).birdname;
   ename = SwitchStruct.bird(bb).exptnum(ee).exptname;
   lt_plot_text(allswitch_dat_traindur(i,indscol), allswitch_dat(i,indscol), [bname '-' ename '-sw' num2str(sw)], 'r', 9);
end
lt_plot_zeroline;



% =================== 1) 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same type learning');
ylabel('learning(z)');
indscol = 2;
y = allswitch_dat(:,indscol);
indsrows = ~isnan(y);
y = y(indsrows,:);
lt_plot_bar(1, mean(y), {'Errors', lt_sem(y)})
lt_plot_MultDist({y}, 1, 0, 'r', 1);
lt_plot_zeroline;
[~, p] = ttest(y);
lt_plot_pvalue(p, 'ttest');


% =================== 1) 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('targ vs. same');
xlabel('TARG - SAME');
ylabel('learning(z)');
indscol = [1 2];

y = allswitch_dat(:,indscol);
indsrows = ~any(isnan(y'));
y = y(indsrows,:);

lt_plot_bar(indscol, mean(y), {'Errors', lt_sem(y)});
plot(indscol, y', '-r');

[~, p] = ttest(y(:,1), y(:,2));
lt_plot_pvalue(p, 'ttestvs');
lt_plot_zeroline;


% =================== 1) 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('targ vs. diff');
xlabel('TARG - DIFF');
ylabel('learning(z)');
indscol = [1 3];

y = allswitch_dat(:,indscol);
indsrows = ~any(isnan(y'));
y = y(indsrows,:);

lt_plot_bar(indscol, mean(y), {'Errors', lt_sem(y)});
plot(indscol, y', '-r');

[~, p] = ttest(y(:,1), y(:,2));
lt_plot_pvalue(p, 'ttestvs');
lt_plot_zeroline;



%% ================ ENTIRE DISTRIBUTION OF NONTARG (DAT = SYL)
% ============= SAME
indsthis = OUTSTRUCT.issame==1 & OUTSTRUCT.istarg==0;
y = OUTSTRUCT.LearnTargDir_Z(indsthis);
y = unique(y);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('dat = syl [SAME]');
lt_plot_MultDist({y}, 1);
lt_plot_zeroline;
ylim([-2 2]);
ylabel('z, rel base');

% ============= DIFF
indsthis = OUTSTRUCT.issame==0 & OUTSTRUCT.istarg==0;
y = OUTSTRUCT.LearnTargDir_Z(indsthis);
y = unique(y);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('dat = syl [DIFF]');
lt_plot_MultDist({y}, 1);
lt_plot_zeroline;
ylim([-2 2]);
ylabel('z, rel base');



