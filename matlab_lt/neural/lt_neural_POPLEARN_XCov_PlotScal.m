function lt_neural_POPLEARN_XCov_PlotScal(Yscalar, OUTSTRUCT, SwitchCohStruct, SwitchStruct, ...
    MOTIFSTATS_Compiled, PARAMS, clim, onlygoodexpt, expttype, plotlevel, OUTSTRUCT_LFP, ...
    useoldInds)
%% 2/1/19 - Plots each scalar, each experiement...

% clim = [-0.1 0.1];
useglobalmotifname = 0; % then all expt across a bird will be aligned... motif names will not necessarily be correct for each experiement though;
minDiffN = 3; % will only noramlize to diff ytpe if diff type has N this or larger.
takeindsWNsecondhalf = 0; % then essentially takes last quartile to compute learning.
% LEAVE AT 0 - if 1, then inds will not necessariyl be exactly matched
% between xcov, coh, and laerning.


%% ======== filter experiments
if onlygoodexpt==1
    % ===== filter outstruct
    [OUTSTRUCT, indstokeep] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT, SwitchStruct, expttype);
    % ===== filter the scalar array
    Yscalar = Yscalar(indstokeep);
end

% ====== SANITY CHECK.
assert(length(Yscalar) == length(OUTSTRUCT.bnum));


%% ============= 1) GET UNIQUE INDS FOR EACH SWITCH, AND EACH CHAN PAIR
[indsgrp_switch, indsgrp_switch_unique] = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch});
[indsgrp_chanpair, indsgrp_chanpair_unique] = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch, ...
    OUTSTRUCT.chanpair});


%% ############ ONE PLOT FOR EACH SWITCH, SHOW ALL MOTIFS.
figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% ========== TO COLLECT FOR MEAN PLOTS
Yall = []; % targ, same, diff, nontarg (means across chans, motifs)
% Ydiffall = {};
Yall_allsyls = cell(length(indsgrp_switch_unique), 4); % same as above, but do not average withni syl type.
All_bname ={};
All_bnum = [];
All_ename = {};
All_swnum = [];
All_learndir = [];
All_enum = [];
Allchanpair_Y = [];
Allchanpair_idx = [];

All_Nmotifs = [];

% ========= GO THRU ALL SWITCHES AND CHAN PAIRS.
for i=1:length(indsgrp_switch_unique)
    swgrpthis = indsgrp_switch_unique(i);
    bnum = unique(OUTSTRUCT.bnum(indsgrp_switch==swgrpthis));
    enum = unique(OUTSTRUCT.enum(indsgrp_switch==swgrpthis));
    swnum = unique(OUTSTRUCT.switch(indsgrp_switch==swgrpthis));
    bname = MOTIFSTATS_Compiled.birds(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    
    % ============= ONE PLOT FOR EACH SWITCH.
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(swnum)]);
    ylabel(['WN - base']);
    
    % ########################### 1) plot each channel pair its own line
    Ytmp = nan(1,4); % targ, same, diff, nontarg
    for chanpair = indsgrp_chanpair_unique'
        indsthis = indsgrp_switch==swgrpthis & indsgrp_chanpair==chanpair;
        if ~any(indsthis)
            continue
        end
        %         OUTSTRUCT.motifnum(indsthis)
        if useglobalmotifname==1
            motifID = lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)); % ---- get positions within global motif
        else
            motifID = OUTSTRUCT.motifnum(indsthis);
        end
        
        if isfield(OUTSTRUCT, 'chanpair_actual')
            chnums = OUTSTRUCT.chanpair_actual(indsthis,:);
        else
            chnums = OUTSTRUCT.chanpair(indsthis,:);
        end
        assert(length(unique(chnums)')==2);
        chnums = chnums(1,:);
        
        % ----- might also have neuron ID
        if isfield(OUTSTRUCT, 'neurpair')
            % ---- then has neuron id
            neurpair = OUTSTRUCT.neurpair(indsthis,:);
            assert(length(unique(neurpair)')==2);
            neurpair = neurpair(1,:);
        end
        
        
        cohscal = Yscalar(indsthis);
        assert(length(unique(motifID)) == length(motifID), 'each motif shoudl have at most 1 datapoint');
        
        % ---- sort in order of motifs
        [~, indsort] = sort(motifID);
        motifID = motifID(indsort);
        cohscal = cohscal(indsort);
        
        % =================== PLOT FOR THIS CHAN PAIR.
        plot(motifID, cohscal, 'o-k');
        if isempty(neurpair)
            tmp = num2str(chnums);
        else
            tmp = [num2str(chnums) ' [nID:' num2str(neurpair) ']'];
        end
        lt_plot_text(motifID(end)+0.3, cohscal(end), tmp, 'm', 8);
        
        % ================================ SAVE OUTPUT
        
        %         Allchanpair_Y = [Allchanpair_Y; ];
    end
    lt_plot_zeroline;
    
    
    % ############################# MEAN ACROSS ALL CHANNEL PAIRS
    indsthis = indsgrp_switch==swgrpthis;
    
    istarg = OUTSTRUCT.istarg(indsthis);
    issame = OUTSTRUCT.issame(indsthis);
    cohscal = Yscalar(indsthis);
    learndir = unique(OUTSTRUCT.learndirTarg(indsthis));
    assert(length(learndir)==1);
    
    if useglobalmotifname==1
        motifID = lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)); % ---- get positions within global motif
    else
        motifID = OUTSTRUCT.motifnum(indsthis); % ---- get positions within global motif
        motifnames = OUTSTRUCT.motifname(indsthis);
    end
    
    
    
    % ================ SAVE INDIVBIDUAL SYLS EBFORE PERFORM MEAN
       Yall_allsyls{i, 1} = cohscal(istarg==1);
       Yall_allsyls{i, 2} = cohscal(istarg==0 & issame==1);
       Yall_allsyls{i, 3} = cohscal(istarg==0 & issame==0);

    
    % ============= GET MEAN ACROSS CHANNEL PAIRS.
    [ymean, ysem] = grpstats(cohscal, motifID, {'mean', 'sem'});
    [istarg, Ntarg] = grpstats(istarg, motifID, {'mean', 'numel'});
    istarg = logical(istarg);
    [issame, Nsame] = grpstats(issame, motifID, {'mean', 'numel'});
    issame = logical(issame);
    x = unique(motifID);
    
    % ================= PLOT MEAN
    disp(ymean);
    disp(x);
    if length(x)==length(ymean)
        lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
    end
    
    % ===================== ANNOTATE PLOT.
    % -------- NOTE DOWN POSITION OF TARGET SYSL
    %     indsthis = indsgrp_switch==swgrpthis & OUTSTRUCT.istarg==1;
    %     xtarg = unique(lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)));
    %     plot(xtarg, clim(1)+0.02, '^r');
    if any(istarg)
    plot(x(istarg), clim(1)+0.02, '^r');
    end
    
    % ------- NOTE POSITION OF SAME_TYPES
    %     indsthis = indsgrp_switch==swgrpthis & OUTSTRUCT.istarg==0 & OUTSTRUCT.issame==1;
    %     xtarg = unique(lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)));
    %     if ~isempty(xtarg)
    %         plot(xtarg, clim(1)+0.02, '^b');
    %     end
    if any(issame)
        plot(x(issame), clim(1)+0.02, '^b');
    end
    
    % ------- note down global mean [of diff types]
    indtmp = istarg==0 & issame==0;
    yglob = mean(ymean(indtmp));
    yglob = mean(ymean);
    line(xlim, [yglob yglob], 'Color', 'r');
    
    
    % ----- labels
    if useglobalmotifname==1
        [~, motiflabels] = lt_neural_QUICK_MotifID(bname);
    else
        motiflabels = {};
        for j=1:max(x)
            if ~any(x==j)
                motiftmp = ' ';
            else
                indstmp = indsgrp_switch==swgrpthis & OUTSTRUCT.motifnum==j;
                motiftmp = unique(OUTSTRUCT.motifname(indstmp));
                motiftmp = motiftmp{1};
            end
            motiflabels = [motiflabels motiftmp];
        end
    end
    set(gca, 'XTick', 1:length(motiflabels));
    set(gca, 'XTickLabel', motiflabels);
    rotateXLabels(gca, 90);
    ylim(clim);
    
    
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
%     Ydiffall = [Ydiffall; {ymean(istarg==0 & issame==0)}];
    All_bname = [All_bname; bname];
    All_ename = [All_ename; ename];
    All_swnum = [All_swnum; swnum];
    All_bnum = [All_bnum; bnum];
    All_enum = [All_enum; enum];
    All_learndir = [All_learndir; learndir];
    
    % =========== SAVE STUFF FOR MOTIFS.
    All_Nmotifs = [All_Nmotifs; [sum(istarg) sum(issame) sum(istarg==0 & issame==0)]];
    
    
end


%% ======== ONE DATAPOINT PER CHANNEL -

Yall_scal = {};
Yall_bnum = [];
Yall_enum = [];
Yall_chanID = [];
Yall_istarg = {};
Yall_issame = {};
Yall_sw = [];
for i=1:length(indsgrp_chanpair_unique)
    
    indsthis = indsgrp_chanpair == indsgrp_chanpair_unique(i);
    
    bnum = unique(OUTSTRUCT.bnum(indsthis));
    enum = unique(OUTSTRUCT.enum(indsthis));
    swnum = unique(OUTSTRUCT.switch(indsthis));
    
    bname = MOTIFSTATS_Compiled.birds(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    
    % =========================== DATA
    cohscal = Yscalar(indsthis);
    istarg = OUTSTRUCT.istarg(indsthis);
    issame = OUTSTRUCT.issame(indsthis);
    motifID = lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)); % ---- get positions within global motif
    chnums = OUTSTRUCT.chanpair_actual(indsthis,:);
    neurpair = OUTSTRUCT.neurpair(indsthis,:);
    
    % ================================ SAVE OUTPUT
    Yall_scal = [Yall_scal; cohscal];
    Yall_bnum = [Yall_bnum ; bnum];
    Yall_enum = [Yall_enum ; enum];
    Yall_chanID = [Yall_chanID; indsgrp_chanpair_unique(i)];
    Yall_istarg = [Yall_istarg; istarg];
    Yall_issame = [Yall_issame; issame];
    Yall_sw = [Yall_sw; swnum];
end

lt_figure; hold on;
% === 1) sort in order of number of 1syls
Nall =cellfun(@length, Yall_scal);
% [~, indsort] = sort(Nall);
[~, ~, indrank] = lt_tools_sort(Nall);


% ################################# NOT CENTERED AT DIFF TYPE
lt_subplot(1,2,1); hold on;
xlabel('xcov (minus base)');
ylabel('chan');
Ylabstr = {};
pcol_expt = lt_make_plot_colors(max(Yall_enum+Yall_bnum), 0, 0);
XMIN = 1.1*min(Yscalar);
for j=1:length(indrank)
    
    y = indrank(j);
    dat = Yall_scal{j};
    istarg = Yall_istarg{j};
    issame = Yall_issame{j};
    
    % ---
    indtmp = istarg==1;
    plot(dat(indtmp), y, 'or');
    %    lt_plot(dat(indtmp), y, {'Color', 'r'});
    
    % --
    indtmp = istarg==0 & issame==1;
    if any(indtmp)
        plot(dat(indtmp), y, 'xb');
    end
    
    % ---
    indtmp = istarg==0 & issame==0;
    if any(indtmp)
        plot(dat(indtmp), y, 'xk');
    end
    % ---
    
    % --- mark to help distinguish different expts
    pcoltmp = pcol_expt{Yall_enum(j)+Yall_bnum(j)};
    lt_plot(XMIN, y, {'Color', pcoltmp});
    
    
    str = [num2str(Yall_bnum(j)) '-' num2str(Yall_enum(j)) '-' num2str(Yall_sw(j))];
    %    lt_plot_text(1.05*max(dat), y, str, 'm', 8);
    
    Ylabstr = [Ylabstr; str];
end
[~, indtmp] = sort(indrank);
Ylabstr = Ylabstr(indtmp);

set(gca, 'YTick', 1:length(Ylabstr), 'YTickLabel', Ylabstr);

% ################################# CENTERED AT DIFF TYPE
lt_subplot(1,2,2); hold on;
xlabel('xcov (minus base)(centered using diff type)');
ylabel('chan');
Ylabstr = {};
pcol_expt = lt_make_plot_colors(max(Yall_enum+Yall_bnum), 0, 0);
XMIN = 1.1*min(Yscalar);
for j=1:length(indrank)
    
    y = indrank(j);
    dat = Yall_scal{j};
    istarg = Yall_istarg{j};
    issame = Yall_issame{j};
    
    % -- center using mean of diff type
    dat = dat-mean(dat(istarg==0 & issame==0));
    
    % ---
    indtmp = istarg==1;
    plot(dat(indtmp), y, 'or');
    %    lt_plot(dat(indtmp), y, {'Color', 'r'});
    
    % --
    indtmp = istarg==0 & issame==1;
    if any(indtmp)
        plot(dat(indtmp), y, 'xb');
    end
    
    % ---
    indtmp = istarg==0 & issame==0;
    if any(indtmp)
        plot(dat(indtmp), y, 'xk');
    end
    % ---
    
    % --- mark to help distinguish different expts
    pcoltmp = pcol_expt{Yall_enum(j)+Yall_bnum(j)};
    lt_plot(XMIN, y, {'Color', pcoltmp});
    
    
    str = [num2str(Yall_bnum(j)) '-' num2str(Yall_enum(j)) '-' num2str(Yall_sw(j))];
    %    lt_plot_text(1.05*max(dat), y, str, 'm', 8);
    
    Ylabstr = [Ylabstr; str];
end
[~, indtmp] = sort(indrank);
Ylabstr = Ylabstr(indtmp);

set(gca, 'YTick', 1:length(Ylabstr), 'YTickLabel', Ylabstr);
lt_plot_zeroline_vert;





%% ======== ONE DATAPOINT PER CHANNEL -

Yall_scal = {};
Yall_bnum = [];
Yall_enum = [];
Yall_chanID = [];
Yall_istarg = {};
Yall_issame = {};
Yall_sw = [];
for i=1:length(indsgrp_chanpair_unique)
    
    indsthis = indsgrp_chanpair == indsgrp_chanpair_unique(i);
    
    bnum = unique(OUTSTRUCT.bnum(indsthis));
    enum = unique(OUTSTRUCT.enum(indsthis));
    swnum = unique(OUTSTRUCT.switch(indsthis));
    
    bname = MOTIFSTATS_Compiled.birds(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    
    % =========================== DATA
    cohscal = Yscalar(indsthis);
    istarg = OUTSTRUCT.istarg(indsthis);
    issame = OUTSTRUCT.issame(indsthis);
    motifID = lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)); % ---- get positions within global motif
    chnums = OUTSTRUCT.chanpair_actual(indsthis,:);
    neurpair = OUTSTRUCT.neurpair(indsthis,:);
    
    % ================================ SAVE OUTPUT
    Yall_scal = [Yall_scal; cohscal];
    Yall_bnum = [Yall_bnum ; bnum];
    Yall_enum = [Yall_enum ; enum];
    Yall_chanID = [Yall_chanID; indsgrp_chanpair_unique(i)];
    Yall_istarg = [Yall_istarg; istarg];
    Yall_issame = [Yall_issame; issame];
    Yall_sw = [Yall_sw; swnum];
end

lt_figure; hold on;
% === 1) sort in order of number of 1syls
Nall =cellfun(@length, Yall_scal);
% [~, indsort] = sort(Nall);
[~, ~, indrank] = lt_tools_sort(Nall);


%% ======= COPY THINGS FROM OLD OUTSTRUCT TO NEW OUTSTRUCT
disp('NOTE: using the exact same trial inds to compute the xcov and learning and cohscalar');

All_learn_targdir_z = nan(length(OUTSTRUCT.bnum), 1);
All_cohdiff_z = nan(length(OUTSTRUCT.bnum), 1);
for i=1:length(OUTSTRUCT.bnum)
    
    % === FIND THE MATCHING BIRD, EXPT, AND MOTIF
    bnum = OUTSTRUCT.bnum(i);
    enum = OUTSTRUCT.enum(i);
    swnum = OUTSTRUCT.switch(i);
    mm = OUTSTRUCT.motifnum(i);
    
    inds_old = OUTSTRUCT_LFP.bnum==bnum & OUTSTRUCT_LFP.enum==enum & ...
        OUTSTRUCT_LFP.switch==swnum & OUTSTRUCT_LFP.motifnum==mm;
    assert(sum(inds_old)>0, 'why cant find?');
    assert(length(unique(cellfun('length', OUTSTRUCT_LFP.tvals(inds_old))))==1, 'why different?');
    indthis = find(inds_old);
    indthis = indthis(1);
    
    % ================ EXTRACT THINGS
    % ======= 1) LEARNING
    ffvals = OUTSTRUCT_LFP.ffvals{indthis};
    if useoldInds==1
        indsbase = OUTSTRUCT_LFP.indsbase_epoch{indthis};
        indswn = OUTSTRUCT_LFP.indsWN_epoch{indthis};
    else
       indsbase = OUTSTRUCT.inds_base_epoch{i};
       indswn = OUTSTRUCT.inds_WN_epoch{i};
    end
    % ------ DIRECTION OF LEARNING?
    learndir = OUTSTRUCT_LFP.learndirTarg(indthis);
    
    if takeindsWNsecondhalf==1
        indswn = indswn(round(length(indswn)/2):end);
    end
    
    learn_targdir_z = (mean(ffvals(indswn))-mean(ffvals(indsbase)))./std(ffvals(indsbase)); % zscore
    learn_targdir_z = learn_targdir_z*learndir;
    
    All_learn_targdir_z(i) = learn_targdir_z;
    
    % ========= COHERENCE CHANGE (z-score) (take average over all chan
    % pairs
    inds_old = find(inds_old);
    tmp = [];
    for ii=1:length(inds_old)
           cohscal = OUTSTRUCT_LFP.cohscal{inds_old(ii)};
           cohdiff_z = (mean(cohscal(indswn)) - mean(cohscal(indsbase)))./std(cohscal(indsbase));
%            cohdiff_z = (mean(cohscal(indswn)) - mean(cohscal(indsbase)));
           tmp = [tmp cohdiff_z]; 
    end
    All_cohdiff_z(i) = mean(tmp);
    
    
    % sanity check...
%     cohscal = OUTSTRUCT_LFP.cohscal{indthis}
%     mean(cohscal(OUTSTRUCT_LFP.indsbase_epoch{indthis})) - mean(cohscal(OUTSTRUCT_LFP.indsWN_epoch{indthis}))
%     OUTSTRUCT_LFP.cohscal_diff(indthis)
end

OUTSTRUCT.All_learn_targdir_z = All_learn_targdir_z;
OUTSTRUCT.All_cohdiff_z = All_cohdiff_z;


%% ======= EXTRACT OTHER THINGS (This is better code, but confirmed that aligns perfectly with previous

% ====================== GET BOTH WINDOWS
% ========= SCALAR 1
windthis = 1;
OUTSTRUCT.yscal1 = cellfun(@(x)(x(windthis,2)-x(windthis,1)), OUTSTRUCT.Xcovscal_window_BaseWN);

% ========= SCALAR 2
if size(OUTSTRUCT.Xcovscal_window_BaseWN{1},1)==2
windthis = 2;
OUTSTRUCT.yscal2 = cellfun(@(x)(x(windthis,2)-x(windthis,1)), OUTSTRUCT.Xcovscal_window_BaseWN);
else
    OUTSTRUCT.yscal2 = OUTSTRUCT.yscal1;
    disp('NOTE: window 2 scalar doesnt exist, so copied exactly from scalar 1 [PAUSED]');
    pause;
end
% ========= LAERNING
if strcmp(plotlevel, 'chanpair')
    [ ~, ~, ~, All_yscal1, ~, ~, ~, ~] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT, 'yscal1');
    [ ~, ~, ~, All_yscal2, ~, ~, ~, ~] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT, 'yscal2');
    [ ~, ~, ~, All_learn_targdir_z, ~, ~, ~, ~] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT, 'All_learn_targdir_z');
    [ ~, ~, ~, All_cohdiff_z, ~, ~, ~, ~] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT, 'All_cohdiff_z');
elseif strcmp(plotlevel, 'switch')
    [ ~, ~, ~, ~, ~, ~, ~, All_yscal1] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT, 'yscal1');
    [ ~, ~, ~, ~, ~, ~, ~, All_yscal2] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT, 'yscal2');
    [ ~, ~, ~, ~, ~, ~, ~, All_learn_targdir_z] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT, 'All_learn_targdir_z');
    [ ~, ~, ~, ~, ~, ~, ~, All_cohdiff_z] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT, 'All_cohdiff_z');
end

    

%% === sanity check, plot coherenec change

lt_figure; hold on;
y = squeeze(All_cohdiff_z(1,1,1,:) - nanmean(All_cohdiff_z(1,1,[2:3],:), 3)); % targ minus mean(same,diff), coherence change from baseline
plot(y, 'ok');
title('change in coherence?');
ylabel('coh change (z) (targ - [same,diff]))');

%% ============= CORRELATION BETWEEN THE TWO SCALARS? BETWEEN SCALARS AN COHERENCE?
lt_figure; hold on;

lt_subplot(3,2,1); hold on;
title('[TARG] WN-base, corr scalar');
indthis = 1;

xlabel('first window');
ylabel('second window');
x = squeeze(All_yscal1(1,1,indthis, :));
y = squeeze(All_yscal2(1,1,indthis, :));
lt_regress(y, x, 1, 0, 1, 1, 'k');
% plot(x,y, 'ok');

lt_subplot(3,2,2); hold on;
title('TARG');
xlabel('Change in spike corr (window 1)');
ylabel('Change in coherence (z rel base)');
x = squeeze(All_yscal1(1,1,1,:));
y = squeeze(All_cohdiff_z(1,1,1,:));
lt_regress(y, x, 1, 0, 1, 1, 'k');


lt_subplot(3,2,3); hold on;
title('TARG');
xlabel('Change in spike corr (window 2)');
ylabel('Change in coherence (z rel base)');
x = squeeze(All_yscal2(1,1,1,:));
y = squeeze(All_cohdiff_z(1,1,1,:));
lt_regress(y, x, 1, 0, 1, 1, 'k');


lt_subplot(3,2,4); hold on;
title('TARG');
xlabel('Change in spike corr (mean of wind 1 and 2)');
ylabel('Change in coherence (z rel base)');
x = mean([squeeze(All_yscal1(1,1,1,:)) squeeze(All_yscal2(1,1,1,:))], 2);
y = squeeze(All_cohdiff_z(1,1,1,:));
lt_regress(y, x, 1, 0, 1, 1, 'k');


%% ============= CORRELATION BETWEEN THE TWO SCALARS? BETWEEN SCALARS AN COHERENCE?
% TARGET MINUS DIFF TYPE
lt_figure; hold on;

lt_subplot(3,2,1); hold on;
title('[TARG-DIFF] WN-base, corr scalar');
indthis = 1;

xlabel('first window');
ylabel('second window');
x = squeeze(All_yscal1(1,1,indthis, :) - All_yscal1(1,1,3, :));
y = squeeze(All_yscal2(1,1,indthis, :) - All_yscal2(1,1,3, :));
lt_regress(y, x, 1, 0, 1, 1, 'k');
% plot(x,y, 'ok');
    
lt_subplot(3,2,2); hold on;
title('TARG-DIFF');
xlabel('Change in spike corr (window 1)');
ylabel('Change in coherence (z rel base)');
x = squeeze(All_yscal1(1,1,1,:)-All_yscal1(1,1,3,:));
y = squeeze(All_cohdiff_z(1,1,1,:) - All_cohdiff_z(1,1,3,:));
lt_regress(y, x, 1, 0, 1, 1, 'k');


lt_subplot(3,2,3); hold on;
title('TARG-DIFF');
xlabel('Change in spike corr (window 2)');
ylabel('Change in coherence (z rel base)');
x = squeeze(All_yscal2(1,1,1,:) - All_yscal2(1,1,3,:));
y = squeeze(All_cohdiff_z(1,1,1,:) - All_cohdiff_z(1,1,3,:));
lt_regress(y, x, 1, 0, 1, 1, 'k');


lt_subplot(3,2,4); hold on;
title('TARG-DIFF');
xlabel('Change in spike corr (mean of wind 1 and 2)');
ylabel('Change in coherence (z rel base)');
x1 = mean([squeeze(All_yscal1(1,1,1,:)) squeeze(All_yscal2(1,1,1,:))], 2);
x2 = mean([squeeze(All_yscal1(1,1,3,:)) squeeze(All_yscal2(1,1,3,:))], 2);
x = x1-x2;
y = squeeze(All_cohdiff_z(1,1,1,:) - All_cohdiff_z(1,1,3,:));
lt_regress(y, x, 1, 0, 1, 1, 'k');


%% ==========


if strcmp(plotlevel, 'chanpair')
    % ============== 1) for each channel, collapse across syl types.
    OUTSTRUCT.Yscalar = Yscalar;
    [All_bnum, allenum, All_swnum, allDat, ~, ~, ~, ~, ...
        All_Nmotifs, ~] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT, 'Yscalar');
    %     [~,~,~,~, All_bnum, allenum, All_swnum, allDat] = ...
    %         lt_neural_LFP_GrpStats(OUTSTRUCT, 'Yscalar');
    
    % ======= convert to correct variables.
    Yall = squeeze(allDat)';
    All_bname = {SwitchStruct.bird(All_bnum).birdname}';
    
    % ==== collect expt names and laern dir
    All_ename = {};
    All_learndir = [];
    for i=1:length(allenum)
        
        indstmp = OUTSTRUCT.bnum==All_bnum(i) & OUTSTRUCT.enum==allenum(i) ...
            & OUTSTRUCT.switch==All_swnum(i);
        
        learndir = unique(OUTSTRUCT.learndirTarg(indstmp));
        ename = SwitchStruct.bird(All_bnum(i)).exptnum(allenum(i)).exptname;
        
        % ====
        All_ename = [All_ename; ename];
        All_learndir = [All_learndir; learndir];
    end
else
    OUTSTRUCT.Yscalar = Yscalar;
    [~, ~, ~, ~, ~, ~, ~, ~, ...
        ~, All_Nmotifs] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT, 'Yscalar');
end

All_Nmotifs = squeeze(All_Nmotifs)';
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
xthis = colsthis;
bnumthis = All_bnum(indsthis);
% ------ 1) PLOT ALL CHANS
lt_plot(colsthis+0.1, mean(ythis,1), {'Errors', lt_sem(ythis), 'Color', 'r'});
% ------ 2) comaprisons (stats)
[p] = signrank(ythis(:,1), ythis(:,2));
% [~, p] = ttest(ythis(:,1), ythis(:,2));
lt_plot_pvalue(p, 'srank (1vs2)', 1);
for i=1:size(ythis,1)
    plot(xthis, ythis(i,:), '-o', 'Color', pcols{bnumthis(i)});
end
[p] = signrank(ythis(:,2), ythis(:,3));
% [~, p] = ttest(ythis(:,2), ythis(:,3));
lt_plot_pvalue(p, 'srank (2vs3)', 1);
for i=1:size(ythis,1)
    plot(xthis, ythis(i,:), '-o', 'Color', pcols{bnumthis(i)});
end
% ------ 2) OVERLAY EACH BIRD MEAN
[ymean, ysem] = grpstats(ythis, bnumthis, {'mean', 'sem'});
bnumtmp = unique(bnumthis);
for i=1:size(ymean,1)
    lt_plot(xthis+0.7*rand-0.3, ymean(i,:), {'Errors', ysem(i,:), 'LineStyle', '-'});
    % -- within bird p val
%     [~, p] = ttest(ythis(bnumthis==bnumtmp(i), 1), ythis(bnumthis==bnumtmp(i), 2));
    [p] = signrank(ythis(bnumthis==bnumtmp(i), 1), ythis(bnumthis==bnumtmp(i), 2));
    
    % --- difference from zero
    [p1] = signrank(ythis(bnumthis==bnumtmp(i), 1));
   [p2] = signrank(ythis(bnumthis==bnumtmp(i), 2));
    % -- plot bird info
    lt_plot_text(xthis(end)+1, ymean(i, end), ['bnum ' num2str(bnumtmp(i)) ',p1=' num2str(p1) 'p2=' num2str(p2) ',p(v2)=' num2str(p)]);
end
xlim([0 5]);

% ========== TARG VS DIFF
lt_subplot(3,2,2); hold on;
xlabel('TARG - DIFF');
ylabel('coh (WN - base)');
colsthis = [1 3];

indsthis = all(~isnan(Yall(:,colsthis))');
ythis = Yall(indsthis,colsthis);
xthis = colsthis;
bnumthis = All_bnum(indsthis);
% ------ 1) PLOT ALL CHANS
lt_plot(colsthis+0.1, mean(ythis,1), {'Errors', lt_sem(ythis), 'Color', 'r'});
% ------ 2) comaprisons (stats)
[p] = signrank(ythis(:,1), ythis(:,2));
% [~, p] = ttest(ythis(:,1), ythis(:,2));
lt_plot_pvalue(p, 'srank (1vs2)', 1);
for i=1:size(ythis,1)
    plot(xthis, ythis(i,:), '-o', 'Color', pcols{bnumthis(i)});
end
% ------ 2) OVERLAY EACH BIRD MEAN
[ymean, ysem] = grpstats(ythis, bnumthis, {'mean', 'sem'});
bnumtmp = unique(bnumthis);
for i=1:size(ymean,1)
    lt_plot(xthis+0.7*rand-0.3, ymean(i,:), {'Errors', ysem(i,:), 'LineStyle', '-'});
    % -- within bird p val
%     [~, p] = ttest(ythis(bnumthis==bnumtmp(i), 1), ythis(bnumthis==bnumtmp(i), 2));
    [p] = signrank(ythis(bnumthis==bnumtmp(i), 1), ythis(bnumthis==bnumtmp(i), 2));
    [p1] = signrank(ythis(bnumthis==bnumtmp(i), 1));
   [p2] = signrank(ythis(bnumthis==bnumtmp(i), 2));
    % -- plot bird info
    lt_plot_text(xthis(end)+1, ymean(i, end), ['bnum ' num2str(bnumtmp(i)) ',p1=' num2str(p1) 'p2=' num2str(p2) ',p(v2)=' num2str(p)]);
end
xlim([0 5]);


% % ################################## SAME PLOT, PLOT ALL DATA
% lt_subplot(3,2,3); hold on;
% xlabel('SAME - TARG - DIFF');
% 
% % ================================== 1) PLOT THOSE THAT HAVE ALL 3
% colsthis = 1:3;
% xthis = [2 1 3];
% 
% indsthis = all(~isnan(Yall(:,colsthis))');
% ythis = Yall(indsthis,colsthis);
% bnumthis = All_bnum(indsthis);
% 
% % -- sort so that is in order of x
% [~, indstmp] = sort(xthis);
% xthis = xthis(indstmp);
% ythis = ythis(:, indstmp);
% 
% 
% % ----- 3) overlay all chans
% for i=1:size(ythis,1)
% %     plot(xthis, ythis(i,:), '-o', 'Color', pcols{bnumthis(i)});
%     plot(xthis, ythis(i,:), '-ok');
% end
% 
% 
% % ================================== 1) PLOT THOSE THAT ONLY HAVE 2
% xthis = [2 3];
% 
% indsthis = ~all(~isnan(Yall(:,colsthis))');
% ythis = Yall(indsthis, [1 3]);
% bnumthis = All_bnum(indsthis);
% 
% 
% % ----- 3) overlay all chans
% for i=1:size(ythis,1)
% %     plot(xthis, ythis(i,:), '-o', 'Color', pcols{bnumthis(i)});
%     plot(xthis, ythis(i,:), '-ok');
% end
% 
% % ------ 1) PLOT MEAN
% % lt_plot(xthis+0.1, mean(ythis,1), {'Errors', lt_sem(ythis), 'Color', 'r'});
% lt_plot_bar(xthis, mean(ythis,1), {'Errors', lt_sem(ythis), 'Color', 'r'});
% 
% % ------ 2) OVERLAY EACH BIRD MEAN
% [ymean, ysem] = grpstats(ythis, bnumthis, {'mean', 'sem'});
% bnumtmp = unique(bnumthis);
% for i=1:size(ymean,1)
%     lt_plot(xthis+0.7*rand-0.3, ymean(i,:), {'Errors', ysem(i,:), 'LineStyle', '-'});
%     % -- within bird p val
% %     [~, p] = ttest(ythis(bnumthis==bnumtmp(i), 1), ythis(bnumthis==bnumtmp(i), 2));
%     [p] = signrank(ythis(bnumthis==bnumtmp(i), 1), ythis(bnumthis==bnumtmp(i), 2));
%     
%     % --- difference from zero
%     [p1] = signrank(ythis(bnumthis==bnumtmp(i), 1));
%    [p2] = signrank(ythis(bnumthis==bnumtmp(i), 2));
%     % -- plot bird info
%     lt_plot_text(xthis(end)+1, ymean(i, end), ['bnum ' num2str(bnumtmp(i)) ',p1=' num2str(p1) 'p2=' num2str(p2) ',p(vs same)=' num2str(p)]);
% end
% 


% 
% % ------ 2) comaprisons (stats)
% [p] = signrank(ythis(:,1), ythis(:,2));
% % [~, p] = ttest(ythis(:,1), ythis(:,2));
% lt_plot_pvalue(p, 'srank (targ vs. same)', 1);
% [p] = signrank(ythis(:,2), ythis(:,3));
% % [~, p] = ttest(ythis(:,2), ythis(:,3));
% lt_plot_pvalue(p, 'srank (same vs. diff)', 1);







% % ===== only those with 2 (targ diff)
% lt_subplot(3,2,4); hold on;
% xlabel('TARG - DIFF');
% ylabel('coh (WN - base)');
% 
% colsthis = [1 3];
% 
% indsthis = all(~isnan(Yall(:,colsthis))');
% ythis = Yall(indsthis,colsthis);
% xthis = colsthis;
% plot(xthis, ythis', 'o-k');
% xlim([0 5]);
% lt_plot(colsthis+0.1, mean(ythis,1), {'Errors', lt_sem(ythis), 'Color', 'r'});
% if length(colsthis)==2
%     [p] = signrank(ythis(:,1), ythis(:,2));
%     lt_plot_pvalue(p, 'srank', 1);
% end
% bnumthis = All_bnum(indsthis);
% for i=1:size(ythis,1)
%     plot(xthis, ythis(i,:), '-o', 'Color', pcols{bnumthis(i)});
% end


% =================== RELATIVE TO NON TARG
if size(Yall,2)>3
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
    
    ythis = Yall(indsthis,colsthis);
    xthis = colsthis + XSHIFT;
    plot(xthis, ythis', 'o-k');
    xlim([0 5+XSHIFT]);
    if size(ythis,1)>1po
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
    lt_plot_text(xthis(end)+1, ymean(i, end), ['bnum ' num2str(bnumtmp(i))]);
end
xlim([0 5]);

% =====
lt_subplot(3,2,6); hold on;
xlabel('TRAIN UP/DOWN [lines: TARG--DIFF]');
ylabel('coh (WN - base)');

colsthis = [1 3];

% ==== UP
indsthis = all(~isnan(Yall(:,colsthis))')' & All_learndir==1;
XSHIFT = 0;

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

%% ====== [PLOTS MORE]
figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

binsize = 0.035;

% ===== only those with 2 (targ diff)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('Targ, corr (WN - base)');

y = Yall(:,1);
tmp = round(100*max(abs(y)))/100;
tmp = tmp+(binsize-mod(tmp, binsize));
% tmp = tmp+binsize/2;
x = -tmp:binsize:tmp;
lt_plot_histogram(y, x, 1, 1, '', 1, 'k');
[p] = signrank(y);
lt_plot_pvalue(p, 'srank', 1);
lt_plot_text(0, 0, ['n=' num2str(length(y))], 'm');
YLIM = ylim;
% --- overlay each birdf
for i=1:max(All_bnum)
    y = Yall(All_bnum==i, 1);
    if isempty(y)
        continue
    end
    ymean = mean(y);
    ysem = lt_sem(y);
    p = signrank(y);
    
    lt_plot(ymean, 0.95*YLIM(2), {'Xerrors', ysem});
    lt_plot_text(ymean, YLIM(2), ['b' num2str(i) ', p=' num2str(p)], 'r');
    
end
xlim([-0.3 0.3]);
lt_plot_zeroline_vert;


% ===== DIFF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('DIFF, corr (WN - base)');

y = cell2mat(Yall_allsyls(:,3));
% y = Yall(:,1);
tmp = round(100*max(abs(y)))/100;
tmp = tmp+(binsize-mod(tmp, binsize));
% tmp = tmp+binsize/2;
x = -tmp:binsize:tmp;
lt_plot_histogram(y, x, 1, 1, '', 1, 'r');
[p] = signrank(y);
lt_plot_pvalue(p, 'srank', 1);
lt_plot_text(0, 0, ['n=' num2str(length(y))], 'm');
YLIM = ylim;
% --- overlay each birdf
for i=1:max(All_bnum)
    y = Yall(All_bnum==i, 1);
    if isempty(y)
        continue
    end
    ymean = mean(y);
    ysem = lt_sem(y);
    p = signrank(y);
    
    lt_plot(ymean, 0.95*YLIM(2), {'Xerrors', ysem});
    lt_plot_text(ymean, YLIM(2), ['b' num2str(i) ', p=' num2str(p)], 'r');
    
end
xlim([-0.3 0.3]);
lt_plot_zeroline_vert;


% ===== only those with 2 (targ diff)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title(['only if Ndiff >= ' num2str(minDiffN)]);
xlabel('TARG - DIFF');
ylabel('coh (WN - base)');

colsthis = [1 3];
indsthis = all(~isnan(Yall(:,colsthis))') & All_Nmotifs(:,3)'>=minDiffN';
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


% ==================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title(['only if Ndiff >= ' num2str(minDiffN)]);
xlabel('TARG - DIFF');
ylabel('coh (WN - base)');
colsthis = [1 3];

indsthis = all(~isnan(Yall(:,colsthis))') & All_Nmotifs(:,3)'>=minDiffN';
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


% ========== TARG VS DIFF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title(['only if Ndiff >= ' num2str(minDiffN)]);
xlabel('TARG - DIFF');
ylabel('coh (WN - base)');
colsthis = [1 3];

indsthis = all(~isnan(Yall(:,colsthis))') & All_Nmotifs(:,3)'>=minDiffN';
% indsthis = all(~isnan(Yall(:,colsthis))');
ythis = Yall(indsthis,colsthis);
xthis = colsthis;
bnumthis = All_bnum(indsthis);
% ------ 1) PLOT ALL CHANS
lt_plot(colsthis+0.1, mean(ythis,1), {'Errors', lt_sem(ythis), 'Color', 'r'});
% ------ 2) comaprisons (stats)
[p] = signrank(ythis(:,1), ythis(:,2));
% [~, p] = ttest(ythis(:,1), ythis(:,2));
lt_plot_pvalue(p, 'srank (1vs2)', 1);
for i=1:size(ythis,1)
    plot(xthis, ythis(i,:), '-o', 'Color', pcols{bnumthis(i)});
end
% ------ 2) OVERLAY EACH BIRD MEAN
[ymean, ysem] = grpstats(ythis, bnumthis, {'mean', 'sem'});
bnumtmp = unique(bnumthis);
for i=1:size(ymean,1)
    lt_plot(xthis+0.7*rand-0.3, ymean(i,:), {'Errors', ysem(i,:), 'LineStyle', '-'});
    % -- within bird p val
%     [~, p] = ttest(ythis(bnumthis==bnumtmp(i), 1), ythis(bnumthis==bnumtmp(i), 2));
    [p] = signrank(ythis(bnumthis==bnumtmp(i), 1), ythis(bnumthis==bnumtmp(i), 2));
    [p1] = signrank(ythis(bnumthis==bnumtmp(i), 1));
   [p2] = signrank(ythis(bnumthis==bnumtmp(i), 2));
    % -- plot bird info
    lt_plot_text(xthis(end)+1, ymean(i, end), ['bnum ' num2str(bnumtmp(i)) ',p1=' num2str(p1) 'p2=' num2str(p2) ',p(v2)=' num2str(p)]);
end
xlim([0 5]);





%% ############## [HISTOGRAMS, lowest level data, chan x syl]
figcount=1;
subplotrows=1;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


% ############################# 1) TARG -- DIFF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('Targ -- DIFF');
ylabel('xcov, change');
title('dat = chan x syl');
x = [1 3];

lt_neural_POPLEARN_XCov_PlotScal_sub2;


% ############################# 1) TARG -- SAME
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('Targ -- SAME -- DIFF');
ylabel('xcov, change');
title('dat = chan x syl');
x = [1 2 3];

lt_neural_POPLEARN_XCov_PlotScal_sub2;

%% ########## [MORE PLOTS] NOT DIRECTLY REALTED TO YSCALAR.
% ==== CORREALTING SCALARS WITH BEHAVIOR

figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];




% ================== RELATED TO ELARNIG?
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('TARG');
indthis = 1; % targ, same,. diff
xlabel('learning (z, targ dir)');
ylabel('WN-base, corr scalar(window 1)');

x = squeeze(All_learn_targdir_z(1,1,indthis,:));
y = squeeze(All_yscal1(1,1,indthis,:));
lt_regress(y, x, 1, 0, 1, 1);
% lt_regress(y, c, 1, 0, 1, 1)

% ================== RELATED TO ELARNIG?
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('TARG');
indthis = 1; % targ, same,. diff
xlabel('learning (z, targ dir)');
ylabel('WN-base, corr scalar(window 2)');

x = squeeze(All_learn_targdir_z(1,1,indthis,:));
y = squeeze(All_yscal2(1,1,indthis,:));
lt_regress(y, x, 1, 0, 1, 1);


% ================== RELATED TO ELARNIG?
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('TARG');
indthis = 1; % targ, same,. diff
xlabel('learning (z, targ dir)');
ylabel('WN-base, corr scalar(mean of wind 1 and 2)');

x = squeeze(All_learn_targdir_z(1,1,indthis,:));
y = mean([squeeze(All_yscal1(1,1,indthis,:)) squeeze(All_yscal2(1,1,indthis,:))], 2);
lt_regress(y, x, 1, 0, 1, 1);



% ================== RELATED TO ELARNIG?
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('TARG-DIFF');
indthis = 1; % targ, same,. diff
xlabel('learning (z, targ dir)');
ylabel('WN-base, corr scalar(window 1)');

x = squeeze(All_learn_targdir_z(1,1,indthis,:) - All_learn_targdir_z(1,1,3,:));
y = squeeze(All_yscal1(1,1,indthis,:) - All_yscal1(1,1,3,:));
lt_regress(y, x, 1, 0, 1, 1);
% lt_regress(y, c, 1, 0, 1, 1)

% ================== RELATED TO ELARNIG?
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('TARG-DIFF');
indthis = 1; % targ, same,. diff
xlabel('learning (z, targ dir)');
ylabel('WN-base, corr scalar(window 2)');

x = squeeze(All_learn_targdir_z(1,1,indthis,:) - All_learn_targdir_z(1,1,3,:));
y = squeeze(All_yscal2(1,1,indthis,:) - All_yscal2(1,1,3,:));
lt_regress(y, x, 1, 0, 1, 1);


% ================== RELATED TO ELARNIG?
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('TARG-DIFF');
indthis = 1; % targ, same,. diff
xlabel('learning (z, targ dir)');
ylabel('WN-base, corr scalar(mean of wind 1 and 2)');

x = squeeze(All_learn_targdir_z(1,1,indthis,:));
y1 = mean([squeeze(All_yscal1(1,1,indthis,:)) squeeze(All_yscal2(1,1,indthis,:))], 2);
y2 = mean([squeeze(All_yscal1(1,1,3,:)) squeeze(All_yscal2(1,1,3,:))], 2);
y = y1-y2;
lt_regress(y, x, 1, 0, 1, 1);



%% ====== change in coherence correlate with learning?

figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];




% ================== RELATED TO ELARNIG?
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('TARG');
indthis = 1; % targ, same,. diff
xlabel('learning (z, targ dir)');
ylabel('WN-base, coh (z)');

x = squeeze(All_learn_targdir_z(1,1,indthis,:));
y = squeeze(All_cohdiff_z(1,1,indthis,:));
lt_regress(y, x, 1, 0, 1, 1);
% lt_regress(y, c, 1, 0, 1, 1)


% ================== RELATED TO ELARNIG?
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('TARG-DIFF');
indthis = 1; % targ, same,. diff
xlabel('learning (z, targ dir)');
ylabel('WN-base, coh (z)');

x = squeeze(All_learn_targdir_z(1,1,indthis,:) - All_learn_targdir_z(1,1,3,:));
y = squeeze(All_cohdiff_z(1,1,indthis,:) - All_cohdiff_z(1,1,3,:));
lt_regress(y, x, 1, 0, 1, 1);
% lt_regress(y, c, 1, 0, 1, 1)
