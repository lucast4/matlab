function DATSTRUCT_POP = lt_neural_POPLEARN_SylLocked_CoorLearn(DATSTRUCT_SPK, DATSTRUCT_LFP, ...
    PARAMS, OUTSTRUCT_XCOV, SwitchStruct)
%% lt 2/25/19 - spk-spk coorelation (WITHIN AREA) CHANGE DURING ELARNIGN?
% USING SMOOTHED FR (LIKE HAMISH STUFF)

%% ===
disp('NOTE: any cases with nan for correlation I replace by median acros tirlas (for that neuron pair x syllable)');

% =========== PARAMS FOR COHERENCE
ntapers = []; % leave [] to set as default.
movingwin = [0.1 0.01]; % leave exmpty for default. [applies for both welches and multiutaper]
tw = [];

fwind = [22 36]; % window to get scalar;

% ==== match movingwin to size of data
wintmp = size(DATSTRUCT_LFP.LFP_dat{1},1)/1500;
disp(['changin cohernece window from: ' num2str(movingwin) ', to: ' num2str([wintmp movingwin(2)])]);
movingwin(1) = wintmp;

%% ============= for each experiment, get coordination in spiking
[indsgrp, indsgrpU] = lt_tools_grp2idx({DATSTRUCT_SPK.bnum, DATSTRUCT_SPK.enum, ...
    DATSTRUCT_SPK.switch, DATSTRUCT_SPK.motifnum});

DATSTRUCT_POP = struct;
DATSTRUCT_POP.bnum = [];
DATSTRUCT_POP.enum = [];
DATSTRUCT_POP.switch = [];
DATSTRUCT_POP.motifnum = [];
DATSTRUCT_POP.bregion = {};
DATSTRUCT_POP.fr_RhoPairwise = {};
% DATSTRUCT_POP.lfp_Pow = {};
DATSTRUCT_POP.fr_MeanFr = {};

DATSTRUCT_POP.indsbase_epoch = {};
DATSTRUCT_POP.indsWN_epoch = {};

% DATSTRUCT_POP.coher_cohgram_all = {};

for i=1:length(indsgrpU)
    disp(i);
    
    % ######################################### LMAN
    bregionthis = 'LMAN';
    
    indsthis_spk = find(indsgrp==indsgrpU(i) & ...
        strcmp(DATSTRUCT_SPK.spike_bregions, bregionthis));
    
    if length(indsthis_spk)>1
        % need at least two spiking units to do this analysis.
        
        bnum = unique(DATSTRUCT_SPK.bnum(indsthis_spk));
        enum = unique(DATSTRUCT_SPK.enum(indsthis_spk));
        sw = unique(DATSTRUCT_SPK.switch(indsthis_spk));
        mm = unique(DATSTRUCT_SPK.motifnum(indsthis_spk));
        
        % --------------- COLLECT DATA
        indsbase_epoch = DATSTRUCT_SPK.inds_base_epoch{indsthis_spk(1)};
        indsWN_epoch = DATSTRUCT_SPK.inds_WN_epoch{indsthis_spk(1)};
        
        
        % ============ 1) compute all pairwise correlations
        frmat = DATSTRUCT_SPK.frmat(indsthis_spk);
        rhoall = lt_neural_POPLEARN_SylLocked_paircorr(frmat);
        
        
        % ============ 1b) collect all mean firing rates
        frmean = cellfun(@(x)mean(x,1), frmat, 'UniformOutput', 0);
        frmean = cell2mat(frmean)';
        
        
        % ================== SAVE OUTPUT
        DATSTRUCT_POP.bnum = [DATSTRUCT_POP.bnum; bnum];
        DATSTRUCT_POP.enum = [DATSTRUCT_POP.enum; enum];
        DATSTRUCT_POP.switch = [DATSTRUCT_POP.switch; sw];
        DATSTRUCT_POP.motifnum = [DATSTRUCT_POP.motifnum; mm];
        DATSTRUCT_POP.bregion = [DATSTRUCT_POP.bregion; bregionthis];
        DATSTRUCT_POP.fr_RhoPairwise = [DATSTRUCT_POP.fr_RhoPairwise; rhoall];
        DATSTRUCT_POP.fr_MeanFr = [DATSTRUCT_POP.fr_MeanFr; frmean];
        DATSTRUCT_POP.indsbase_epoch = [DATSTRUCT_POP.indsbase_epoch; indsbase_epoch];
        DATSTRUCT_POP.indsWN_epoch = [DATSTRUCT_POP.indsWN_epoch; indsWN_epoch];
        
    end
    
    
    % ######################################### LMAN
    bregionthis = 'RA';
    
    indsthis_spk = find(indsgrp==indsgrpU(i) & ...
        strcmp(DATSTRUCT_SPK.spike_bregions, bregionthis));
    
    if length(indsthis_spk)>1
        % need at least two spiking units to do this analysis.
        
        bnum = unique(DATSTRUCT_SPK.bnum(indsthis_spk));
        enum = unique(DATSTRUCT_SPK.enum(indsthis_spk));
        sw = unique(DATSTRUCT_SPK.switch(indsthis_spk));
        mm = unique(DATSTRUCT_SPK.motifnum(indsthis_spk));
        
        % --------------- COLLECT DATA
        indsbase_epoch = DATSTRUCT_SPK.inds_base_epoch{indsthis_spk(1)};
        indsWN_epoch = DATSTRUCT_SPK.inds_WN_epoch{indsthis_spk(1)};
        
        
        % ============ 1) compute all pairwise correlations
        frmat = DATSTRUCT_SPK.frmat(indsthis_spk);
        rhoall = lt_neural_POPLEARN_SylLocked_paircorr(frmat);
        
        
        % ============ 1b) collect all mean firing rates
        frmean = cellfun(@(x)mean(x,1), frmat, 'UniformOutput', 0);
        frmean = cell2mat(frmean)';
        
        
        % ================== SAVE OUTPUT
        DATSTRUCT_POP.bnum = [DATSTRUCT_POP.bnum; bnum];
        DATSTRUCT_POP.enum = [DATSTRUCT_POP.enum; enum];
        DATSTRUCT_POP.switch = [DATSTRUCT_POP.switch; sw];
        DATSTRUCT_POP.motifnum = [DATSTRUCT_POP.motifnum; mm];
        DATSTRUCT_POP.bregion = [DATSTRUCT_POP.bregion; bregionthis];
        DATSTRUCT_POP.fr_RhoPairwise = [DATSTRUCT_POP.fr_RhoPairwise; rhoall];
        DATSTRUCT_POP.fr_MeanFr = [DATSTRUCT_POP.fr_MeanFr; frmean];
        DATSTRUCT_POP.indsbase_epoch = [DATSTRUCT_POP.indsbase_epoch; indsbase_epoch];
        DATSTRUCT_POP.indsWN_epoch = [DATSTRUCT_POP.indsWN_epoch; indsWN_epoch];
        
    end
    
end



%% ================= [PROCESS] get mean correlation during base and WN

Rho_BaseWN_All = cell(length(DATSTRUCT_POP.bnum),1);
istarg_all = nan(length(DATSTRUCT_POP.bnum),1);
issame_all = nan(length(DATSTRUCT_POP.bnum),1);

for i=1:length(DATSTRUCT_POP.bnum)
    
    indsbase = DATSTRUCT_POP.indsbase_epoch{i};
    indsWN = DATSTRUCT_POP.indsWN_epoch{i};
    
    rhoall = DATSTRUCT_POP.fr_RhoPairwise{i};
    
    % ==== get baseline means
    rhobase = mean(rhoall(indsbase, :),1);
    rhoWN = mean(rhoall(indsWN, :),1);
    
    rho_BaseWn = [rhobase' rhoWN'];
    Rho_BaseWN_All{i} = rho_BaseWn;
    
    
    % ==== FIGURE OUT SYL TYPE
    bnum = DATSTRUCT_POP.bnum(i);
    enum = DATSTRUCT_POP.enum(i);
    sw = DATSTRUCT_POP.switch(i);
    mm = DATSTRUCT_POP.motifnum(i);
    indsout = find(OUTSTRUCT_XCOV.bnum==bnum & OUTSTRUCT_XCOV.enum==enum ...
        & OUTSTRUCT_XCOV.switch==sw & OUTSTRUCT_XCOV.motifnum==mm);
    istarg = unique(OUTSTRUCT_XCOV.istarg(indsout));
    issame = unique(OUTSTRUCT_XCOV.issame(indsout));
    
    istarg_all(i) = istarg;
    issame_all(i) = issame;
end

DATSTRUCT_POP.Rho_BaseWN_All = Rho_BaseWN_All;
DATSTRUCT_POP.istarg_all = istarg_all;
DATSTRUCT_POP.issame_all = issame_all;

%% ================= [PLOT] EACH SWITCH SEPARATELY
[indsgrp, indsgrpU] = lt_tools_grp2idx({DATSTRUCT_POP.bnum, DATSTRUCT_POP.enum, ...
    DATSTRUCT_POP.switch});

bregionthis = 'LMAN';


% ====================== PLOT
figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

if strcmp(bregionthis, 'LMAN')
    pcol = [0.2 0.5 0.2];
elseif strcmp(bregionthis, 'RA')
    pcol = 'r';
end

    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(bregionthis);

Yall_base = nan(length(indsgrpU), 3); % switch x [targ, same, diff]
Yall_wn = nan(length(indsgrpU), 3); % switch x [targ, same, diff]

for i=1:length(indsgrpU)
    
    indsthis= find(indsgrp==indsgrpU(i));
    
    bnum = unique(DATSTRUCT_POP.bnum(indsthis));
    enum = unique(DATSTRUCT_POP.enum(indsthis));
    sw = unique(DATSTRUCT_POP.switch(indsthis));
    
    bname = SwitchStruct.bird(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(sw)]);
    xlabel('Base [unit-unit corr]');
    ylabel('WN [dat = motif x chanpair]');
    
    indsthis= find(indsgrp==indsgrpU(i) & strcmp(DATSTRUCT_POP.bregion, bregionthis));
    if isempty(indsthis)
        lt_plot_text(0, 0.5, 'no dat');
        continue
    end
    
    % ========== TARG SYL
    pcolsyl = 'k';
    indsthis= find(indsgrp==indsgrpU(i) & DATSTRUCT_POP.istarg_all==1 ...
        & strcmp(DATSTRUCT_POP.bregion, bregionthis));
    
    rhodat = cell2mat(DATSTRUCT_POP.Rho_BaseWN_All(indsthis));
    
    plot(rhodat(:,1), rhodat(:,2), 'o', 'Color', pcolsyl);
    
    % --- save output
    Yall_base(i, 1) = mean(rhodat(:,1));
    Yall_wn(i, 1) = mean(rhodat(:,2));
    
    % ========== SAME SYL
    pcolsyl = 'b';
    indsthis= find(indsgrp==indsgrpU(i) & DATSTRUCT_POP.istarg_all==0 & ...
        DATSTRUCT_POP.issame_all==1 & strcmp(DATSTRUCT_POP.bregion, bregionthis));
    if ~isempty(indsthis)
        rhodat = cell2mat(DATSTRUCT_POP.Rho_BaseWN_All(indsthis));
        
        plot(rhodat(:,1), rhodat(:,2), 'o', 'Color', pcolsyl);
        
        % --- save output
        Yall_base(i, 2) = mean(rhodat(:,1));
        Yall_wn(i, 2) = mean(rhodat(:,2));
    end
    
    % ========== DIFF SYL
    pcolsyl = 'r';
    indsthis= find(indsgrp==indsgrpU(i) & DATSTRUCT_POP.istarg_all==0 & ...
        DATSTRUCT_POP.issame_all==0 & strcmp(DATSTRUCT_POP.bregion, bregionthis));
    
    rhodat = cell2mat(DATSTRUCT_POP.Rho_BaseWN_All(indsthis));
    
    plot(rhodat(:,1), rhodat(:,2), 'o', 'Color', pcolsyl);
    % --- save output
    Yall_base(i, 3) = mean(rhodat(:,1));
    Yall_wn(i, 3) = mean(rhodat(:,2));
    
    % ========== format plot
    lt_plot_makesquare_plot45line(gca, 'k', -0.2);
end


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('all exp[eriments (each one line)');
ylabel('(rho, all unit pairs)');
xlabel('TARG -- SAME -- DIFF');

indsthis = ~any(isnan(Yall_base)');
ybase = Yall_base(indsthis,:);
ywn = Yall_wn(indsthis,:);
x = [1 2 3];

plot(x-0.2, ybase', '-ok');
plot(x+0.2, ywn', '-or');
lt_plot_zeroline;



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('all exp[eriments (each one line)');
ylabel('(rho, all unit pairs)');
xlabel('TARG -- SAME -- DIFF [BASE-WN]');

indsthis = ~any(isnan(Yall_base)');
ybase = Yall_base(indsthis,:);
ywn = Yall_wn(indsthis,:);

x = [1 2 3];

for i=1:length(x)
    xthis = [x(i)-0.2 x(i)+0.2];
    ythis = [ybase(:,i) ywn(:,i)];
    plot(xthis, ythis, '-ko');
    if size(ythis,1)>1
    lt_plot(xthis+0.2, mean(ythis,1), {'Errors', lt_sem(ythis)});
    end
end
lt_plot_zeroline;



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('all exp[eriments (each one line)');
ylabel('(rho, all unit pairs)');
xlabel('TARG -- SAME -- DIFF');
colsget = [1 3];
indsthis = ~any(isnan(Yall_base(:, colsget)'));
ybase = Yall_base(indsthis,colsget);
ywn = Yall_wn(indsthis,colsget);
x = colsget;

plot(x-0.2, ybase', '-ok');
plot(x+0.2, ywn', '-or');
lt_plot_zeroline;



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('all exp[eriments (each one line)');
ylabel('(rho, all unit pairs)');
xlabel('TARG -- SAME -- DIFF [BASE-WN]');

indsthis = ~any(isnan(Yall_base(:, colsget)'));
ybase = Yall_base(indsthis,colsget);
ywn = Yall_wn(indsthis,colsget);
x = colsget;

for i=1:length(x)
    xthis = [x(i)-0.2 x(i)+0.2];
    ythis = [ybase(:,i) ywn(:,i)];
    plot(xthis, ythis, '-ko');
    lt_plot(xthis+0.2, mean(ythis,1), {'Errors', lt_sem(ythis)});
end
lt_plot_zeroline;


%% ================= [PLOT] EACH SWITCH SEPARATELY
[indsgrp, indsgrpU] = lt_tools_grp2idx({DATSTRUCT_POP.bnum, DATSTRUCT_POP.enum, ...
    DATSTRUCT_POP.switch});

bregionthis = 'RA';


% ====================== PLOT
figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

if strcmp(bregionthis, 'LMAN')
    pcol = [0.2 0.5 0.2];
elseif strcmp(bregionthis, 'RA')
    pcol = 'r';
end


    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(bregionthis);

Yall_base = nan(length(indsgrpU), 3); % switch x [targ, same, diff]
Yall_wn = nan(length(indsgrpU), 3); % switch x [targ, same, diff]

for i=1:length(indsgrpU)
    
    indsthis= find(indsgrp==indsgrpU(i));
    
    bnum = unique(DATSTRUCT_POP.bnum(indsthis));
    enum = unique(DATSTRUCT_POP.enum(indsthis));
    sw = unique(DATSTRUCT_POP.switch(indsthis));
    
    bname = SwitchStruct.bird(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(sw)]);
    xlabel('Base [unit-unit corr]');
    ylabel('WN [dat = motif x chanpair]');
    
    indsthis= find(indsgrp==indsgrpU(i) & strcmp(DATSTRUCT_POP.bregion, bregionthis));
    if isempty(indsthis)
        lt_plot_text(0, 0.5, 'no dat');
        continue
    end
    
    % ========== TARG SYL
    pcolsyl = 'k';
    indsthis= find(indsgrp==indsgrpU(i) & DATSTRUCT_POP.istarg_all==1 ...
        & strcmp(DATSTRUCT_POP.bregion, bregionthis));
    
    rhodat = cell2mat(DATSTRUCT_POP.Rho_BaseWN_All(indsthis));
    
%     plot(rhodat(:,1), rhodat(:,2), 'o', 'Color', pcolsyl);
    lt_plot(rhodat(:,1), rhodat(:,2));
    
    % --- save output
    Yall_base(i, 1) = mean(rhodat(:,1));
    Yall_wn(i, 1) = mean(rhodat(:,2));
    
    % ========== SAME SYL
    pcolsyl = 'b';
    indsthis= find(indsgrp==indsgrpU(i) & DATSTRUCT_POP.istarg_all==0 & ...
        DATSTRUCT_POP.issame_all==1 & strcmp(DATSTRUCT_POP.bregion, bregionthis));
    if ~isempty(indsthis)
        rhodat = cell2mat(DATSTRUCT_POP.Rho_BaseWN_All(indsthis));
        
        plot(rhodat(:,1), rhodat(:,2), 'o', 'Color', pcolsyl);
        
        % --- save output
        Yall_base(i, 2) = mean(rhodat(:,1));
        Yall_wn(i, 2) = mean(rhodat(:,2));
    end
    
    % ========== DIFF SYL
    pcolsyl = 'r';
    indsthis= find(indsgrp==indsgrpU(i) & DATSTRUCT_POP.istarg_all==0 & ...
        DATSTRUCT_POP.issame_all==0 & strcmp(DATSTRUCT_POP.bregion, bregionthis));
    
    rhodat = cell2mat(DATSTRUCT_POP.Rho_BaseWN_All(indsthis));
    
    plot(rhodat(:,1), rhodat(:,2), 'o', 'Color', pcolsyl);
    % --- save output
    Yall_base(i, 3) = mean(rhodat(:,1));
    Yall_wn(i, 3) = mean(rhodat(:,2));
    
    % ========== format plot
    lt_plot_makesquare_plot45line(gca, 'k', -0.2);
end


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('all exp[eriments (each one line)');
ylabel('(rho, all unit pairs)');
xlabel('TARG -- SAME -- DIFF');

indsthis = ~any(isnan(Yall_base)');
ybase = Yall_base(indsthis,:);
ywn = Yall_wn(indsthis,:);
x = [1 2 3];

plot(x-0.2, ybase', '-ok');
plot(x+0.2, ywn', '-or');
lt_plot_zeroline;



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('all exp[eriments (each one line)');
ylabel('(rho, all unit pairs)');
xlabel('TARG -- SAME -- DIFF [BASE-WN]');

indsthis = ~any(isnan(Yall_base)');
ybase = Yall_base(indsthis,:);
ywn = Yall_wn(indsthis,:);

x = [1 2 3];

for i=1:length(x)
    xthis = [x(i)-0.2 x(i)+0.2];
    ythis = [ybase(:,i) ywn(:,i)];
    plot(xthis, ythis, '-ko');
    lt_plot(xthis+0.2, mean(ythis,1), {'Errors', lt_sem(ythis)});
end
lt_plot_zeroline;



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('all exp[eriments (each one line)');
ylabel('(rho, all unit pairs)');
xlabel('TARG -- SAME -- DIFF');
colsget = [1 3];
indsthis = ~any(isnan(Yall_base(:, colsget)'));
ybase = Yall_base(indsthis,colsget);
ywn = Yall_wn(indsthis,colsget);
x = colsget;

plot(x-0.2, ybase', '-ok');
plot(x+0.2, ywn', '-or');
lt_plot_zeroline;



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('all exp[eriments (each one line)');
ylabel('(rho, all unit pairs)');
xlabel('TARG -- SAME -- DIFF [BASE-WN]');

indsthis = ~any(isnan(Yall_base(:, colsget)'));
ybase = Yall_base(indsthis,colsget);
ywn = Yall_wn(indsthis,colsget);
x = colsget;

for i=1:length(x)
    xthis = [x(i)-0.2 x(i)+0.2];
    ythis = [ybase(:,i) ywn(:,i)];
    plot(xthis, ythis, '-ko');
    lt_plot(xthis+0.2, mean(ythis,1), {'Errors', lt_sem(ythis)});
end
lt_plot_zeroline;

