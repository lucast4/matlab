function lt_neural_POPLEARN_SylLocked_Coord(DATSTRUCT_SPK, DATSTRUCT_LFP, ...
    PARAMS, plotRawRand)
%% lt 2/25/19 - spk-spk coorelation related to LFP?

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


DATSTRUCT_POP.bnum = [];
DATSTRUCT_POP.enum = [];
DATSTRUCT_POP.switch = [];
DATSTRUCT_POP.motifnum = [];
DATSTRUCT_POP.bregion = {};
DATSTRUCT_POP.fr_RhoPairwise = {};
DATSTRUCT_POP.lfp_Pow = {};
DATSTRUCT_POP.fr_MeanFr = {};

DATSTRUCT_POP.coher_cohgram_all = {};

for i=1:length(indsgrpU)
    
    % ######################################### LMAN
    bregionthis = 'LMAN';
    
    indsthis_spk = find(indsgrp==indsgrpU(i) & ...
        strcmp(DATSTRUCT_SPK.spike_bregions, bregionthis));
    
    if length(indsthis_spk)<2
        % need at least two spiking units to do this analysis.
        continue
    end
    
    bnum = unique(DATSTRUCT_SPK.bnum(indsthis_spk));
    enum = unique(DATSTRUCT_SPK.enum(indsthis_spk));
    sw = unique(DATSTRUCT_SPK.switch(indsthis_spk));
    mm = unique(DATSTRUCT_SPK.motifnum(indsthis_spk));
    
    
    indsthis_lfp = find(DATSTRUCT_LFP.bnum==bnum & DATSTRUCT_LFP.enum==enum ...
        & DATSTRUCT_LFP.switch==sw & DATSTRUCT_LFP.motifnum==mm & ...
        strcmp(DATSTRUCT_LFP.LFP_bregions, bregionthis));
    
    assert(length(indsthis_spk)==length(indsthis_lfp));
    
    % --------------- COLLECT DATA
    indstrials = DATSTRUCT_SPK.inds_base{indsthis_spk(1)};
    frmat = DATSTRUCT_SPK.frmat(indsthis_spk);
    spkdat = DATSTRUCT_SPK.spike_dat(indsthis_spk);
    lfpdat = DATSTRUCT_LFP.LFP_dat(indsthis_lfp);
    assert(size(frmat{1},2) == size(lfpdat{1},2));
    
    % ---------- GET TRIALS OF INTEREST
    frmat = cellfun(@(x)x(:, indstrials), frmat, 'UniformOutput', 0);
    spkdat = cellfun(@(x)x(indstrials), spkdat, 'UniformOutput', 0);
    lfpdat = cellfun(@(x)x(:, indstrials), lfpdat, 'UniformOutput', 0);
    
    % -- save for later
    lfpdat_LMAN = lfpdat;
    spkdat_LMAN = spkdat;
    
    % ============ 1) compute all pairwise correlations
    rhoall = lt_neural_POPLEARN_SylLocked_paircorr(frmat);
    % --- take average across channel pairs
    rhoall = mean(rhoall,2);
    
    
    % ============ 1b) collect all mean firing rates
    frmean = cellfun(@(x)mean(x,1), frmat, 'UniformOutput', 0);
    frmean = cell2mat(frmean);
    frmean = mean(frmean,1)';
    
    % ============ 2) LFP POWER
    % - if multiple channels, get power for each and then take average
    ntrials = size(lfpdat{1}, 2);
    nchans = length(lfpdat);
    lfppow_all = nan(ntrials, nchans);
    for j=1:length(lfpdat)
        lfppow = sqrt(mean(lfpdat{j}.^2, 1));
        %         lfppow = std(lfpdat{j}, [], 1);
        lfppow_all(:, j) = lfppow;
    end
    % -- take average across channels
    lfppow_all = mean(lfppow_all,2);
    
    
    % ================== SAVE OUTPUT
    DATSTRUCT_POP.bnum = [DATSTRUCT_POP.bnum; bnum];
    DATSTRUCT_POP.enum = [DATSTRUCT_POP.enum; enum];
    DATSTRUCT_POP.switch = [DATSTRUCT_POP.switch; sw];
    DATSTRUCT_POP.motifnum = [DATSTRUCT_POP.motifnum; mm];
    DATSTRUCT_POP.bregion = [DATSTRUCT_POP.bregion; bregionthis];
    DATSTRUCT_POP.fr_RhoPairwise = [DATSTRUCT_POP.fr_RhoPairwise; rhoall];
    DATSTRUCT_POP.fr_MeanFr = [DATSTRUCT_POP.fr_MeanFr; frmean];
    DATSTRUCT_POP.lfp_Pow = [DATSTRUCT_POP.lfp_Pow; lfppow_all];
    DATSTRUCT_POP.coher_cohgram_all = [DATSTRUCT_POP.coher_cohgram_all; {[]}];
    
    
%     if bnum==4 & enum==5 & sw==1 & mm==9
%         keyboard
%     end
    
    
    % ####################### GET COHERENCE BETWEEN ALL PAIRS OF LFP
    
    % plot random case, low, med, and high correlation cases?
    if plotRawRand ==1
        if rand<0.15
            lt_neural_POPLEARN_SylLocked_Coord_sub1;
        end
    end
    
    
    
    % ######################################### RA
    bregionthis = 'RA';
    
    indsthis_spk = find(indsgrp==indsgrpU(i) & ...
        strcmp(DATSTRUCT_SPK.spike_bregions, bregionthis));
    
    if length(indsthis_spk)<2
        % need at least two spiking units to do this analysis.
        continue
    end
    
    bnum = unique(DATSTRUCT_SPK.bnum(indsthis_spk));
    enum = unique(DATSTRUCT_SPK.enum(indsthis_spk));
    sw = unique(DATSTRUCT_SPK.switch(indsthis_spk));
    mm = unique(DATSTRUCT_SPK.motifnum(indsthis_spk));
    
    indsthis_lfp = find(DATSTRUCT_LFP.bnum==bnum & DATSTRUCT_LFP.enum==enum ...
        & DATSTRUCT_LFP.switch==sw & DATSTRUCT_LFP.motifnum==mm & ...
        strcmp(DATSTRUCT_LFP.LFP_bregions, bregionthis));
    
    assert(length(indsthis_spk)==length(indsthis_lfp));
    
    % --------------- COLLECT DATA
    indstrials = DATSTRUCT_SPK.inds_base{indsthis_spk(1)};
    frmat = DATSTRUCT_SPK.frmat(indsthis_spk);
    spkdat = DATSTRUCT_SPK.spike_dat(indsthis_spk);
    lfpdat = DATSTRUCT_LFP.LFP_dat(indsthis_lfp);
    assert(size(frmat{1},2) == size(lfpdat{1},2));
    
    % ---------- GET TRIALS OF INTEREST
    frmat = cellfun(@(x)x(:, indstrials), frmat, 'UniformOutput', 0);
    spkdat = cellfun(@(x)x(indstrials), spkdat, 'UniformOutput', 0);
    lfpdat = cellfun(@(x)x(:, indstrials), lfpdat, 'UniformOutput', 0);
    
    % -- save for later
    lfpdat_RA = lfpdat;
    spkdat_RA = spkdat;
    
    % ============ 1) compute all pairwise correlations
    rhoall = lt_neural_POPLEARN_SylLocked_paircorr(frmat);
    % --- take average across channel pairs
    rhoall = mean(rhoall,2);
    
    
    % ============ 1b) collect all mean firing rates
    frmean = cellfun(@(x)mean(x,1), frmat, 'UniformOutput', 0);
    frmean = cell2mat(frmean);
    frmean = mean(frmean,1)';
    
    % ============ 2) LFP POWER
    % - if multiple channels, get power for each and then take average
    ntrials = size(lfpdat{1}, 2);
    nchans = length(lfpdat);
    lfppow_all = nan(ntrials, nchans);
    for j=1:length(lfpdat)
        lfppow = sqrt(mean(lfpdat{j}.^2, 1));
        %         lfppow = std(lfpdat{j}, [], 1);
        lfppow_all(:, j) = lfppow;
    end
    % -- take average across channels
    lfppow_all = mean(lfppow_all,2);
    
    
    % ================== SAVE OUTPUT
    DATSTRUCT_POP.bnum = [DATSTRUCT_POP.bnum; bnum];
    DATSTRUCT_POP.enum = [DATSTRUCT_POP.enum; enum];
    DATSTRUCT_POP.switch = [DATSTRUCT_POP.switch; sw];
    DATSTRUCT_POP.motifnum = [DATSTRUCT_POP.motifnum; mm];
    DATSTRUCT_POP.bregion = [DATSTRUCT_POP.bregion; bregionthis];
    DATSTRUCT_POP.fr_RhoPairwise = [DATSTRUCT_POP.fr_RhoPairwise; rhoall];
    DATSTRUCT_POP.fr_MeanFr = [DATSTRUCT_POP.fr_MeanFr; frmean];
    DATSTRUCT_POP.lfp_Pow = [DATSTRUCT_POP.lfp_Pow; lfppow_all];
    DATSTRUCT_POP.coher_cohgram_all = [DATSTRUCT_POP.coher_cohgram_all; {[]}];
    
    % plot random case, low, med, and high correlation cases?
    if plotRawRand ==1
        if rand<0.15
            lt_neural_POPLEARN_SylLocked_Coord_sub1;
        end
    end
    
    
    %% ======= get coherence between LMAN and RA LFP channels
    
    ntrials = size(lfpdat_LMAN{1},2);
    nchan_LMAN = length(lfpdat_LMAN);
    nchan_RA = length(lfpdat_RA);
    t_LFP = PARAMS.THIS.lfpx;
    
    % --- go thru all pairs of LMAN and RA channels, once for each trials
    cohgram_all = cell(1, nchan_LMAN*nchan_RA);
    cc=1;
    for nn=1:nchan_LMAN
        for nnn=1:nchan_RA
            
            lfpL = lfpdat_LMAN{nn};
            lfpR = lfpdat_RA{nnn};
            [~, t, f, ~,~,~,~, Ctrials, ~] = ...
                lt_neural_Coher_BatchCoher(lfpL, lfpR, 'mtaper_trials', ...
                t_LFP, ntapers, movingwin, tw, []);
            
            % ------------ SAVE COHEROGRAM
            cohgram_all{cc} = Ctrials;
            cc=cc+1;
        end
    end
    % ---- take mean across channel pairs
    cohgram_all = lt_neural_Coher_Cell2Mat(cohgram_all);
    cohgram_all = mean(cohgram_all, 3);
    assert(all(size(cohgram_all) == [length(f) ntrials]), ' probably multiple time bins?');
    
    PARAMS.THIS.coher_f = f;
    DATSTRUCT_POP.bnum = [DATSTRUCT_POP.bnum; bnum];
    DATSTRUCT_POP.enum = [DATSTRUCT_POP.enum; enum];
    DATSTRUCT_POP.switch = [DATSTRUCT_POP.switch; sw];
    DATSTRUCT_POP.motifnum = [DATSTRUCT_POP.motifnum; mm];
    DATSTRUCT_POP.bregion = [DATSTRUCT_POP.bregion; 'LMAN-RA'];
    
    DATSTRUCT_POP.fr_RhoPairwise = [DATSTRUCT_POP.fr_RhoPairwise; {[]}];
    DATSTRUCT_POP.fr_MeanFr = [DATSTRUCT_POP.fr_MeanFr; {[]}];
    DATSTRUCT_POP.lfp_Pow = [DATSTRUCT_POP.lfp_Pow; {[]}];
    
    DATSTRUCT_POP.coher_cohgram_all = [DATSTRUCT_POP.coher_cohgram_all; cohgram_all];
    
    xbins = f>=fwind(1) & f<=fwind(2);
    cohscal_tmp = mean(cohgram_all(xbins,:),1);

    if plotRawRand ==1
        if rand<0.1
            lt_neural_POPLEARN_SylLocked_Coord_sub2;
        end
    end
end


%% ========= [PROCESS] CONVERT COHERENCE TO SCALAR
xbins = PARAMS.THIS.coher_f>=fwind(1) & PARAMS.THIS.coher_f<=fwind(2);
cohscal_all = cell(length(DATSTRUCT_POP.bnum),1);
for i=1:length(DATSTRUCT_POP.bnum)
    cohgram = DATSTRUCT_POP.coher_cohgram_all{i};
    if isempty(cohgram)
        continue
    end
    
    cohscal_all{i} = mean(cohgram(xbins,:),1);
end

DATSTRUCT_POP.coher_cohscal = cohscal_all;

%% ======= [PROCESS] - FOR EACH CASE GET CORRELATION BETWEEN SPIKE CORR AND LFP POW
% corrtype = 'Pearson';
corrtype = 'Spearman';
corrall_frcorrVslfppow = [];
corrall_frcorrVslfppow_pval = [];
corrall_frmeanVslfppow = [];
corrall_frmeanVslfppow_pval = [];
regress_b_IntMeanCorr = {};
regress_bInt_IntMeanCorr = {};

for i=1:length(DATSTRUCT_POP.bnum)
    
    spkcorr = DATSTRUCT_POP.fr_RhoPairwise{i};
    lfppow = DATSTRUCT_POP.lfp_Pow{i};
    spkmean = DATSTRUCT_POP.fr_MeanFr{i};
    
    if isempty(spkcorr)
    corrall_frcorrVslfppow = [corrall_frcorrVslfppow; nan];
    corrall_frcorrVslfppow_pval = [corrall_frcorrVslfppow_pval; nan];
    corrall_frmeanVslfppow = [corrall_frmeanVslfppow; nan];
    corrall_frmeanVslfppow_pval = [corrall_frmeanVslfppow_pval; nan];
    regress_b_IntMeanCorr = [regress_b_IntMeanCorr; {[]}];
    regress_bInt_IntMeanCorr = [regress_bInt_IntMeanCorr; {[]}];
        continue
    end
    
    % -- 2) zscore spk mean and spk corr
    spkmean = (spkmean-mean(spkmean))./std(spkmean);
    spkcorr = (spkcorr-mean(spkcorr))./std(spkcorr);
    
    % ==== get correlation
    [rho, p] = corr(spkcorr, lfppow, 'type', corrtype);
    
    corrall_frcorrVslfppow = [corrall_frcorrVslfppow; rho];
    corrall_frcorrVslfppow_pval = [corrall_frcorrVslfppow_pval; p];
    
    
    % ==== get correlation
    [rho, p] = corr(spkmean, lfppow, 'type', corrtype);
    
    corrall_frmeanVslfppow = [corrall_frmeanVslfppow; rho];
    corrall_frmeanVslfppow_pval = [corrall_frmeanVslfppow_pval; p];
    
    
    % ===== multiple regression
    % -- 1) conver lfp power to rank order\
    if (1)
        [~, indsort] = sort(lfppow);
        lfppow_rank = nan(size(lfppow));
        lfppow_rank(indsort) = 1:length(lfppow);
        % --- convert to percentile
        lfppow_rank = lfppow_rank./length(lfppow_rank); 
    else
        lfppow_rank = lfppow;
    end
    
    % --- 2) do mulitp[le regression
    X = [ones(size(spkmean)) spkmean spkcorr];
    y = lfppow_rank;
    [b, bint] = regress(y, X);
    
    regress_b_IntMeanCorr = [regress_b_IntMeanCorr; b];
    regress_bInt_IntMeanCorr = [regress_bInt_IntMeanCorr; bint];
    
%     if b(3)<-0.048 & b(3)>-0.05
%         keyboard
%     end
        
end

DATSTRUCT_POP.Rho_frcorrVslfppow = corrall_frcorrVslfppow;
DATSTRUCT_POP.Rho_frcorrVslfppow_pval = corrall_frcorrVslfppow_pval;

DATSTRUCT_POP.Rho_frmeanVslfppow = corrall_frmeanVslfppow;
DATSTRUCT_POP.Rho_frmeanVslfppow_pval = corrall_frmeanVslfppow_pval;


DATSTRUCT_POP.regress_b_IntMeanCorr = regress_b_IntMeanCorr;
DATSTRUCT_POP.regress_bInt_IntMeanCorr = regress_bInt_IntMeanCorr;


%% ============== [PROCESS] correlations between spiking and coherence
[indsgrp, indsgrpU] = lt_tools_grp2idx({DATSTRUCT_POP.bnum, DATSTRUCT_POP.enum, ...
    DATSTRUCT_POP.switch, DATSTRUCT_POP.motifnum});

corrtype = 'Spearman';
Rho_frcorr_vs_coherence_LMAN = [];
Rho_frcorr_vs_coherence_RA = [];
Rho_frcorr_vs_coherence_LmanRaMean = [];

Rho_lfpow_vs_coherence_LMAN = [];
Rho_lfpow_vs_coherence_RA = [];
Rho_lfpow_vs_coherence_LmanRaMean = [];


for i=1:length(indsgrpU)
    
    % =============== COLLECT THINGS
    % ------ COHERENCE
    indsthis = find(indsgrp==indsgrpU(i) & strcmp(DATSTRUCT_POP.bregion, 'LMAN-RA'));
    assert(length(indsthis)==1);
    
    cohscalar = DATSTRUCT_POP.coher_cohscal{indsthis};
    
    % ------ SPIKE CORRELATIONS (LMAN)
    indsthis = find(indsgrp==indsgrpU(i) & strcmp(DATSTRUCT_POP.bregion, 'LMAN'));
    assert(length(indsthis)==1);
    
    frcorr_LMAN = DATSTRUCT_POP.fr_RhoPairwise{indsthis};
    
    % ------ SPIKE CORRELATIONS (RA)
    indsthis = find(indsgrp==indsgrpU(i) & strcmp(DATSTRUCT_POP.bregion, 'RA'));
    assert(length(indsthis)==1);
    
    frcorr_RA = DATSTRUCT_POP.fr_RhoPairwise{indsthis};
    
    
    % ========= take zscore
    frcorr_LMAN = (frcorr_LMAN-mean(frcorr_LMAN))./std(frcorr_LMAN);
    frcorr_RA= (frcorr_RA-mean(frcorr_RA))./std(frcorr_RA);
    frcorr_LMAN_RA_mean = mean([frcorr_LMAN frcorr_RA],2);
    
    % ================= GET CORRELATIONS
    % -- 1) LMAN vs. coherence
    [rho, p] = corr(frcorr_LMAN, cohscalar', 'type', corrtype);
    Rho_frcorr_vs_coherence_LMAN = [Rho_frcorr_vs_coherence_LMAN; rho];
    
    % -- 1) RA vs. coherence
    [rho, p] = corr(frcorr_RA, cohscalar', 'type', corrtype);
    Rho_frcorr_vs_coherence_RA = [Rho_frcorr_vs_coherence_RA; rho];
    
    % -- 1) mean(LMAN, RA) vs. coherence
    [rho, p] = corr(frcorr_LMAN_RA_mean, cohscalar', 'type', corrtype);
    Rho_frcorr_vs_coherence_LmanRaMean = [Rho_frcorr_vs_coherence_LmanRaMean; rho];
    
    
    
    % ################################ RELATIONSHIP BETWEEN LFP POWER AND
    % COHERENCE?
    
    % ------ LFP POWER (LMAN)
    indsthis = find(indsgrp==indsgrpU(i) & strcmp(DATSTRUCT_POP.bregion, 'LMAN'));
    lfppow_LMAN = DATSTRUCT_POP.lfp_Pow{indsthis};
    
    % ------ LFP POWER (RA)
    indsthis = find(indsgrp==indsgrpU(i) & strcmp(DATSTRUCT_POP.bregion, 'RA'));
    lfppow_RA = DATSTRUCT_POP.lfp_Pow{indsthis};
    
    % ========= take zscore
    lfppow_LMAN = (lfppow_LMAN-mean(lfppow_LMAN))./std(lfppow_LMAN);
    lfppow_RA= (lfppow_RA-mean(lfppow_RA))./std(lfppow_RA);
    lfppow_LmanRa = mean([lfppow_LMAN lfppow_RA],2);
    
    
    % ------ DO CORERLATIONS
    % -- 1) LMAN vs. coherence
    [rho, p] = corr(lfppow_LMAN, cohscalar', 'type', corrtype);
    Rho_lfpow_vs_coherence_LMAN = [Rho_lfpow_vs_coherence_LMAN; rho];
    
    % -- 1) RA vs. coherence
    [rho, p] = corr(lfppow_RA, cohscalar', 'type', corrtype);
    Rho_lfpow_vs_coherence_RA = [Rho_lfpow_vs_coherence_RA; rho];
    
    % -- 1) LMAN/RA vs. coherence
    [rho, p] = corr(lfppow_LmanRa, cohscalar', 'type', corrtype);
    Rho_lfpow_vs_coherence_LmanRaMean = [Rho_lfpow_vs_coherence_LmanRaMean; rho];
end


%% ============== [PLOT] - correlation between fr correaltion and coherence?
lt_figure; hold on;

lt_subplot(3,2,1); hold on;
xlabel('LMAN');
ylabel('RA');
title('frcorr vs. coherence(lfp)');
x = Rho_frcorr_vs_coherence_LMAN;
y = Rho_frcorr_vs_coherence_RA;
plot(x,y, 'ok');
lt_plot_makesquare_plot45line(gca, 'r');
p = signrank(x, y);
lt_plot_pvalue(p, 'srank');

lt_subplot(3,2,2); hold on;
xlabel('LMAN - RA - mean(LMAN,RA)');
ylabel('rho (frcorr vs. coherence)');
x = [1 2 3];
% Y = {Rho_frcorr_vs_coherence_LMAN, Rho_frcorr_vs_coherence_RA, Rho_frcorr_vs_coherence_LmanRaMean};
Y = [Rho_frcorr_vs_coherence_LMAN, Rho_frcorr_vs_coherence_RA, Rho_frcorr_vs_coherence_LmanRaMean];
plot(x, Y', '-', 'Color', [0.7 0.7 0.7]);
lt_plot(x+0.2, mean(Y,1), {'Errors', lt_sem(Y), 'Color', 'r'});
xlim([0 4]);
lt_plot_zeroline;
for j=1:size(Y,2)
    p = signrank(Y(:,j));
    lt_plot_text(j, max(Y(:,j)), ['p=' num2str(p)], 'm');
end



lt_subplot(3,2,3); hold on;
xlabel('LMAN');
ylabel('RA');
title('lfp pow vs. coherence(lfp)');
x = Rho_lfpow_vs_coherence_LMAN;
y = Rho_lfpow_vs_coherence_RA;
plot(x,y, 'ok');
lt_plot_makesquare_plot45line(gca, 'r');
p = signrank(x, y);
lt_plot_pvalue(p, 'srank');

lt_subplot(3,2,4); hold on;
xlabel('LMAN - RA - mean(LMAN,RA)');
ylabel('rho (lfp pow vs. coherence)');
x = [1 2 3];
% Y = {Rho_frcorr_vs_coherence_LMAN, Rho_frcorr_vs_coherence_RA, Rho_frcorr_vs_coherence_LmanRaMean};
Y = [Rho_lfpow_vs_coherence_LMAN, Rho_lfpow_vs_coherence_RA, Rho_lfpow_vs_coherence_LmanRaMean];
plot(x, Y', '-', 'Color', [0.7 0.7 0.7]);
lt_plot(x+0.2, mean(Y,1), {'Errors', lt_sem(Y), 'Color', 'r'});
xlim([0 4]);
lt_plot_zeroline;
for j=1:size(Y,2)
    p = signrank(Y(:,j));
    lt_plot_text(j, max(Y(:,j)), ['p=' num2str(p)], 'm');
end


%% ============== [PLOT] - correlation between spk corr and LFP power?
lt_figure;

% ===
lt_subplot(3,2,1); hold on;
xlabel('pval(log10)');
ylabel('rho (frcorr vs. lfp pow)');
title('LMAN');
indsthis = strcmp(DATSTRUCT_POP.bregion, 'LMAN');
x = log10(DATSTRUCT_POP.Rho_frcorrVslfppow_pval(indsthis));
y = DATSTRUCT_POP.Rho_frcorrVslfppow(indsthis);
plot(x,y,'ok');
lt_plot_zeroline;
line(log10([0.05 0.05]), ylim);


% ===
lt_subplot(3,2,2); hold on;
xlabel('rho (frcorr vs. lfp pow)');
indsthis = strcmp(DATSTRUCT_POP.bregion, 'LMAN');
x = DATSTRUCT_POP.Rho_frcorrVslfppow(indsthis);
lt_plot_histogram(x);
lt_plot_zeroline_vert;
p = signrank(x);
lt_plot_pvalue(p, 'vs0');


% ===
lt_subplot(3,2,3); hold on;
xlabel('pval(log10)');
ylabel('rho (frMean vs. lfp pow)');
indsthis = strcmp(DATSTRUCT_POP.bregion, 'LMAN');
x = log10(DATSTRUCT_POP.Rho_frmeanVslfppow_pval(indsthis));
y = DATSTRUCT_POP.Rho_frmeanVslfppow(indsthis);
plot(x,y,'ok');
lt_plot_zeroline;
line(log10([0.05 0.05]), ylim);


% ===
lt_subplot(3,2,4); hold on;
xlabel('rho (frMean vs. lfp pow)');
indsthis = strcmp(DATSTRUCT_POP.bregion, 'LMAN');
x = DATSTRUCT_POP.Rho_frmeanVslfppow(indsthis);
lt_plot_histogram(x);
lt_plot_zeroline_vert;
p = signrank(x);
lt_plot_pvalue(p, 'vs0');

% ====
lt_subplot(3,2,5); hold on;
title('multipel regression');
xlabel('effects [meanFR, FRcorr]');
ylabel('b');
indsthis = strcmp(DATSTRUCT_POP.bregion, 'LMAN');
b_meanfr = cellfun(@(x)x(2), DATSTRUCT_POP.regress_b_IntMeanCorr(indsthis));
b_corrfr = cellfun(@(x)x(3), DATSTRUCT_POP.regress_b_IntMeanCorr(indsthis));
Y = {b_meanfr, b_corrfr};
x = [1 2];
lt_plot_MultDist(Y, x);
for j=1:length(Y)
    p = signrank(Y{j});
    lt_plot_text(j, max(Y{j}), ['p=' num2str(p)], 'r');
end
xlim([0 3]);
lt_plot_zeroline;


lt_figure;

% ===
lt_subplot(3,2,1); hold on;
xlabel('pval(log10)');
ylabel('rho (frcorr vs. lfp pow)');
title('RA');
indsthis = strcmp(DATSTRUCT_POP.bregion, 'RA');
x = log10(DATSTRUCT_POP.Rho_frcorrVslfppow_pval(indsthis));
y = DATSTRUCT_POP.Rho_frcorrVslfppow(indsthis);
plot(x,y,'ok');
lt_plot_zeroline;
line(log10([0.05 0.05]), ylim);


% ===
lt_subplot(3,2,2); hold on;
xlabel('rho (frcorr vs. lfp pow)');
indsthis = strcmp(DATSTRUCT_POP.bregion, 'RA');
x = DATSTRUCT_POP.Rho_frcorrVslfppow(indsthis);
lt_plot_histogram(x);
lt_plot_zeroline_vert;
p = signrank(x);
lt_plot_pvalue(p, 'vs0');


% ===
lt_subplot(3,2,3); hold on;
xlabel('pval(log10)');
ylabel('rho (frMean vs. lfp pow)');
indsthis = strcmp(DATSTRUCT_POP.bregion, 'RA');
x = log10(DATSTRUCT_POP.Rho_frmeanVslfppow_pval(indsthis));
y = DATSTRUCT_POP.Rho_frmeanVslfppow(indsthis);
plot(x,y,'ok');
lt_plot_zeroline;
line(log10([0.05 0.05]), ylim);


% ===
lt_subplot(3,2,4); hold on;
xlabel('rho (frMean vs. lfp pow)');
indsthis = strcmp(DATSTRUCT_POP.bregion, 'RA');
x = DATSTRUCT_POP.Rho_frmeanVslfppow(indsthis);
lt_plot_histogram(x);
lt_plot_zeroline_vert;
p = signrank(x);
lt_plot_pvalue(p, 'vs0');

% ====
lt_subplot(3,2,5); hold on;
title('multipel regression');
xlabel('effects [meanFR, FRcorr]');
ylabel('b');
indsthis = strcmp(DATSTRUCT_POP.bregion, 'RA');
b_meanfr = cellfun(@(x)x(2), DATSTRUCT_POP.regress_b_IntMeanCorr(indsthis));
b_corrfr = cellfun(@(x)x(3), DATSTRUCT_POP.regress_b_IntMeanCorr(indsthis));
Y = {b_meanfr, b_corrfr};
x = [1 2];
lt_plot_MultDist(Y, x);
for j=1:length(Y)
    p = signrank(Y{j});
    lt_plot_text(j, max(Y{j}), ['p=' num2str(p)], 'r');
end
xlim([0 3]);
lt_plot_zeroline;


%% ====== sanity check
if (0)
   

    % ==== print bnum, enum, sw, and motif, for a given index
    indtmp = 16;
    DATSTRUCT_POP.bnum(indtmp)
    DATSTRUCT_POP.enum(indtmp)
    DATSTRUCT_POP.switch(indtmp)
    DATSTRUCT_POP.motifnum(indtmp)
    DATSTRUCT_POP.bregion(indtmp)
end

%% ======= compare LMAN vs. RA
lt_figure; hold on;

% --- first make sure that each case has LMAN and RA
tmp = strcmp(DATSTRUCT_POP.bregion, 'LMAN');
tmp2 = strcmp(DATSTRUCT_POP.bregion, 'RA');
assert(all(tmp(1:3:end)==1));
assert(all(tmp2(2:3:end)==1));
assert(mod(length(tmp), 3)==0);
% note; this just checks that lman and ra alternate, almost certainly means
% they are matched././/


% =====
lt_subplot(3,2,1); hold on;
xlabel('LMAN [b(corr) - b(mean)]');
ylabel('RA [b(corr) - b(mean)]');
title('mult regression b');

indsthis = strcmp(DATSTRUCT_POP.bregion, 'LMAN');
y1 = DATSTRUCT_POP.regress_b_IntMeanCorr(indsthis);
y1 = cellfun(@(x)(x(3)-x(2)), y1); % - get diff of b values

indsthis = strcmp(DATSTRUCT_POP.bregion, 'RA');
y2 = DATSTRUCT_POP.regress_b_IntMeanCorr(indsthis);
y2 = cellfun(@(x)(x(3)-x(2)), y2); % - get diff of b values

plot(y1, y2, 'ok');
lt_plot_zeroline;
lt_plot_zeroline_vert;

y1mean = mean(y1);
y1sem = lt_sem(y1);
y2mean = mean(y2);
y2sem = lt_sem(y2);
lt_plot(y1mean, y2mean, {'Errors', y2sem, 'Xerrors', y1sem});


% =====
lt_subplot(3,2,2); hold on;
xlabel('LMAN - RA');
title('[b(corr) - b(mean)]');
ylabel('multiple regression');
indsthis = strcmp(DATSTRUCT_POP.bregion, 'LMAN');
y1 = DATSTRUCT_POP.regress_b_IntMeanCorr(indsthis);
y1 = cellfun(@(x)(x(3)-x(2)), y1); % - get diff of b values

indsthis = strcmp(DATSTRUCT_POP.bregion, 'RA');
y2 = DATSTRUCT_POP.regress_b_IntMeanCorr(indsthis);
y2 = cellfun(@(x)(x(3)-x(2)), y2); % - get diff of b values

x = [1 2];
Y = [y1 y2];
plot(x, Y, '-', 'Color', [0.7 0.7 0.7]);

lt_plot(x+0.15, mean(Y), {'Errors', lt_sem(Y), 'Color', 'r'});
xlim([0 3]);
lt_plot_zeroline;




% =====
lt_subplot(3,2,3); hold on;
xlabel('LMAN - RA');
title('[b(corr) - b(mean)]');
ylabel('separate correlations');
indsthis = strcmp(DATSTRUCT_POP.bregion, 'LMAN');
y1 = DATSTRUCT_POP.Rho_frcorrVslfppow(indsthis) - DATSTRUCT_POP.Rho_frmeanVslfppow(indsthis);

indsthis = strcmp(DATSTRUCT_POP.bregion, 'RA');
y2 = DATSTRUCT_POP.Rho_frcorrVslfppow(indsthis) - DATSTRUCT_POP.Rho_frmeanVslfppow(indsthis);

x = [1 2];
Y = [y1 y2];
plot(x, Y, '-', 'Color', [0.7 0.7 0.7]);

lt_plot(x+0.15, mean(Y), {'Errors', lt_sem(Y), 'Color', 'r'});
xlim([0 3]);
lt_plot_zeroline;



% =====
lt_subplot(3,2,4); hold on;
xlabel('LMAN - RA');
title('multiple regression');
ylabel('b(corr)');

indsthis = strcmp(DATSTRUCT_POP.bregion, 'LMAN');
y1 = DATSTRUCT_POP.regress_b_IntMeanCorr(indsthis);
y1 = cellfun(@(x)(x(3)), y1); % - get diff of b values

indsthis = strcmp(DATSTRUCT_POP.bregion, 'RA');
y2 = DATSTRUCT_POP.regress_b_IntMeanCorr(indsthis);
y2 = cellfun(@(x)(x(3)), y2); % - get diff of b values

x = [1 2];
Y = [y1 y2];
plot(x, Y, '-', 'Color', [0.7 0.7 0.7]);

lt_plot(x+0.15, mean(Y), {'Errors', lt_sem(Y), 'Color', 'r'});
xlim([0 3]);
lt_plot_zeroline;

p = signrank(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'srank(vs)');
for i=1:size(Y,2)
    p = signrank(Y(:,i));
    lt_plot_text(i, max(Y(:,i)), ['p=' num2str(p)],'m');
end


% =====
lt_subplot(3,2,5); hold on;
xlabel('LMAN - RA');
title('multiple regression');
ylabel('b(frmean)');

indsthis = strcmp(DATSTRUCT_POP.bregion, 'LMAN');
y1 = DATSTRUCT_POP.regress_b_IntMeanCorr(indsthis);
y1 = cellfun(@(x)(x(2)), y1); % - get diff of b values

indsthis = strcmp(DATSTRUCT_POP.bregion, 'RA');
y2 = DATSTRUCT_POP.regress_b_IntMeanCorr(indsthis);
y2 = cellfun(@(x)(x(2)), y2); % - get diff of b values

x = [1 2];
Y = [y1 y2];
plot(x, Y, '-', 'Color', [0.7 0.7 0.7]);

lt_plot(x+0.15, mean(Y), {'Errors', lt_sem(Y), 'Color', 'r'});
xlim([0 3]);
lt_plot_zeroline;

p = signrank(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'srank(vs)');
for i=1:size(Y,2)
    p = signrank(Y(:,i));
    lt_plot_text(i, max(Y(:,i)), ['p=' num2str(p)],'m');
end


