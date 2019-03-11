function lt_neural_POPLEARN_SylLocked_Spectra(DATSTRUCT_SPK, DATSTRUCT_LFP, PARAMS, ...
    xtoplot, SwitchCohStruct, SwitchStruct)
%% lt 3/4/19 - trial by trial power spectra for spiking and LFP.
% Also spike-LFP coherence

ntapers = 1; % for spike field coh will always use 5 tapers.
powscale = 'db_prop'; % first take db, then subtract min, then take proportion.
% powscale = 'db'; % just db

% window size in sec
T = (PARAMS.THIS.frmat_x(end)-PARAMS.THIS.frmat_x(1))+(PARAMS.THIS.frmat_x(2)-PARAMS.THIS.frmat_x(1)); 

%% 
lt_switch_chronux(1);

params = struct;

params.fpass = [1/T 150];

% tw = 3;
% w = tw/movingwin(1); % in hz, for desired frequency resolution of tapers. % note, t is set to movingwin(1)
% % tw = movingwin(1)*w;
% if isempty(ntapers)
%     k = 2*tw-1;
% else
%     k = ntapers;
% end

TW = (ntapers+1)/2;

params.tapers = [TW ntapers];
params.Fs = round(1/(PARAMS.THIS.lfpx(2)-PARAMS.THIS.lfpx(1))); % hard coded fs for LFP;
params.trialave = 0;

disp(params);
disp('ARE THESE PARAMS OK?')
% pause;

% ------- params for spiking
params_spk = params;
params_spk.Fs = 1000;
params_spk.Fs = 1500;
fscorr = 1;
mintime = xtoplot(1);
maxtime = xtoplot(2);
dt=1/params_spk.Fs; % sampling time
t_spk=mintime-dt:dt:maxtime+dt; % time grid for prolates

% --- params for coh
k = 5;
params_coh = struct;
params_coh.fpass = [1/T 150];
TW = (k+1)/2;

params_coh.tapers = [TW k];
params_coh.Fs = round(1/(PARAMS.THIS.lfpx(2)-PARAMS.THIS.lfpx(1))); % hard coded fs for LFP;
params_coh.trialave = 0;

%% ========= 1) EXTRACTION [LFP] SPECTRA

SpecAll = {};

for i=1:length(DATSTRUCT_LFP.bnum)
    disp(i);
    
%     chanthis = DATSTRUCT_LFP.LFP_chan(i);
    lfpdat = DATSTRUCT_LFP.LFP_dat{i};
    
    [S, f] =  mtspectrumc(lfpdat, params);
    
    % --- SCALE POWER...
    if strcmp(powscale, 'db')
       S = 10*log10(S); 
    elseif strcmp(powscale, 'db_prop')
        S = 10*log10(S);
        Smin = min(S, [], 1);
        S = S-repmat(Smin, size(S,1), 1);
        S = S./sum(S, 1);
%     elseif strcmp(powscale', 'db_prop_max1')
%         S = 10*log10(S);
%         Smin = min(S, [], 1);
%         S = S-repmat(Smin, size(S,1), 1);
%         tmp = max(S,[],1);
%         S = S./tmp;
    end
    

    SpecAll = [SpecAll; S];
   
end

% ======= SAVE OUTPUT
DATSTRUCT_LFP.SpecAll = SpecAll;
PARAMS.THIS.spec_lfp_f = f;

%% ========= 1) EXTRACTION [SPIKE] SPECTRA

% SpecAll_LFP= {};
SpecAll_SPK = {};
CohAll = {};

for i=1:length(DATSTRUCT_SPK.bnum)
    disp(i);
    
    % ####################### SPIKE DATA
    spkdat = DATSTRUCT_SPK.spike_dat{i};
    
    % --- convert to structure arrauy
    spkdat2 = struct;
    for j=1:length(spkdat)
        spkdat2(j).spktimes = double(spkdat{j})';
    end
    
%     [Sspk, fspk, Rspk] = mtspectrumpt(spkdat2(1).spktimes, params_spk, fscorr, t_spk);
    [S, fspk, Rspk] = mtspectrumpt(spkdat2, params_spk, fscorr, t_spk);
    
    
    % ============== RESCALE
    if strcmp(powscale, 'db')
       S = 10*log10(S); 
    elseif strcmp(powscale, 'db_prop')
        S = 10*log10(S);
        Smin = min(S, [], 1);
        S = S-repmat(Smin, size(S,1), 1);
        S = S./sum(S, 1);
    end
    
    % ============ SAVE
    SpecAll_SPK= [SpecAll_SPK; S];
    
    
    
    
    % ########################### GET LFP DATA
    inds_lfp = DATSTRUCT_SPK.LFP_indthis(i);
    lfpdat = DATSTRUCT_LFP.LFP_dat{inds_lfp};
    
    assert(length(spkdat)==size(lfpdat,2));
    
    % ========== GET COHERENCE
    t = PARAMS.THIS.lfpx;
    [C,phi,S12,S1,S2,fcoh]=coherencycpt(lfpdat, spkdat2, params_coh,fscorr,t);   
    
    CohAll = [CohAll; C];
    
end

% ======= SAVE OUTPUT
DATSTRUCT_SPK.SpecAll = SpecAll_SPK;
PARAMS.THIS.spec_spk_f = fspk;

DATSTRUCT_SPK.CohAll = CohAll;
PARAMS.THIS.coh_f = fcoh;




%% ========= PLOT SPECTRA FOR ALL SYLS
% --- 1 plot for each motif x brain region

[indsgrp, indsgrpU] = lt_tools_grp2idx({DATSTRUCT_SPK.bnum, DATSTRUCT_SPK.enum, DATSTRUCT_SPK.switch, ...
    DATSTRUCT_SPK.motifnum, DATSTRUCT_SPK.spike_bregions});

figcount=1;
subplotrows=6;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


for i=1:length(indsgrpU)
   
    inds_spk = find(indsgrp==indsgrpU(i));
    spec_spk = DATSTRUCT_SPK.SpecAll(inds_spk);
    
    inds_lfp = DATSTRUCT_SPK.LFP_indthis(inds_spk);
   
    spec_lfp = DATSTRUCT_LFP.SpecAll(inds_lfp);
    
    coh = DATSTRUCT_SPK.CohAll(inds_spk);
    
    % ============= metadat
    bnum = unique(DATSTRUCT_SPK.bnum(inds_spk));
    enum = unique(DATSTRUCT_SPK.enum(inds_spk));
    sw = unique(DATSTRUCT_SPK.switch(inds_spk));
    bregion = unique(DATSTRUCT_SPK.spike_bregions(inds_spk)); bregion = bregion{1};
    mm = unique(DATSTRUCT_SPK.motifnum(inds_spk));

    % ------
    assert(all(DATSTRUCT_LFP.bnum(inds_lfp)==bnum));
    assert(all(DATSTRUCT_LFP.enum(inds_lfp)==enum));
    assert(all(DATSTRUCT_LFP.switch(inds_lfp)==sw));
    assert(all(DATSTRUCT_LFP.motifnum(inds_lfp)==mm));
    assert(all(strcmp(DATSTRUCT_LFP.LFP_bregions(inds_lfp), bregion)));
    
    % --------
    bname = SwitchStruct.bird(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    motifname = SwitchCohStruct.bird(bnum).exptnum(enum).switchlist(sw).motifnum(mm).motifname;
    
    inds_base = DATSTRUCT_SPK.inds_base{inds_spk(1)};
    
    if strcmp(bregion, 'LMAN')
        pcol = [0.2 0.7 0.2];
    elseif strcmp(bregion, 'RA')
        pcol = 'r';
    end
    
    % ################ SPIKES
    % ============== Initiate figure
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['[SPIKES]'])
    ylabel([bname '-' ename '-sw' num2str(sw) ',' motifname]);
    
    % ===== GET MEANS OVER TRIALS
    Y_spk = cellfun(@(x)mean(x(:, inds_base),2), spec_spk, 'UniformOutput', 0);
    Y_spk = squeeze(lt_neural_Coher_Cell2Mat(Y_spk));
    
    % ---- PLOT
    f = PARAMS.THIS.spec_spk_f;
    plot(f, Y_spk, 'Color', pcol);
    plot(f, mean(Y_spk,2), 'Color', pcol, 'LineWidth', 2)
    axis tight;
    
    % ################ LFP
    % ============== Initiate figure
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['[LFP]']);
        
    % ===== GET MEANS OVER TRIALS
    Y_lfp = cellfun(@(x)mean(x(:, inds_base),2), spec_lfp, 'UniformOutput', 0);
    Y_lfp = squeeze(lt_neural_Coher_Cell2Mat(Y_lfp));
    
    % ---- PLOT
    f = PARAMS.THIS.spec_lfp_f;
    plot(f, Y_lfp, 'Color', pcol);
    plot(f, mean(Y_lfp,2), 'Color', pcol, 'LineWidth', 2)
    axis tight;
    
    
    % ################ SPIKE-LFP COHERENCE
    % ============== Initiate figure
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['[SPIKE-LFP]']);
        
    % ===== GET MEANS OVER TRIALS
    Y_coh = cellfun(@(x)mean(x(:, inds_base),2), coh, 'UniformOutput', 0);
    Y_coh = squeeze(lt_neural_Coher_Cell2Mat(Y_coh));
    
    % ---- PLOT
    f = PARAMS.THIS.coh_f;
    plot(f, Y_coh, 'Color', pcol);
    plot(f, mean(Y_coh,2), 'Color', pcol, 'LineWidth', 2)
    axis tight;
    
    
    
   % ################## OVELRYA SPIKES AND LFP 
    % ============== Initiate figure
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['[SPIKES (k), LFP(b), COH(r)']);
    
    % ---- PLOT
    % - LFP
    % -- normalize so max is 1
    Y_lfp = Y_lfp./max(Y_lfp,[], 1);
    y = mean(Y_lfp,2);
    if size(Y_lfp,2)>1
    ysem = lt_sem(Y_lfp');
        shadedErrorBar(PARAMS.THIS.spec_lfp_f, y, ysem, {'Color', 'b'}, 1);
    else
        plot(PARAMS.THIS.spec_lfp_f, y, '-b');
    end
    
    % - SPK
    Y_spk = Y_spk./max(Y_spk,[], 1);
    y = mean(Y_spk,2);
if size(Y_spk,2)>1
    ysem = lt_sem(Y_spk');
    shadedErrorBar(PARAMS.THIS.spec_spk_f, y, ysem, {'Color', 'k'}, 1);
    else
        plot(PARAMS.THIS.spec_spk_f, y, '-k');
end

    % - coh
    % ---- normalize everything to 
    Y_coh = Y_coh./max(Y_coh,[], 1);
    y = mean(Y_coh,2);
    if size(Y_coh,2)>1
        ysem = lt_sem(Y_coh');
        shadedErrorBar(PARAMS.THIS.coh_f, y, ysem, {'Color', 'r'}, 1);
    else
        plot(PARAMS.THIS.coh_f, y, '-r');
    end

axis tight;
end





%% ================ [PLOT] PLTO MEAN, AND MEAN OVER BIRDS
% -- each dat = unit x motif
figcount=1;
subplotrows=5;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];



% ======== all birds
indstoplot = find(strcmp(DATSTRUCT_SPK.spike_bregions, 'LMAN'));
pcol = [0.2 0.6 0.2];
ptit = 'LMAN';

Yall  = [];
Yall_lfp = [];
Yall_coh = [];
for i=1:length(indstoplot)
    indthis = indstoplot(i);
    
    indsbase = DATSTRUCT_SPK.inds_base{indthis};

    % ====== SPK
    y = DATSTRUCT_SPK.SpecAll{indthis};
    y = mean(y(:, indsbase), 2);

    Yall = [Yall; y'];
    
    % ===== LFP
    indlfp = DATSTRUCT_SPK.LFP_indthis(indthis);
    y = DATSTRUCT_LFP.SpecAll{indlfp};
    y = mean(y(:, indsbase), 2);

    Yall_lfp = [Yall_lfp; y'];
    
    % ===== cohernece
    y = DATSTRUCT_SPK.CohAll{indlfp};
    y = mean(y(:, indsbase), 2);

    Yall_coh = [Yall_coh; y'];
end

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title([ptit '(LFP)']);
x = PARAMS.THIS.spec_lfp_f;
plot(x, Yall_lfp', '-');

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('(LFP)');
x = PARAMS.THIS.spec_lfp_f;
ymean = nanmean(Yall_lfp, 1);
ysem = lt_sem(Yall_lfp);
shadedErrorBar(x, ymean, ysem, {'Color', pcol}, 1);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('(SPK)');
x = PARAMS.THIS.spec_spk_f;
plot(x, Yall', '-k');


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('(SPK)');
x = PARAMS.THIS.spec_spk_f;
ymean = nanmean(Yall, 1);
ysem = lt_sem(Yall);
shadedErrorBar(x, ymean, ysem, {'Color', pcol}, 1);


% == cojhernece
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('(COH)');
x = PARAMS.THIS.coh_f;
plot(x, Yall_coh', '-k');

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('(COH)');
x = PARAMS.THIS.coh_f;
ymean = nanmean(Yall_coh, 1);
ysem = lt_sem(Yall_coh);
shadedErrorBar(x, ymean, ysem, {'Color', pcol}, 1);




% ======== all birds
indstoplot = find(strcmp(DATSTRUCT_SPK.spike_bregions, 'RA'));
pcol = 'r';
ptit = 'RA';

Yall  = [];
Yall_lfp = [];
Yall_coh = [];
for i=1:length(indstoplot)
    indthis = indstoplot(i);
    
    indsbase = DATSTRUCT_SPK.inds_base{indthis};

    % ====== SPK
    y = DATSTRUCT_SPK.SpecAll{indthis};
    y = mean(y(:, indsbase), 2);

    Yall = [Yall; y'];
    
    % ===== LFP
    indlfp = DATSTRUCT_SPK.LFP_indthis(indthis);
    y = DATSTRUCT_LFP.SpecAll{indlfp};
    y = mean(y(:, indsbase), 2);

    Yall_lfp = [Yall_lfp; y'];
    
    % ===== cohernece
    y = DATSTRUCT_SPK.CohAll{indlfp};
    y = mean(y(:, indsbase), 2);

    Yall_coh = [Yall_coh; y'];
end

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title([ptit '(LFP)']);
x = PARAMS.THIS.spec_lfp_f;
plot(x, Yall_lfp', '-');

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('(LFP)');
x = PARAMS.THIS.spec_lfp_f;
ymean = nanmean(Yall_lfp, 1);
ysem = lt_sem(Yall_lfp);
shadedErrorBar(x, ymean, ysem, {'Color', pcol}, 1);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('(SPK)');
x = PARAMS.THIS.spec_spk_f;
plot(x, Yall', '-k');


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('(SPK)');
x = PARAMS.THIS.spec_spk_f;
ymean = nanmean(Yall, 1);
ysem = lt_sem(Yall);
shadedErrorBar(x, ymean, ysem, {'Color', pcol}, 1);


% == cojhernece
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('(COH)');
x = PARAMS.THIS.coh_f;
plot(x, Yall_coh', '-k');

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('(COH)');
x = PARAMS.THIS.coh_f;
ymean = nanmean(Yall_coh, 1);
ysem = lt_sem(Yall_coh);
shadedErrorBar(x, ymean, ysem, {'Color', pcol}, 1);
