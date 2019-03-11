function lt_neural_LFP_PlotEgRaw_v2(DatAll, t_onoff, fs, bregionlist, chanlist_toget, ...
    bnum, enum, swplot, motifnum, filt_low, filt_hi, SwitchCohStruct, ...
    SwitchStruct, OUTSTRUCT, MOTIFSTATS_pop, PARAMS, ...
    SummaryStruct, OUTSTRUCT_CohMatOnly, plotUnfiltered, trialOrder, noiseband, ...
    cohwind_f)
%%  lt 3/2/19 - plots all trials, heat map of neural activity

lt_switch_chronux(1);

%% COLLECT DATA


i=bnum;
ii=enum;
mm = motifnum;
filt_fs = fs;

bname = SwitchStruct.bird(i).birdname;
ename = SwitchStruct.bird(i).exptnum(ii).exptname;

% ==== FROM OUTSTRUCT
indsthis = find(OUTSTRUCT.bnum==i & OUTSTRUCT.enum==ii & OUTSTRUCT.switch==swplot & ...
    OUTSTRUCT.motifnum==mm);
triallist = [find(OUTSTRUCT.indsbase{indsthis(1)}) find(OUTSTRUCT.indsWN{indsthis(1)})];
% triallist_orig = triallist;

trialWNon = sum(OUTSTRUCT.indsbase{indsthis(1)})+1;
chanpairs_OUT = OUTSTRUCT.chanpair(indsthis,:);
cohscal_OUT = OUTSTRUCT.cohscal(indsthis);

if strcmp(trialOrder, 'coh')
    cohscal_OUT_trials = cellfun(@(x)x(triallist), OUTSTRUCT.cohscal(indsthis), 'UniformOutput', 0);
    cohscalTMP = mean(cell2mat(cohscal_OUT_trials),1); % take average over all chan pairs
    
    [~, indsort] = sort(cohscalTMP);
    triallist = indsort;
end

%%
datcohstruct = SwitchCohStruct.bird(i).exptnum(ii).switchlist(swplot).motifnum(mm);

% -----
% indsbase = datcohstruct.indsbase_epoch;
% indsWN = datcohstruct.indsWN_epoch;
assert(size(datcohstruct.lfpall,2) == size(DatAll,1)) % n chans
assert(size(datcohstruct.lfpall,1) == size(DatAll{1},2)) % n trials

LFPdat = datcohstruct.lfpall;
t_lfp = datcohstruct.t_lfp;


%% ====== COLLECT SPIKING DATA

neurset = SwitchCohStruct.bird(bnum).exptnum(enum).switchlist(swplot).motifnum(motifnum).neursetused;
segAll = MOTIFSTATS_pop.birds(bnum).exptnum(enum).DAT.setnum(neurset).motif(motifnum).SegExtr_neurfakeID;

assert(length(segAll(1).SegmentsExtract)==size(DatAll{1},2));

pcols_chans = lt_make_plot_colors(size(DatAll,1), 0, 0);

%% ==== collect LFP data
lfpdat = SwitchCohStruct.bird(bnum).exptnum(enum).switchlist(swplot).motifnum(motifnum).lfpall;
lfpchans = SwitchCohStruct.bird(bnum).exptnum(enum).switchlist(swplot).motifnum(motifnum).lfpall_chans;
lfp_t = SwitchCohStruct.bird(bnum).exptnum(enum).switchlist(swplot).motifnum(motifnum).t_lfp;

datthis  = SwitchCohStruct.bird(bnum).exptnum(enum).switchlist(swplot).motifnum(motifnum);

filename = [datthis.fileprefix '/phi' datthis.filesuffix];
pairstoget = datthis.chanpairstokeep;
phimat = lt_neural_LFP_loadProcessDat(filename, pairstoget);
phi_chanpair = datthis.chanpair;

%% WHAT trial to plot?

% for i = 1:length(trialtoplot)
%
%     j = trialtoplot(i);
%
%     % ============================ PLOT FOR THIS TRIAL.
%
%     ons = segAll(1).SegmentsExtract(j).motifsylOnsets - PARAMS.motif_predur;
%     offs = segAll(1).SegmentsExtract(j).motifsylOffsets - PARAMS.motif_predur;
%
%     % ============================ PLOT ALL RAW NEURAL
%     datfiltall = [];
%     pcols_all = {};

figcount=1;
subplotrows=5;
if ~isempty(noiseband)
    subplotcols=4;
else
    subplotcols=3;
end
fignums_alreadyused=[];
hfigs=[];
hsplots = [];



for cc = 1:length(chanlist_toget)
    
    chanthis = chanlist_toget(cc);
    bregionthis = bregionlist{cc};
    
    % =============== RAW DATA
    datthis = DatAll{cc}(:, triallist);
    x = linspace(t_onoff{cc}(1), t_onoff{cc}(2), size(datthis,1));
    f = 1:length(triallist); % trial num
    
    % ---- filter (jhigh pass);
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(['(spikeband)']);
    ylabel(['chan' num2str(chanthis) ',[' bregionthis ']']);
    datthis_high = lt_neural_filter(double(datthis), fs);
    clim = prctile(datthis_high(:), [2.5 97.5]);
    lt_neural_Coher_Plot(datthis_high, x, f, 1, '', clim);
    line([-0.1 -0.1], ylim, 'Color', 'w', 'LineStyle', '-');
    line([-0 -0], ylim, 'Color', 'w', 'LineStyle', '-');
    xlim([-0.15 0.05]);
    
    %     % --- filter (high pass) - and square (power)
    %     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    %     datthis_high_sq = datthis_high.^2;
    %     clim = prctile(datthis_high_sq(:), [5 95]);
    %     lt_neural_Coher_Plot(datthis_high_sq, x, f, 1, '', clim);
    
    % --- filter (low pass)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('8-350Hz');
    [filt_b,filt_a]=butter(2,[8*2/fs 350*2/fs]);
    datthis_low =filtfilt(filt_b,filt_a,double(datthis));
    clim = prctile(datthis_low(:), [2.5 97.5]);
    lt_neural_Coher_Plot(datthis_low, x, f, 1, '', clim);
    line([-0.1 -0.1], ylim, 'Color', 'w', 'LineStyle', '-');
    line([-0 -0], ylim, 'Color', 'w', 'LineStyle', '-');
    xlim([-0.15 0.05]);
    
    % --- filter (band pass)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title([num2str(filt_low) '-' num2str(filt_hi) 'Hz']);
    [filt_b,filt_a]=butter(2,[filt_low*2/fs filt_hi*2/fs]);
    datthis_band =filtfilt(filt_b,filt_a,double(datthis));
    clim = prctile(datthis_band(:), [2.5 97.5]);
    lt_neural_Coher_Plot(datthis_band, x, f, 1, '', clim);
    line([-0.1 -0.1], ylim, 'Color', 'w', 'LineStyle', '-');
    line([-0 -0], ylim, 'Color', 'w', 'LineStyle', '-');
    xlim([-0.15 0.05]);
    
    
    % --- filter (noise band)
    if ~isempty(noiseband)
        % -- params
        movingwin = [0.1 0.01];
        params = struct;
        params.fpass = [1/movingwin(1) noiseband(2)+50];
        w = 30; % in hz, for desired frequency resolution of tapers. % note, t is set to movingwin(1)
        % w = 20; % in hz, for desired frequency resolution of tapers. % note, t is set to movingwin(1)
        tw = movingwin(1)*w;
        params.tapers = [tw 2*tw-1];
        params.Fs = fs; % hard coded fs for LFP;
        
        [S, tSp, fSp] = mtspecgramc(datthis,movingwin,params);
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        
        S = 10*log10(S);
        tSp = linspace(x(1)+movingwin(1)/2, x(end)-movingwin(1)/2, size(S,1));
        
        % -- get noise band
        ftoget = fSp>noiseband(1) & fSp<noiseband(2);
        S2 = squeeze(mean(S(:, ftoget, :),2));
        
        clim = prctile(S2(:), [2.5 97.5]);
        lt_neural_Coher_Plot(S2, tSp, f, 1, '', clim);
        line([-0.1 -0.1], ylim, 'Color', 'w', 'LineStyle', '-');
        line([-0 -0], ylim, 'Color', 'w', 'LineStyle', '-');
    end
    
end


%% ======= plot Phi
for i=1:size(phi_chanpair,1)
    cpthis = phi_chanpair(i,:);
    phithis = phimat(:,:,:, i);
    
    % --- get tri9als
    phithis = phithis(:,:, triallist);
    % -- get scalar
    assert(length(PARAMS.ffbins)==size(phithis,2));
    assert(length(PARAMS.tbins)==size(phithis,1));
    indtmp = PARAMS.ffbins>cohwind_f(1) & PARAMS.ffbins<cohwind_f(2);
    
    phiscal = squeeze(circ_mean(phithis(:, indtmp, :), [], 2));
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(['phi, ' num2str(cpthis)]);
    clim = [-pi pi];
    lt_neural_Coher_Plot(phiscal, PARAMS.tbins, 1:size(phiscal,2), 1, '', clim);
    line([-0.1 -0.1], ylim, 'Color', 'w', 'LineStyle', '-');
    line([-0 -0], ylim, 'Color', 'w', 'LineStyle', '-');
    xlim([-0.15 0.05]);
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(['phi, ' num2str(cpthis)]);
    clim = [0 pi];
    lt_neural_Coher_Plot(abs(phiscal), PARAMS.tbins, 1:size(phiscal,2), 1, '', clim);
    lt_plot_colormap('pval')
    line([-0.1 -0.1], ylim, 'Color', 'w', 'LineStyle', '-');
    line([-0 -0], ylim, 'Color', 'w', 'LineStyle', '-');
    xlim([-0.15 0.05]);
end


% ========= plot COHERENCE FOR ALL PAIRS FOR THIS TRIAL
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('coherence');
xlabel('coh');
hsplots = [hsplots hsplot];
pcols = lt_make_plot_colors(length(cohscal_OUT), 0,0);
for cc=1:length(cohscal_OUT)
    cohthis = cohscal_OUT{cc}(triallist);
    chanpairthis = chanpairs_OUT(cc,:);
    ythis = 1:size(datthis_high,2);
    %     plot(cohthis, ythis, '-', 'Color', pcols{cc});
    plot(cohthis, ythis, '.', 'Color', pcols{cc});
    lt_plot_text(cohthis(end), ythis(end)+1, num2str(chanpairthis), pcols{cc}, 8);
end
axis tight;


% ========= plot phi FOR ALL PAIRS FOR THIS TRIAL
% -- get phi scalar for all
    indT = PARAMS.tbins>PARAMS.cohscal_twind(1) & PARAMS.tbins<PARAMS.cohscal_twind(2);
    indF = PARAMS.ffbins>cohwind_f(1) & PARAMS.ffbins<cohwind_f(2);
    
    phiscal = squeeze(circ_mean(circ_mean(phimat(indT, indF, triallist, :), [], 2), [], 1));

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('phi');
xlabel('phi');
hsplots = [hsplots hsplot];
pcols = lt_make_plot_colors(size(phi_chanpair,1), 0,0);
for cc=1:size(phi_chanpair,1)
    phithis = phiscal(:,cc);
    chanpairthis = phi_chanpair(cc,:);
    ythis = 1:length(phithis);
    %     plot(cohthis, ythis, '-', 'Color', pcols{cc});
    plot(phithis, ythis, '.', 'Color', pcols{cc});
    lt_plot_text(phithis(end), ythis(end)+1, num2str(chanpairthis), pcols{cc}, 8);
end
axis tight;
lt_plot_zeroline_vert;


% ============= plot actual trial number
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('trial num');
xlabel('actual trial num');
ylabel('sorted trial num');
hsplots = [hsplots hsplot];
plot(triallist, 1:length(triallist), '.');
axis tight;
line([trialWNon+0.5 trialWNon+0.5], ylim);





% ========= fromating
linkaxes(hsplots, 'y');



lt_switch_chronux(0);
