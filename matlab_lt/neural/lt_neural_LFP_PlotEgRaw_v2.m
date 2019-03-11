function lt_neural_LFP_PlotEgRaw_v2(DatAll, t_onoff, fs, bregionlist, chanlist_toget, ...
    bnum, enum, swplot, motifnum, trialtoplot, filt_low, filt_hi, SwitchCohStruct, ...
    SwitchStruct, OUTSTRUCT, MOTIFSTATS_pop, PARAMS, ...
    SummaryStruct, OUTSTRUCT_CohMatOnly, plotUnfiltered)
%% lt 11/14/18 - plots output ofr lt_neural_LFP_PlotEgRaw_Extract
% NOTE: needs to run that function right before runnign this. run this as
% many times as desired for plotting examples. Can also save if desired.


%%
% if ~exist('plotCohScalExtremes', 'var')
%     plotCohScalExtremes=0;
% end

% chanpairs, n x 2, n mnumber of pairs (array)
% cohscal_allpairs, cell nx1, each is 1xtrials.


i=bnum;
ii=enum;
mm = motifnum;
filt_fs = fs;

bname = SwitchStruct.bird(i).birdname;
ename = SwitchStruct.bird(i).exptnum(ii).exptname;

%%
datcohstruct = SwitchCohStruct.bird(i).exptnum(ii).switchlist(swplot).motifnum(mm);

% -----
indsbase = datcohstruct.indsbase_epoch;
indsWN = datcohstruct.indsWN_epoch;
assert(size(datcohstruct.lfpall,2) == size(DatAll,1)) % n chans
assert(size(datcohstruct.lfpall,1) == size(DatAll{1},2)) % n trials

LFPdat = datcohstruct.lfpall;
t_lfp = datcohstruct.t_lfp;



%% ====== COLLECT SPIKING DATA

neurset = SwitchCohStruct.bird(bnum).exptnum(enum).switchlist(swplot).motifnum(motifnum).neursetused;
segAll = MOTIFSTATS_pop.birds(bnum).exptnum(enum).DAT.setnum(neurset).motif(motifnum).SegExtr_neurfakeID;

assert(length(segAll(1).SegmentsExtract)==size(DatAll{1},2));

pcols_chans = lt_make_plot_colors(size(DatAll,1), 0, 0);

%% WHAT trial to plot?
        lt_switch_chronux(1);

for i = 1:length(trialtoplot)
    
    j = trialtoplot(i);
    
    % ============================ PLOT FOR THIS TRIAL.
    figcount=1;
    subplotrows=7;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    % -- get ons and offwets of sound
    ons = segAll(1).SegmentsExtract(j).motifsylOnsets - PARAMS.motif_predur;
    offs = segAll(1).SegmentsExtract(j).motifsylOffsets - PARAMS.motif_predur;
    
    % ============================ PLOT ALL RAW NEURAL
    datfiltall = [];
    pcols_all = {};
    for cc = 1:length(chanlist_toget)
        
        % == what color to plot
        if strcmp(bregionlist{cc}, 'LMAN')
            pcol = [0.2 0.6 0.2];
        elseif strcmp(bregionlist{cc}, 'RA')
            pcol = 'r';
        end
        
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(['trial' num2str(j) ',chan' num2str(chanlist_toget(cc)) ',' bregionlist{cc}]);
        ylabel([bname '-' ename '-sw' num2str(swplot)]);
        %     lt_plot_annotation(1, ['coh=' num2str(cohscaltmp(j))], 'r');
        chanthis = chanlist_toget(cc);
        
        % =============== RAW DATA
        datmat = DatAll{cc};
        y = datmat(:, j);
        x = linspace(t_onoff{cc}(1), t_onoff{cc}(2), length(y));
        plot(x,y, 'Color', [0.6 0.6 0.6]);
        
        % ==================== HIGH PASS
        if plotUnfiltered==0
            [filt_b,filt_a]=butter(2,[250*2/fs], 'high');
            y=filtfilt(filt_b,filt_a,double(y));
            plot(x,y, 'Color', 'k');
        end
        
        %     lt_plot_text(x(end), y(end), ['ch' num2str(chanlist_toget(cc))], pcols_chans{cc});
        lt_plot_text(x(end), y(end), ['ch' num2str(chanlist_toget(cc))], pcol);
        
        % ==================== SPIKING
        % --- find units corresponding to this channel
        indstoget_n = find([SummaryStruct.birds(bnum).neurons([segAll.neurID_orig]).channel]  == chanthis);
        for i=1:length(indstoget_n)
            segthis = segAll(indstoget_n(i));
            spkthis = segthis.SegmentsExtract(j).spk_Times;
            spkthis = spkthis - PARAMS.motif_predur;
            % === plot raster
            lt_neural_PLOT_rasterline(spkthis, 0+20*i, 'r', 0, 50);
        end
        
        axis tight;
        
        line([-0.1 -0.1], ylim, 'Color', 'm', 'LineStyle', '-');
        line([-0 -0], ylim, 'Color', 'm', 'LineStyle', '-');
        
        
        
        % ############################ POWER SPECTRUM FOR RAW NEURAL
        y = datmat(:, j);
        % -- params
        movingwin = [0.1 0.01]; 
        params = struct;
        params.fpass = [1/movingwin(1) 4000];
        w = 30; % in hz, for desired frequency resolution of tapers. % note, t is set to movingwin(1)
        % w = 20; % in hz, for desired frequency resolution of tapers. % note, t is set to movingwin(1)
        tw = movingwin(1)*w;
        params.tapers = [tw 2*tw-1];
        params.Fs = fs; % hard coded fs for LFP;
        
        [S, tSp, fSp] = mtspecgramc(y,movingwin,params);
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        S = 10*log10(S);
%         S = (S-mean(S,1))./std(S,[],1);
        clim = prctile(S(:), [2.5 97.5]);
        tSp = linspace(x(1)+movingwin(1)/2, x(end)-movingwin(1)/2, size(S,1));
%         tSp = tSp - PARAMS.motif_predur;
        lt_neural_Coher_Plot(S, tSp, fSp, 1, '', clim);
        lt_plot_colormap('pval');
        
        
        % ##################################3
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(['trial' num2str(j) ',chan' num2str(chanlist_toget(cc)) ',' bregionlist{cc}]);
        
        % =================== LOW PASS
        y = datmat(:, j);
        [filt_b,filt_a]=butter(2,[8*2/fs 350*2/fs]);
        y=filtfilt(filt_b,filt_a,double(y));
        plot(x,y, 'Color', pcols_chans{cc});
        plot(x,y, 'Color', pcol);
        
        
        %     % ================= LFP
        %     datlfp = LFPdat{j,cc};
        %     %        plot(t_lfp, datlfp, 'Color', pcols_chans{cc}, 'LineWidth', 2);
        %     plot(t_lfp, datlfp, 'Color', 'k', 'LineWidth', 2);
        
        % =========== lfp (filtere)
        datfilt = lt_neural_filter(double(y), filt_fs, 0, filt_low, filt_hi);
        %        plot(x, datfilt, 'Color', pcols_chans{cc});
        plot(x, datfilt, 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
        
        axis tight;
        
        datfiltall = [datfiltall; datfilt'];
        pcols_all = [pcols_all; pcol];
        line([-0.1 -0.1], ylim, 'Color', 'm', 'LineStyle', '-');
        line([-0 -0], ylim, 'Color', 'm', 'LineStyle', '-');
    end
    
    
    % ========= FILTER ALL AND OVERLAY
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    for cc=1:size(DatAll,1)
        %     plot(x, datfiltall(cc,:), 'Color', pcols_chans{cc}, 'LineWidth', 2);
        plot(x, datfiltall(cc,:), 'Color', pcols_all{cc}, 'LineWidth', 2);
    end
    
    axis tight;
    line([-0.1 -0.1], ylim, 'Color', 'm', 'LineStyle', '-');
    line([-0 -0], ylim, 'Color', 'm', 'LineStyle', '-');
    YLIM = ylim;
    XLIM = xlim;
    lt_neural_QUICK_PlotSylPatches(ons, offs, [YLIM(end)-1 YLIM(end)]);
    xlim(XLIM);
    
    % ========= plot COHERENCE FOR ALL PAIRS FOR THIS TRIAL
    indsOUT = find(OUTSTRUCT.bnum==bnum & OUTSTRUCT.enum==enum & OUTSTRUCT.switch==swplot ...
        & OUTSTRUCT.motifnum==motifnum);
    chanpairs_OUT = OUTSTRUCT.chanpair(indsOUT,:);
    cohscal_OUT = OUTSTRUCT.cohscal(indsOUT);
    assert(all(cellfun('length', cohscal_OUT)==length(segAll(1).SegmentsExtract)));
    cohscal_OUT = cellfun(@(x)x(j), cohscal_OUT); % get just this trial.
    
    cohmat_all = OUTSTRUCT_CohMatOnly(indsOUT);
    assert(all(cellfun(@(x)size(x,3), cohmat_all)==length(segAll(1).SegmentsExtract)));
    cohmat_all = cellfun(@(x)x(:,:,j), cohmat_all, 'UniformOutput', 0); % get this trial;
    
    %  ---- PLOT
    for i=1:size(chanpairs_OUT,1)
        chanpairthis = chanpairs_OUT(i,:);
        cohscalthis = cohscal_OUT(i);
        cohmatthis = cohmat_all{i};
        
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        clim = [0.2 1];
        lt_neural_Coher_Plot(cohmatthis, PARAMS.tbins, PARAMS.ffbins, 1, '', clim);
        title(['chanpair: ' num2str(chanpairthis) ',coh=' num2str(cohscalthis)]);
        axis tight;
        ylim([10 80]);
    end
    
    
    linkaxes(hsplots, 'x');
    xlim([-0.15 0.05]);
    
    
    % [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    % title('coh scalar');
    % xlabel('chan pair');
    % cohthis = cellfun(@(x)x(j), cohscal_allpairs);
    % lt_plot_bar(1:length(cohthis), cohthis);
    % set(gca, 'XTick', 1:length(cohthis), 'XTickLabel', xlabel_pairs);
    % ylim([0 1]);
    % xlim([0 length(cohthis)+1]);
    %
    %
    % axis tight;
    % linkaxes(hsplots, 'xy');
end
        lt_switch_chronux(0);
