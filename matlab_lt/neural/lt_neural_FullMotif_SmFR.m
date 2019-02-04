function lt_neural_FullMotif_SmFR(MOTIFSTATS_Compiled, birdtoplot, ...
    motiftoplot, plotcv, bregionsToPlot, plotRaw, doshuffle, whichhalf)


lt_figure; hold on;
lt_plot_text(0,1, 'NOTE: make sure edges do not influence xcov...');
lt_plot_text(0,1, 'NOTE: currently using entire time window');
%% lt - plots smoothed fr for a given motif, across all channels (diff brain regions)

% === NOTE: will assert that has been timewarped (all neurons simultaneously) ...

%%
% ------------ OPTIONS
i = find(strcmp({MOTIFSTATS_Compiled.birds.birdname}, birdtoplot));
% i=2; % birdnum
mm = find(strcmp(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.params.motif_regexpr_str, motiftoplot));
% mm = 2; % motif num
% plotcv = 0; % if 1, then plots running cv. if 0 then plots fr and running mean.
% bregionsToPlot = {'RA', 'LMAN'};

assert(length(i)==1, 'asdfsd');
assert(length(mm)==1, 'sdafsd');

%%

doztransform = 1; % fr trace (performs on the mean);
if plotcv==1
    doztransform=0;
end


%% ------------- RUNS


%% ------------- COLLECT AND PLOT EACH NEURON
plotcols = lt_make_plot_colors(length(bregionsToPlot), 0,0);
motifthis = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.params.motif_regexpr_str{mm};

if plotRaw==1
    % ------ one figure for each motif
    figcount=1;
    subplotrows=6;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
end

numneurons = length(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons);

% ------- to collect
AllBregion = {};
AllFRsm = {};
AllFRsmActual = {}; % will be used if cv is taking up AllFRsm
AllOnsets = [];
AllOffsets =[];
for ii=1:numneurons
    
    bregionthis = MOTIFSTATS_Compiled.birds(i).SummaryStruct.birds(1).neurons(ii).NOTE_Location;
    pcol = plotcols{strcmp(bregionsToPlot, bregionthis)};
    
    % ============= PLOT ALL TRIALS FOR THIS NEURON
    segextract = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(mm).SegmentsExtract;
    
    if isempty(segextract)
        continue
    end
%     
%     if doshuffle==1
%         for tt=1:length(segextract)
%             
%             segextract(tt).spk_Times 
%             
%     end
    
    % -- make sure has been timewarped
    assert(isfield(segextract, 'LinTWSeg_OriginalSpkTimes'), 'need to time warp!!');
    
    segextract = lt_neural_SmoothFR(segextract);
    
    x = segextract(1).FRsmooth_xbin_CommonTrialDur;
    frmat = [segextract.FRsmooth_rate_CommonTrialDur];
    if plotcv==0
        y = mean(frmat,2);
    elseif plotcv==1
        y = std(frmat, 0, 2);
        y = y./mean(frmat,2);
    end
    
    if plotRaw==1
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(['n' num2str(ii) '-' motifthis '-'  bregionthis]);
        
        if plotcv==0
            plot(x, frmat, '-', 'Color', [0.7 0.7 0.7]);
            plot(x, y, '-', 'LineWidth', 3, 'Color', pcol);
        elseif plotcv==1
            
            plot(x, y, '-', 'LineWidth', 3, 'Color', pcol);
        end
        axis tight
        lt_plot_zeroline;
        YLIM = ylim;
        
    end
    
    
    AllFRsmActual = [AllFRsmActual; mean(frmat,2)];
    
    % --- overlay syl onsets/offsets
    onsets = segextract(1).motifsylOnsets;
    offsets = segextract(1).motifsylOffsets;
    
    if plotRaw==1
        for j=1:length(onsets)
            line([onsets(j) offsets(j)], [3 3], 'Color', 'k', 'LineWidth', 4);
            line([onsets(j) offsets(j)], [YLIM(2)-4 YLIM(2)-4], 'Color', 'k', 'LineWidth', 4);
            line([onsets(j) onsets(j)], ylim, 'Color', 'r');
            line([offsets(j) offsets(j)], ylim, 'Color', 'r');
        end
    end
    
    % =========================== COLLECT FOR PLOTTING
    AllBregion = [AllBregion;  bregionthis];
    AllFRsm = [AllFRsm; y];
    AllOnsets = [AllOnsets; onsets];
    AllOffsets =[AllOffsets; offsets];
end

if plotRaw==1
    linkaxes(hsplots, 'x');
end


%% ===== maximum timepoint shared by all trial

% XMAX = 
%% ================== SPLIT DATA BY HALVES?

if whichhalf==1
    
    AllFRsm
    
elseif whichhalf==2
    
elseif isempty(whichhalf)
   % then keep all data.
  
end

%% =================== each brain region own plot, overlay all neurons

lt_figure; hold on;

% ======================
for i=1:length(bregionsToPlot)
    bregion = bregionsToPlot{i};
    
    lt_subplot(2,1,i); hold on;
    title(bregion);
    
    indstoplot = find(strcmp(AllBregion, bregion));
    
    tmp = cellfun('length', AllFRsm(indstoplot));
    xmax = min(tmp);
    
    for k = 1:length(indstoplot)
        j = indstoplot(k);
        
        pcol = 0.2+0.6*[rand rand rand];
        
        y = AllFRsm{j}(1:xmax);
        x = 0:0.001:0.001*(length(y)-1);
        
        % -------- z-transform fr
        if doztransform==1
            y = zscore(y);
        end
        
        % ------- plot
        plot(x,y+3*(k-1), '-', 'Color', pcol, 'LineWidth', 2);
        % plot(x,y+3*(k-1), '-', 'Color', 'k', 'LineWidth', 2);
        
    end
    
    % === overlay onsets/offsets
    onsets = median(AllOnsets,1);
    offsets = median(AllOffsets,1);
    axis tight
    YLIM = ylim;
    
    for j=1:length(onsets)
        line([onsets(j) offsets(j)], [YLIM(1)+1 YLIM(1)+1], 'Color', 'k', 'LineWidth', 4);
        line([onsets(j) offsets(j)], [YLIM(2)-1 YLIM(2)-1], 'Color', 'k', 'LineWidth', 4);
        line([onsets(j) onsets(j)], ylim, 'Color', 'k');
        line([offsets(j) offsets(j)], ylim, 'Color', 'k');
    end
    
end

%%
if doztransform==1
    AllFRsm = cellfun(@zscore, AllFRsm, 'UniformOutput', 0);
    AllFRsmActual = cellfun(@zscore, AllFRsmActual, 'UniformOutput', 0);
end

%% ================= [HEAT MAP, OVER ALL NEURONS]

lt_figure; hold on;
xmax = min(cellfun('length', AllFRsm));

% ======================
for i=1:length(bregionsToPlot)
    bregion = bregionsToPlot{i};
    
    lt_subplot(2,1,i); hold on;
    title(bregion);
    
    indstoplot = find(strcmp(AllBregion, bregion));
    
    %     tmp = cellfun('length', AllFRsm(indstoplot));
    
    y = cellfun(@(x)x(1:xmax)', AllFRsm(indstoplot), 'UniformOutput', 0);
    y = cell2mat(y);
    x = 0:0.001:0.001*(size(y,2)-1);
    xx = 1:size(y,1);
    imagesc(x, xx, y, [-3 3]);
    lt_plot_colormap('centered');
    colorbar('EastOutside');
    
    % === overlay onsets/offsets
    onsets = median(AllOnsets,1);
    offsets = median(AllOffsets,1);
    axis tight
    YLIM = ylim;
    
    for j=1:length(onsets)
        line([onsets(j) offsets(j)], [YLIM(1)+1 YLIM(1)+1], 'Color', 'k', 'LineWidth', 4);
        line([onsets(j) offsets(j)], [YLIM(2)-1 YLIM(2)-1], 'Color', 'k', 'LineWidth', 4);
        line([onsets(j) onsets(j)], ylim, 'Color', 'k');
        line([offsets(j) offsets(j)], ylim, 'Color', 'k');
    end
    
end

%%  ======================== PLOT OVELRAY
lt_figure; hold on;

% ====================== 1) overlay all trials
lt_subplot(3,2,1); hold on;

tmp = cellfun('length', AllFRsm);
xmax = min(tmp);
for j=1:length(AllFRsm)
    
    pcol = plotcols{strcmp(bregionsToPlot, AllBregion{j})};
    
    y = AllFRsm{j}(1:xmax);
    x = 0:0.001:0.001*(length(y)-1);
    
    % -------- z-transform fr
    if doztransform==1
        y = zscore(y);
    end
    
    % ------- plot
    plot(x,y, '-', 'Color', pcol);
    
end

% ====================== 2) plot each with mean
splotcount = 1;
for i=1:length(bregionsToPlot)
    splotcount = splotcount+1;
    indstoplot = find(strcmp(AllBregion, bregionsToPlot{i}));
    lt_subplot(3,2,splotcount); hold on;
    title(bregionsToPlot{i});
    
    Yall = [];
    YallFR = [];
    for j=indstoplot'
        
        pcol = plotcols{strcmp(bregionsToPlot, AllBregion{j})};
        
        y = AllFRsm{j}(1:xmax);
        x = 0:0.001:0.001*(length(y)-1);
        
        % -------- z-transform fr
        if doztransform==1
            y = zscore(y);
        end
        
        % ------- plot
        plot(x,y, '-', 'Color', pcol);
        
        Yall = [Yall; y'];
        
        % ----- collect mean FR, even if doing cv
        YallFR = [YallFR; zscore(AllFRsmActual{j}(1:xmax))'];
        
    end
    
    % -- plot mean
    plot(x, mean(Yall,1), '-', 'LineWidth', 3, 'Color', pcol);
    
    
    if plotcv==1
        % -- then overlay fr mean
        plot(x, mean(YallFR,1), '-k');
        
    end
    
    % --- overlay syllables
    onsets = median(AllOnsets,1);
    offsets = median(AllOffsets,1);
    axis tight
    YLIM = ylim;
    
    for j=1:length(onsets);
        line([onsets(j) offsets(j)], [YLIM(1)+1 YLIM(1)+1], 'Color', 'k', 'LineWidth', 4);
        line([onsets(j) offsets(j)], [YLIM(2)-1 YLIM(2)-1], 'Color', 'k', 'LineWidth', 4);
        line([onsets(j) onsets(j)], ylim, 'Color', 'r');
        line([offsets(j) offsets(j)], ylim, 'Color', 'r');
    end
    
end


% ====================== 2) overlay each bregion mean
lt_subplot(3,2,splotcount+1); hold on;

for i=1:length(bregionsToPlot)
    
    indstoplot = find(strcmp(AllBregion, bregionsToPlot{i}));
    pcol = plotcols{i};
    Yall = [];
    YallFR = [];
    for j=indstoplot'
        
        
        y = AllFRsm{j}(1:xmax);
        x = 0:0.001:0.001*(length(y)-1);
        
        % -------- z-transform fr
        if doztransform==1
            y = zscore(y);
        end
        
        Yall = [Yall; y'];
        
        % ----- collect mean FR, even if doing cv
        YallFR = [YallFR; zscore(AllFRsmActual{j}(1:xmax))'];
        
    end
    
    % -- plot mean
    shadedErrorBar(x, mean(Yall,1), lt_sem(Yall), {'Color', pcol},1);
    
end
% --- overlay syllables
onsets = median(AllOnsets,1);
offsets = median(AllOffsets,1);
axis tight
YLIM = ylim;

for j=1:length(onsets);
    line([onsets(j) offsets(j)], [YLIM(1)+1 YLIM(1)+1], 'Color', 'k', 'LineWidth', 4);
    line([onsets(j) offsets(j)], [YLIM(2)-1 YLIM(2)-1], 'Color', 'k', 'LineWidth', 4);
    line([onsets(j) onsets(j)], ylim, 'Color', 'r');
    line([offsets(j) offsets(j)], ylim, 'Color', 'r');
end




%% ====================== CROSS CORRELATION
xcovwindmax = 0.1;
xbin = 0.001;
if plotcv==0
    
    % --- only do if is for FR
    
    Ybyregion = [];
    % ============================== LMAN
    indstoplot = find(strcmp(AllBregion, 'LMAN'));
    Yall = [];
    for j=indstoplot'
        
        y = AllFRsm{j}(1:xmax);
        %         x = 0:0.001:0.001*(length(y)-1);
        
        % -------- z-transform fr
        if doztransform==1
            y = zscore(y);
        end
        
        Yall = [Yall; y'];
    end
    Ybyregion = [Ybyregion; mean(Yall,1)];
    
    
    % ============================== RA
    indstoplot = find(strcmp(AllBregion, 'RA'));
    Yall = [];
    for j=indstoplot'
        
        y = AllFRsm{j}(1:xmax);
        %         x = 0:0.001:0.001*(length(y)-1);
        
        % -------- z-transform fr
        if doztransform==1
            y = zscore(y);
        end
        
        Yall = [Yall; y'];
    end
    Ybyregion = [Ybyregion; mean(Yall,1)];
    
    
    % ================ DO XCORR
    [cc, lags] = xcov(Ybyregion(1,:), Ybyregion(2,:), xcovwindmax/0.001, 'coeff');
    %     [cc, lags] = xcov(Ybyregion(1,:), Ybyregion(2,:), xcovwindmax/0.001, 'unbiased');
    %     [cc, lags] = xcov(Ybyregion(1,:), Ybyregion(2,:), xcovwindmax/0.001);
    lt_figure; hold on ;
    lt_subplot(2,2,1); hold on;
    xlabel('LMAN ---- RA');
    ylabel('xcov');
    plot(lags*0.001, cc, '-k');
    axis tight;
    lt_plot_zeroline_vert;
    
    
    % ================ GET COHERENCE
    % ======== CHORNUX PARAMS
    lt_switch_chronux(1);
    
    params = struct;
    params.Fs = 1/xbin;
    %     tw = 5;
    %     params.tapers = [tw 2*tw-1];
    params.tapers = [5 xmax/params.Fs 1];
    params.fpass = [5 60];
    params.err = [2 0.05];
    
    [C,phi,S12,S1,S2,f, confC,phistd,Cerr] = coherencyc(Ybyregion(1,:), Ybyregion(2,:), ...;
        params);
    
    %     lt_figure; hold on;
    lt_subplot(2,2,2); hold on;
    %     plot(f, C, 'g');
    ylabel('coherence');
    xlabel('f');
    shadedErrorBar(f, C, Cerr, {'Color', 'k'},1);
    axis tight;
    
    
    % ================= SPEC
    % ------- LMAN
    [S,f_spec]=mtspectrumc(Ybyregion(1,:),params);
    lt_subplot(2,2,3); hold on;
    title('SPEC, LMAN');
    xlabel('f');
    ylabel('spec power');
    plot(f_spec, S);
    axis tight;
    
    % ------- RA
    [S,f_spec]=mtspectrumc(Ybyregion(2,:),params);
    lt_subplot(2,2,4); hold on;
    title('SPEC, RA');
    xlabel('f');
    ylabel('spec power');
    plot(f_spec, S);
    axis tight;
    lt_switch_chronux(0);
end


%% =================== ztramnsform all data if wanted
if doztransform==1
    AllFRsm = cellfun(@zscore, AllFRsm, 'UniformOutput', 0);
    AllFRsmActual = cellfun(@zscore, AllFRsmActual, 'UniformOutput', 0);
end
%% ================= XCOV AND COHERENCE ACROSS ALL PAIRS.
xcovwindmax = 0.1;
xbin = 0.001;
if plotcv==0
    
    xcorr_all = [];
    coh_all = [];
    phi_all = [];
    spec_all = [];
    
    % ======== CHORNUX PARAMS
    lt_switch_chronux(1);
    
    params = struct;
    params.Fs = 1/xbin;
    %     tw = 5;
    %     params.tapers = [tw 2*tw-1];
    params.tapers = [5 xmax/params.Fs 1];
    params.fpass = [5 60];
    params.err = [2 0.05];
    
    
    % ============= RUN
    ncase = length(AllBregion);
    for i=1:ncase
        disp(i);
        for ii=i+1:ncase
            
            
            if strcmp(AllBregion{i}, 'LMAN') & strcmp(AllBregion{ii}, 'RA')
                dat1 = AllFRsm{i}(1:xmax);
                dat2 = AllFRsm{ii}(1:xmax);
            elseif strcmp(AllBregion{i}, 'RA') & strcmp(AllBregion{ii}, 'LMAN')
                dat1 = AllFRsm{ii}(1:xmax);
                dat2 = AllFRsm{i}(1:xmax);
            else
                continue
            end
            
            % ===================== 1) XCORR
            [cc, lags] = xcov(dat1, dat2, xcovwindmax/xbin, 'coeff');
            xcorr_all = [xcorr_all; cc'];
            
            % ===================== 2) COHERENCE
            [C,phi,S12,S1,S2,f, confC,phistd,Cerr] = coherencyc(dat1, dat2, ...;
                params);
            coh_all = [coh_all; C'];
            phi_all = [phi_all; phi'];
            
            
            
            
        end
        
        % ===================== 3) SPECTRUM
        datthis = AllFRsm{i}(1:xmax);
        [S,f_spec]=mtspectrumc(datthis,params);
        spec_all = [spec_all; S'];
    end
    
    
    % ############ PLOT
    lt_switch_chronux(0);
    
    lt_figure; hold on;
    
    % =====
    lt_subplot(2,2,1); hold on;
    title('xcorr (coeff)');
    xlabel('lag(s) [LMAN leads ---- RA leads]');
    plot(lags*xbin, xcorr_all', '-', 'Color', [0.7 0.7 0.7]);
    ymean = mean(xcorr_all, 1);
    ysem = lt_sem(xcorr_all);
    shadedErrorBar(lags*xbin, ymean, ysem, {'Color', 'k'},1);
    axis tight;
    lt_plot_zeroline_vert;
    
    
    % =====
    lt_subplot(2,2,2); hold on;
    title('coherence, all chan pairs');
    xlabel('f');
    ylabel('coh');
    plot(f, coh_all', '-', 'Color', [0.7 0.7 0.7]);
    ymean = mean(coh_all, 1);
    ysem = lt_sem(coh_all);
    shadedErrorBar(f, ymean, ysem, {'Color', 'k'},1);
    axis tight;
    %     lt_plot_zeroline_vert;
    
    % =====
    lt_subplot(2,2,3); hold on;
    title('SPEC, LMAN');
    xlabel('f');
    ylabel('spec power');
    indstmp = strcmp(AllBregion, 'LMAN');
    %     y = 10*log10(spec_all(indstmp,:));
    y = spec_all(indstmp,:);
    plot(f_spec, y', '-', 'Color', [0.7 0.7 0.7]);
    ymean = mean(y, 1);
    ysem = lt_sem(y);
    shadedErrorBar(f_spec, ymean, ysem, {'Color', 'k'},1);
    axis tight;
    %     lt_plot_zeroline_vert;
    
    
    % =====
    lt_subplot(2,2,4); hold on;
    title('SPEC, RA');
    xlabel('f');
    ylabel('spec power');
    indstmp = strcmp(AllBregion, 'RA');
    %     y = 10*log10(spec_all(indstmp,:));
    y = spec_all(indstmp,:);
    plot(f_spec, y', '-', 'Color', [0.7 0.7 0.7]);
    ymean = mean(y, 1);
    ysem = lt_sem(y);
    shadedErrorBar(f_spec, ymean, ysem, {'Color', 'k'},1);
    axis tight;
    %     lt_plot_zeroline_vert;
    
    % ============ COHERENCE PHI
    lt_figure; hold on;
    lt_subplot(2,2,1); hold on;
    xlabel('f');
    ylabel('phi (pos = RA lead) [std]');
%     plot(f, phi_all', '-', 'Color', [0.7 0.7 0.7]);
%     circ_stats(phi_all(1,:), []); 
    ymean = circ_mean(phi_all, [], 1);
    ystd = circ_std(phi_all, [], [], 1);
    shadedErrorBar(f, ymean, ystd, {'Color', 'k'}, 1);
    axis tight;
    lt_plot_zeroline;
end






