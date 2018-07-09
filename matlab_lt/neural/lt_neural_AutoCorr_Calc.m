function lt_neural_AutoCorr_Calc(MOTIFSTATS_Compiled)

normmethod = 'unbiased'; % for xcorr
xcovwindmax = 0.1; % sec

if ~exist('binsize_spk', 'var')
    binsize_spk = 0.001; % default, 5ms bins for cross corr
end

%% =========== PREPROCESSING RECOMMNENDED:

% See steps in lt_neural_v2_FullMotifActivity;
% want to do linear time warping so that can get across trial correlations
% (i.e. shift predictor). I have only done this on full motif acrtivity,
% but should also work for partial motif.

% === NOTE: will assert that has been timewarped (all neurons simultaneously) ...

%%

OUTSTRUCT.AllBirdnum = [];
OUTSTRUCT.AllMotifnum = [];
OUTSTRUCT.AllNeurnum = [];
OUTSTRUCT.AllBregion = {};
OUTSTRUCT.AllCCraw = {};
OUTSTRUCT.AllCCshift = {};

% ------------ OPTIONS
Numbirds = length(MOTIFSTATS_Compiled.birds);
for i=1:Numbirds
    
    Motiflist = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.params.motif_regexpr_str;
    Nummotifs = length(Motiflist);
    for mm=1:Nummotifs
        
        motifthis = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.params.motif_regexpr_str{mm};
        
        % ===================== COLLECT for each NEURON
        numneurons = length(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons);
        
        for ii=1:numneurons
            
            bregionthis = MOTIFSTATS_Compiled.birds(i).SummaryStruct.birds(1).neurons(ii).NOTE_Location;
            
            % ============= PLOT ALL TRIALS FOR THIS NEURON
            segextract = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(mm).SegmentsExtract;
            
            if isempty(segextract)
                continue
            end
            
            % -- make sure has been timewarped
            assert(isfield(segextract, 'LinTWSeg_OriginalSpkTimes'), 'need to time warp!!');
            
            
            
            % ########################## CALCULATE AUTOCORR (using binned spikes)
            % ------------ what is maximum time for timebins? since must be
            % time warped, use single trial maximum time
            maxtime = segextract(1).motifsylOffsets(end);
            segextract = lt_neural_QUICK_SpkBinned(segextract, maxtime, ...
                binsize_spk, 1);
            
            
            spkmat = [segextract.spk_Binned];
            spkmat_t = segextract(1).spk_Binned_x;
            
            
            % ====================== prune time windows
            % x window to take, from onset of token (i.e. don't go
            % as far back as onset of data, since lin time warp so
            % some trials could be shorter than this...) to offset of
            % last syl
            
            tmin = MOTIFSTATS_Compiled.birds(1).Params_regexp.motif_predur; % this is where time warping aligned to (i.e. sure that data not clipped here)
            indstoremove = spkmat_t<tmin;
            spkmat(indstoremove,:) = [];
            spkmat_t(indstoremove) = [];
            
            % =================== CC raw
            nrends = size(spkmat,2);
            CCraw = [];
            CCshift = [];
            for j=1:nrends
                
                
                % ============================ CC RAW
                y1 = spkmat(:,j);
                y2 = spkmat(:,j);
                
                [cc, lags] = xcorr(y1, y2, xcovwindmax/binsize_spk, normmethod);
                
                CCraw = [CCraw; cc'];
                
                
                % =========================== CC SHIFT
                % ---- take next trial for y2. if end of trials,
                % then take 3rd to last trial
                if j<nrends
                    y2 = spkmat(:,j+1);
                elseif j==nrends
                    y2 = spkmat(:, j-2);
                end
                
                [cc, lags] = xcorr(y1, y2, xcovwindmax/binsize_spk, normmethod);
                
                CCshift = [CCshift; cc'];
            end
            
            % ################################## OUTPUT
            OUTSTRUCT.AllBirdnum = [OUTSTRUCT.AllBirdnum; i];
            OUTSTRUCT.AllMotifnum = [OUTSTRUCT.AllMotifnum; mm];
            OUTSTRUCT.AllNeurnum = [OUTSTRUCT.AllNeurnum; ii];
            OUTSTRUCT.AllBregion = [OUTSTRUCT.AllBregion; bregionthis];
            OUTSTRUCT.AllCCraw = [OUTSTRUCT.AllCCraw; CCraw];
            OUTSTRUCT.AllCCshift = [OUTSTRUCT.AllCCshift; CCshift];
            
            % ========================================= SANITY CHECK
            if (0)
                lt_figure; hold on;
                cc = mean(CCraw,1) - mean(CCshift,1);
                plot(lags*binsize_spk, cc, '-k');
            end
            
            
        end
    end
end

OUTPARAMS.lags = lags;
OUTPARAMS.binsize = binsize_spk;

%%

Maxbirds = max(OUTSTRUCT.AllBirdnum);
Maxmotifs = max(OUTSTRUCT.AllMotifnum);
Maxneur = max(OUTSTRUCT.AllNeurnum);
Bregionstoplot = unique(OUTSTRUCT.AllBregion)';

%% ============ [PLOT] Summary

pcolors = lt_make_plot_colors(length(Bregionstoplot), 0,0);
Xlags = OUTPARAMS.lags*OUTPARAMS.binsize;

figcount=1;
subplotrows=5;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for i=1:Maxbirds
    
    for ii=1:Maxmotifs
        
        inds = OUTSTRUCT.AllBirdnum==i & OUTSTRUCT.AllMotifnum==ii;
        if ~any(inds)
            continue
        end
        
        for bregion = Bregionstoplot
            
            indstoget = find(OUTSTRUCT.AllBirdnum==i & OUTSTRUCT.AllMotifnum==ii ...
                & strcmp(OUTSTRUCT.AllBregion, bregion))';
     
            if isempty(indstoget)
                continue
            end
            
            pcol = pcolors{strcmp(Bregionstoplot, bregion)};
            
            % ############################### 1) PLOT RAW AND SHIFTED XCORR
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots; hsplot];
            bname = MOTIFSTATS_Compiled.birds(i).birdname;
            motifname = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.params.motif_regexpr_str{ii};
            title([bname '-' motifname '-' bregion{1}], 'Color', pcol);
            xlabel('lag');
            ylabel('raw xcorr (raw:bk; shift:rd)');
            
            
            % =============== for each neuron plot xcov
            for j=indstoget
                
                ccraw = OUTSTRUCT.AllCCraw{j};
                ccshift = OUTSTRUCT.AllCCshift{j};
                
                % ---- take means across trials
                ccraw_mean = mean(ccraw,1);
                ccraw_std = std(ccraw, [],1);
                
                ccshift_mean = mean(ccshift,1);               
                ccshift_std = std(ccshift,[],1);               
                
                if (0)
                shadedErrorBar(OUTPARAMS.lags*OUTPARAMS.binsize, ccraw_mean, ...
                    ccraw_std, {'Color' ,'k'}, 1);
                shadedErrorBar(OUTPARAMS.lags*OUTPARAMS.binsize, ccshift_mean, ...
                    ccshift_std, {'Color' ,'r'}, 1);
                else
                plot(Xlags, ccraw_mean, '-', 'Color' ,'k');
                    plot(Xlags, ccshift_mean, '-', 'Color' ,'r');
                end
            end
            
            
            
            % ############################ 1) PLOT DIFFERENCE (XCOVARIOGRAM)
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots; hsplot];
            title([bname '-' motifname '-' bregion{1}], 'Color', pcol);
            xlabel('lag');
            ylabel('xcovariogram');
            
            % =============== for each neuron plot xcov
            cov_all = [];
            for j=indstoget
                
                ccraw = OUTSTRUCT.AllCCraw{j};
                ccshift = OUTSTRUCT.AllCCshift{j};
                
                cov = mean(ccraw,1) - mean(ccshift,1);
                
                plot(Xlags, cov, '-', 'Color', pcol);
                
                % ======== collect
                cov_all = [cov_all; cov];
            end
            
            % ================= PLOT MEANS
            cov_mean = mean(cov_all,1);
            cov_sem = lt_sem(cov_all);
            plot(Xlags, cov_mean, '-k', 'LineWidth', 3);

            % ================
            lt_plot_zeroline;
            lt_plot_zeroline_vert;
            
        end
    end
end



linkaxes(hsplots, 'x');




%% ========= [PLOT] for each brain region, show all motifs, all animals

makesymmetric = 1; % duplicates ccshift and flips.

pcolors = lt_make_plot_colors(length(Bregionstoplot), 0,0);
Xlags = OUTPARAMS.lags*OUTPARAMS.binsize;

figcount=1;
subplotrows=5;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for i=1:Maxbirds
    
    for bregion = Bregionstoplot
        
        pcol = pcolors{strcmp(Bregionstoplot, bregion)};
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots; hsplot];
        bname = MOTIFSTATS_Compiled.birds(i).birdname;
        %         motifname = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.params.motif_regexpr_str{ii};
        title([bname '-' bregion{1}], 'Color', pcol);
        xlabel('lag');
        ylabel('raw xcorr (raw:bk; shift:rd)');
        
        cov_all = [];
        for ii=1:Maxmotifs
            
            indstoget = find(OUTSTRUCT.AllBirdnum==i & OUTSTRUCT.AllMotifnum==ii ...
                & strcmp(OUTSTRUCT.AllBregion, bregion))';
            
            if isempty(indstoget)
                continue
            end
            
            % ############################### 1) PLOT RAW AND SHIFTED XCORR
            
            % =============== for each neuron plot xcov
            for j=indstoget % iterate thru neurons
                
                ccraw = OUTSTRUCT.AllCCraw{j};
                ccshift = OUTSTRUCT.AllCCshift{j};
                
                % ====================== make symmetric by flipping ccshift
                if makesymmetric==1
                    ccshift = [ccshift; fliplr(ccshift)];
                end
                
                % ================== get covariogram
                cov = mean(ccraw,1) - mean(ccshift,1);
                                
                % ==================== NORMALIZE TO 0 LAG
                indzero = Xlags==0; assert(sum(indzero)==1,'sadf');
                cov = cov./cov(indzero);
                
                % ==================== PLOT
                if makesymmetric==1
                    x = Xlags(Xlags>=0);
                    y = cov(Xlags>=0);
                else
                    x = Xlags;
                    y = cov;
                end                
                plot(x, y, '-', 'Color', pcol);
                
                % ======== collect
                cov_all = [cov_all; cov];
            end
            
        end
        
        % ========= plot mean
        cov_mean = mean(cov_all,1);
        cov_sem = lt_sem(cov_all);
        if makesymmetric==1
           x = Xlags(Xlags>=0);
           y = cov_mean(Xlags>=0);
           ysem = cov_sem(Xlags>=0);
        else
           x = Xlags;
           y = cov_mean;
           ysem = cov_sem;
        end
        plot(x, y, 'Color', 'k', 'LineWidth', 3);
        
        % ================
        lt_plot_zeroline;
        lt_plot_zeroline_vert;
        
    end
end



linkaxes(hsplots, 'x');





















