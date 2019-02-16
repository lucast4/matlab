function lt_neural_POPLEARN_PairRast(BirdExptPairsToPlot, SwitchToPlot, ...
    neurpair_globID, motiftoplot, OUTSTRUCT, OUTSTRUCT_XCOV, SwitchStruct, ...
    SwitchCohStruct, MOTIFSTATS_Compiled, MOTIFSTATS_pop, SummaryStruct, ...
    PARAMS, xgramxlim, xgramylim, clim)


%% ==================== PARAMS
% BirdExptPairsToPlot = {'pu69wh78', 'RALMANOvernightLearn1'}; % ordered p[aors {bird, expt, bird, expt } ....
% SwitchToPlot = [1]; % array of swiches.
% neurpair_globID = [17 19]; % only one pair allowed: e.g. [9 17], actual global ID.
% TypeOfPairToPlot = {'LMAN-RA'}; % e.g. 'LMAN-RA' (in alphabetical order)
% motiftoplot = {}; % if empty, then plots target(s).

% ========= PLOTTING PARAMS
Nrand = 20; % number of trials from base and WN to extract
% clim = [-0.125 0.125]; % for cov gram;

% defgault
corrwind = [-0.115 0.015]; % rel syl onset, to get pairwise trial by trial correlations.
% corrwind = [-0.1 0.0]; % rel syl onset, to get pairwise trial by trial correlations.
% 
% corrwind = [-0.05 0.02]; % rel syl onset, to get pairwise trial by trial correlations.
% corrwind = [-0.12 -0.05]; % rel syl onset, to get pairwise trial by trial correlations.

% =========== SMOOTHING PARAMS
kernelSD = 0.005;
binsize = 0.001;

% ========= NOTE;
% uses inds (trials) stores in OUTSTRUCT (i.e. OUTSTRUCT_XCOV doesn't
% store, but uses the same ones as in OUTSTRUCT0.
disp('NOTE: only using random subset for the rasters - rest use all trials in epoch');


neurpair_globID = sort(neurpair_globID); % must be in order, for code to work below.
%%
assert(all(size(neurpair_globID)==[1 2]));

%%
numbirds = length(SwitchStruct.bird);
for i=1:numbirds
    numexpts = length(SwitchStruct.bird(i).exptnum);
    birdname = SwitchStruct.bird(i).birdname;
    for ii=1:numexpts
        numswitches = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
        
        % ----------------- ONLY PLOT SPECIFIC BIRD?
        ind1 = find(strcmp(BirdExptPairsToPlot, birdname));
        ind2 = find(strcmp(BirdExptPairsToPlot, exptname));
        
        if ~any(ind1+1 == ind2)
            disp(['SKIPPED ' birdname '-' exptname]);
            continue
        end
        
        % ----------------- GO THRU ALL SWITCHES
        for iii=1:numswitches
            
            if ~isempty(SwitchToPlot)
                if ~any(SwitchToPlot == iii)
                    continue
                end
            end
            
            % ======= figure out which dataset to use
            neurset = SwitchCohStruct.bird(i).exptnum(ii).switchlist(iii).motifnum(1).neursetused;
            
            % ========= get segextract for the neuron pair (using the correct neuron set)
            
            DAT = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(neurset);
            neuridx_list = find(ismember(MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons{neurset}, neurpair_globID));
            assert(length(neuridx_list)==2, 'entered wrong neuron pair, must be length 2 array');
            assert(all(MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons{neurset}(neuridx_list) == neurpair_globID), 'not in correct order - fix code [or sort first intoorder..');
            bregion_list = {SummaryStruct.birds(i).neurons(neurpair_globID).NOTE_Location};
            
            
            
            % ########################################## GO THRU ALL MOTIFS
            % ====================== 1) get liust of motifst ot plot
            if isempty(motiftoplot)
                % ===== then find out what are valid targ motifs
                indstmp = OUTSTRUCT_XCOV.bnum==i & OUTSTRUCT_XCOV.enum==ii ...
                    & OUTSTRUCT_XCOV.switch==iii & OUTSTRUCT_XCOV.istarg==1;
            else
                % === get those desired.
                indstmp = OUTSTRUCT_XCOV.bnum==i & OUTSTRUCT_XCOV.enum==ii ...
                    & OUTSTRUCT_XCOV.switch==iii & ismember(OUTSTRUCT_XCOV.motifname, motiftoplot);
            end
            motiflist = unique(OUTSTRUCT_XCOV.motifnum(indstmp));
            
            % ################################ PLOT THE DESIRED MOTIFS.
            for m=1:length(motiflist)
                mm = motiflist(m);
                
                % ================= COLLECT DAT
                % --- 1) neural dat
                seg1 = DAT.motif(mm).SegExtr_neurfakeID(neuridx_list(1)).SegmentsExtract;
                assert(DAT.motif(mm).SegExtr_neurfakeID(neuridx_list(1)).neurID_orig==neurpair_globID(1));
                seg2 = DAT.motif(mm).SegExtr_neurfakeID(neuridx_list(2)).SegmentsExtract;
                assert(DAT.motif(mm).SegExtr_neurfakeID(neuridx_list(2)).neurID_orig==neurpair_globID(2));
                seg_global = DAT.motif(mm).SegExtr_neurfakeID(1).SegmentsExtract;
                
                % --- 2) other stuff
                motifname = DAT.motif(mm).regexpstr;
                motifpredur = MOTIFSTATS_Compiled.birds(1).MOTIFSTATS.params.motif_predur;
                
                %                indstmp = OUTSTRUCT_XCOV.bnum==i & OUTSTRUCT_XCOV.enum==ii ...
                %                    & OUTSTRUCT_XCOV.switch==iii & OUTSTRUCT_XCOV.motifnum==mm ...
                %                    & all(OUTSTRUCT_XCOV.neurpair'==neurpair_globID')';
                %                assert(sum(indstmp)==1, 'if not 1, then means crucial error somewjhere...');
                
                indstmp = find(OUTSTRUCT.bnum==i & OUTSTRUCT.enum==ii & OUTSTRUCT.switch==iii ...
                    & OUTSTRUCT.motifnum==mm);
                indstmp=indstmp(1); % assumes that multip[el cases are multip[el chans..
                assert(length(indstmp)==1);
                
                if (0)
                    inds_base = find(OUTSTRUCT.indsbase{indstmp});
                else
                    inds_base = OUTSTRUCT.indsbase_epoch{indstmp};
                end
                if (0) % entire WN
                    inds_wn = find(OUTSTRUCT.indsWN{indstmp});
                else % WN epoch that was used for analsyis.
                    inds_wn = OUTSTRUCT.indsWN_epoch{indstmp};
                end
                
                
                %% =================== PLOT RASTERS
                lt_figure;
                hsplots = [];
                
                hsplot =lt_subplot(8,4,[1:3 5:7 9:11]);
                hsplots = [hsplots hsplot];
                
                if strcmp(bregion_list{1}, 'RA')
                %                 pcol1 = 'r';
                pcol1 = [0.8 0.2 0.2];
                %                 pcol2 = [0.1 0.2 0.1];
                pcol2 = [0.1 0.2 0.1];
                else
                pcol2 = [0.8 0.2 0.2];
                pcol1 = [0.1 0.2 0.1];
                end
                
                % === plot random subset of inds from base and WN
                inds_base_rand = inds_base(sort(randperm(length(inds_base), Nrand)));
                inds_wn_rand = inds_wn(sort(randperm(length(inds_wn), Nrand)));
                
                indstoplot = [inds_base_rand inds_wn_rand];
                ind_wnon = length(inds_base_rand)+1;
                
                % ==== region 1
                Y = {seg1(indstoplot).spk_Times};
                Y = cellfun(@(x)x-motifpredur, Y, 'UniformOutput', 0);
                pcol = pcol1;
                yoffset = -0.1;
                
                for kk=1:length(Y)
                    lt_neural_PLOT_rasterline(Y{kk}, kk+yoffset, pcol, 0, 0.4);
                end
                lt_plot_text(0, kk+1, [bregion_list{1} '-' num2str(neurpair_globID(1))], pcol);
                
                
                % ==== region 2
                Y = {seg2(indstoplot).spk_Times};
                Y = cellfun(@(x)x-motifpredur, Y, 'UniformOutput', 0);
                pcol = pcol2;
                yoffset = 0.11;
                
                for kk=1:length(Y)
                    lt_neural_PLOT_rasterline(Y{kk}, kk+yoffset, pcol, 0, 0.4);
                end
                lt_plot_text(0, kk+2, [bregion_list{2} '-' num2str(neurpair_globID(2))], pcol);
                
                % ============ ANNOTATE
                axis tight;
                line(xlim, [ind_wnon ind_wnon]-0.5)
                xlabel('t rel syl');
                ylabel('trial (random subset, in order) [line=wnon]');
                title([birdname '-' exptname '-sw' num2str(iii) '-' motifname]);
                lt_plot_zeroline_vert;
                
                
                
                %% ============== PLOT SMOOTHED FR.
                % ====== 1) extract smoothed FR
                seg1 = lt_neural_SmoothFR(seg1, [], kernelSD, binsize, 0, seg_global);
                seg2 = lt_neural_SmoothFR(seg2, [], kernelSD, binsize, 0, seg_global);
                
                
                % ##################################### UNIT 1
                hsplot = lt_subplot(8,4,[13:15]); hold on;
                hsplots = [hsplots hsplot];
                title([bregion_list{1} '-' num2str(neurpair_globID(1))]);
                segthis = seg1;
                pcolthis = pcol1;
                ylabel('BASE(solid), WN(dash)');
                
                % ==== BASELINE
                indsthis = inds_base;
                linstyle = '-';
                
                if (0)
                    Ycell = {segthis(indsthis).spk_Times};
                    [xbin, ymean, ysem, ymean_hz, ysem_hz, ystd, ystd_hz] = ...
                        lt_neural_plotRastMean(Ycell, 0.01, 0.001, 0, []);
                    shadedErrorBar(xbin, ymean_hz, ysem_hz, {}, 1);
                else
                    frx = segthis(indsthis).FRsmooth_xbin_CommonTrialDur;
                    frx = frx-motifpredur;
                    
                    frmat = [segthis(indsthis).FRsmooth_rate_CommonTrialDur];
                    ymean = mean(frmat,2);
                    ysem = lt_sem(frmat');
                    shadedErrorBar(frx, ymean, ysem, ...
                        {'Color', pcolthis, 'LineStyle', linstyle},1);
                end
                
                % ==== WN
                indsthis = inds_wn;
                linstyle = '--';
                
                frx = segthis(indsthis).FRsmooth_xbin_CommonTrialDur;
                frx = frx-motifpredur;
                
                frmat = [segthis(indsthis).FRsmooth_rate_CommonTrialDur];
                ymean = mean(frmat,2);
                ysem = lt_sem(frmat');
                shadedErrorBar(frx, ymean, ysem, ...
                    {'Color', pcolthis, 'LineStyle', linstyle},1);
                axis tight;
                lt_plot_zeroline_vert;
                
                % ############################## UNIT 2
                hsplot = lt_subplot(8,4,[17:19]); hold on;
                hsplots = [hsplots hsplot];
                
                title([bregion_list{2} '-' num2str(neurpair_globID(2))]);
                segthis = seg2;
                pcolthis = pcol2;
                ylabel('BASE(solid), WN(dash)');
                
                % ==== BASELINE
                indsthis = inds_base;
                linstyle = '-';
                
                frx = segthis(indsthis).FRsmooth_xbin_CommonTrialDur;
                frx = frx-motifpredur;
                
                frmat = [segthis(indsthis).FRsmooth_rate_CommonTrialDur];
                ymean = mean(frmat,2);
                ysem = lt_sem(frmat');
                shadedErrorBar(frx, ymean, ysem, ...
                    {'Color', pcolthis, 'LineStyle', linstyle},1);
                
                % ==== WN
                indsthis = inds_wn;
                linstyle = '--';
                
                frx = segthis(indsthis).FRsmooth_xbin_CommonTrialDur;
                frx = frx-motifpredur;
                
                frmat = [segthis(indsthis).FRsmooth_rate_CommonTrialDur];
                ymean = mean(frmat,2);
                ysem = lt_sem(frmat');
                shadedErrorBar(frx, ymean, ysem, ...
                    {'Color', pcolthis, 'LineStyle', linstyle},1);
                axis tight;
                lt_plot_zeroline_vert;
                
                %% ============== OVERLAY Z-SCORED SMOOTHED FR
                hsplot = lt_subplot(8,4,[21:23]); hold on;
                hsplots = [hsplots hsplot];
                
                % =========================== 1) UNIT 1
                segthis = seg1;
                pcolthis = pcol1;
                
                % ==== BASELINE
                indsthis = inds_base;
                linstyle = '-';
                
                frmat = [segthis(indsthis).FRsmooth_rate_CommonTrialDur];
                ymean = mean(frmat,2);
                ymean = (ymean-mean(ymean(10:end-9)))./std(ymean(10:end-9));
                plot(frx, ymean, 'Color', pcolthis, 'LineStyle', linstyle, 'LineWidth', 2);
                
                % ==== WN
                indsthis = inds_wn;
                linstyle = '--';
                
                frmat = [segthis(indsthis).FRsmooth_rate_CommonTrialDur];
                ymean = mean(frmat,2);
                ymean = (ymean-mean(ymean(10:end-9)))./std(ymean(10:end-9));
                plot(frx, ymean, 'Color', pcolthis, 'LineStyle', linstyle, 'LineWidth', 2);
                
                
                % =========================== 1) UNIT 1
                segthis = seg2;
                pcolthis = pcol2;
                
                % ==== BASELINE
                indsthis = inds_base;
                linstyle = '-';
                
                frmat = [segthis(indsthis).FRsmooth_rate_CommonTrialDur];
                ymean = mean(frmat,2);
                ymean = (ymean-mean(ymean(10:end-9)))./std(ymean(10:end-9));
                plot(frx, ymean, 'Color', pcolthis, 'LineStyle', linstyle, 'LineWidth', 2);
                
                % ==== WN
                indsthis = inds_wn;
                linstyle = '--';
                
                frmat = [segthis(indsthis).FRsmooth_rate_CommonTrialDur];
                ymean = mean(frmat,2);
                ymean = (ymean-mean(ymean(10:end-9)))./std(ymean(10:end-9));
                plot(frx, ymean, 'Color', pcolthis, 'LineStyle', linstyle, 'LineWidth', 2);
                
                
                % =======================
                ylim([-3 3]);
                lt_plot_zeroline_vert;
                %% ============== PLOT XCOV GRAMS
                indthis = OUTSTRUCT_XCOV.bnum==i & OUTSTRUCT_XCOV.enum==ii ...
                    & OUTSTRUCT_XCOV.switch==iii & OUTSTRUCT_XCOV.motifnum==mm ...
                    & (all(OUTSTRUCT_XCOV.neurpair'==neurpair_globID')' | all(OUTSTRUCT_XCOV.neurpair'==fliplr(neurpair_globID)')');
                assert(sum(indthis)==1, 'if not 1, then means crucial error somewjhere...');
                
                % ==== important: determine the order of neurons
                neurpair_ordered = OUTSTRUCT_XCOV.neurpair(indthis,:);
                bregion_ordered = {SummaryStruct.birds(i).neurons(neurpair_ordered).NOTE_Location};
                
                % === pare down to relevant region
                x = PARAMS.xcenters_gram;
                y = PARAMS.Xcov_ccLags;
                indx = x>=xgramxlim(1) & x<=xgramxlim(2);
                indy = y>=xgramylim(1) & y<=xgramylim(2);
                x = x(indx);
                y = y(indy);
                

                % ==================== BASELINE
                hsplot = lt_subplot(8,4,[25:27]); hold on;
                hsplots = [hsplots hsplot];
                xcovgram = OUTSTRUCT_XCOV.XcovgramBase{indthis};
                xcovgram = xcovgram(indx, indy);
                title('BASE');
                ylabel([bregion_ordered{1} '<-->' bregion_ordered{2}]);
                lt_neural_Coher_Plot(xcovgram, x, y, ...
                    1, '', clim, '', '', '', 'East');
                lt_plot_zeroline;
                
                % ==================== WN
                hsplot = lt_subplot(8,4,[29:31]); hold on;
                hsplots = [hsplots hsplot];
                xcovgram = OUTSTRUCT_XCOV.XcovgramWN{indthis};
                xcovgram = xcovgram(indx, indy);
                title('WN');
                ylabel([bregion_ordered{1} '<-->' bregion_ordered{2}]);
                lt_neural_Coher_Plot(xcovgram, x, y, ...
                    1, '', clim, '', '', '', 'East');
                lt_plot_zeroline;
                
                
                %% ======= format plot
                linkaxes(hsplots, 'x');
                
                
                %% ############################# PLOT INDIVIDUAL TRIALS
                hsplots = [];
                
                % ##################################### UNIT 1
                segthis = seg1;
                pcolthis = pcol1;
                
                % ==== BASELINE
                hsplot = lt_subplot(8,4,4); hold on;
                hsplots = [hsplots hsplot];
                ylabel('BASE');
                indsthis = inds_base_rand;
                linstyle = '-';
                
                frmat = [segthis(indsthis).FRsmooth_rate_CommonTrialDur];
                plot(frx, frmat, 'Color', pcolthis, 'LineStyle', linstyle);
                axis tight;
                line([corrwind(1) corrwind(1)], ylim);
                line([corrwind(2) corrwind(2)], ylim);
                lt_plot_zeroline_vert;
                
                
                % ==== WN
                hsplot = lt_subplot(8,4,8); hold on;
                hsplots = [hsplots hsplot];
                ylabel('WN');
                indsthis = inds_wn_rand;
                linstyle = '-';
                
                frmat = [segthis(indsthis).FRsmooth_rate_CommonTrialDur];
                plot(frx, frmat, 'Color', pcolthis, 'LineStyle', linstyle);
                xlim(corrwind);
                axis tight;
                line([corrwind(1) corrwind(1)], ylim);
                line([corrwind(2) corrwind(2)], ylim);
                lt_plot_zeroline_vert;
                
                
                % ##################################### UNIT 1
                segthis = seg2;
                pcolthis = pcol2;
                
                % ==== BASELINE
                hsplot = lt_subplot(8,4,12); hold on;
                hsplots = [hsplots hsplot];
                ylabel('BASE');
                indsthis = inds_base_rand;
                linstyle = '-';
                
                frmat = [segthis(indsthis).FRsmooth_rate_CommonTrialDur];
                plot(frx, frmat, 'Color', pcolthis, 'LineStyle', linstyle);
                xlim(corrwind);
                axis tight;
                line([corrwind(1) corrwind(1)], ylim);
                line([corrwind(2) corrwind(2)], ylim);
                lt_plot_zeroline_vert;
                
                
                % ==== WN
                hsplot = lt_subplot(8,4,16); hold on;
                hsplots = [hsplots hsplot];
                ylabel('WN');
                indsthis = inds_wn_rand;
                linstyle = '-';
                
                frmat = [segthis(indsthis).FRsmooth_rate_CommonTrialDur];
                plot(frx, frmat, 'Color', pcolthis, 'LineStyle', linstyle);
                xlim(corrwind);
                axis tight;
                line([corrwind(1) corrwind(1)], ylim);
                line([corrwind(2) corrwind(2)], ylim);
                lt_plot_zeroline_vert;
                
                
                %% ############################## SOME STATISTICS
                %% ======== CORRELATIONS OF SMOOTHED FR (ONLY ADJACENT TRIALS)
                Yall = cell(1,4); % unit1(base, wn), 2(base,wn)
                XlabAll = {};
                % ##################################### UNIT 1
                segthis = seg1;
                
                % ================= BASELINE
                indsthis = inds_base;
                yidx = 1;
                xlab = ['N' num2str(neurpair_globID(1)) '-B'];
                
                % ----------------- COLLECT
                frmat = [segthis(indsthis).FRsmooth_rate_CommonTrialDur];
                frx = segthis(indsthis(1)).FRsmooth_xbin_CommonTrialDur - motifpredur;
                % --- only keep within time bin
                indx = frx>=corrwind(1) & frx<=corrwind(2);
                frmat = frmat(indx,:);
                % --- GET ALL PAIRWISE CORRELATIONS (n vs. n+1);
                rhoall = [];
                for kk=1:size(frmat,2)-1
                    fr1 = frmat(:,kk);
                    fr2 = frmat(:,kk+1);
                    %                     rho =
                    rhoall = [rhoall; corr(fr1, fr2)];
                end
                Yall{yidx} = rhoall;
                XlabAll{yidx} = xlab;
                
                % ================= WN
                indsthis = inds_wn;
                yidx = 2;
                xlab = ['N' num2str(neurpair_globID(1)) '-W'];
                
                % ----------------- COLLECT
                frmat = [segthis(indsthis).FRsmooth_rate_CommonTrialDur];
                frx = segthis(indsthis(1)).FRsmooth_xbin_CommonTrialDur - motifpredur;
                % --- only keep within time bin
                indx = frx>=corrwind(1) & frx<=corrwind(2);
                frmat = frmat(indx,:);
                % --- GET ALL PAIRWISE CORRELATIONS (n vs. n+1);
                rhoall = [];
                for kk=1:size(frmat,2)-1
                    fr1 = frmat(:,kk);
                    fr2 = frmat(:,kk+1);
                    %                     rho =
                    rhoall = [rhoall; corr(fr1, fr2)];
                end
                Yall{yidx} = rhoall;
                XlabAll{yidx} = xlab;
                
                % ##################################### UNIT 2
                segthis = seg2;
                
                % ================= BASELINE
                indsthis = inds_base;
                yidx = 3;
                xlab = ['N' num2str(neurpair_globID(2)) '-B'];
                
                % ----------------- COLLECT
                frmat = [segthis(indsthis).FRsmooth_rate_CommonTrialDur];
                frx = segthis(indsthis(1)).FRsmooth_xbin_CommonTrialDur - motifpredur;
                % --- only keep within time bin
                indx = frx>=corrwind(1) & frx<=corrwind(2);
                frmat = frmat(indx,:);
                % --- GET ALL PAIRWISE CORRELATIONS (n vs. n+1);
                rhoall = [];
                for kk=1:size(frmat,2)-1
                    fr1 = frmat(:,kk);
                    fr2 = frmat(:,kk+1);
                    %                     rho =
                    rhoall = [rhoall; corr(fr1, fr2)];
                end
                Yall{yidx} = rhoall;
                XlabAll{yidx} = xlab;
                
                % ================= WN
                indsthis = inds_wn;
                yidx = 4;
                xlab = ['N' num2str(neurpair_globID(2)) '-W'];
                
                % ----------------- COLLECT
                frmat = [segthis(indsthis).FRsmooth_rate_CommonTrialDur];
                frx = segthis(indsthis(1)).FRsmooth_xbin_CommonTrialDur - motifpredur;
                % --- only keep within time bin
                indx = frx>=corrwind(1) & frx<=corrwind(2);
                frmat = frmat(indx,:);
                % --- GET ALL PAIRWISE CORRELATIONS (n vs. n+1);
                rhoall = [];
                for kk=1:size(frmat,2)-1
                    fr1 = frmat(:,kk);
                    fr2 = frmat(:,kk+1);
                    %                     rho =
                    rhoall = [rhoall; corr(fr1, fr2)];
                end
                Yall{yidx} = rhoall;
                XlabAll{yidx} = xlab;
                
                % ##################################### PLOT
                lt_subplot(8,4, [20 24]); hold on;
                lt_plot_MultDist(Yall, [1:4], 0);
                set(gca, 'XTick', 1:4, 'XTickLabel', XlabAll);
                rotateXLabels(gca, 90);
                axis tight;
                lt_plot_zeroline;
                
                %% ================ LMAN-RA CORRELATION (NOT XCORR)
                % ===== ON EACH TRIAL, COMPUTE CORRELATION OF SMOOTHED FR
                
                Yall = cell(1,2); % unit1(base, wn), 2(base,wn)
                XlabAll = {};

                % =================================== BASELINE
                indsthis = inds_base;
                yidx = 1;

                % ----------------- COLLECT
                frmat1 = [seg1(indsthis).FRsmooth_rate_CommonTrialDur];
                frmat2 = [seg2(indsthis).FRsmooth_rate_CommonTrialDur];
                frx = seg1(indsthis(1)).FRsmooth_xbin_CommonTrialDur - motifpredur;
                
                % --- only keep within time bin
                indx = frx>=corrwind(1) & frx<=corrwind(2);
                frmat1 = frmat1(indx,:);
                frmat2 = frmat2(indx,:);
                
                % ======== one rho for each trial
                rho = corr(frmat1, frmat2);
                rho = diag(rho);
                
                % ==== output
                Yall{yidx} = rho;
                
                
                % =================================== WN
                indsthis = inds_wn;
                yidx = 2;

                % ----------------- COLLECT
                frmat1 = [seg1(indsthis).FRsmooth_rate_CommonTrialDur];
                frmat2 = [seg2(indsthis).FRsmooth_rate_CommonTrialDur];
                frx = seg1(indsthis(1)).FRsmooth_xbin_CommonTrialDur - motifpredur;
                
                % --- only keep within time bin
                indx = frx>=corrwind(1) & frx<=corrwind(2);
                frmat1 = frmat1(indx,:);
                frmat2 = frmat2(indx,:);
                
                % ======== one rho for each trial
                rho = corr(frmat1, frmat2);
                rho = diag(rho);
                
                % ==== output
                Yall{yidx} = rho;
                
                
                % ###################################### PLOT
                lt_subplot(8,4, [28 32]); hold on;
                ylabel('within trial rho');
                xlabel('BASE - WN');
                lt_plot_MultDist(Yall, [1 2], 1);
                axis tight;
                lt_plot_zeroline;
                
            end
            
            
            
        end
        
        
    end
end
