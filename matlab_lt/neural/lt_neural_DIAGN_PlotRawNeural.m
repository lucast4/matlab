function lt_neural_DIAGN_PlotRawNeural(SummaryStruct, BirdToPlot, NeurToPlot, motiflist, ...
    motifpredur, motifpostdur, PlotDirSong, preAndPostDurRelSameTimept, saveON, ...
    Nmax, savedirmain)

% Nmax = 20; % num trials to plot, if more, then takes random
    tstamp = lt_get_timestamp(0);

if ~exist('savedirmain', 'var')
    savedirmain = ['/bluejay5/lucas/analyses/neural/FIGS/DIAGN_PlotRawNeural/' tstamp];
elseif isempty(savedirmain)
    savedirmain = ['/bluejay5/lucas/analyses/neural/FIGS/DIAGN_PlotRawNeural/' tstamp];
end
%% === initiate savedir
if saveON==1
    mkdir(savedirmain);
%                 mkdir(savedir);
                cd(savedirmain);
    save('SummaryStruct', 'SummaryStruct');
end

%% lt 4/5/18 - plots multiple trials of raw neural activity

numbirds = length(SummaryStruct.birds);

for i=1:numbirds
    birdname = SummaryStruct.birds(i).birdname;
    
    if ~strcmp(birdname, BirdToPlot)
        continue
    end
    
    numneurons = length(SummaryStruct.birds(i).neurons);
    
    for ii=1:numneurons
        
        if ~isempty(NeurToPlot)
            if all(ii~=NeurToPlot)
                continue
            end
        end
       
        bregion = SummaryStruct.birds(i).neurons(ii).NOTE_Location;
        isSU = SummaryStruct.birds(i).neurons(ii).NOTE_is_single_unit;
        chan = SummaryStruct.birds(i).neurons(ii).channel;
        
        % ==================== PLOT
        for j = 1:length(motiflist)
            motiftoplot = motiflist{j};
            
            
            % ================ EXTRACT METADAT
            [SongDat, NeurDat, Params] = lt_neural_ExtractDat2(SummaryStruct, i, ii);
            
            if isempty(SongDat.AllLabels)
                continue
            end
            
            collectWNhit = 0;
            %             preAndPostDurRelSameTimept = 1;
            RemoveIfTooLongGapDur = 1;
            FFparams.collectFF=0;
            clustnum = SummaryStruct.birds(i).neurons(ii).clustnum;
            keepRawNeuralDat =1;
            [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
                motiftoplot, motifpredur, motifpostdur, 1, '', FFparams, ...
                0, 1, collectWNhit, 0, 0, preAndPostDurRelSameTimept, RemoveIfTooLongGapDur, ...
                clustnum, PlotDirSong, keepRawNeuralDat);
            
            if isempty(SegmentsExtract)
                continue
            end
            
            % =============================== DIR OR UNDIR?
            if PlotDirSong==0
                % then remove any dir
                SegmentsExtract([SegmentsExtract.DirSong]==1) = [];
            elseif PlotDirSong==1
                % then remove any undir
                SegmentsExtract([SegmentsExtract.DirSong]==0) = [];
            end
            
            
            %% =============== overlay neural + rasters
            figcount=1;
            subplotrows=10;
            subplotcols=2;
            fignums_alreadyused=[];
            hfigs=[];
            hsplots = [];
            
            
            if length(SegmentsExtract) > Nmax
                disp('getting random subset of trials');
                indstoplot = sort(randperm(length(SegmentsExtract), Nmax));
            else
                indstoplot = 1:length(SegmentsExtract);
            end
            
            for ind = indstoplot
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title([birdname ',n' num2str(ii) ',' motiftoplot ',tr' num2str(ind)]);
                ylabel([bregion '-SU' num2str(isSU) '-ch' num2str(chan)]);
                hsplots = [hsplots hsplot];
                datneur = SegmentsExtract(ind).neurdatseg_filt;
                tneur = (1:length(datneur))/NeurDat.metaDat(1).fs;
                spktimes = SegmentsExtract(ind).spk_Times;
                
                plot(tneur, datneur, '-k');
                lt_neural_PLOT_rasterline(spktimes, 0, 'r');
%                 plot(spktimes, 0, 'or');
                line([motifpredur motifpredur], ylim, 'Color', 'b');
                
                if preAndPostDurRelSameTimept==0
                   % -- put line for all syls in motifs
                   for k=1:length(SegmentsExtract(ind).motifsylOnsets)
                       ons = SegmentsExtract(ind).motifsylOnsets(k);
                       off = SegmentsExtract(ind).motifsylOffsets(k);
                       line([ons off], [200 200], 'Color', 'm', 'LineWidth', 2);
                       line([ons ons], [-200 200], 'Color', 'm', 'LineStyle', '--');
                       line([off off], [-200 200], 'Color', 'm', 'LineStyle', '--');
                   end
                end
            end
            
            linkaxes(hsplots, 'xy');
            if preAndPostDurRelSameTimept == 1
                xlim([0 motifpredur+motifpostdur]);
            else
                
                axis tight
            end
            
            
            %%  SAVE
            
            if saveON ==1
%                 savedir = '
                savedir = [savedirmain '/' birdname '-neur' num2str(ii) '-' motiftoplot];
                mkdir(savedir);
                cd(savedir);
                disp(savedir);
                lt_save_all_figs;
                close all; fclose all;
            end
            
            
            %% 
            %             % ------------- 2) PLOT SMOOTHED FR
            %             SegmentsExtract = lt_neural_SmoothFR(SegmentsExtract, '');
            %             FRmat = [SegmentsExtract.FRsmooth_rate_CommonTrialDur];
            %             X = SegmentsExtract(1).FRsmooth_xbin_CommonTrialDur;
            %             % 1) each trial
            %             if plotIndivRaster==1
            %                 hsplot = lt_subplot(6,1,5); hold on;
            %                 hsplots = [hsplots hsplot];
            %                 for tt =1:length(SegmentsExtract)
            %                     plot(X, FRmat(:,tt), 'Color', [0.7 0.7 0.7]);
            %                 end
            %                 line([motifpredur motifpredur], ylim);
            %             end
            %
            %             % 2) mean
            %             if length(SegmentsExtract)>2
            %                 FRmean = mean(FRmat,2);
            %                 FRsem = lt_sem(FRmat');
            %
            %                 if plotIndivRaster==1
            %                     hsplot = lt_subplot(6,1,6); hold on;
            %                     hsplots = [hsplots hsplot];
            %                     shadedErrorBar(X, FRmean, FRsem, {'Color','r'},1);
            %                     line([motifpredur motifpredur], ylim);
            %                 end
            %
            %                 % ############################### COLLECT SMOOTHED FR
            %                 AllSmoothFR = [AllSmoothFR FRmean];
            %                 AllSmoothFR_sem = [AllSmoothFR_sem FRsem];
            %                 AllSmoothFR_x = [AllSmoothFR_x X];
            %             end
            %
            %
            %             % ------------- 3) PLOT syl onset/offsets
            %             if plotIndivRaster==1
            %                 hsplot = lt_subplot(6,1,1); hold on;
            %                 %         hsplots = [hsplots hsplot];
            %                 numsyltrialstoplot = 15;
            %                 maxdur = motifpostdur+motifpredur -0.005;
            %                 [SylContours, x] = lt_neural_v2_ANALY_GetSylContours(SegmentsExtract, ...
            %                     numsyltrialstoplot, maxdur);
            %                 spy(SylContours);
            %                 %         plot(x, SylContours, '+');
            %                 set(gca, 'XTick', 1:50:size(SylContours,2), 'XTickLabel', 100*x(1:50:end))
            %                 line([motifpredur motifpredur], ylim);
            %
            %                 % ---
            %                 linkaxes(hsplots, 'x');
            %             end
            %
            %             % ###################### collect stuff
            %             AllMotif = [AllMotif motiftoplot];
            %             AllNeurNum = [AllNeurNum ii];
        end
    end
end