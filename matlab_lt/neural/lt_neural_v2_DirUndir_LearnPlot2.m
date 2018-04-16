function lt_neural_v2_DirUndir_LearnPlot2(MOTIFSTATS_Compiled, SwitchStruct, BirdExptPairsToPlot, ...
    plotIndTrial, SummaryStruct, motiftoplot, xwindplot, onlyPlotIfBothPrePostTrials)

%% given syl, plots learning trajectory + UDNIR and DIR mean FR overlayed


%%


numbirds = length(SwitchStruct.bird);
for i=1:numbirds
    
    numexpts = length(SwitchStruct.bird(i).exptnum);
    birdname = SwitchStruct.bird(i).birdname;
    
    for ii=1:numexpts
        
        exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
        
        % ------------------------------- PLOT?
        if ~isempty(BirdExptPairsToPlot)
            
            ind1 = find(strcmp(BirdExptPairsToPlot, birdname));
            ind2 = find(strcmp(BirdExptPairsToPlot, exptname));
            
            if ~any(ind1+1 == ind2)
                disp(['SKIPPED ' birdname '-' exptname]);
                continue
            end
        end
        
%         dayfirst = datestr(floor(SwitchStruct.bird(i).exptnum(ii).switchlist(1).switchdnum), 'ddmmmyyyy');
        numswitches = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        neuronsThisExpt = find(strcmp({SummaryStruct.birds(i).neurons.exptID}, exptname));
        
        % ======= get list of motifs for this expt
        % -- first find neurons for this expt
        motiflist = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(neuronsThisExpt(1)).motif_regexpr_str;
        % -- check all neurons make sure motiflist identical
        for nn=neuronsThisExpt
           
            motiflistthis = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(nn).motif_regexpr_str;
            
            assert(all(strcmp(motiflist, motiflistthis)), 'sadfasd');
        end
        
        % -- make sure they are identical
%         motiflist = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.params.motif_regexpr_str;
        motifpredur = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.params.motif_predur;
        xwindplot_fromzero = motifpredur+xwindplot;
        
        
        % ===============================================
        % -- for each motif plot each neuron for each switch.
        for mm=1:length(motiflist)
            if ~any(strcmp(motiflist{mm}, motiftoplot))
                continue
            end
            
            %% =================== FIRST PLOT LEARNING TRAJECTORY
            lt_figure; hold on;
            batchesdone = {};
            for jj=1:length(neuronsThisExpt)
                nn = neuronsThisExpt(jj);
                
                batchthis = SummaryStruct.birds(i).neurons(nn).batchfilename;
                if any(strcmp(batchesdone, batchthis))
                    continue
                end
                batchesdone = [batchesdone batchthis];
%                 bregionthis = MOTIFSTATS_Compiled.birds(i).SummaryStruct.birds(1).neurons(nn).NOTE_Location;
%                 
%                 % ===== does this neuron have pre or post data for this
%                 % switch?
%                 songtimes = MOTIFSTATS_Compiled.birds(i).SummaryStruct.birds(1).neurons(nn).Filedatenum_unsorted;
%                 %                 disp(songtimes)
%                 inds_pre = find(songtimes>swthis.switchdnum_previous & songtimes<swthis.switchdnum);
%                 inds_post = find(songtimes>swthis.switchdnum & songtimes<swthis.switchdnum_next);
%                 
%                 if isempty(inds_pre) & isempty(inds_post)
%                     continue
%                 end
%                 
%                 if onlyPlotIfBothPrePostTrials ==1
%                     if isempty(inds_pre) | isempty(inds_post)
%                         continue
%                     end
%                 end
%                 disp(['analyzing: ' birdname '-' exptname '-sw' num2str(iii) '-neurset' num2str(ss)]);
%                 
%                 
%                 
                % ======================
                segextract = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(nn).motif(mm).SegmentsExtract;
                
                if isempty(segextract)
                    continue
                end
                
                % --- get pre and post inds
                tvals = [segextract.song_datenum];
                ffvals = [segextract.FF_val];
                indsdir = [segextract.DirSong];
                
                % ---- UNDIR
                plot(tvals(indsdir==0), ffvals(indsdir==0), 'ok');
                
                % ---- DIR
                plot(tvals(indsdir==1), ffvals(indsdir==1), 'ob', 'MarkerFaceColor', 'b');
                
            end
            
                % =========== PUT LINES FOR SWITCHES
                for iii=1:numswitches
                    swthis = SwitchStruct.bird(i).exptnum(ii).switchlist(iii);
                    
                    line([swthis.switchdnum swthis.switchdnum], ylim);
%                     earliest = max([swthis.switchdnum_previous, floor(min(tvals))]);
%                     latest = min([swthis.switchdnum_next, ceil(min(tvals))]);
%                     
%                     line([earliest earliest], ylim, 'LineStyle', '--');
%                     line([latest latest], ylim, 'LineStyle', '--');
                    tmp = swthis.learningContingencies;
                    for jjj=1:length(tmp)/2
                       lt_plot_text(swthis.switchdnum, max(ffvals), [num2str(tmp{2*jjj-1}) ':' num2str(tmp{2*jjj})], 'r'); 
                       lt_plot_text(swthis.switchdnum, 1.01*max(ffvals), ['sw' num2str(iii)], 'r'); 
                    end
                    axis tight;
                    datetick('x', 'ddmmm-HHMM');
%                     rotateXLabels(gca, 45)
                end
            
            
            %% ========= SECOND, plot neural 
            figcount=1;
            subplotrows=4;
            subplotcols=8;
            fignums_alreadyused=[];
            hfigs=[];
            hsplots =[];
            figureopened = 1;
            
            motifstr = motiflist{mm};
            
            % =================== GO THRU EACH SWITCH. FOR EACH SWITCH PLOT
            % ALL THE NEURONS THAT HAVE DATA EITHER IN PRE OR POST PERIOD
            % FOR THIS SWITCH
            
            for iii=1:numswitches
                swthis = SwitchStruct.bird(i).exptnum(ii).switchlist(iii);
                
                titleplotted = 0; % plot once to indicate new switch
                for nn = neuronsThisExpt
                    bregionthis = MOTIFSTATS_Compiled.birds(i).SummaryStruct.birds(1).neurons(nn).NOTE_Location;
                    
                    % ===== does this neuron have pre or post data for this
                    % switch?
                    songtimes = MOTIFSTATS_Compiled.birds(i).SummaryStruct.birds(1).neurons(nn).Filedatenum_unsorted;
                    %                 disp(songtimes)
                    inds_pre = find(songtimes>swthis.switchdnum_previous & songtimes<swthis.switchdnum);
                    inds_post = find(songtimes>swthis.switchdnum & songtimes<swthis.switchdnum_next);
                    
                    if isempty(inds_pre) & isempty(inds_post)
                        continue
                    end
                    
                    if onlyPlotIfBothPrePostTrials ==1
                        if isempty(inds_pre) | isempty(inds_post)
                            continue
                        end
                    end
                    disp(['analyzing: ' birdname '-' exptname '-sw' num2str(iii)]);
                    
                    % ======================
                    segextract = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(nn).motif(mm).SegmentsExtract;
                    
                    
                    if isempty(segextract)
                        continue
                    end
                    
                    % --- get pre and post inds
                    tvals = [segextract.song_datenum];
                    
                    indspre = tvals>swthis.switchdnum_previous ...
                        & tvals < swthis.switchdnum;
                    indspost = tvals>swthis.switchdnum ...
                        & tvals<swthis.switchdnum_next;

                    % ---- get smoothed FR
                    segextract = lt_neural_SmoothFR(segextract);
                    
                    %% === INITIATE COMBINED SUBPLOT
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    hsplots = [hsplots hsplot];
                    title(motifstr);
                    xlabel({['PRE(k=U;b=D)'], ['POST(r=U;m=D)']});
                        
                        title([motifstr '-neur' num2str(nn) '-' bregionthis]);
                        if titleplotted==0
                            ylabel([birdname '-' exptname '-sw' num2str(iii)]);
                            titleplotted=1;
                        end
                    
                    %% ##################### PLOT PRE
                    indsEpoch = indspre;
                    
                    if any(indsEpoch)
                        % -------------- UNDIR
                        inds = [segextract.DirSong]==0 & indsEpoch;
                        plottitle = [motifstr '-UNDIR'];
                        
                        frmat = [segextract(inds).FRsmooth_rate_CommonTrialDur];
                        x = segextract(1).FRsmooth_xbin_CommonTrialDur;
                        y = mean(frmat,2);
                        if plotIndTrial==1
                            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                            hsplots = [hsplots hsplot];
                            title(plottitle);
                            plot(x, frmat, '-', 'Color', [0.6 0.6 0.6]);
                            plot(x, y, '-k', 'LineWidth', 3);
                            line([motifpredur motifpredur], ylim, 'Color','r');
                        end
                        
                        x1 = x;
                        y1 = y;
                        y1sem = lt_sem(frmat');
                        
                        
                        % -------------- DIR
                        inds = [segextract.DirSong]==1 & indsEpoch;
                        plottitle = [motifstr '-DIR'];

                        frmat = [segextract(inds).FRsmooth_rate_CommonTrialDur];
                        x = segextract(1).FRsmooth_xbin_CommonTrialDur;
                        y = mean(frmat,2);
                        if plotIndTrial==1
                            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                            hsplots = [hsplots hsplot];
                            title(plottitle);
                            plot(x, frmat, '-', 'Color', [0.6 0.6 0.6]);
                            plot(x, y, '-k', 'LineWidth', 3);
                            line([motifpredur motifpredur], ylim, 'Color','r');
                        end
                        x2 = x;
                        y2 = y;
                        y2sem = lt_sem(frmat');
                        
%                         if iii==1 & any(inds)
%                             keyboard
%                         end
                        
                        % -------------- COMBINED
                        
                        shadedErrorBar(x1, y1, y1sem, {'Color', 'k'}, 1);
                        if size(frmat,2)>1
                            shadedErrorBar(x2, y2, y2sem, {'Color', 'b'}, 1);
                        end
                        line([motifpredur motifpredur], ylim, 'Color','r');
                    end
                    
                  %% PLOT POST (overlay)
                  indsEpoch = indspost;
                  
                  if any(indsEpoch)
                      
                      % -------------- UNDIR
                      inds = [segextract.DirSong]==0 & indsEpoch;
                      plottitle = [motifstr '-UNDIR'];
                      
                      frmat = [segextract(inds).FRsmooth_rate_CommonTrialDur];
                      x = segextract(1).FRsmooth_xbin_CommonTrialDur;
                      y = mean(frmat,2);
                      if plotIndTrial==1
                          [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                          hsplots = [hsplots hsplot];
                          title(plottitle);
                          plot(x, frmat, '-', 'Color', [0.6 0.6 0.6]);
                          plot(x, y, '-k', 'LineWidth', 3);
                          line([motifpredur motifpredur], ylim, 'Color','r');
                      end
                      
                      x1 = x;
                      y1 = y;
                      y1sem = lt_sem(frmat');
                      
                      
                      % -------------- DIR
                      inds = [segextract.DirSong]==1 & indsEpoch;
                      plottitle = [motifstr '-DIR'];
                      
                      
                      frmat = [segextract(inds).FRsmooth_rate_CommonTrialDur];
                      x = segextract(1).FRsmooth_xbin_CommonTrialDur;
                      y = mean(frmat,2);
                      if plotIndTrial==1
                          [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                          hsplots = [hsplots hsplot];
                          title(plottitle);
                          plot(x, frmat, '-', 'Color', [0.6 0.6 0.6]);
                          plot(x, y, '-k', 'LineWidth', 3);
                          line([motifpredur motifpredur], ylim, 'Color','r');
                      end
                      x2 = x;
                      y2 = y;
                      y2sem = lt_sem(frmat');
                      
                      % -------------- COMBINED
                      %                     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                      %                     hsplots = [hsplots hsplot];
                      %                     title(motifstr);
                      %                     xlabel(EpochLab);
                      %
                      %                     title([motifstr '-neur' num2str(nn) '-' bregionthis]);
                      %                     if titleplotted==0
                      %                         ylabel([birdname '-' exptname '-sw' num2str(iii)]);
                      %                         titleplotted=1;
                      %                     end
                      
                      shadedErrorBar(x1, y1, y1sem, {'Color', 'r'}, 1);
                      if size(frmat,2)>1
                          shadedErrorBar(x2, y2, y2sem, {'Color', 'm'}, 1);
                      end
                      %                     line([motifpredur motifpredur], ylim, 'Color','r');
                  end
                     
                end
                
                
            end
            
              if ~isempty(hsplots)
                    linkaxes(hsplots, 'xy');
                    % --- zoom in on x axis
                    xlim([xwindplot_fromzero]);
                end

            
            
        end
    end
end
