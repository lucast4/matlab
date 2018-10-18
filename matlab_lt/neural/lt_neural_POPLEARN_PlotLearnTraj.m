function lt_neural_POPLEARN_PlotLearnTraj(MOTIFSTATS_pop, ...
    SwitchStruct, SummaryStruct, BirdExptPairsToPlot, motiftoplot)


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
        
            
        dayfirst = datestr(floor(SwitchStruct.bird(i).exptnum(ii).switchlist(1).switchdnum), 'ddmmmyyyy');
        numswitches = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        
        % ===============================================
        % -- for each motif plot each neuron for each switch.
        
        %% ====================== PLOT EACH NEURON SET
        lt_figure; hold on;
        DAT = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT;
        numsets = length(MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons);
        plotcols = lt_make_plot_colors(numsets, 0,0);
        
        hsplot1 = lt_subplot(3,1,1:2); hold on;
        title([birdname '-' exptname '-' motiftoplot])
        hsplot2 = lt_subplot(3,1,3); hold on;
        title('neuron sets');
        
        for ss =1:numsets
           
            if isempty(motiftoplot)
                    motifnum = find(~isempty({DAT.setnum(ss).motif.SegExtr_neurfakeID}), 1, 'first');
            else                
            motifnum = find(strcmp({DAT.setnum(ss).motif.regexpstr}, ...
                motiftoplot));
            end
            
            assert(~isempty(motifnum), 'this motif doesnt exist?');
            
            ff = [DAT.setnum(ss).motif(motifnum).SegExtr_neurfakeID(1).SegmentsExtract.FF_val];
            t = [DAT.setnum(ss).motif(motifnum).SegExtr_neurfakeID(1).SegmentsExtract.song_datenum];
            isdir = [DAT.setnum(ss).motif(motifnum).SegExtr_neurfakeID(1).SegmentsExtract.DirSong];
            
            % === print list of motifs
            disp(['SET' num2str(ss) ': ' {DAT.setnum(ss).motif.regexpstr}])
            
            % --- convert t to days
            %            tmp = lt_convert_EventTimes_to_RelTimes(dayfirst, t);
            %            t = tmp.FinalValue;
            
            %  ############################## PLOT FF
            lt_subplot(3,1,1:2); hold on;
            % ============ PLOT UNDIR
            plot(t(isdir==0), ff(isdir==0), 'o', 'Color', plotcols{ss});
            
            % ============ PLOT DUR
            plot(t(isdir==1), ff(isdir==1), 'o', 'Color', plotcols{ss}, ...
                'MarkerFaceColor', 'k')
            
            % ##################### indicate extent of files for this set
            lt_subplot(3,1,3); hold on;
            line([min(t) max(t)], [ss ss], 'Color', plotcols{ss}, 'LineWidth', 2);
            ylim([0 ss+1]);
            
            % ---- indicate what units are in this set
            neurthis = MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons{ss};
            locthis = [SummaryStruct.birds(i).neurons(neurthis).NOTE_Location];
            lt_plot_text(min(t), ss+0.2,  ['nn: ' num2str(neurthis) ',' locthis], plotcols{ss});
        end
        %
        datetick('x', 'ddmmm-HHMM');
        linkaxes([hsplot1 hsplot2], 'x');
        
        % ================= OVERLAY SWITCH TIMES AND CONTINGENCIES
        lt_subplot(3,1,1:2);
        for iii=1:numswitches
            swthis = SwitchStruct.bird(i).exptnum(ii).switchlist(iii);
            
            line([swthis.switchdnum swthis.switchdnum], ylim);
            %                     earliest = max([swthis.switchdnum_previous, floor(min(tvals))]);
            %                     latest = min([swthis.switchdnum_next, ceil(min(tvals))]);
            %
            %                     line([earliest earliest], ylim, 'LineStyle', '--');
            %                     line([latest latest], ylim, 'LineStyle', '--');
            tmp = swthis.learningContingencies;
            lt_plot_text(swthis.switchdnum, 1.01*max(ff), ['sw' num2str(iii)], 'k');
            for jjj=1:length(tmp)/2
                lt_plot_text(swthis.switchdnum, (1.01+jjj*0.01)*max(ff), [num2str(tmp{2*jjj-1}) ':' num2str(tmp{2*jjj})], 'k');
            end
            axis tight;
            datetick('x', 'ddmmm-HHMM');
            %                     rotateXLabels(gca, 45)
        end
        
        
    end
end
