function lt_neural_v2_DirUndir_LearnSummary1(MOTIFSTATS_Compiled, SwitchStruct, BirdExptPairsToPlot, ...
    SummaryStruct, xwindplot)
Nmin = 3; % min num dir songs, really only applies for DIR. checks on a case by case (neuron x motif) basis.

%% 3/21/18 - lt, summarize effect of DIR during learning expts

numbirds = length(SwitchStruct.bird);
DATSTRUCT = struct;
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
        
        neuronsThisExpt = find(strcmp({SummaryStruct.birds(i).neurons.exptID}, exptname));
        
        % ======= get list of motifs for this expt
        % -- first find neurons for this expt
        motiflist = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(neuronsThisExpt(1)).motif_regexpr_str;
        % -- check all neurons make sure motiflist identical
        for nn=neuronsThisExpt
            
            motiflistthis = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(nn).motif_regexpr_str;
            disp(length(motiflistthis));
            assert(all(strcmp(motiflist, motiflistthis)), 'sadfasd');
        end
        
        %         motiflist = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.params.motif_regexpr_str;
        motifpredur = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.params.motif_predur;
        xwindplot_fromzero = motifpredur+xwindplot;
        
        
        % ===============================================
        % -- for each motif plot each neuron for each switch.
        DatCell = cell(length(motiflist), 2, 2, numswitches); % motifnum x pre/post x UNDIR/DIR x switchnum
        DatCell_neuronID = cell(length(motiflist), 2, 2, numswitches);
        DatCell_xtimes = cell(length(motiflist), 2, 2, numswitches);
        DatCell_bregion = cell(length(motiflist), 2, 2, numswitches);
        
        for mm=1:length(motiflist)
            
            motifstr = motiflist{mm};
            
            % =================== GO THRU EACH SWITCH. FOR EACH SWITCH PLOT
            % ALL THE NEURONS THAT HAVE DATA EITHER IN PRE OR POST PERIOD
            % FOR THIS SWITCH
            for iii=1:numswitches
                swthis = SwitchStruct.bird(i).exptnum(ii).switchlist(iii);
                
                titleplotted = 0;
                
                for jj = 1:length(neuronsThisExpt)
                    nn = neuronsThisExpt(jj);
                    bregionthis = SummaryStruct.birds(i).neurons(nn).NOTE_Location;
                    
                    % ===== does this neuron have pre or post data for this
                    % switch?
                    songtimes = MOTIFSTATS_Compiled.birds(i).SummaryStruct.birds(1).neurons(nn).Filedatenum_unsorted;
                    %                 disp(songtimes)
                    inds_pre = find(songtimes>swthis.switchdnum_previous & songtimes<swthis.switchdnum);
                    inds_post = find(songtimes>swthis.switchdnum & songtimes<swthis.switchdnum_next);
                    
                    if isempty(inds_pre) & isempty(inds_post)
                        continue
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
                    xtmp = segextract(1).FRsmooth_xbin_CommonTrialDur;
                    xinds = xtmp >= xwindplot_fromzero(1) ...
                        & xtmp <=xwindplot_fromzero(2); % common trial duration
                    
                    %% ##################### COLLECT PRE
                    indsEpoch = indspre;
                    
                    if ~any(indsEpoch)
                        continue
                    end
                    
                    % =========================== UNDIR
                    inds = [segextract.DirSong]==0 & indsEpoch;
                    %                     plottitle = [motifstr '-UNDIR'];
                    
                    frmat = [segextract(inds).FRsmooth_rate_CommonTrialDur];
                    x = segextract(1).FRsmooth_xbin_CommonTrialDur;
                    y = mean(frmat,2);
                    % --- shorten to common duration
                    y = y(xinds);
                    x = x(xinds);
                    
                    % ---- put into output
                    DatCell{mm, 1, 1, iii} = [DatCell{mm, 1, 1, iii}; y'];
                    DatCell_xtimes{mm, 1, 1, iii} = x;
                    DatCell_neuronID{mm, 1, 1, iii} = [DatCell_neuronID{mm, 1, 1, iii}; nn];
                    
                    % ============================== DIR
                    inds = [segextract.DirSong]==1 & indsEpoch;
                    ndirsongs = length(unique({segextract([segextract.DirSong]==1).song_filename}));
                    
                    if ndirsongs>=Nmin & any(inds)
                        frmat = [segextract(inds).FRsmooth_rate_CommonTrialDur];
                        x = segextract(1).FRsmooth_xbin_CommonTrialDur;
                        y = mean(frmat,2);
                        % --- shorten to common duration
                        y = y(xinds);
                        x = x(xinds);
                        
                        % ---- put into output
                        DatCell{mm, 1, 2, iii} = [DatCell{mm, 1, 2, iii}; y'];
                        DatCell_xtimes{mm, 1, 2, iii} = x;
                        DatCell_neuronID{mm, 1, 2, iii} = [DatCell_neuronID{mm, 1, 2, iii}; nn];
                    end
                    
                    %% PLOT POST (overlay)
                    indsEpoch = indspost;
                    
                    if ~any(indsEpoch)
                        continue
                    end
                    % ========================== UNDIR
                    inds = [segextract.DirSong]==0 & indsEpoch;
                    
                    frmat = [segextract(inds).FRsmooth_rate_CommonTrialDur];
                    x = segextract(1).FRsmooth_xbin_CommonTrialDur;
                    y = mean(frmat,2);
                    % --- shorten to common duration
                    y = y(xinds);
                    x = x(xinds);
                    
                    % ---- put into output
                    DatCell{mm, 2, 1, iii} = [ DatCell{mm, 2, 1, iii} ; y'];
                    DatCell_xtimes{mm, 2,1, iii} = x;
                    DatCell_neuronID{mm, 2, 1, iii} = [DatCell_neuronID{mm, 2, 1, iii}; nn];
                    
                    % ========================== DIR
                    inds = [segextract.DirSong]==1 & indsEpoch;
                    ndirsongs = length(unique({segextract([segextract.DirSong]==1).song_filename}));
                    
                    if ndirsongs>=Nmin & any(inds)
                        frmat = [segextract(inds).FRsmooth_rate_CommonTrialDur];
                        x = segextract(1).FRsmooth_xbin_CommonTrialDur;
                        y = mean(frmat,2);
                        % --- shorten to common duration
                        y = y(xinds);
                        x = x(xinds);
                        
                        % ---- put into output
                        DatCell{mm, 2, 2, iii} = [DatCell{mm, 2, 2, iii}; y'];
                        DatCell_xtimes{mm, 2,2, iii} = x;
                        DatCell_neuronID{mm, 2, 2, iii} = [DatCell_neuronID{mm, 2, 2, iii}; nn];
                    end
                    
                    % ============== COLLECT XTIMES
                end
            end
            
        end
        
        % ---- put cell for this expt into output struct
        DATSTRUCT.bird(i).expt(ii).DatCell = DatCell;
        DATSTRUCT.bird(i).expt(ii).DatCell_dimensions = 'motifnum x pre/post x UNDIR/DIR x switchnum';
        DATSTRUCT.bird(i).expt(ii).DatCell_xtimes = DatCell_xtimes;
        DATSTRUCT.bird(i).expt(ii).DatCell_neuronID = DatCell_neuronID;
        
    end
end

%% ======= plot
numbirds = length(DATSTRUCT.bird);

% plotstat = 'corr';
plotstat = 'prctdiff';

for i=1:numbirds
    birdname = SwitchStruct.bird(i).birdname;
    
    numexpts = length(DATSTRUCT.bird(i).expt);
    for ii=1:numexpts
        exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
        
        % ============= get motiflist
        motiflist = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(neuronsThisExpt(1)).motif_regexpr_str;
        % -- check all neurons make sure motiflist identical
        for nn=neuronsThisExpt
            
            motiflistthis = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(nn).motif_regexpr_str;
            disp(length(motiflistthis));
            assert(all(strcmp(motiflist, motiflistthis)), 'sadfasd');
        end
        
        % =============
        DatCell = DATSTRUCT.bird(i).expt(ii).DatCell;
        DatCell_dimensions = DATSTRUCT.bird(i).expt(ii).DatCell_dimensions;
        DatCell_xtimes= DATSTRUCT.bird(i).expt(ii).DatCell_xtimes;
        DatCell_neuronID= DATSTRUCT.bird(i).expt(ii).DatCell_neuronID;
        
        numswitches = size(DatCell, 4);
        
        % ------- prepare figures
        figcount=1;
        subplotrows=6;
        subplotcols=1;
        fignums_alreadyused=[];
        hfigs=[];
        hsplots = [];
        
        for ss =1:numswitches
            
            % ===== FOR EACH NEURON AND MOTIF, GET CORRELATION BETWEEN DIR
            % AND UNDIR
            
            
            % ######################## PRE
            preind = 1;
            dat = squeeze(DatCell(:, preind, :, ss));
            dat_neurID = squeeze(DatCell_neuronID(:, preind, :, ss));
            
            nummotifs = size(dat,1);
            assert(length(motiflist) == nummotifs, 'asdafds');
            
            CorrCell = cell(nummotifs, 2); % motif x pre/post
            DiffPrctCell = cell(nummotifs, 2);
            BrainLocationCell = cell(nummotifs, 2);
            
            % --- go thru each motif
            for mm=1:nummotifs
                
                if isempty(dat{mm, 1}) | isempty(dat{mm,2})
                    % then either dir or undir is empty, skip...
                    continue
                end
                
                y1 = dat{mm,1}; % neuron x timbin
                y2 = dat{mm,2};
                
                % ================== make sure same neurons are matched up
                y1_neurID = dat_neurID{mm, 1};
                y2_neurID = dat_neurID{mm, 2};
                
                % ---- find overlapping neurons
                [~, indneur1, indneur2] = intersect(y1_neurID, y2_neurID);
                y1 = y1(indneur1, :);
                y2 = y2(indneur2, :);
                
                % --- get braing regions for these neurons
                BrainLocationCell{mm, preind} = ...
                    {SummaryStruct.birds(i).neurons(y1_neurID(indneur1)).NOTE_Location};
                
                % ================= 1) correlation
                rho = corr(y1', y2');
                rho = diag(rho); % these are corr within same neuron ...
                
                CorrCell{mm, preind} = rho;
                
                % ================= 2) mean % change in FR
                tmp = 100*(y2-y1)./y1;
                DiffPrctCell{mm, preind} = mean(tmp,2);
            end
            
            % ######################## POST
            preind = 2;
            dat = squeeze(DatCell(:, preind, :, ss));
            dat_neurID = squeeze(DatCell_neuronID(:, preind, :, ss));
            
            nummotifs = size(dat,1);
            assert(length(motiflist) == nummotifs, 'asdafds');
            
            % --- go thru each motif
            for mm=1:nummotifs
                
                if isempty(dat{mm, 1}) | isempty(dat{mm,2})
                    % then either dir or undir is empty, skip...
                    continue
                end
                
                y1 = dat{mm,1}; % neuron x timbin
                y2 = dat{mm,2};
                
                % ================== make sure same neurons are matched up
                y1_neurID = dat_neurID{mm, 1};
                y2_neurID = dat_neurID{mm, 2};
                
                % ---- find overlapping neurons
                [~, indneur1, indneur2] = intersect(y1_neurID, y2_neurID);
                y1 = y1(indneur1, :);
                y2 = y2(indneur2, :);
                
                % --- get braing regions for these neurons
                BrainLocationCell{mm, preind} = ...
                    {SummaryStruct.birds(i).neurons(y1_neurID(indneur1)).NOTE_Location};
                
                % ================= 1) correlation
                rho = corr(y1', y2');
                rho = diag(rho); % these are corr within same neuron ...
                
                CorrCell{mm, preind} = rho;
                
                % ================= 2) mean % change in FR
                tmp = 100*(y2-y1)./y1;
                DiffPrctCell{mm, preind} = mean(tmp,2);
            end
            
            % ################################## PLOT FOR THIS SWITCH
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title([birdname ',' exptname ',sw' num2str(ss)]);
            xlabel('motif');
            ylabel('corr (dir vs. undir)');
            
            if strcmp(plotstat, 'prctdiff')
                Stattoplot = DiffPrctCell;
            elseif strcmp(plotstat, 'corr')
                Stattoplot = CorrCell;
                
            end
            for mm=1:nummotifs
                
                % --------------- pre
                y1 = Stattoplot{mm, 1};
                breglist = BrainLocationCell{mm, 1};
                if ~isempty(y1)
                    % -- LMAN
                    indtmp = strcmp(breglist, 'LMAN');
                    plot(mm-0.2, y1(indtmp), 'og');
                    % -- RA
                    indtmp = strcmp(breglist, 'RA');
                    plot(mm-0.2, y1(indtmp), 'or')
                    % -- mean
                    lt_plot(mm-0.1, mean(y1), {'Errors', lt_sem(y1), 'Color', 'k'});
                end
                
                % --------------- post
                y2 = Stattoplot{mm, 2};
                breglist = BrainLocationCell{mm, 2};
                if ~isempty(y2)
                    % -- LMAN
                    indtmp = strcmp(breglist, 'LMAN');
                    plot(mm+0.2, y1(indtmp), 'og');
                    % -- RA
                    indtmp = strcmp(breglist, 'RA');
                    plot(mm+0.2, y1(indtmp), 'or')
                    % -- mean
                    lt_plot(mm+0.3, mean(y2), {'Errors', lt_sem(y2), 'Color', 'k'});
                end
                
                % ------------- line connecting
                if ~isempty(y1) & ~isempty(y2)
                    plot([mm-0.2 mm+0.2], [y1'; y2], '-');
                end
                
            end
            
            % ====================== PLOT LINE THRU ALL MOTIFS
            % ------ PRE
            preind = 1;
            
            ymat = [Stattoplot{:,preind}];
            if ~isempty(ymat)
                x = find(~cellfun(@isempty, Stattoplot(:,preind)));
                plot(x-0.2, ymat, '-', 'Color', [0.7 0.7 0.7]);
            end
            
            % ------ POST
            preind = 2;
            
            ymat = [Stattoplot{:,preind}];
            if ~isempty(ymat)
                x = find(~cellfun(@isempty, Stattoplot(:,preind)));
                plot(x-0.2, ymat, '-', 'Color', [0.7 0.7 0.7]);
            end
            
            % ----------------------
            lt_plot_zeroline;
            set(gca, 'XTickLabel', motiflist);
            
            if (0)
                % to plot raw fr traces [currently only written for pre...]
                
                % ==== plot mean FR
                lt_figure; hold on;
                
                % ==== pre (undir)
                preind = 1;
                dirind = 1;
                plotcol = 'k';
                
                dat = DatCell(:, preind, dirind, ss);
                datx = DatCell_xtimes(:, preind, dirind, ss);
                
                nummotifs = size(dat,1);
                count = 1;
                assert(length(motiflist) == nummotifs, 'asdafds');
                for mm=1:nummotifs
                    
                    y = dat{mm}; % each row one neuron
                    x = datx{mm};
                    
                    numneurons = size(y,1);
                    
                    for nn=1:numneurons
                        
                        lt_subplot(nummotifs, numneurons, count); hold on;
                        count = count+1;
                        title(['mot' num2str(mm) '-n' num2str(nn)]);
                        plot(x, y(nn,:), '-', 'Color', plotcol);
                        
                    end
                end
                
                % ==== pre (dir)
                preind = 1;
                dirind = 2;
                plotcol = 'b';
                
                dat = DatCell(:, preind, dirind, ss);
                datx = DatCell_xtimes(:, preind, dirind, ss);
                
                nummotifs = size(dat,1);
                count = 1;
                for mm=1:nummotifs
                    
                    y = dat{mm}; % each row one neuron
                    x = datx{mm};
                    
                    numneurons = size(y,1);
                    
                    for nn=1:numneurons
                        
                        lt_subplot(nummotifs, numneurons, count); hold on;
                        count = count+1;
                        title(['mot' num2str(mm) '-n' num2str(nn)]);
                        plot(x, y(nn,:), '-', 'Color', plotcol);
                        
                    end
                end
                
            end
            
            
        end
        
        linkaxes(hsplots, 'xy');
    end
end













