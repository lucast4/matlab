function lt_neural_NGRAMS_PlotEgAll(OUTSTRUCT, SummaryStruct, Params, ...
    birdtoplot, neurtoplot, ngramstoplot, pairtypetoplot)
%% lt 3/2019 - for a given neuron, plot all cases for each pairtype

%%
savedir = '/bluejay5/lucas/analyses/neural/NGRAMS';
savedir = [savedir '/' Params.dirname];

%%
motifpredur = Params.regexpr.motifpredur;
windowx = motifpredur+Params.window_prem;

%% ========== go thru all birds and extract
numbirds = length(SummaryStruct.birds);
for i=1:numbirds
    
    birdname = SummaryStruct.birds(i).birdname;
    if ~strcmp(birdname, birdtoplot)
        continue
    end
    
    tmp = load([savedir '/bird' num2str(i) '.mat']);
    birdstruct = tmp.birdstruct;
    % ============================== calculate things across all neurons
    numneurons = length(birdstruct.neuron);
    
    for nn=1:numneurons
        if nn ~= neurtoplot
            continue
        end
        
        disp(' ############################################################### ');
        disp([birdname ', ' num2str(i) '-' num2str(nn)]);
        ngramlist = birdstruct.neuron(nn).ngramlist;
        
        % -----------
        if isempty(ngramstoplot)
            
            if ~isempty(pairtypetoplot)
                ptypeind = find(strcmp(OUTSTRUCT.PairTypesInOrder, pairtypetoplot));
                indtmp = find(OUTSTRUCT.All_birdnum==i & OUTSTRUCT.All_neurnum==nn ...
                    & OUTSTRUCT.All_diffsyl_PairType==ptypeind);
            else
                indtmp = find(OUTSTRUCT.All_birdnum==i & OUTSTRUCT.All_neurnum==nn);
            end
            % -- pikc a random pair
            indOutStr = indtmp(randperm(length(indtmp), 1));
            tmp = OUTSTRUCT.All_ngrampair_inorder(indOutStr, :);
            j = tmp(1);
            jj = tmp(2);
        else
            asdfasdfasdf;
        end
        % ----------- all pairwise ngrams
        
        % =================== 1) FOR EACH PAIR TY
        % ---- syl 1
        motif1 = birdstruct.neuron(nn).ngramnum(j).regexprstr;
        motif2 = birdstruct.neuron(nn).ngramnum(jj).regexprstr;
        
        % ---- sample size
        N1 = size(birdstruct.neuron(nn).ngramnum(j).DAT.frmat, 2);
        N2 = size(birdstruct.neuron(nn).ngramnum(jj).DAT.frmat, 2);
        
        % ---- pairtypes
        pairtype = OUTSTRUCT.All_diffsyl_PairType(indOutStr);
        
        % ---- osnets and offests
        sylon1 = median(birdstruct.neuron(nn).ngramnum(j).DAT.motifsylOn,1);
        syloff1 = median(birdstruct.neuron(nn).ngramnum(j).DAT.motifsylOff,1);
        
        sylon2 = median(birdstruct.neuron(nn).ngramnum(jj).DAT.motifsylOn,1);
        syloff2 = median(birdstruct.neuron(nn).ngramnum(jj).DAT.motifsylOff,1);
        
        % ---- convert to time relative to onset of premotor window
        sylon1 = sylon1-motifpredur-Params.window_prem(1);
        syloff1 = syloff1-motifpredur-Params.window_prem(1);
        sylon2 = sylon2-motifpredur-Params.window_prem(1);
        syloff2 = syloff2-motifpredur-Params.window_prem(1);
        
        %         % --- get the first offset that is after onset
        %         syloff1 = min(syloff1(syloff1>motifpredur)) - motifpredur;
        %         syloff2 = min(syloff2(syloff2>motifpredur)) - motifpredur;
        
        
        % ==================
        
        %% ################################## PLOT
        lt_figure; hold on;
        pcolors = lt_make_plot_colors(3,0, 1);
        hsplots = [];
        
        % ================ FR STUFF
        frmat1 = birdstruct.neuron(nn).ngramnum(j).DAT.frmat;
        frmat2 = birdstruct.neuron(nn).ngramnum(jj).DAT.frmat;
        
        x = lt_neural_QUICK_XfromFRmat(frmat1);
        x = x-motifpredur-Params.window_prem(1);
        
        % ========================= ACTUAL DAT
        hsplot = lt_subplot(2,2,1); hold on;
        hsplots = [hsplots hsplot];
        title(['[' motif1 '-' motif2 '] actual dat']);
        plot(x, mean(sqrt(frmat1),2), '-', 'Color', pcolors{1});
        plot(x, mean(sqrt(frmat2),2), '-', 'Color', pcolors{2});
        
        % -------------- line for syl onset and offsets
        YLIM = ylim;
        for k=1:length(sylon1)
            on = sylon1(k);
            off = syloff1(k);
            pcol = pcolors{1};
            h = patch([on off off on], [YLIM(1) YLIM(1) YLIM(2) YLIM(2)], pcol);
            set(h, 'FaceAlpha', 0.1);
            set(h, 'EdgeColor', 'none');
            
            on = sylon2(k);
            off = syloff2(k);
            pcol = pcolors{2};
            h = patch([on off off on], [YLIM(1) YLIM(1) YLIM(2) YLIM(2)], pcol);
            set(h, 'FaceAlpha', 0.1);
            set(h, 'EdgeColor', 'none');
        end
        
        
        % ============================ FR window
        nshufftmp=2;
        Nmin = Params.Nmin;
        downfactor= 1;
        [FRdiffDAT, FRdiffShuff, FRdiff_Z, Nboth, FRmat, TrialInds] = ...
            lt_neural_NGRAMS_QUICK_FRdiff(frmat1, frmat2, downfactor, windowx, ...
            Nmin, nshufftmp, Params, Params.dodecode);
        
        
        % --------------- plot
        x = lt_neural_QUICK_XfromFRmat(FRmat) + 0.001;
        
        % ============================= SUBPLOT 1 - ACTUAL DAT
        
        Yall = {};
        % --- motif 1
        indmotif = 1;
        plot(x, mean(FRmat(:, TrialInds{indmotif}),2), '-', 'Color', pcolors{indmotif});
        Yall{indmotif} = mean(FRmat(:, TrialInds{indmotif}),2);
        
        % --- motif 2
        indmotif = 2;
        plot(x, mean(FRmat(:, TrialInds{indmotif}),2), '-', 'Color', pcolors{indmotif});
        Yall{indmotif} = mean(FRmat(:, TrialInds{indmotif}),2);
        
        yabsdiff = abs(Yall{2} - Yall{1});
        ylower = min([Yall{1} Yall{2}]');
        
        for k = 1:length(x)
            xx = x(k);
            line([xx xx], [ylower(k) ylower(k)+yabsdiff(k)], 'Color', [0.6 0.6 0.6]);
        end
        lt_plot_text(0, min(min(cell2mat(Yall)))+1, ['FRdiff=' num2str(FRdiffDAT)], 'r');
        assert(FRdiffDAT == mean(abs(Yall{2} - Yall{1})), 'asfasd');
        
        
        %% =================== same thing, but shuffle frmats first
        frmat1 = birdstruct.neuron(nn).ngramnum(j).DAT.frmat;
        frmat2 = birdstruct.neuron(nn).ngramnum(jj).DAT.frmat;
        
        FRmat = [frmat1 frmat2];
        FRmat = FRmat(:, randperm(size(FRmat,2)));
        frmat1 = FRmat(:, TrialInds{1});
        frmat2 = FRmat(:, TrialInds{2});
        
        x = lt_neural_QUICK_XfromFRmat(frmat1);
        x = x-motifpredur-Params.window_prem(1);
        
        % ========================= ACTUAL DAT
        hsplot = lt_subplot(2,2,2); hold on;
        hsplots = [hsplots hsplot];
        
        title(['[' motif1 '-' motif2 '] shuffle']);
        plot(x, mean(sqrt(frmat1),2), '-', 'Color', pcolors{1});
        plot(x, mean(sqrt(frmat2),2), '-', 'Color', pcolors{2});
        
        % -------------- line for syl onset and offsets
        YLIM = ylim;
        for k=1:length(sylon1)
            on = sylon1(k);
            off = syloff1(k);
            pcol = pcolors{1};
            h = patch([on off off on], [YLIM(1) YLIM(1) YLIM(2) YLIM(2)], pcol);
            set(h, 'FaceAlpha', 0.1);
            set(h, 'EdgeColor', 'none');
            
            on = sylon2(k);
            off = syloff2(k);
            pcol = pcolors{2};
            h = patch([on off off on], [YLIM(1) YLIM(1) YLIM(2) YLIM(2)], pcol);
            set(h, 'FaceAlpha', 0.1);
            set(h, 'EdgeColor', 'none');
        end
        
        
        % ============================ FR window
        nshufftmp=2;
        Nmin = Params.Nmin;
        downfactor= 1;
        [FRdiffDAT, FRdiffShuff, FRdiff_Z, Nboth, FRmat, TrialInds] = ...
            lt_neural_NGRAMS_QUICK_FRdiff(frmat1, frmat2, downfactor, windowx, ...
            Nmin, nshufftmp, Params, Params.dodecode);
        
        
        % --------------- plot
        x = lt_neural_QUICK_XfromFRmat(FRmat) + 0.001;
        
        % ============================= SUBPLOT 1 - ACTUAL DAT
        
        Yall = {};
        % --- motif 1
        indmotif = 1;
        plot(x, mean(FRmat(:, TrialInds{indmotif}),2), '-', 'Color', pcolors{indmotif});
        Yall{indmotif} = mean(FRmat(:, TrialInds{indmotif}),2);
        
        % --- motif 2
        indmotif = 2;
        plot(x, mean(FRmat(:, TrialInds{indmotif}),2), '-', 'Color', pcolors{indmotif});
        Yall{indmotif} = mean(FRmat(:, TrialInds{indmotif}),2);
        
        yabsdiff = abs(Yall{2} - Yall{1});
        ylower = min([Yall{1} Yall{2}]');
        
        for k = 1:length(x)
            xx = x(k);
            line([xx xx], [ylower(k) ylower(k)+yabsdiff(k)], 'Color', [0.6 0.6 0.6]);
        end
        lt_plot_text(0, min(min(cell2mat(Yall)))+1, ['FRdiff=' num2str(FRdiffDAT)], 'r');
        assert(FRdiffDAT == mean(abs(Yall{2} - Yall{1})), 'asfasd');
        
        
        % ============
        linkaxes(hsplots, 'xy');
    end
end



if (0) % OLD VERSION, only plots within premotor window ...
    
    %% ========== go thru all birds and extract
    numbirds = length(SummaryStruct.birds);
    for i=1:numbirds
        
        birdname = SummaryStruct.birds(i).birdname;
        if ~strcmp(birdname, birdtoplot)
            continue
        end
        
        tmp = load([savedir '/bird' num2str(i) '.mat']);
        birdstruct = tmp.birdstruct;
        % ============================== calculate things across all neurons
        numneurons = length(birdstruct.neuron);
        
        for nn=1:numneurons
            if nn ~= neurtoplot
                continue
            end
            
            disp(' ############################################################### ');
            disp([birdname ', ' num2str(i) '-' num2str(nn)]);
            ngramlist = birdstruct.neuron(nn).ngramlist;
            
            % -----------
            if isempty(ngramstoplot)
                
                if ~isempty(pairtypetoplot)
                    ptypeind = find(strcmp(OUTSTRUCT.PairTypesInOrder, pairtypetoplot));
                    indtmp = find(OUTSTRUCT.All_birdnum==i & OUTSTRUCT.All_neurnum==nn ...
                        & OUTSTRUCT.All_diffsyl_PairType==ptypeind);
                else
                    indtmp = find(OUTSTRUCT.All_birdnum==i & OUTSTRUCT.All_neurnum==nn);
                end
                % -- pikc a random pair
                indOutStr = indtmp(randperm(length(indtmp), 1));
                tmp = OUTSTRUCT.All_ngrampair_inorder(indOutStr, :);
                j = tmp(1);
                jj = tmp(2);
            else
                asdfasdfasdf;
            end
            % ----------- all pairwise ngrams
            
            % =================== 1) FOR EACH PAIR TY
            % ---- syl 1
            motif1 = birdstruct.neuron(nn).ngramnum(j).regexprstr;
            motif2 = birdstruct.neuron(nn).ngramnum(jj).regexprstr;
            
            % ---- sample size
            N1 = size(birdstruct.neuron(nn).ngramnum(j).DAT.frmat, 2);
            N2 = size(birdstruct.neuron(nn).ngramnum(jj).DAT.frmat, 2);
            
            % ---- pairtypes
            pairtype = OUTSTRUCT.All_diffsyl_PairType(indOutStr);
            
            % ---- osnets and offests
            %         sylon1 = median(birdstruct.neuron(nn).ngramnum(j).DAT.motifsylOn,1);
            syloff1 = median(birdstruct.neuron(nn).ngramnum(j).DAT.motifsylOff,1);
            
            %         sylon2 = median(birdstruct.neuron(nn).ngramnum(jj).DAT.motifsylOn,1);
            syloff2 = median(birdstruct.neuron(nn).ngramnum(jj).DAT.motifsylOff,1);
            
            % --- get the first offset that is after onset
            syloff1 = min(syloff1(syloff1>motifpredur)) - motifpredur;
            syloff2 = min(syloff2(syloff2>motifpredur)) - motifpredur;
            
            
            % ================ FR STUFF
            frmat1 = birdstruct.neuron(nn).ngramnum(j).DAT.frmat;
            frmat2 = birdstruct.neuron(nn).ngramnum(jj).DAT.frmat;
            
            
            % ================ FR WITHIN PREMOTOR WINDOW
            nshufftmp=2;
            Nmin = Params.Nmin;
            downfactor= 1;
            [FRdiffDAT, FRdiffShuff, FRdiff_Z, Nboth, FRmat, TrialInds] = ...
                lt_neural_NGRAMS_QUICK_FRdiff(frmat1, frmat2, downfactor, windowx, ...
                Nmin, nshufftmp);
            
            
            
            % ########################################################## plot
            lt_figure; hold on;
            x = lt_neural_QUICK_XfromFRmat(FRmat);
            pcolors = lt_make_plot_colors(3,0, 1);
            
            
            
            % ============================= SUBPLOT 1 - ACTUAL DAT
            lt_subplot(2,2,1); hold on;
            title(['[' motif1 '-' motif2 '] actual dat']);
            
            Yall = {};
            % --- motif 1
            indmotif = 1;
            plot(x, mean(FRmat(:, TrialInds{indmotif}),2), '-', 'Color', pcolors{indmotif});
            Yall{indmotif} = mean(FRmat(:, TrialInds{indmotif}),2);
            
            % --- motif 2
            indmotif = 2;
            plot(x, mean(FRmat(:, TrialInds{indmotif}),2), '-', 'Color', pcolors{indmotif});
            Yall{indmotif} = mean(FRmat(:, TrialInds{indmotif}),2);
            
            yabsdiff = abs(Yall{2} - Yall{1});
            ylower = min([Yall{1} Yall{2}]');
            
            for k = 1:length(x)
                xx = x(k);
                line([xx xx], [ylower(k) ylower(k)+yabsdiff(k)], 'Color', [0.6 0.6 0.6]);
            end
            lt_plot_text(0, min(min(cell2mat(Yall)))+1, ['FRdiff=' num2str(FRdiffDAT)], 'r');
            assert(FRdiffDAT == mean(abs(Yall{2} - Yall{1})), 'asfasd');
            
            % --- indicate syl onsets/offset
            onset = -Params.window_prem(1);
            line([onset onset], ylim, 'Color', 'g');
            offset = syloff1+onset;
            line([offset offset], ylim, 'Color', 'g');
            
            
            
            % ============================ SUBPLOT 2 - SINGLE SHUFF TRIAL
            lt_subplot(2,2,2); hold on;
            title('single shufffle');
            
            % ------- first, shuffle TrialInds
            FRmatShuff = FRmat(:, randperm(TrialInds{end}(end)));
            
            Yall = {};
            % --- motif 1
            indmotif = 1;
            plot(x, mean(FRmatShuff(:, TrialInds{indmotif}),2), '-', 'Color', pcolors{indmotif});
            Yall{indmotif} = mean(FRmatShuff(:, TrialInds{indmotif}),2);
            
            % --- motif 2
            indmotif = 2;
            plot(x, mean(FRmatShuff(:, TrialInds{indmotif}),2), '-', 'Color', pcolors{indmotif});
            Yall{indmotif} = mean(FRmatShuff(:, TrialInds{indmotif}),2);
            
            yabsdiff = abs(Yall{2} - Yall{1});
            ylower = min([Yall{1} Yall{2}]');
            
            for k = 1:length(x)
                xx = x(k);
                line([xx xx], [ylower(k) ylower(k)+yabsdiff(k)], 'Color', [0.6 0.6 0.6]);
            end
            lt_plot_text(0, min(min(cell2mat(Yall)))+1, ...
                ['FRdiff=' num2str(mean(abs(Yall{2} - Yall{1})))], 'r');
            
            % --- indicate syl onsets/offset
            onset = -Params.window_prem(1);
            line([onset onset], ylim, 'Color', 'g');
            offset = syloff1+onset;
            line([offset offset], ylim, 'Color', 'g');
            
        end
    end
    
    
end









