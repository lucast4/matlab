function lt_neural_NGRAMS_PlotEgPair(OUTSTRUCT, SummaryStruct, Params, ...
    birdtoplot, neurtoplot, pairtypes, plotsqrt)

%%
savedir = '/bluejay5/lucas/analyses/neural/NGRAMS';
savedir = [savedir '/' Params.dirname];

%%
motifpredur = Params.regexpr.motifpredur;
windowx = motifpredur+Params.window_prem;

%% ========== go thru all birds and extract
% numbirds = length(SummaryStruct.birds);

birdname = birdtoplot;
i = find(strcmp({SummaryStruct.birds.birdname}, birdname));
nn = neurtoplot;

tmp = load([savedir '/bird' num2str(i) '.mat']);
birdstruct = tmp.birdstruct;


% ============================== calculate things across all neurons
disp(' ############################################################### ');
disp([birdname ', ' num2str(i) '-' num2str(nn)]);
ngramlist = birdstruct.neuron(nn).ngramlist;


%% ====================
bnum =  find(strcmp({SummaryStruct.birds.birdname}, birdname));
neur = nn;
bregion = SummaryStruct.birds(bnum).neurons(nn).NOTE_Location;
%% what pair types to plot?

for i=1:length(pairtypes)
    
    figcount=1;
    subplotrows=4;
    subplotcols=8;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    ptypeind = find(strcmp(OUTSTRUCT.PairTypesInOrder, pairtypes{i}));
    
    indtmp = find(OUTSTRUCT.All_birdnum==bnum & OUTSTRUCT.All_neurnum==nn ...
        & OUTSTRUCT.All_diffsyl_PairType==ptypeind);
    
    
    ngrampairs = OUTSTRUCT.All_ngrampair_inorder(indtmp, :);
    ngrampairs_str = OUTSTRUCT.All_ngramstring_inorder(indtmp, :);
    
    
    % ================ plot random subset if too many
    if size(ngrampairs,1)>32
        indstmp = randperm(size(ngrampairs, 1), 32);
        ngrampairs = ngrampairs(indstmp,:);
        ngrampairs_str = ngrampairs_str(indstmp, :);
    end
    
    for ii=1:size(ngrampairs,1)
        
        % -- 
        j = ngrampairs(ii,1); % first motif
        jj = ngrampairs(ii,2); % second motif
        
        
        % =================== 1) FOR EACH PAIR TY
        motif1 = birdstruct.neuron(nn).ngramnum(j).regexprstr;
        motif2 = birdstruct.neuron(nn).ngramnum(jj).regexprstr;
        assert(strcmp(motif1(regexp(motif1, '[a-z]')), ngrampairs_str{ii, 1}));
        assert(strcmp(motif2(regexp(motif2, '[a-z]')), ngrampairs_str{ii, 2}));
        
        % ---- sample size
        %     N1 = size(birdstruct.neuron(nn).ngramnum(j).DAT.frmat, 2);
        %     N2 = size(birdstruct.neuron(nn).ngramnum(jj).DAT.frmat, 2);
        
        % ---- pairtypes
        %     pairtype = OUTSTRUCT.All_diffsyl_PairType(indOutStr);
        
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
        
        
        
        % ############################ PLOT THIS PAIR
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(['[' motif1 '-' motif2 ']']);
        ylabel(pairtypes{i});
        if ii==1
            xlabel([birdname '-neur' num2str(neur) '-' bregion]);
        end
        pcolors = lt_make_plot_colors(3,0, 1);
        
        % ================ FR STUFF
        frmat1 = birdstruct.neuron(nn).ngramnum(j).DAT.frmat;
        frmat2 = birdstruct.neuron(nn).ngramnum(jj).DAT.frmat;
        x = lt_neural_QUICK_XfromFRmat(frmat1);
        x = x-motifpredur-Params.window_prem(1);
        
        
        % ######################## ACTUAL DAT
        if plotsqrt==1
        shadedErrorBar(x, mean(sqrt(frmat1),2), lt_sem(sqrt(frmat1)'), {'Color', pcolors{1}, ...
            'LineWidth', 2}, 1)
        shadedErrorBar(x, mean(sqrt(frmat2),2), lt_sem(sqrt(frmat2)'), {'Color', pcolors{2}, ...
            'LineWidth', 2}, 1)
        else
        shadedErrorBar(x, mean(frmat1,2), lt_sem(frmat1'), {'Color', pcolors{1}, ...
            'LineWidth', 2}, 1)
        shadedErrorBar(x, mean(frmat2,2), lt_sem(frmat2'), {'Color', pcolors{2}, ...
            'LineWidth', 2}, 1)
        end


        % ################################### OVERLAY SHUFFLE
        % ---------- LINE FOR PREMOTOR WINDOW
        [~, ~, ~, ~, ~, TrialInds] = ...
            lt_neural_NGRAMS_QUICK_FRdiff(frmat1, frmat2, 1, windowx, ...
            Params.Nmin, 2, Params, Params.dodecode);

        FRmat = [frmat1 frmat2];
        FRmat = FRmat(:, randperm(size(FRmat,2)));
        frmat1 = FRmat(:, TrialInds{1});
        frmat2 = FRmat(:, TrialInds{2});
        
        
        if plotsqrt==1
%         shadedErrorBar(x, mean(sqrt(frmat1),2), lt_sem(sqrt(frmat1)'), {'Color', pcolors{1}}, 1)
%         shadedErrorBar(x, mean(sqrt(frmat2),2), lt_sem(sqrt(frmat2)'), {'Color', pcolors{2}}, 1)
        plot(x, mean(sqrt(frmat1),2), '-', 'Color', pcolors{1})
        plot(x, mean(sqrt(frmat2),2), '-', 'Color', pcolors{2})
        else
        plot(x, mean(frmat1,2), '-', 'Color', pcolors{1})
        plot(x, mean(frmat2,2), '-', 'Color', pcolors{2})
        end

        
        
        xlim([-0.05 0.1]);
        if plotsqrt==1
        ylim([0 30]);
        else
            ylim([0 500]);
        end
        % -------------- line for syl onset and offsets
        YLIM = ylim;
        %         lt_neural_QUICK_PlotSylPatches(sylon1, syloff1, YLIM, 0, pcolors{1})
        lt_neural_QUICK_PlotSylPatches(sylon1, syloff1, YLIM(2), 1, pcolors{1})
        lt_neural_QUICK_PlotSylPatches(sylon2, syloff2, 0.9*YLIM(2), 1, pcolors{2})
        

        % ---------- LINE FOR PREMOTOR WINDOW
        [FRdiffDAT, FRdiffShuff, FRdiff_Z, Nboth, FRmat, TrialInds] = ...
            lt_neural_NGRAMS_QUICK_FRdiff(frmat1, frmat2, 1, windowx, ...
            Params.Nmin, 2, Params, Params.dodecode);
        
        % --------------- plot
        x = lt_neural_QUICK_XfromFRmat(FRmat) + 0.001;
        line([x(1) x(1)], ylim, 'Color', 'k');
        line([x(end) x(end)], ylim, 'Color', 'k');
        
        
        % ============
        linkaxes(hsplots, 'xy');
        
        
        
    end
end
% ==================

