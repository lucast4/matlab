function lt_neural_POPLEARN_SummaryPlot1(OUTSTRUCT, MOTIFSTATS_pop, SwitchStruct, ...
    birdnum, exptnum)

%% lt 3/20/18 - plot all CC across
i=birdnum;

numexpts = length(OUTSTRUCT.bird(i).expt);

for ii=1:numexpts
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT
    
    numswitches = length(OUTSTRUCT.bird(i).expt(ii).swnum);
    nummotifs = length(OUTSTRUCT.bird(i).expt(ii).swnum(1).PRE.motif);
    
    % --- collect motif list
    motiflist = {MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(1).motif.regexpstr};
    for j=2:length(MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum)
        assert(all(strcmp(motiflist, {MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(j).motif.regexpstr})), 'saf');
    end
    
    assert(nummotifs == length(motiflist),  'asfsd');
    
    % ============ COLLECT FOR PRE AND POST SEPARATELY.
    epoch = 'PRE';
    DATTMP = struct;
    for m=1:nummotifs
        
        SwnumAll = [];
        CCAll = [];
        XlagsAll = [];
        for ss = 1:numswitches
            
            DAT = OUTSTRUCT.bird(i).expt(ii).swnum(ss);
            
            if length(DAT.(epoch).motif)<m
                continue
            end
            if isempty(DAT.(epoch).motif(m).neurset)
                continue
            end
            
            if isempty(DAT.(epoch))
                continue
            end
            
            numsets = length(DAT.(epoch).motif(m).neurset);
                        
            % --- COLLECT ALL Xcov (even across diff sets of neurons pairs)
            CCthis = [DAT.(epoch).motif(m).neurset.CCallpairs];
            xlagstmp = [DAT.(epoch).motif(m).neurset.xlags_sec];
            if ~isempty(xlagstmp)
                xlags = xlagstmp(:,1);
            end
            
            SwnumAll = [SwnumAll, ones(1, size(CCthis,2))*ss];
            CCAll = [CCAll CCthis];
            
        end
        DATTMP.(epoch).motif(m).SwnumAll = SwnumAll;
        DATTMP.(epoch).motif(m).CCAll = CCAll;
        DATTMP.(epoch).motif(m).xlags = xlags;
    end
    
    epoch = 'POST';
    for m=1:nummotifs
        
        SwnumAll = [];
        CCAll = [];
        XlagsAll = [];
        for ss = 1:numswitches
            
            DAT = OUTSTRUCT.bird(i).expt(ii).swnum(ss);
            
            if isempty(DAT.(epoch))
                continue
            end
            
            numsets = length(DAT.(epoch).motif(m).neurset);
            
            % --- COLLECT ALL Xcov (even across diff sets of neurons pairs)
            CCthis = [DAT.(epoch).motif(m).neurset.CCallpairs];
            xlags = [DAT.(epoch).motif(m).neurset.xlags_sec];
            xlagstmp = [DAT.(epoch).motif(m).neurset.xlags_sec];
            if ~isempty(xlagstmp)
                xlags = xlagstmp(:,1);
            end
            
            SwnumAll = [SwnumAll, ones(1, size(CCthis,2))*ss];
            CCAll = [CCAll CCthis];
        end
        DATTMP.(epoch).motif(m).SwnumAll = SwnumAll;
        DATTMP.(epoch).motif(m).CCAll = CCAll;
        DATTMP.(epoch).motif(m).xlags = xlags;
    end
    
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT ACROSS SWITCHES FOR EACH  MOTIF
    numswitches = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
    
    figcount=1;
    subplotrows=6;
    subplotcols= numswitches*2;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    for m=1:nummotifs
        
        % ---------- INDIVIDUAL SUBPLOT FOR EACH SWITCH
        for ss = 1:numswitches
            
            
            % ======= PRE
            epoch = 'PRE';
            plotcol = 'k';
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title([motiflist{m} ',sw' num2str(ss) '[' epoch ']']);
            
            SwnumAll = DATTMP.(epoch).motif(m).SwnumAll;
            CCall = DATTMP.(epoch).motif(m).CCAll;
            xlags = DATTMP.(epoch).motif(m).xlags;
            ccthis = CCall(:, SwnumAll==ss);
            
            if ~isempty(ccthis)
                
                plot(xlags, ccthis, '-', 'Color', plotcol);
                ccmean = mean(ccthis,2);
                ccsem = lt_sem(ccthis');
                shadedErrorBar(xlags, ccmean, ccsem, {'Color', plotcol}, 1);
                lt_plot_zeroline;
                
            end
            
            % ======= POST
            epoch = 'POST';
            plotcol = 'b';
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title([motiflist{m} ',sw' num2str(ss) '[' epoch ']']);
            
            SwnumAll = DATTMP.(epoch).motif(m).SwnumAll;
            CCall = DATTMP.(epoch).motif(m).CCAll;
            xlags = DATTMP.(epoch).motif(m).xlags;
            ccthis = CCall(:, SwnumAll==ss);
            
            if ~isempty(ccthis)
                
                plot(xlags, ccthis, '-', 'Color', plotcol);
                ccmean = mean(ccthis,2);
                ccsem = lt_sem(ccthis');
                shadedErrorBar(xlags, ccmean, ccsem, {'Color', plotcol}, 1);
                lt_plot_zeroline;
                
            end
            
        end
        
        
    end
    
    linkaxes(hsplots, 'xy');
    
    %% ============= PLOT ALL SWITCHES, PAIRED (WITHIN NEURON) ACROSS SWITCH
    numswitches = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
    
    figcount=1;
    subplotrows=6;
    subplotcols= numswitches*2;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    for m=1:nummotifs
        
        % ---------- INDIVIDUAL SUBPLOT FOR EACH SWITCH
        for ss = 1:numswitches
            
            % =================== only plot if this switch has both pre and
            % post
            tmp1 = isempty(DATTMP.PRE.motif(m).CCAll);
            tmp2 = isempty(DATTMP.POST.motif(m).CCAll);
            if any([tmp1 tmp2] ==1)
                continue
            end
            
            
            % ======= PRE
            epoch = 'PRE';
            
            SwnumAll = DATTMP.(epoch).motif(m).SwnumAll;
            CCall = DATTMP.(epoch).motif(m).CCAll;
            xlags = DATTMP.(epoch).motif(m).xlags;
            ccthisPRE = CCall(:, SwnumAll==ss);
            
            
            % ======= POST
            epoch = 'POST';
            SwnumAll = DATTMP.(epoch).motif(m).SwnumAll;
            CCall = DATTMP.(epoch).motif(m).CCAll;
            xlags = DATTMP.(epoch).motif(m).xlags;
            ccthisPOST = CCall(:, SwnumAll==ss);
            
            % ======= GET CHANGE (POST - PRE);
            assert(all(size(ccthisPRE) == size(ccthisPOST)), 'asdf');
            
            % --- get diff
            ccDIFF = ccthisPOST - ccthisPRE;
            
            
            % ======================== PLOT
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title([motiflist{m} ',sw' num2str(ss) '[POST-PRE]']);
            plot(xlags, ccDIFF, '-', 'Color', 'r');
            ccmean = mean(ccDIFF,2);
            ccsem = lt_sem(ccDIFF');
            shadedErrorBar(xlags, ccmean, ccsem, {'Color', 'r'}, 1);
            lt_plot_zeroline;
            
        end
        
        
    end
    
    linkaxes(hsplots, 'xy');
    
end

