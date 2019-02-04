function lt_neural_POPLEARN_SummaryPlot3(OUTSTRUCT, MOTIFSTATS_pop, SwitchStruct, ...
    birdnum, windowmean, SkipIfTargsDiffSyls)

%% -- BREAK UP INTO TARG, NONTARG CONTEXT, NONTARG SYL
numbirds = length(OUTSTRUCT.bird);
for i=1:numbirds
    
    numexpts = length(OUTSTRUCT.bird(i).expt);
    birdname = SwitchStruct.bird(i).birdname;
    for ii=1:numexpts
        
        % --- initiate figure
        figcount=1;
        subplotrows=3;
        subplotcols=5;
        fignums_alreadyused=[];
        hfigs=[];
        hsplots = [];
        
        % =========== go thru switches
        exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
        
        numswitches = length(OUTSTRUCT.bird(i).expt(ii).swnum);
        assert(numswitches == length(SwitchStruct.bird(i).exptnum(ii).switchlist), 'asdf');
        
        % =========== COLLECT MOTIF LIST
        
        % --- collect motif list
        motiflist = {MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(1).motif.regexpstr};
        for j=2:length(MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum)
            assert(all(strcmp(motiflist, {MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(j).motif.regexpstr})), 'saf');
        end
        
        
        % ===========
        for ss =1:numswitches
            swthis = SwitchStruct.bird(i).exptnum(ii).switchlist(ss);
            targsyls = swthis.learningContingencies(1:2:end);
            if swthis.targsAreSameSyl==0
                continue
            end
            
            if ~any([swthis.neurFakeID.haspresongs]) & ...
                    ~any([swthis.neurFakeID.haspostsongs]);
                %  then no songs, so skip
                continue
            end
            
            % ##################################### PRE
            epochthis = 'PRE';
            nummotifs = length(OUTSTRUCT.bird(i).expt(ii).swnum(ss).(epochthis).motif);
            assert(nummotifs == length(motiflist),  'asfsd');
            
            % ======================== collect all CC
            DAT = OUTSTRUCT.bird(i).expt(ii).swnum(ss).(epochthis);
            
            AllCC = [];
            AllIsTarg = [];
            AllIsSame = [];
            AllMotifNum = [];
            for mm=1:nummotifs
                
                % ==== what class is this motif?
                sylthis = motiflist{mm};
                
                % -------------- TARG?
                if any(strcmp(sylthis, targsyls))==1
                    istarg =1;
                else
                    istarg=0;
                end
                
                % -------------- SAME?
                issameall = [];
                for kk=1:length(targsyls)
                    [issame, syl1, syl2] = lt_neural_QUICK_SylPairType(sylthis, targsyls{kk});
                    issameall = [issameall issame];
                end
                issame = unique(issameall);
                assert(length(issame)==1, 'asfdas');
                
                
                % === go thru each set of neurons and collect CC
                numsets = length(DAT.motif(mm).neurset);
                ccMeans = [];
                for tt =1:numsets
                    CCall = DAT.motif(mm).neurset(tt).CCallpairs;
                    xlags = DAT.motif(mm).neurset(tt).xlags_sec;
                    
                    if isempty(CCall)
                        continue
                    end
                    % ==== get sum within a lag window
                    indtmp = xlags >= windowmean(1) & xlags <= windowmean(2);
                    ccMeans = mean(CCall(indtmp,:),1);
                end
                
                % ======================== COLLECT DATA
                AllCC = [AllCC; ccMeans'];
                AllIsTarg = [AllIsTarg; istarg*ones(size(ccMeans'))];
                AllIsSame = [AllIsSame; issame*ones(size(ccMeans'))];
                AllMotifNum = [AllMotifNum; mm*ones(size(ccMeans'))];
            end
            
            
            % ================== PLOT HISTOGRAM
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([birdname ',' exptname ',sw' num2str(ss) ',[' epochthis ']']);
            xlabel('xcov');
            
            if ~isempty(AllCC)
                % === annotate targets and switch status
                lt_plot_annotation(1, swthis.learningContingencies);
                
                % --- get xbins
                [Ybinned, Xcenters, hbar]=lt_plot_histogram(AllCC, ...
                    [], 0, 0, 0.5, 0, [0.6 0.6 0.6]);
                
                % --- diff syls
                indstmp = AllIsTarg==0 & AllIsSame==0;
                lt_plot_histogram(AllCC(indstmp), ...
                    Xcenters, 1, 0.4, 0.5, 1, [0.7 0.7 0.7]);
                
                % --- histogram of target
                indstmp = AllIsTarg==1;
                targmotifs = unique(AllMotifNum(indstmp));
                pcols = lt_make_plot_colors(length(targmotifs), 0,0);
                for mmm=1:length(targmotifs)
                    indstmp = AllIsTarg==1 & AllMotifNum==targmotifs(mmm);
                    pcol = pcols{mmm};
                    lt_plot_histogram(AllCC(indstmp), ...
                        Xcenters, 1, 0.4, 0.5, 1, pcol);
                    lt_plot_annotation(2, ['targ ' motiflist{targmotifs(mmm)}], pcol);
                end
                
                % --- histogram of same types
                indstmp = AllIsSame==1 & AllIsTarg==0;
                lt_plot_histogram(AllCC(indstmp), ...
                    Xcenters, 1, 0.4, 0.5, 1, [0.3 0.3 0.7]);
            end
            
            
            
            % ##################################### POST
            epochthis = 'POST';
            nummotifs = length(OUTSTRUCT.bird(i).expt(ii).swnum(ss).(epochthis).motif);
            assert(nummotifs == length(motiflist),  'asfsd');
            
            % ======================== collect all CC
            DAT = OUTSTRUCT.bird(i).expt(ii).swnum(ss).(epochthis);
            
            AllCC = [];
            AllIsTarg = [];
            AllIsSame = [];
            AllMotifNum = [];
            for mm=1:nummotifs
                
                % ==== what class is this motif?
                sylthis = motiflist{mm};
                
                % -------------- TARG?
                if any(strcmp(sylthis, targsyls))==1
                    istarg =1;
                else
                    istarg=0;
                end
                
                % -------------- SAME?
                issameall = [];
                for kk=1:length(targsyls)
                    [issame, syl1, syl2] = lt_neural_QUICK_SylPairType(sylthis, targsyls{kk});
                    issameall = [issameall issame];
                end
                issame = unique(issameall);
                assert(length(issame)==1, 'asfdas');
                
                
                % === go thru each set of neurons and collect CC
                numsets = length(DAT.motif(mm).neurset);
                ccMeans = [];
                for tt =1:numsets
                    CCall = DAT.motif(mm).neurset(tt).CCallpairs;
                    xlags = DAT.motif(mm).neurset(tt).xlags_sec;
                    
                    if isempty(CCall)
                        continue
                    end
                    % ==== get sum within a lag window
                    indtmp = xlags >= windowmean(1) & xlags <= windowmean(2);
                    ccMeans = mean(CCall(indtmp,:),1);
                end
                
                % ======================== COLLECT DATA
                AllCC = [AllCC; ccMeans'];
                AllIsTarg = [AllIsTarg; istarg*ones(size(ccMeans'))];
                AllIsSame = [AllIsSame; issame*ones(size(ccMeans'))];
                AllMotifNum = [AllMotifNum; mm*ones(size(ccMeans'))];
            end
            
            
            % ================== PLOT HISTOGRAM
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([birdname ',' exptname ',sw' num2str(ss) ',[' epochthis ']']);
            
            if ~isempty(AllCC)
                % === annotate targets and switch status
                lt_plot_annotation(1, swthis.learningContingencies);
                
                % --- get xbins
                [Ybinned, Xcenters, hbar]=lt_plot_histogram(AllCC, ...
                    [], 0, 0, 0.5, 0, [0.6 0.6 0.6]);
                
                % --- diff syls
                indstmp = AllIsTarg==0 & AllIsSame==0;
                lt_plot_histogram(AllCC(indstmp), ...
                    Xcenters, 1, 0.4, 0.5, 1, [0.7 0.7 0.7]);
                
                % --- histogram of target
                indstmp = AllIsTarg==1;
                targmotifs = unique(AllMotifNum(indstmp));
                pcols = lt_make_plot_colors(length(targmotifs), 0,0);
                for mmm=1:length(targmotifs)
                    indstmp = AllIsTarg==1 & AllMotifNum==targmotifs(mmm);
                    pcol = pcols{mmm};
                    lt_plot_histogram(AllCC(indstmp), ...
                        Xcenters, 1, 0.4, 0.5, 1, pcol);
                    lt_plot_annotation(2, ['targ ' motiflist{targmotifs(mmm)}], pcol);
                end
                
                % --- histogram of same types
                indstmp = AllIsSame==1 & AllIsTarg==0;
                lt_plot_histogram(AllCC(indstmp), ...
                    Xcenters, 1, 0.4, 0.5, 1, [0.3 0.3 0.7]);
            end
            
            
        end
    end
end
