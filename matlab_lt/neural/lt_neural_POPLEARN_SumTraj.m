function lt_neural_POPLEARN_SumTraj(MOTIFSTATS_pop, SwitchStruct, ...
    metadatstruct, bregionwanted)
%%

xwindow_mean = [-0.021 0.01]; % lags to take average over. (negative means LMAN leads)
% not inclusive


%% INPUTS

% e.g. metadatcell{1} = {birdname, exptname, targsyl, [baseline sets], ...
% {{setnum, syltoignore}}, [trainingsets], {{setnum, syltoignore}}};

% e.g. 
% metadatcell{1} = {'wh44wh39', 'RALMANlearn3', 'nh(h)', [3], {{}}, ...
%     [4 5 6 7 8], {{}, {}, {}, {}, {}}};


% changed to metadatstruct format. e.g.:
% ind = 1;
% metadatstruct(ind).bname = 'wh44wh39';
% metadatstruct(ind).exptname = 'RALMANlearn3';
% metadatstruct(ind).targsyl = 'nh(h)';
% metadatstruct(ind).base_set = [3];
% metadatstruct(ind).base_sylsignore = {{}};
% metadatstruct(ind).train_set = [4 5 6 7 8];
% metadatstruct(ind).train_sylsignore = {{}, {}, {}, {}, {}};
% 
%% lt 3/31/18 - hand code trajectory of a given syllable over the course of an expeirment
% useful ebcause of many switches etc.


% ----------- to collect across all experiments/birds
% [one ind for each neuron pair]
AllAll_CCtarg_z = [];
AllAll_CCtarg_z_mean = [];
AllAll_TargLearnZ = [];

for j=1:length(metadatstruct)
    
%     birdthis = metadatcell{j}{1};
%     exptthis = metadatcell{j}{2};
%     targsylthis = metadatcell{j}{3};
%     sets_base = metadatcell{j}{4};
%     sylstoignore_base = metadatcell{j}{5};
%     sets_train = metadatcell{j}{6};
%     sylstoignore_train = metadatcell{j}{7};
 
    birdthis = metadatstruct(j).bname;
    exptthis = metadatstruct(j).exptname;
    targsylthis = metadatstruct(j).targsyl;
    sets_base = metadatstruct(j).base_set;
    switches_base = metadatstruct(j).base_switches;
    sylstoignore_base = metadatstruct(j).base_sylsignore;
    sets_train = metadatstruct(j).train_set;
    switches_train = metadatstruct(j).train_switches;
    sylstoignore_train = metadatstruct(j).train_sylsignore;
    targlearndir = metadatstruct(j).targlearndir;
    
    i = find(strcmp({MOTIFSTATS_pop.birds.birdname}, birdthis));
    ii = find(strcmp({MOTIFSTATS_pop.birds(i).exptnum.exptname}, exptthis));
    
    DAT = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT;
    
    % ===== global motif list
    motiflist = {MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(1).motif.regexpstr};
    nummotifs = length(motiflist);
    
    
    
    %% =========== get distribution of FF at baseline
    switchconting_this = switches_base;
    bb = 1; % ITERATE, base set num
    mm = find(strcmp(motiflist, targsylthis));
    ss = sets_base(bb);
    ff = [DAT.setnum(ss).motif(mm).SegExtr_neurfakeID(1).SegmentsExtract.FF_val];

    % ================= GET CORRECT TRIALS
    % i) remove DIR if exist; ii) if baseline, then make sure if
    % pre WN onset
    dirinds = [DAT.setnum(ss).motif(mm).SegExtr_neurfakeID(1).SegmentsExtract.DirSong];
    goodtrials = ~dirinds;
    
    % ii) --- if desired, make sure is on one or other side of a
    % switch.
    if ~isempty(switchconting_this{bb})
        % then need to take data on either before or after a switch
        swnumthis = str2num(switchconting_this{bb}(1));
        swtiming = switchconting_this{bb}(2);
        tvals = [DAT.setnum(ss).motif(mm).SegExtr_neurfakeID(1).SegmentsExtract.song_datenum];
        
        if swtiming=='l'
            % -- then get data before switch
            t1 = SwitchStruct.bird(i).exptnum(ii).switchlist(swnumthis).switchdnum_previous;
            t2 = SwitchStruct.bird(i).exptnum(ii).switchlist(swnumthis).switchdnum;
        elseif swtiming =='r'
            % -- then get data after switch
            t1 = SwitchStruct.bird(i).exptnum(ii).switchlist(swnumthis).switchdnum;
            t2 = SwitchStruct.bird(i).exptnum(ii).switchlist(swnumthis).switchdnum_next;
        else
            % then PROBLEM
            asdfasdfasdfsdf;
        end
        
        indtiming = tvals>t1 & tvals<t2;
        assert(any(indtiming), 'mistake in hand entering swithc, no data ...');
        goodtrials = goodtrials & indtiming;
    end
    
%     tvals = tvals(goodtrials);
    ff_base = ff(goodtrials);
    
    %% ################## TRAINING
    % ###### go thru all desired sets and collect CC, separated by motif classes
    setsthis = sets_train;
    motifstoignore = sylstoignore_train;
    assert(length(motifstoignore) == length(setsthis), 'asdfasdf');
    switchconting_this = switches_train;
    assert(length(switchconting_this) == length(setsthis), 'asdfasdf');
    
    for ss = setsthis
        assert(all(strcmp(motiflist, ...
            {MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(ss).motif.regexpstr})), 'asfda');
    end
    
    AllCC = [];
    AllMotifID = [];
    AllSetnum = [];
    AllNeurPairID = [];
    AllIsTarg = [];
    AllIsSame = [];
    AllTargLearning = []; % z score, adaptive direction
    % ======== go thru all sets
    for k = 1:length(setsthis)
        ss = setsthis(k);
        motifignore = motifstoignore{k};
        
        datsizeall = []; % sanity check, make sure (neur1, neur2) all same sizes
        for mm=1:nummotifs
            
            % ================ skip this motif?
            if any(strcmp(motifignore, motiflist{mm}))
                disp(['ignored ' motiflist{mm}]);
                continue
            end
            
            % =================== COLLECT DATA
            XcovDat = DAT.setnum(ss).motif(mm).XCov_neurpair;
            if isempty(XcovDat)
                continue
            end
            
            datsizeall = [datsizeall; size(XcovDat)];
            
            % ====== extract xcov for all pairs (across trials)
            % automatically flips to match desired order
            [ccRealAll, ccShiftAll, xlags_sec, indspairID] = ...
                lt_neural_POPLEARN_FlipXcovDat(XcovDat, bregionwanted);
            
            
            % ================= GET CORRECT TRIALS
            % i) remove DIR if exist; ii) if baseline, then make sure if
            % pre WN onset
            dirinds = [DAT.setnum(ss).motif(mm).SegExtr_neurfakeID(1).SegmentsExtract.DirSong];
            assert(length(dirinds) == size(ccRealAll{1},1), 'asdfd');
            goodtrials = ~dirinds;
            
            % ii) --- if desired, make sure is on one or other side of a
            % switch.
            if ~isempty(switchconting_this{k})
                % then need to take data on either before or after a switch
               swnumthis = str2num(switchconting_this{k}(1));
               swtiming = switchconting_this{k}(2);
               tvals = [DAT.setnum(ss).motif(mm).SegExtr_neurfakeID(1).SegmentsExtract.song_datenum];
               
               if swtiming=='l'
                   % -- then get data before switch
                   t1 = SwitchStruct.bird(i).exptnum(ii).switchlist(swnumthis).switchdnum_previous;
                   t2 = SwitchStruct.bird(i).exptnum(ii).switchlist(swnumthis).switchdnum;
               elseif swtiming =='r'
                   % -- then get data after switch
                   t1 = SwitchStruct.bird(i).exptnum(ii).switchlist(swnumthis).switchdnum;
                   t2 = SwitchStruct.bird(i).exptnum(ii).switchlist(swnumthis).switchdnum_next;
               else
                   % then PROBLEM
                   asdfasdfasdfsdf;
               end
               
               indtiming = tvals>t1 & tvals<t2;
               assert(any(indtiming), 'mistake in hand entering swithc, no data ...');
               goodtrials = goodtrials & indtiming;
            end
            
            
            
            
            % ======= get one mean, shuffle corrected, xcov for each pair
            npairs = length(ccRealAll);
            CCall = [];
            for np =1:npairs
                cc = mean(ccRealAll{np}(goodtrials,:) - ccShiftAll{np}(goodtrials,:) ,1);
                CCall = [CCall; cc];
            end
            
            
            % ======== GET LEARNING FOR TARGET SYL
            indtmp = find(strcmp(motiflist, targsylthis));
            ff_targ = [DAT.setnum(ss).motif(indtmp).SegExtr_neurfakeID(1).SegmentsExtract.FF_val];
            tvals = [DAT.setnum(ss).motif(indtmp).SegExtr_neurfakeID(1).SegmentsExtract.song_datenum];
            isdir = [DAT.setnum(ss).motif(indtmp).SegExtr_neurfakeID(1).SegmentsExtract.DirSong];
            
            % -- restrict to good trials
            if ~isempty(switchconting_this{k})
                indtrialtmp = isdir==0 & tvals>t1 & tvals<t2;
            else
                indtrialtmp = isdir==0;
            end
            
            ff_targ = ff_targ(indtrialtmp);
            % -- convert to z relative to baseline
            ffz = (mean(ff_targ) - mean(ff_base))./std(ff_base);
            ffz = ffz*targlearndir; % in adaptive direction.

            
            % ################## COLLECT, ACROSS ALL MOTIFS
            AllCC = [AllCC; CCall];
            AllMotifID = [AllMotifID; mm*ones(size(CCall,1),1)];
            AllSetnum = [AllSetnum; ss*ones(size(CCall,1),1)];
            AllNeurPairID = [AllNeurPairID; indspairID'];
            
            % --- what is the class of this motif?
            istarg = strcmp(motiflist{mm}, targsylthis);
            issame = lt_neural_QUICK_SylPairType(motiflist{mm}, targsylthis);
            
            AllIsTarg = [AllIsTarg; istarg*ones(size(CCall,1),1)];
            AllIsSame = [AllIsSame; issame*ones(size(CCall,1),1)];
            AllTargLearning = [AllTargLearning; ffz*ones(size(CCall,1),1)];
            
        end
        tmp = diff(datsizeall);
        assert(all(tmp(:) ==0), 'neur pair id mght not be correct ...');
    end
    
    %% == note: for baseline same thing, but make sure not during WN
    
    %% ############### PLOT
    lt_figure; hold on;
    hsplots = [];
    
    % ======= TARGET
    hsplot = lt_subplot(3,2,1); hold on;
    hsplots = [hsplots hsplot];
    title('targ');
    motifind= find(strcmp(motiflist, targsylthis));
    
    % ----
    indstoplot = AllMotifID == motifind;
    ccmat = AllCC(indstoplot, :);
    plot(xlags_sec, ccmat, '-', 'Color', [0.7 0.7 0.7]);
    plot(xlags_sec, mean(ccmat,1), '-k', 'LineWidth', 3);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    % --- indicate magnitude of learning
    lt_plot_text(xlags_sec(1), 0.8, ...
        ['targ learn mean: ' num2str(mean(AllTargLearning(indstoplot)))], 'r');
    
    % ======= ALL NONTARG
    hsplot= lt_subplot(3,2,2); hold on;
        hsplots = [hsplots hsplot];
    title('all nontarg');
    motifind= find(strcmp(motiflist, targsylthis));
    % ----
    indstoplot = AllMotifID ~= motifind;
    ccmat = AllCC(indstoplot, :);
    plot(xlags_sec, ccmat, '-', 'Color', [0.7 0.7 0.7]);
    plot(xlags_sec, mean(ccmat,1), '-k', 'LineWidth', 3);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    
    % ======= SAME TYPE
    hsplot = lt_subplot(3,2,3); hold on;
        hsplots = [hsplots hsplot];
    title('same type');
    [sametypes] = lt_neural_QUICK_ExtractSameType(motiflist, targsylthis);
    motifinds = find(ismember(motiflist, sametypes));
    
    % ----
    indstoplot = ismember(AllMotifID, motifinds);
    ccmat = AllCC(indstoplot, :);
    if ~isempty(ccmat)
    plot(xlags_sec, ccmat, '-', 'Color', [0.7 0.7 0.7]);
    plot(xlags_sec, mean(ccmat,1), '-k', 'LineWidth', 3);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    end
    lt_plot_text(xlags_sec(1), 0.8, ...
        [birdthis '-' exptthis '-' targsylthis], 'r');
    
    % ======= DIFF TYPE
    hsplot = lt_subplot(3,2,4); hold on;
        hsplots = [hsplots hsplot];
    title('diff type');
    
    motifinds = find(~strcmp(motiflist, targsylthis) & ~ismember(motiflist, sametypes));
    
    % ----
    indstoplot = ismember(AllMotifID, motifinds);
    ccmat = AllCC(indstoplot, :);
    plot(xlags_sec, ccmat, '-', 'Color', [0.7 0.7 0.7]);
    plot(xlags_sec, mean(ccmat,1), '-k', 'LineWidth', 3);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    % ==================
    linkaxes(hsplots, 'xy');
    
   %% ====================== PLOT, SEPARATELY FOR EACH NEURAL SET
    
   % ===== for each set of neurons (individually for each pair of neurons
   % within that set) subtract diff-type xcov from the target xcov.
   lt_figure; hold on;
   maxsets = max(AllSetnum);
   maxneurpairs = max(AllNeurPairID);
   
   for ss = 1:maxsets
      for nn = 1:maxneurpairs
          
          % --- find all the motifs for this neurpair
          inds = AllSetnum==ss & AllNeurPairID==nn;
          if ~any(inds)
              continue
          end
          % -- sanity check - must be one and only one targ
          assert(sum(AllIsTarg(inds))==1, 'why not exactly one targ?');
          
          % ========== get mean CC across all nontargs
          inds = AllSetnum==ss & AllNeurPairID==nn & AllIsTarg==0 & AllIsSame==0;
          ccmean_nontarg = mean(AllCC(inds,:));
          ccstd_nontarg = std(AllCC(inds,:));
          
          % ====== convert the targ cc to zscore
          inds = AllSetnum==ss & AllNeurPairID==nn & AllIsTarg==1;
          cctarg = AllCC(inds,:);
          cctarg_z = (cctarg - ccmean_nontarg)./(ccstd_nontarg);
          
          plot(xlags_sec, cctarg_z, '-k');
          
          % =============== COLLECT ACROSS ALL EXPERIMENTS
          AllAll_CCtarg_z = [AllAll_CCtarg_z; cctarg_z];
          
          % =============== GET MEAN Z-SCORE IN NEGATIVE TIME WINDOW.
          indstmp = xlags_sec>xwindow_mean(1) & xlags_sec<xwindow_mean(2);
          ccmean = mean(cctarg_z(indstmp));
          AllAll_CCtarg_z_mean = [AllAll_CCtarg_z_mean; ccmean];
          
          % ============== learning at target
          inds = AllSetnum==ss & AllNeurPairID==nn;
          targlearn = unique(AllTargLearning(inds));
          assert(length(targlearn)==1, 'all shoudl have same targ ...');
          
          AllAll_TargLearnZ = [AllAll_TargLearnZ; targlearn];

      end       
   end    
    
    
    
end

%% ############### PLOTS, ACROSS ALL EXPERIMENTS
lt_figure; hold on;
title('mean across all neuron pairs');
ylabel('zscore, rel to nontarg syls');
plot(xlags_sec, AllAll_CCtarg_z, '-', 'Color', [0.8 0.8 0.8]);
ccmean = mean(AllAll_CCtarg_z,1);
ccsem = lt_sem(AllAll_CCtarg_z);
shadedErrorBar(xlags_sec, ccmean, ccsem, {'Color', 'r'},1);


%% ################ PLOTS
% ====== across neural sets, scatter of z(targ) vs. learning(targ)
lt_figure; hold on;
xlabel('learn at targ (z, adaptive)');
ylabel('mean xcov, LMANlead');
plot(AllAll_TargLearnZ, AllAll_CCtarg_z_mean, 'or');
















