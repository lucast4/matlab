%% ==== 1) EXTRACT MOTIFSTATS POP

% ======================== EXTRACT SEGMENTS FOR POPULATIONS

close all; clear MOTIFSTATS_Compiled;
collectWNhit=0;
onlyCollectTargSyl=0;
LearnKeepOnlyBase = 0;
saveOn = 1;
OrganizeByExpt =0;
collectFF=1;
Params_regexp.motif_predur = [];
Params_regexp.motif_postdur = [];
Params_regexp.preAndPostDurRelSameTimept = 1;
Params_regexp.RemoveIfTooLongGapDur = [];
Params_regexp.extractDirSong = 1;

MOTIFSTATS_Compiled = lt_neural_v2_ANALY_MultExtractMotif(SummaryStruct, ...
    collectWNhit, LearnKeepOnlyBase, saveOn, onlyCollectTargSyl, OrganizeByExpt,...
    collectFF, [], Params_regexp);


%% ==== MAKE SURE SUMMARY STRCUT IS THE CORRECT ONE
SummaryStruct = MOTIFSTATS_Compiled.SummaryStruct;

%% ==== REMOVE DIR SONG
MOTIFSTATS_Compiled = lt_neural_QUICK_MotCom_RemoveDIR(MOTIFSTATS_Compiled);

%% ================ extraction continued
close all;
MOTIFSTATS_pop = lt_neural_v2_POP_ExtractMotifs(MOTIFSTATS_Compiled, SummaryStruct);
% clear MOTIFSTATS_Compiled;


%% #######################################################################
%% ############################# OLD STUFF
%% ================ PLOT [CORRELATION WITH FF]
close all;
xcov_dattotake = [-0.01 0.05];
xcov_dattotake = [-0.08 0.04];
xcov_dattotake = [-0.075 0.025];
xcovwindmax = 0.04;
binsize_spk = 0.0025;

MOTIFSTATS_pop = lt_neural_POP_ExtractXCov(MOTIFSTATS_pop, SummaryStruct, ...
    xcov_dattotake, xcovwindmax, binsize_spk);


%% #######################################################################
%% ############################# LEARNING STUFF;
%% ==== 2) EXTRACT LEARNING SWITCH STRUCT

SwitchStruct = lt_neural_LEARN_getswitch(SummaryStruct);


%% =========== SUMMARIZE LEARNING TRAJECTORY (PLUS NEURON SETS)
close all;

% BirdExptPairsToPlot = {'wh44wh39', 'RALMANlearn1'};
% motiftoplot = 'c(b)';
BirdExptPairsToPlot = {'pu69wh78', 'RALMANOvernightLearn1'};
motiftoplot = 'aa(b)';

lt_neural_POPLEARN_PlotLearnTraj(MOTIFSTATS_pop ,SwitchStruct, ...
    SummaryStruct, BirdExptPairsToPlot, motiftoplot);


%% #######################################################################
%% ############################# CROSS CORR CHANGE DURING LAERNING
%% ================ PLOT CROSS CORR WRT TO LEARNING
close all;
BirdExptPairsToPlot = {'wh44wh39', 'RALMANlearn1'};
% BirdExptPairsToPlot = {'wh44wh39', 'RALMANlearn2'};
SwitchToPlot = [2];
BregionWantedList = {{'LMAN', 'RA'}};
onlyPlotIfBothPrePostTrials = 0;
lt_neural_POPLEARN_Plot(MOTIFSTATS_pop, SwitchStruct, BirdExptPairsToPlot, ...
    SwitchToPlot, BregionWantedList, onlyPlotIfBothPrePostTrials);


%% ================ SUMMARIZE CROSS CORRELATION OVER COURSE OF EXPERIMENT
% over multiple switches
close all;
exptnum = [1];
birdnum = [2];
BregionWantedList = {{'LMAN', 'RA'}};
[OUTSTRUCT, birdnum] = lt_neural_POPLEARN_Summary(MOTIFSTATS_pop, SwitchStruct, ...
    birdnum, exptnum, BregionWantedList);


%% ======= PLOT SUMMARY OF XCOV traces
close all;
lt_neural_POPLEARN_SummaryPlot1(OUTSTRUCT, MOTIFSTATS_pop, SwitchStruct, ...
    birdnum, exptnum)

%% ======= PLOT MEAN XCOV
close all;
windowmean = [-0.05 0.02]; % window, in ms, relative to lag = 0;
lt_neural_POPLEARN_SummaryPlot2(OUTSTRUCT, MOTIFSTATS_pop, SwitchStruct, ...
    birdnum, windowmean)



%% ======== SUMMARIZE OVER ALL EXPERIMENTS
windowmean = [-0.05 0.02]; % window, in s, relative to lag = 0;
SkipIfTargsDiffSyls = 1; % skips switches where targets are different syl types.

lt_neural_POPLEARN_SummaryPlot3(OUTSTRUCT, MOTIFSTATS_pop, SwitchStruct, ...
    birdnum, windowmean, SkipIfTargsDiffSyls)




%% [GOOD] ################ FOR EACH BIRD, SUMMARIZE ACROSS ALL EXPTS
% DONE: removed dir inds
% TO DO: 1) baseline, 2) use early or late period in epoch

% INDICATE BY HAND WHICH NEURON SETS ARE APPROPRIATE
lt_neural_POPLEARN_SumTraj_Input; % go in here and select which dataset

% --- RUN
close all;
bregionwanted = {'LMAN', 'RA'};
lt_neural_POPLEARN_SumTraj(MOTIFSTATS_pop, SwitchStruct, ...
    metadatstruct, bregionwanted);

%% ================ PLOT PAIRED RASTERS WRT TO LEARNING
% [IN PROGRESS!!!]
lt_neural_POPLEARN_PairRast




%% ###################################################################
%% ############################################ COHERENCE
% OVERALL: for each motif, look at coherence of raw data, aligned to syl
% onset. Does that change during learning?

% NOTE: to see some old progress, see:
lt_neural_MasterScript_Pop;



%% ================ EXTRACT COHSTRUCT

COHSTRUCT = lt_neural_Coher_ExtrCohStruct(MOTIFSTATS_pop, SummaryStruct);


% ============== SAVE COHERENCE
if (0)
    marker = '03Oct2018_0040';
    fname = ['/bluejay5/lucas/analyses/neural/COHERENCE/COHSTRUCT_' marker '.mat'];
    save(fname, 'COHSTRUCT');
end



%% ############### SUMMARY PLOT OF COHERANCE
% PLOT COHGRAM AND FF BANDS FOR EVERY CHANNEL PAIR AND MOTIF

close all;
lt_neural_Coher_PlotRawSum(COHSTRUCT, MOTIFSTATS_pop);

%% ################ NOTE DOWN BRAIN REGION PAIRS

numbirds = length(COHSTRUCT.bird);
for i=1:numbirds
    numexpt = length(COHSTRUCT.bird(i).experiment);
    for ii=1:numexpt
        numsets = length(COHSTRUCT.bird(i).experiment(ii).setnum);
        for ss=1:numsets
            
            nummotifs = length(COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif);
            
            % ====== get list of neurons, bregions, and chans for this
            % dataset
            Neurlist = MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons{ss};
            Chanlist = [SummaryStruct.birds(i).neurons(Neurlist).channel];
            Bregionlist = {SummaryStruct.birds(i).neurons(Neurlist).NOTE_Location};
            
            for mm=1:nummotifs
                Chanpairs = COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif(mm).Chanpairs;
                numchanpairs = size(Chanpairs, 1);
                
                bregionpairs_sorted = {};
                for cc=1:numchanpairs
                    chansthis = Chanpairs(cc,:);
                    bregionsthis = {Bregionlist{Chanlist==chansthis(1)}, ...
                        Bregionlist{Chanlist==chansthis(2)}}; assert(length(bregionsthis)==2);
                    bregionsthis = sort(bregionsthis);
                    
                    bregionpairs_sorted = [bregionpairs_sorted; ...
                        [bregionsthis{1} '-' bregionsthis{2}]];
                end
                
                COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif(mm).bregionpairs_sorted = bregionpairs_sorted;
            end
        end
    end
end



%% ################ PLOT SUMMARY - FIRST EXTRACT MEAN COHGRAMS FOR ALL

All_CohgramMean = {};
All_birdnum = [];
All_enum = [];
All_setnum = [];
All_motifname = {};
All_chanpair = [];
All_bregionpair = {};
All_bregionpair_alphaorder = {};

numbirds = length(COHSTRUCT.bird);
for i=1:numbirds
    numexpt = length(COHSTRUCT.bird(i).experiment);
    birdname = SummaryStruct.birds(i).birdname;
    
    for ii=1:numexpt
        exptid = MOTIFSTATS_pop.birds(i).exptnum(ii).exptname;
        numsets = length(COHSTRUCT.bird(i).experiment(ii).setnum);
        
        for ss=1:numsets
            
            nummotifs = length(COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif);
            
            % ====== get list of neurons, bregions, and chans for this
            % dataset
            Neurlist = MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons{ss};
            Chanlist = [SummaryStruct.birds(i).neurons(Neurlist).channel];
            Bregionlist = {SummaryStruct.birds(i).neurons(Neurlist).NOTE_Location};
            
            for mm=1:nummotifs
                motifname = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(ss).motif(mm).regexpstr;
                
                CohCell = COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif(mm).Coh_ChpairByTrial;
                Chanpairs = COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif(mm).Chanpairs;
                tbins = COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif(mm).t_relons;
                ffbins = COHSTRUCT.bird(i).experiment(ii).setnum(ss).motif(mm).ffbins;
                
                numchanpairs = size(Chanpairs, 1);
                
                for cc=1:numchanpairs
                    chansthis = Chanpairs(cc,:);
                    bregionsthis = {Bregionlist{Chanlist==chansthis(1)}, ...
                        Bregionlist{Chanlist==chansthis(2)}}; assert(length(bregionsthis)==2);
                    
                    % =================== COLLECT COHEROGRAMS ACROSS TRIALS
                    % --------- 1) COLLECT COHEROGRAMS FOR THIS PAIR
                    %                     dim1 = size(CohCell{1},1);
                    %                     dim2 = size(CohCell{1},2);
                    %
                    % --------------- CONVERT CELLS TO A MATRIC
                    % use for loop because soemtime different size..
                    ntrials = size(CohCell,2);
                    cohmat = nan(length(tbins), length(ffbins), ntrials);
                    for tt=1:ntrials
                        if all(size(cohmat(:,:,tt))==size(CohCell{cc,tt}))
                            cohmat(:,:,tt) = CohCell{cc,tt};
                        else
                            disp('skipped! wrong size');
                        end
                    end
                    
                    % =========== get mean coherogram
                    cohmean = nanmean(cohmat, 3);
                    
                    %                     cohmat = reshape(cell2mat(CohCell(cc, :)), dim1, dim2, []);
                    
                    
                    % ========================== SAVE OUTPUT
                    All_CohgramMean = [All_CohgramMean; cohmean];
                    All_birdnum = [All_birdnum; i];
                    All_enum = [All_enum ; ii];
                    All_setnum = [All_setnum ; ss];
                    All_motifname = [All_motifname; motifname];
                    All_chanpair = [All_chanpair; chansthis];
                    All_bregionpair = [All_bregionpair; bregionsthis];
                    
                    bregionsthis_sort = sort(bregionsthis);
                    All_bregionpair_alphaorder = [All_bregionpair_alphaorder; ...
                        [bregionsthis_sort{1} '-' bregionsthis_sort{2}]];
                end
            end
        end
    end
end


%% ======= PLOT SUMMARY
close all;

% ============ 1) Lock to syllable onset, one plot for each brain region
% pair
pairthis = 'LMAN-RA';
lt_figure; hold on;
indsthis = strcmp(All_bregionpair_alphaorder, pairthis);
cohmat = lt_neural_Coher_Cell2Mat(All_CohgramMean(indsthis));

cohmean = nanmean(cohmat,3);
imagesc(tbins, ffbins, cohmean');
figure; plot(cohmean(:,3), 'ok-')



% ============ 2) ONE PLOT FOR EACH MOTIF AND PAIR TYPE
% --------------- INPUTS
ffbinsedges = [15 30 80 150]; % edges, to plot timecourse in frequency bands
pcols = lt_make_plot_colors(length(ffbinsedges)-1, 1, [1 0 0]);
figcount=1;
subplotrows=6;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% ------------- figure out ff bins

maxbird = max(All_birdnum);
bregionlist = unique(All_bregionpair_alphaorder);
for i=1:maxbird
    
    motiflist = unique(All_motifname(All_birdnum==i));
    birdname = SummaryStruct.birds(i).birdname;
    
    for b=1:length(bregionlist)
        bregionthis = bregionlist{b};
        
        % ==== go thru all motifs.
        for j=1:length(motiflist)
            motifthis = motiflist{j};
            
            indsthis = All_birdnum==i & strcmp(All_motifname, motifthis) ...
                & strcmp(All_bregionpair_alphaorder, bregionthis);
            
            if ~any(indsthis)
                continue
            end
            
            % --------- get mean cohgram
            cohmat = lt_neural_Coher_Cell2Mat(All_CohgramMean(indsthis));
            
            % ========== PLOT MEAN COHGRAM
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([motifthis ',' bregionthis]);
            ylabel(birdname);
            
            cohmean = nanmean(cohmat,3);
            imagesc(tbins, ffbins, cohmean');
            axis tight;
            line([0 0], ylim, 'Color', 'k');
            
            % ========= PLOT MEAN TRACE IN MULTIPLE FF BINS
            %             for k=1:size(cohmat,2)
            %                 % ff bins
            %                cohthis = squeeze(cohmat(:, k, :));
            %                cohmean = nanmean(cohthis,2);
            %                cohsem = lt_sem(cohthis');
            %
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([motifthis ',' bregionthis]);
            
            for k=1:length(ffbinsedges)-1
                
                indsff = ffbins>ffbinsedges(k) & ffbins<=ffbinsedges(k+1);
                ffmean = mean(ffbins(indsff));
                % ff bins
                cohthis = squeeze(nanmean(cohmat(:, indsff, :), 2)); % first take mean over the ff bins
                %                cohthis = squeeze(cohmat(:, indsff, :));
                cohmean = nanmean(cohthis,2); % then take mean across trials
                cohsem = lt_sem(cohthis');
                if length(cohsem)==1
                    plot(tbins, cohmean, 'Color', pcols{k});
                else
                    shadedErrorBar(tbins, cohmean, cohsem, {'Color', pcols{k}}, 1);
                end
                lt_plot_text(tbins(end), cohmean(end), [num2str(ffmean)], pcols{k}, 10);
            end
            axis tight;
            ylim([0.2 0.8]);
            line([0 0], ylim, 'Color', 'k');
            
        end
    end
end


%% =================== FOR EACH MOTIF AND BRAIN REGION, ONE PLOT
close all;

MotifLists = {...
    {'pu69wh78', {'(a)ab', 'a(a)b', 'aa(b)', 'aab(h)', 'aabh(h)'}}, ...
    {'pu69wh78', {'(j)jb', 'j(j)b', 'jj(b)', 'jjb(h)', 'jjbh(h)', 'jjbhh(g)'}}, ...
    {'pu69wh78', {'(j)jbhh', 'j(j)bhh', 'jj(b)hh', 'jjb(h)h'}}, ...
    {'pu69wh78', {'h(g)'}}, ...
    {'wh44wh39', {'(j)n', '(n)hh', 'n(h)h', 'nh(h)'}}, ...
    {'wh44wh39', {'(m)d', '(d)kcc', 'd(k)cc', 'dk(c)c', 'dkc(c)', 'c(b)', 'cb(b)'}}, ...
    {'wh44wh39', {'(n)h', 'n(h)'}}, ...
    {'wh44wh39', {'(d)kc', 'd(k)c', 'dk(c)'}}, ...
    };
BregionPair = 'LMAN-RA';

% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

for k=1:length(MotifLists)
    
    bthis = MotifLists{k}{1};
    bnumthis = find(strcmp({MOTIFSTATS_Compiled.birds.birdname}, bthis));
    motiflistthis = MotifLists{k}{2};
    
    
    
    
    figcount=1;
    subplotrows=2;
    subplotcols=length(motiflistthis);
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    % ################### PLOT COHEROGRAMS
    % ================ plot each motif one by one
    for m=1:length(motiflistthis)
        motifthis = motiflistthis{m};
        indthis = All_birdnum==bnumthis & strcmp(All_bregionpair_alphaorder, BregionPair) & ...
            strcmp(All_motifname, motifthis);
        if ~any(indthis)
            continue
        end
        % ============ extract coherence
        cohmat = lt_neural_Coher_Cell2Mat(All_CohgramMean(indthis));
        
        % ################## PLOT
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(motifthis);
        ylabel(bthis);
        
        cohmean = nanmean(cohmat,3);
        imagesc(tbins, ffbins, cohmean');
        axis tight;
        line([0 0], ylim, 'Color', 'k');
        
    end
    
    % ############# PLOT SMOOTHED POWER
    % ================ plot each motif one by one
    for m=1:length(motiflistthis)
        motifthis = motiflistthis{m};
        indthis = All_birdnum==bnumthis & strcmp(All_bregionpair_alphaorder, BregionPair) & ...
            strcmp(All_motifname, motifthis);
        if ~any(indthis)
            continue
        end
        
        % ============ extract coherence
        cohmat = lt_neural_Coher_Cell2Mat(All_CohgramMean(indthis));
        
        % ################## PLOT
        % ========= PLOT MEAN TRACE IN MULTIPLE FF BINS
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(motifthis);
        ylabel(bthis);
        for kk=1:length(ffbinsedges)-1
            
            indsff = ffbins>ffbinsedges(kk) & ffbins<=ffbinsedges(kk+1);
            ffmean = mean(ffbins(indsff));
            % ff bins
            cohthis = squeeze(nanmean(cohmat(:, indsff, :), 2)); % first take mean over the ff bins
            %                cohthis = squeeze(cohmat(:, indsff, :));
            cohmean = nanmean(cohthis,2); % then take mean across trials
            cohsem = lt_sem(cohthis');
            if length(cohsem)==1
                plot(tbins, cohmean, 'Color', pcols{kk});
            else
                shadedErrorBar(tbins, cohmean, cohsem, {'Color', pcols{kk}}, 1);
            end
            lt_plot_text(tbins(end), cohmean(end), [num2str(ffmean)], pcols{kk}, 10);
        end
        axis tight;
        ylim([0.2 0.8]);
        line([0 0], ylim, 'Color', 'k');
    end
end



%% LEARNING [COLLECT, FIR EACH SIWTCH]
SwitchCohStruct = struct;
pairtoget = 'LMAN-RA';
numbirds = length(SwitchStruct.bird);
for i=1:numbirds
    numexpts = length(SwitchStruct.bird(i).exptnum);
    for ii=1:numexpts
        numswitch = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        for ss=1:numswitch
            
            % =========== for this switch, get edge times and switch time
            tstart = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum_previous;
            tswitch = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum;
            tend = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum_next;
            
            % =========== Find channel pairs that have data both pre and post.
            numsets = length(MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons);
            setskept = [];
            for k=1:numsets
                %    COHSTRUCT.bird(i).experiment(ii).setnum(k).motif(1);
                nummotifs = length(MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(k).motif);
                
                
                for mm=1:nummotifs
                    if isempty(MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(k).motif(mm).SegExtr_neurfakeID)
                        continue
                    end
                    segextract = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(k).motif(mm).SegExtr_neurfakeID(1).SegmentsExtract;
                    
                    tvals = [segextract.song_datenum];
                    
                    % ---- ARE THERE BOTH PRE AND POST TVALS?
                    if ~any(tvals>tstart & tvals<tswitch) | ~any(tvals>tswitch & tvals<tend)
                        % then don't keep...
                        continue
                    end
                    
                    %         disp(k);
                    % =============== SAVE DATA FOR THIS MOTIF
                    cohdat = COHSTRUCT.bird(i).experiment(ii).setnum(k).motif(mm);
                    motifname = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(k).motif(mm).regexpstr;
                    
                    indtmp = strcmp(cohdat.bregionpairs_sorted, pairtoget);
                    if ~any(indtmp)
                        continue
                    end
                    cohmat = lt_neural_Coher_Cell2Mat(cohdat.Coh_ChpairByTrial(indtmp,:));
                    bregionpair = cohdat.bregionpairs_sorted(indtmp);
                    chanpair = cohdat.Chanpairs(indtmp,:);
                    tvals = [segextract.song_datenum];
                    ffvals = [segextract.FF_val];
                    assert(length(tvals)==size(cohmat,3));
                    
                    if any(isnan(cohmat(:)))
                        disp('NOTE: currenrtly skipping those with nans in coherence (becasue size of array is rong... should fix')
                        continue
                    end
                    % ====================== SAVE OUTPUT
                    SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).motifname = motifname;
                    SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).cohmat = cohmat;
                    SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).tvals = tvals;
                    SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).ffvals = ffvals;
                    SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).bregionpair = bregionpair;
                    SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).chanpair = chanpair;
                    setskept = [setskept k];
                end
            end
            % --- sanity check, confirm that for any switch at most only
            % one set of neurons is included.
            if ~isempty(setskept)
                assert(length(unique(setskept))==1, 'if this fails, then must combine across sets...');
            end
        end
    end
end

%% ################### LEARNING [PLOT EACH EXPERIMENT]
% ======= 1) EXTRACT
close all;
plotON=0;
OUTSTRUCT = lt_neural_Coher_Learn_Extr(SwitchStruct, SwitchCohStruct, plotON);

% ======== 2) CALCULATE TRAIN MINUS BASE COHERENCE
cohdiff_all = {};
% --- takes 2nd half of training and base
for j=1:length(OUTSTRUCT.bnum)
   
    cohmat = OUTSTRUCT.CohMat{j};
    indsbase = find(OUTSTRUCT.indsbase{j});
    indsWN = find(OUTSTRUCT.indsWN{j});
    
    % -------------- take second half of WN inds
    indsbase = indsbase(round(length(indsbase)/2):end);
    indsWN = indsWN(round(length(indsWN)/2):end);
    
    % ----- 
    cohmean_base = mean(cohmat(:,:, indsbase),3);
    cohmean_WN = mean(cohmat(:,:, indsWN),3);
    cohmean_diff = cohmean_WN - cohmean_base;
    
    cohdiff_all = [cohdiff_all; cohmean_diff];    
end
OUTSTRUCT.CohMean_WNminusBase = cohdiff_all;

%% ======= 2) PLOT
close all;
% ============== PARAMS
useabs = 1; % if 1, then absolute values (wn minus diff)
plotON =0;

% ===============
indsgrp = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch});
indsgrp_unique = unique(indsgrp)';


% ===============
figcount=1;
subplotrows=3;
subplotcols=6;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

clim = [-0.4 0.4];

% ==== to plot mean across all switches
Yallswitch = {};
for j=indsgrp_unique
    
    bname = SwitchStruct.bird(unique(OUTSTRUCT.bnum(indsgrp==j))).birdname;
    ename = SwitchStruct.bird(unique(OUTSTRUCT.bnum(indsgrp==j))).exptnum(unique(OUTSTRUCT.enum(indsgrp==j))).exptname;
    swnum = unique(OUTSTRUCT.switch(indsgrp==j));
    
    % ======================= TRAIN MINUS BASELINE COHERENCE
    Yall = {};
    % --- targets
    indsthis = OUTSTRUCT.istarg==1 & indsgrp==j;
    cohmat = lt_neural_Coher_Cell2Mat(OUTSTRUCT.CohMean_WNminusBase(indsthis));
    if useabs==1
        cohmat = abs(cohmat);
    end
    cohmean = mean(cohmat,3);
    Yall = [Yall cohmean];
    
    % --- nontarg
    indsthis = OUTSTRUCT.istarg==0 & indsgrp==j;
    cohmat = lt_neural_Coher_Cell2Mat(OUTSTRUCT.CohMean_WNminusBase(indsthis));
    if useabs==1
        cohmat = abs(cohmat);
    end
    cohmean = mean(cohmat,3);
    Yall = [Yall cohmean];
    
    % --------------------- targ minus nontarg
    cohdiff = Yall{1} - Yall{2};
    
    % ================= PLOT
    if plotON==1
        % ----- 1) TARG (WN minsu base)
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        ylabel('targ(WN-base)');
        title([bname '-' ename '-sw' num2str(swnum)]);
        lt_neural_Coher_Plot(Yall{1}, tbins, ffbins, 1, '-', clim);
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        ylabel('targ(WN-base)');
        lt_neural_Coher_Plot(Yall{1}, tbins, ffbins, 2, '-', clim);
        lt_plot_zeroline;
        
        % ----- 2) NONTARG (WN minsu base)
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        ylabel('nontarg(WN-base)');
        lt_neural_Coher_Plot(Yall{2}, tbins, ffbins, 1, '-', clim);
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        ylabel('nontarg(WN-base)');
        lt_neural_Coher_Plot(Yall{2}, tbins, ffbins, 2, '-', clim);
        lt_plot_zeroline;
        
        % ----- 3) TARG-NONTARG (WN minsu base)
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        lt_neural_Coher_Plot(cohdiff, tbins, ffbins, 1, '-', clim);
        ylabel('targ(WN-base) - nontarg(WN-base)');
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        lt_neural_Coher_Plot(cohdiff, tbins, ffbins, 2, '-', clim);
        lt_plot_zeroline;
    end
    
    % ============================= COLLECT
    Yallswitch = [Yallswitch; Yall];
end

% ======================= PLOT MEANS
lt_figure; hold on;

% --------- 1) mean heat maps
cohmat = lt_neural_Coher_Cell2Mat(Yallswitch(:,1));

lt_subplot(3,2,1); hold on;
title('mean over switch');
ylabel('TARG(WN-base)');
lt_neural_Coher_Plot(mean(cohmat,3), tbins, ffbins, 1, '-', []);

lt_subplot(3,2,2); hold on;
% title('mean over switch');
ylabel('TARG');
lt_neural_Coher_Plot(mean(cohmat,3), tbins, ffbins, 2, '-', clim);
lt_plot_zeroline;


% --------- 1) mean heat maps
cohmat = lt_neural_Coher_Cell2Mat(Yallswitch(:,2));

lt_subplot(3,2,3); hold on;
title('mean over switch');
ylabel('NONTARG(WN-base)');
lt_neural_Coher_Plot(mean(cohmat,3), tbins, ffbins, 1, '-', []);

lt_subplot(3,2,4); hold on;
% title('mean over switch');
ylabel('NONTARG');
lt_neural_Coher_Plot(mean(cohmat,3), tbins, ffbins, 2, '-', clim);
lt_plot_zeroline;


% ------------
cohmat1 = lt_neural_Coher_Cell2Mat(Yallswitch(:,1)); 
cohmat2 = lt_neural_Coher_Cell2Mat(Yallswitch(:,2));
cohmat = cohmat1 - cohmat2;

lt_subplot(3,2,5); hold on;
title('mean over switch');
ylabel('TARG-NONTARG(WN-base)');
lt_neural_Coher_Plot(mean(cohmat,3), tbins, ffbins, 1, '-', []);

lt_subplot(3,2,6); hold on;
% title('mean over switch');
ylabel('TARG-NONTARG(WN-base)');
lt_neural_Coher_Plot(mean(cohmat,3), tbins, ffbins, 2, '-', clim);
lt_plot_zeroline;

%% ##################
%% ################### LEARNING, FOR EACH BIRD, SUMMARIZE ACROSS ALL EXPTS
% DONE: removed dir inds
% TO DO: 1) baseline, 2) use early or late period in epoch

% INDICATE BY HAND WHICH NEURON SETS ARE APPROPRIATE
lt_neural_POPLEARN_SumTraj_Input; % go in here and select which dataset

% --- RUN
close all;
bregionwanted = {'LMAN', 'RA'};
lt_neural_POPLEARN_SumTraj(MOTIFSTATS_pop, SwitchStruct, ...
    metadatstruct, bregionwanted);



%% plot single cases
if (0)
    lt_figure; hold on;
    
    % === 1) plot multiple trials
    
    
    % === 2) plot average over all trials
    lt_subplot(3,2,2); hold on;
    
    chanpair = 1;
    dim1 = size(CohAllTrials{1},1);
    dim2 = size(CohAllTrials{1},2);
    cohmat = CohAllTrials(chanpair, :);
    
    % -- convert to 3d mat
    cohmat = reshape(cell2mat(cohmat), dim1, dim2, []);
    cohmean = mean(cohmat, 3);
    imagesc(t, ff, 10*log10(cohmean'));
    
    % ==== 3) plot specific ff
    lt_subplot(3,2,3); hold on;
    plot(cohmean(:, 3), 'ok')
    
    
    % ############################### PLOT EACH FREQUENCY BIN SEPARATELY
    figcount=1;
    subplotrows=5;
    subplotcols=3;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    
    for i=1:length(ff)
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        fthis = ff(i);
        title(['ffbin:' num2str(fthis)]);
        hsplots = [hsplots hsplot];
        
        cohthis = squeeze(cohmat(:, i, :));
        ymean = mean(cohthis,2);
        ysem = lt_sem(cohthis');
        x = 1:length(ymean);
        lt_plot(x, ymean, {'Errors', ysem, 'Color', 'k'});
    end
    linkaxes(hsplots, 'xy');
    
    % -- PLOT EACH FF AS A LINE WITH ERROR BARS
end

%% ############### SUMMARY PLOT - FIRST EXTRACT MEAN COHEROGRAM FOR EACH DATAPT

%% ===== convert from segextract to a specific song and time in song

lt_neural_v2_EXTRACT_WithinSongTimings(SummaryStruct, i, neurthis);


