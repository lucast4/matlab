%% ================ EXTRACT COHSTRUCT

COHSTRUCT = lt_neural_Coher_ExtrCohStruct(MOTIFSTATS_pop, SummaryStruct, ...
    MOTIFSTATS_Compiled);


% ============== SAVE COHERENCE
if (0)
    marker = '14Oct2018_2147_NewReExtract';
    fname = ['/bluejay5/lucas/analyses/neural/COHERENCE/COHSTRUCT_' marker '.mat'];
    save(fname, 'COHSTRUCT');
end

%%  FOR LOADING

load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_14Oct2018_2147.mat');
load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_14Oct2018_2147.mat');
load('/bluejay5/lucas/analyses/neural/COHERENCE/COHSTRUCT_14Oct2018_2147_NewReExtract.mat');


%% ################ NOTE DOWN BRAIN REGION PAIRS

COHSTRUCT = lt_neural_Coher_GetBrRegPairs(COHSTRUCT, MOTIFSTATS_pop, SummaryStruct);


%% ================= FIX PROBLEM THAT COHMATS HAVE DIFFERENT SIZES IN COHSTRUCT.
tlengthdesired = 20; % number of time bins desired - ad hoc to make sure all same dimension
fflengthdesired = 19;
COHSTRUCT = lt_neural_Coher_FixCohMatSize(COHSTRUCT, MOTIFSTATS_pop, SummaryStruct, ...
    tlengthdesired, fflengthdesired);

%% ############### SUMMARY PLOT OF COHERANCE
% PLOT COHGRAM AND FF BANDS FOR EVERY CHANNEL PAIR AND MOTIF

close all;
lt_neural_Coher_PlotRawSum(COHSTRUCT, MOTIFSTATS_pop);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUMMARY ANALYSES
%% ================ 1) EXTRACTION (EXTRACT MEAN COHGRAMS FOR ALL)
% tlengthdesired = 20; % number of time bins desired - ad hoc to make sure all same dimension
% fflengthdesired = 19;

tlengthdesired = length(COHSTRUCT.bird(1).experiment(1).setnum(2).motif(1).t_relons);
fflengthdesired = length(COHSTRUCT.bird(1).experiment(1).setnum(2).motif(1).ffbins);

%  Series of summary analyses - here just does extraction
[All_CohgramMean, All_birdnum, All_enum, All_setnum, All_motifname, ...
    All_chanpair, All_bregionpair, All_bregionpair_alphaorder, ...
    All_tbins, All_ffbins] = lt_neural_Coher_Summ_Extr(COHSTRUCT, ...
    MOTIFSTATS_pop, SummaryStruct, tlengthdesired, fflengthdesired, ...
    PARAMS);

%% ======= PLOT SUMMARY
close all;

pairthis = 'LMAN-RA';
ffbinsedges = [10 25 36 80 150]; % edges, to plot timecourse in frequency bands
lt_neural_Coher_Summ_Plot1(pairthis, All_CohgramMean, All_birdnum, All_enum, All_setnum, All_motifname, ...
    All_chanpair, All_bregionpair, All_bregionpair_alphaorder, ffbinsedges, ...
    PARAMS, SummaryStruct);

%% =================== FOR EACH MOTIF AND BRAIN REGION, ONE PLOT
close all;

MotifLists = {...
    {'pu69wh78', {'(a)ab', 'a(a)b', 'aa(b)', 'aab(h)', 'aabh(h)'}}, ...
    {'pu69wh78', {'(j)jb', 'j(j)b', 'jj(b)', 'jjb(h)', 'jjbh(h)', 'jjbhh(g)'}}, ...
    {'wh44wh39', {'(j)n', '(n)hh', 'n(h)h', 'nh(h)'}}, ...
    {'wh44wh39', {'(m)d', '(d)kcc', 'd(k)cc', 'dk(c)c', 'dkc(c)', 'c(b)', 'cb(b)'}}, ...
    };
% MotifLists = {...
%     {'pu69wh78', {'(a)ab', 'a(a)b', 'aa(b)', 'aab(h)', 'aabh(h)'}}, ...
%     {'pu69wh78', {'(j)jb', 'j(j)b', 'jj(b)', 'jjb(h)', 'jjbh(h)', 'jjbhh(g)'}}, ...
%     {'pu69wh78', {'(j)jbhh', 'j(j)bhh', 'jj(b)hh', 'jjb(h)h'}}, ...
%     {'pu69wh78', {'h(g)'}}, ...
%     {'wh44wh39', {'(j)n', '(n)hh', 'n(h)h', 'nh(h)'}}, ...
%     {'wh44wh39', {'(m)d', '(d)kcc', 'd(k)cc', 'dk(c)c', 'dkc(c)', 'c(b)', 'cb(b)'}}, ...
%     {'wh44wh39', {'(n)h', 'n(h)'}}, ...
%     {'wh44wh39', {'(d)kc', 'd(k)c', 'dk(c)'}}, ...
%     };
BregionPair = 'LMAN-RA';
% ffbinsedges = [10 25 36 80]; % edges, to plot timecourse in frequency bands
ffbinsedges = [10 25 32 80]; % edges, to plot timecourse in frequency bands
normtomeancoh=1;  % subtracts within-FF (acorss T) mean (i.e. to see modualtion.
lt_neural_Coher_Summ_Plot2(MotifLists, BregionPair, MOTIFSTATS_Compiled, ...
    All_CohgramMean, All_birdnum, All_enum, All_setnum, All_motifname, All_chanpair, ...
    All_bregionpair, All_bregionpair_alphaorder, PARAMS.tbins, PARAMS.ffbins, ...
    ffbinsedges, normtomeancoh);


%% ============================= NOTE DOWN TIME AND FF BINS

if exist('All_tbins', 'var')
    tbins = All_tbins{1};
    ffbins = All_ffbins{1};
    PARAMS.tbins = tbins;
    PARAMS.ffbins = ffbins;
else
    assert(~isempty(COHSTRUCT.bird(1).experiment(1).setnum(2).motif));
    PARAMS.tbins = COHSTRUCT.bird(1).experiment(1).setnum(2).motif(1).t_relons;
    PARAMS.ffbins = COHSTRUCT.bird(1).experiment(1).setnum(2).motif(1).ffbins;
end

%% #####################################################
%% ################################## LEARNING [COLLECT, FIR EACH SIWTCH]
% NOTE: this is not latest. for latest see one level above.

pairtoget = 'LMAN-RA';
SwitchCohStruct = lt_neural_Coher_LearnExtr(COHSTRUCT, MOTIFSTATS_pop, ...
    SwitchStruct, pairtoget, LFPSTRUCT);


%% =================== PLOT, ONE FIGURE FOR EACH CHANNEL PAIR
% ======= ie plots all cases.
% NOTE: OBSOLETE, SEE RAW PLOTS BELOW.
close all;
lt_neural_Coher_Learn_PlotAll;

%% ======= 1) EXTRACT
close all;
plotON=0; % ONE PLOT FOR EACH MOTIF/CHANPAIR [NOTE: use PlotChanMotif below instead]
averagechanpairs= 0; % for each motif, average over all chan pairs
onlyfirstswitch = 0;
removeBadSyls = 1; % LEAVE AT 1.
OUTSTRUCT = lt_neural_Coher_Learn_Extr(SwitchStruct, SwitchCohStruct, ...
    plotON, averagechanpairs, PARAMS, onlyfirstswitch, removeBadSyls);


%% ======= 2) GET SCALAR MEASURE OF CHANGE IN COHERENCE IN CERTAIN WINDOW

premotor_wind = [-0.07 0.02];
premotor_wind = [-0.07 0];
ff_wind = [12 40];

tbins = PARAMS.tbins;
ffbins = PARAMS.ffbins;

% ###################### MEAN OF POST MINUS PRE
tinds = tbins>premotor_wind(1) & tbins<premotor_wind(2);
finds = ffbins>ff_wind(1) & ffbins<ff_wind(2);
CohMean_WNminusBase_scalar = nan(length(OUTSTRUCT.CohMean_WNminusBase),1);
for i=1:length(OUTSTRUCT.CohMean_WNminusBase)
    tmp = OUTSTRUCT.CohMean_WNminusBase{i}(tinds, finds); % get t,f, window desired
    CohMean_WNminusBase_scalar(i) = nanmean(tmp(:));
end
OUTSTRUCT.CohMean_WNminusBase_scalar = CohMean_WNminusBase_scalar;
PARAMS.scalar_premotor_wind = premotor_wind;
PARAMS.scalar_ff_wind = ff_wind;

% ####################### ALL TRIALS
OUTSTRUCT.CohMat_scalar = cell(length(OUTSTRUCT.CohMean_WNminusBase),1);
for i=1:length(OUTSTRUCT.CohMean_WNminusBase)
    
    cohmat = OUTSTRUCT.CohMat{i};
    
    tmp = squeeze(nanmean(nanmean(cohmat(tinds, finds, :),1),2));
    
    OUTSTRUCT.CohMat_scalar{i} = tmp;
end

%% ######################### RAW PLOT, SINGLE TRIALS, LFP, COH, SPECTROGRAM
% ==== for a given channel pair plot individiual pre and post trials
% ---- what to plot?
close all;
% birdtoplot = 'pu69wh78';
% expttoplot = 'RALMANlearn2';
% swnum = 1;
birdtoplot = 'wh44wh39';
expttoplot = 'RALMANlearn2';
swnum = 1;
motiftoplot = ''; % ------ WILL plot target syl, unless say otherwise.
chanpairtoplot = [];
fs = 1500;

lt_neural_Coher_LFP_plottrials(SwitchCohStruct, OUTSTRUCT, SwitchStruct, ...
    PARAMS, birdtoplot, expttoplot, swnum, motiftoplot, chanpairtoplot, fs);


%% ========== RAW PLOTS

lt_neural_Coher_AllRaw;

%% ======= SUMMARIZE SCALAR RESULT, OVER ALL SWITCHES, MOTIFS, AND CHANNELS
% i.e. for all motifs plot change in coh scalar.
close all;
clim = [-0.2 0.2];

indsgrp_switch = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch});
indsgrp_chanpair = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch, ...
    OUTSTRUCT.chanpair});

figcount=1;
subplotrows=6;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

indsgrp_switch_unique = unique(indsgrp_switch);
indsgrp_chanpair_unique = unique(indsgrp_chanpair);
for i=1:length(indsgrp_switch_unique)
    swgrpthis = indsgrp_switch_unique(i);
    bnum = unique(OUTSTRUCT.bnum(indsgrp_switch==swgrpthis));
    enum = unique(OUTSTRUCT.enum(indsgrp_switch==swgrpthis));
    swnum = unique(OUTSTRUCT.switch(indsgrp_switch==swgrpthis));
    bname = MOTIFSTATS_Compiled.birds(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(swnum)]);
    ylabel(['WN - base, coh']);
    
    % ====== plot each channel pair its own line
    for chanpair = indsgrp_chanpair_unique'
        indsthis = indsgrp_switch==swgrpthis & indsgrp_chanpair==chanpair;
        if ~any(indsthis)
            continue
        end
        motifID = lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)); % ---- get positions within global motif
        cohscal = OUTSTRUCT.CohMean_WNminusBase_scalar(indsthis);
        assert(length(unique(motifID)) == length(motifID), 'each motif shoudl have at most 1 datapoint');
        % ---- sort in order of motifs
        [~, indsort] = sort(motifID);
        motifID = motifID(indsort);
        cohscal = cohscal(indsort);
        plot(motifID, cohscal, 'o-k');
    end
    lt_plot_zeroline;
    
    % ====== overall
    indsthis = indsgrp_switch==swgrpthis;
    
    istarg = OUTSTRUCT.istarg(indsthis);
    motifs = OUTSTRUCT.motifname(indsthis);
    motifID = lt_neural_QUICK_MotifID(bname, motifs); % ---- get positions within global motif
    cohscal = OUTSTRUCT.CohMean_WNminusBase_scalar(indsthis);
    
    [ymean, ysem] = grpstats(cohscal, motifID, {'mean', 'sem'});
    x = unique(motifID);
    lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
    
    % -------- NOTE DOWN POSITION OF TARGET SYSL
    indsthis = indsgrp_switch==swgrpthis & OUTSTRUCT.istarg==1;
    xtarg = unique(lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)));
    plot(xtarg, clim(1)+0.02, '^r');
    
    % ------- NOTE POSITION OF SAME_TYPES
    indsthis = indsgrp_switch==swgrpthis & OUTSTRUCT.istarg==0 & OUTSTRUCT.issame==1;
    xtarg = unique(lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)));
    if ~isempty(xtarg)
        plot(xtarg, clim(1)+0.02, '^b');
    end
    
    % ----- labels
    [~, motiflabels] = lt_neural_QUICK_MotifID(bname);
    set(gca, 'XTick', 1:length(motiflabels));
    set(gca, 'XTickLabel', motiflabels);
    rotateXLabels(gca, 90);
    ylim(clim);
    
end






%% SUMMARY PLOTS (compare diff syl types ...)
close all;
sumplottype = 'switches'; % i.e. what is datapoint?
% switches
% chanpairs
plotAllSwitchRaw = 0;
clim = [-0.05 0.05];
lt_neural_Coher_Learn_PlotSum(OUTSTRUCT, PARAMS, SwitchStruct, sumplottype, ...
    plotAllSwitchRaw, clim);



%%  ######################################### GRAND MEAN FIGURE:
% incorporates functions above with appropriate params.
% see notes within script.
close all;
lt_neural_Coher_Learn_PlotGrandMean;



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