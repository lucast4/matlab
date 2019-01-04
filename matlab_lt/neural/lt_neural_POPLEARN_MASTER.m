%% ==== 1) EXTRACT MOTIFSTATS POP
%% AD HOC CHANGES TO THE SETS OF NEURONS (E.G. TO MAXIMIZE DATASET SIZE)
if (0)
% ===== 1) pu69, combining sets so that can look at learning.
i=1; 
ii=1;
assert(strcmp(SummaryStruct.birds(i).birdname, 'pu69wh78'));
assert(strcmp(SummaryStruct.birds(i).exptnum_pop(ii).exptname, 'RALMANOvernightLearn1'));
assert(all(SummaryStruct.birds(i).exptnum_pop(ii).Sets_neurons{4} == [17 19]));
assert(all(SummaryStruct.birds(i).exptnum_pop(ii).Sets_neurons{5} == [17 18 19]));

% --- get new
newfiles = [SummaryStruct.birds(i).exptnum_pop(ii).Sets_songfiles{4:5}];
newneurons = [17 19];

% -- remove old
SummaryStruct.birds(i).exptnum_pop(ii).Sets_neurons(4:5) = [];
SummaryStruct.birds(i).exptnum_pop(ii).Sets_songfiles(4:5) = [];

% --- add
SummaryStruct.birds(i).exptnum_pop(ii).Sets_neurons = [SummaryStruct.birds(i).exptnum_pop(ii).Sets_neurons ...
    newneurons];
SummaryStruct.birds(i).exptnum_pop(ii).Sets_songfiles = [SummaryStruct.birds(i).exptnum_pop(ii).Sets_songfiles ...
    {newfiles}];
end


if (0)
   % ==== wh44wh39, RALMANlearn3, combining sets
    i=2; 
    ii=3;
    assert(strcmp(SummaryStruct.birds(i).birdname, 'wh44wh39'));
    assert(strcmp(SummaryStruct.birds(i).exptnum_pop(ii).exptname, 'RALMANlearn3'));
%     assert(all(SummaryStruct.birds(i).exptnum_pop(ii).Sets_neurons{2} == [13 14 15 16 19]));
%     assert(all(SummaryStruct.birds(i).exptnum_pop(ii).Sets_neurons{3} == [13 14 15 16 17 19]));
    assert(all(SummaryStruct.birds(i).exptnum_pop(ii).Sets_neurons{2} == [57    58    59    60    63]));
    assert(all(SummaryStruct.birds(i).exptnum_pop(ii).Sets_neurons{3} == [57    58    59    60    61    63]));
    
    % --- get new [i.e. the new combination of neruons]
    newfiles = [SummaryStruct.birds(i).exptnum_pop(ii).Sets_songfiles{2:3}];
%     newneurons = [13 14 15 16 19];
    newneurons = [57    58    59    60    63];
    
    % -- remove old
    SummaryStruct.birds(i).exptnum_pop(ii).Sets_neurons(2:3) = [];
    SummaryStruct.birds(i).exptnum_pop(ii).Sets_songfiles(2:3) = [];

    % --- add
    SummaryStruct.birds(i).exptnum_pop(ii).Sets_neurons = [SummaryStruct.birds(i).exptnum_pop(ii).Sets_neurons ...
        newneurons];
    SummaryStruct.birds(i).exptnum_pop(ii).Sets_songfiles = [SummaryStruct.birds(i).exptnum_pop(ii).Sets_songfiles ...
        {newfiles}];
end

%% ======================== EXTRACT SEGMENTS FOR POPULATIONS

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



%% ==== ALIGN ALL MOTIFS TO COMMON MOTIFS FOR PLOTTING ACROSS EXPERIMENTS

% ==== 1) quickly list all motifs
disp('========================');
numbirds = length(MOTIFSTATS_Compiled.birds);
for i=1:numbirds
    bname = MOTIFSTATS_Compiled.birds(i).birdname;
    MotiflistAll = {};
    numneurons = length(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons);
    for ii=1:numneurons
        
        motiflist = [MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif_regexpr_str];
        MotiflistAll = [MotiflistAll motiflist];
        
    end
    MotiflistAll = unique(MotiflistAll);
    disp(['BIRD: ' bname]);
    disp(MotiflistAll);
end


% ========================== ANNOTATE IN MOTIFSTATS COMPILED

%% ==== GET PARAMS
clear PARAMS;

PARAMS.motif_predur = MOTIFSTATS_Compiled.birds(1).MOTIFSTATS.params.motif_predur;
PARAMS.motif_postdur = MOTIFSTATS_Compiled.birds(1).MOTIFSTATS.params.motif_postdur;
PARAMS.alignbyonset = MOTIFSTATS_Compiled.birds(1).MOTIFSTATS.neurons(1).motif(1).Params.REGEXP.alignByOnset;
assert(~isempty(PARAMS.alignbyonset));
% PARAMS.savemarker = '14Oct2018_2147';
PARAMS.savemarker = input('what is save marker? (e.g. 14Oct2018_2147)? ', 's');

%% ==== REMOVE DIR SONG
MOTIFSTATS_Compiled = lt_neural_QUICK_MotCom_RemoveDIR(MOTIFSTATS_Compiled);

%% ==== get times of WN hits

timewindhit = [-0.15 0.1]; % window, rel syl onset, during which any WN onset or offset will be considered a hit.
% leave expty [] to use pre and postdur from data extraction.

numbirds = length(MOTIFSTATS_Compiled.birds);

for i=1:numbirds
   
    numneurons = length(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons);
    birdname = MOTIFSTATS_Compiled.birds(i).birdname;
    
    for ii=1:numneurons
       
        nummotifs = length(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif);
        exptid = MOTIFSTATS_Compiled.SummaryStruct.birds(i).neurons(ii).exptID;
        
        hasdir = 0;
        for mm=1:nummotifs
           
            
            segextract = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(mm).SegmentsExtract;
            
            if isempty(segextract)
                continue
            end
            
            % ========== get time of WN for all trials
            ntrials = length(segextract);
            for tt = 1:ntrials
               fname = segextract(tt).song_filename;
               disp(fname);
               ons = segextract(tt).WithinSong_TokenOns;
               
               datdur = segextract(tt).actualmotifdur;
               
               if isempty(timewindhit)
               windon = ons-PARAMS.motif_predur;
               windoff = ons + datdur + PARAMS.motif_postdur;
               else
                   windon = ons + timewindhit(1);
                   windoff = ons + timewindhit(2);
               end
               [a] = fileparts(SummaryStruct.birds(i).neurons(ii).dirname);
               tmp = load([a '/' fname '.wntime.mat']);
               tmp.wnstruct;
               
               % ---- FIND WN ONSETS/OFFSETS WITHIN WINDOW (collect if on
               % OR off is within data window)
               indsthis1 = tmp.wnstruct.WNonsets>windon & tmp.wnstruct.WNonsets<windoff;
               onsthis = tmp.wnstruct.WNonsets(indsthis1) - windon;
               
               indsthis2 = tmp.wnstruct.WNoffsets>windon & tmp.wnstruct.WNoffsets<windoff;
               offthis = tmp.wnstruct.WNoffsets(indsthis2) - windon;
               
               wasTrialHit = any([indsthis1 indsthis2]);
               
               segextract(tt).hit_WN = wasTrialHit;
               segextract(tt).WNonset_sec = onsthis;
               segextract(tt).WNoffset_sec = offthis;
            end
            MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif(mm).SegmentsExtract = segextract;
        end
    end
end




%% ================ extraction continued
close all;
MOTIFSTATS_pop = lt_neural_v2_POP_ExtractMotifs(MOTIFSTATS_Compiled, SummaryStruct);
% clear MOTIFSTATS_Compiled;


%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ XCOV (SPIKING)

lt_neural_POPLEARN_Xcov_Master;



%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ COHERENCE
% OVERALL: for each motif, look at coherence of raw data, aligned to syl
% onset. Does that change during learning?

% NOTE: to see some old progress, see:
lt_neural_MasterScript_Pop;


%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% OLD VERSION, EXTRACTING COHERENCE FIRST.

lt_neural_Coher_Script1;


if (0)
% ===== convert from segextract to a specific song and time in song
lt_neural_v2_EXTRACT_WithinSongTimings(SummaryStruct, i, neurthis);
end


%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ################## NEW VERSION, FIRST EXTRACTING LFP, THEN DOING STUFF WTIH LFP
% NOTE: [processing steps summary]
% 1) Extract LFP for each raw song file (save adjacent)
% 2) Extract those segmented data into LFPSTRUCT
% 3) Calculate (and save) coherence using data in LFPSTRUCT

%% ==== EXTRACT LFP AND SPEC STRUCT, WITH ORGANIZATION MATCHED TO COHSTRUCT
close all;
skipifOnlyOneChan = 1; % i.e. if a given dataset/motif not paired, then skip. 
BirdsToPlot = {'pu69wh78', 'wh44wh39', 'wh72pk12'};
% SetsToSkip = {'1-2-2'};
SetsToSkip = {};

LFPSTRUCT = lt_neural_LFP_ExtractStruct(MOTIFSTATS_pop, SummaryStruct, ...
    MOTIFSTATS_Compiled, skipifOnlyOneChan, BirdsToPlot, SetsToSkip);

% ================ save
if (1)
%     marker = '14Oct2018_2147';
    marker = PARAMS.savemarker;
    fname = ['/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_' marker '.mat'];
    save(fname, 'LFPSTRUCT');
    
    save(['/bluejay5/lucas/analyses/neural/LFP/PARAMS_' PARAMS.savemarker '.mat'], 'PARAMS');
end


%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

%% +++++++++++++++++++++++++++++++++ PLOT RAW NEURAL, LFP, LFP(filtred)
if (0)
    % individual trials.
    % 2 VERSIONS:
    % 1) choose which motif to plot.
    % 2) go thru all switches and plot all target sylalbles (if multiple syls, does each separately
    
    %% ############################ VERSION 1- plot chosed motif
    
    % ======== 1) EXTRACT RAW DATA (NOT JUST LFP)
    close all;
    
    birdplot = 'pu69wh78';
    exptplot = 'RALMANlearn1';
    swplot = 1;
    motifplot = []; % [string] leave blank for target
    extrapad = 0.05; % seconds, pre and post...
    [DatAll, t_onoff, fs, bregionlist, chanlist_toget, i, ii, mm] = ...
        lt_neural_LFP_PlotEgRaw_Extract(PARAMS, SwitchStruct, MOTIFSTATS_pop, SummaryStruct, ...
        SwitchCohStruct, birdplot, exptplot, swplot, motifplot, ...
        extrapad);
    
    
    % =========== GET LIST OF HIGH AND LOW COHERENCE TRIALS
    tscalar = -0.045;
    fscalar = 30;
    
    a = SwitchCohStruct.bird(i).exptnum(ii).switchlist(swplot).motifnum(mm).fileprefix;
    b = SwitchCohStruct.bird(i).exptnum(ii).switchlist(swplot).motifnum(mm).filesuffix;
    cohdat = load([a '/Coh' b]);
    cohdat = lt_neural_Coher_Cell2Mat(cohdat.CohAllTrials);
    [~, ind_f] = min(abs(PARAMS.ffbins-fscalar));
    [~, ind_t] = min(abs(PARAMS.tbins-tscalar));
    cohscaltmp = squeeze(cohdat(ind_t, ind_f, :));
    
    
    % ============ 2) FOR BASE AND WN, PLOT N TRIALS OF LFP, NEURAL
    close all;
    ntoplot = 10; % trials pre and post
    filt_low = 25;
    filt_hi = 35;
    % filt_fs = fs;
    
    savedir = ['/bluejay5/lucas/analyses/neural/LFP/FIGS_PlotEgRaw/' PARAMS.savemarker];
    saveON = 0;
    
    plotCohScalExtremes =1; % if 1, then plots example higha and low coherence tirals.
    
    lt_neural_LFP_PlotEgRaw(DatAll, t_onoff, fs, bregionlist, chanlist_toget, ...
        i, ii, swplot, mm, ntoplot, filt_low, filt_hi, SwitchCohStruct, ...
        SwitchStruct, savedir, saveON, plotCohScalExtremes, cohscaltmp);
    
    
    %% ############################ VERSION 2 - GO thru all swithches
    close all;
    for i=1:length(SwitchCohStruct.bird)
        for ii=1:length(SwitchCohStruct.bird(i).exptnum)
            
            for ss=1:length(SwitchCohStruct.bird(i).exptnum(ii).switchlist)
                
                if isempty(SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum)
                    continue
                end
                
                % ========== which motif to plot?
                % get target
                motifplot_list = ...
                    SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningDirs(1:2:end);
                
                % ================ GO THRU ALL MOTIFS
                for mm=1:length(motifplot_list)
                    motifplot = motifplot_list{mm};
                    
                    
                    % ================ RUN
                    birdplot = SwitchStruct.bird(i).birdname;
                    exptplot = SwitchStruct.bird(i).exptnum(ii).exptname;
                    extrapad = 0.05; % seconds, pre and post...
                    [DatAll, t_onoff, fs, bregionlist, chanlist_toget, ~, ~, motifnum] = ...
                        lt_neural_LFP_PlotEgRaw_Extract(PARAMS, SwitchStruct, MOTIFSTATS_pop, SummaryStruct, ...
                        SwitchCohStruct, birdplot, exptplot, ss, motifplot, ...
                        extrapad);
                    
                    % ================== PLTO AND SAVE EXAMPLES
                    close all;
                    ntoplot = 15; % trials pre and post
                    filt_low = 25;
                    filt_hi = 35;
                    % filt_fs = fs;
                    
                    savedir = ['/bluejay5/lucas/analyses/neural/LFP/FIGS_PlotEgRaw/' PARAMS.savemarker];
                    saveON = 1;
                    
                    lt_neural_LFP_PlotEgRaw(DatAll, t_onoff, fs, bregionlist, chanlist_toget, ...
                        i, ii, ss, motifnum, ntoplot, filt_low, filt_hi, SwitchCohStruct, ...
                        SwitchStruct, savedir, saveON);
                    
                end
                
            end
        end
    end
    
    
    %% ================ KEEP LIST OF CHANNELS TO EXCLUDE
    
    {...
        {'pu69wh78', 'RALMANlearn1', 1, 14}, ... % not bad. noise some trials.
        {'pu69wh78', 'RALMANOvernightLearn1', 1, 21}, ... % noise some trials.
        {'wh44wh39', 'RALMANlearn1', 1, 15}, ... % huge slow noise sometimes (rare), leads to power increase in beta
        }; % {bird, expt, swnum, chan}
    
    % ========= GOOD ONES:
    % pu69, RALMANlearn1, 2: rare trials at base see noise in all
    % channels... seems fine.
    
    % pu69, RALMANlearn2, 1: OK
    
    % pu69, RALMANlearn2, 2: OK
    
    % pu69, RALMANlearn2, 3: OK. rare trials see high hz noise -50ms, all
    % chans.

    % pu69, RALMANlearn2, 4: OK

    % pu69, RALMANlearn2, 5: OK

    % wh44wh39-RALMANlearn1-sw2 - OK
    
    % wh44wh39-RALMANlearn1-sw7-motif8 - OK
    
    % wh44wh39-RALMANlearn1-sw7-motif9 - OK

    % wh44wh39-RALMANlearn1-sw9-motif8 - OK

    % wh44wh39-RALMANlearn1-sw9-motif9 - OK
    
    % wh44wh39-RALMANlearn2-sw1-motif8 - OK
    
    % wh44wh39-RALMANlearn3-sw1-motif2 - OK
    
    % wh44wh39-RALMANlearn3-sw1-motif3 - OK
    
    % wh44wh39-RALMANlearn4-sw1-motif7 - OK
    
    % wh44wh39-RALMANlearn4-sw2-motif6 - OK
    
    % SUMMARY: Only noise that seems to affect beta band (25-35hz) is in
    % wh44 rl1, sw1, which is the large deviations. Also sometimes in wh44
    % ral3[? I can't remember] saw some large deviations that were
    % suspiociously synchronized between LMAN and RA, but not obvious
    % qwithin eacjh region that is noise...
    
end


%% =============== [GOOD, CHECK EACH EXPT IN DETAIL]
% 1) timecourse of scalar coherence
% 2) raw neural data for select trials (e.g. high, low coh)
% 3) comparison of multitaper method to across trials.
lt_neural_LFP_SingleExptCheck;


%% ################## [GOOD, RAW PLOTS]
lt_neural_Coher_AllRaw;


%% #######################################################################
%% ====== CALCULATE COHERENCE USING LFP
% GOAL TO GET COHSTRUCT
% NOTE: will save cohstruct also.
close all;
% savemarker = '14Oct2018_2147';
% savemarker = PARAMS.savemarker;
% savemarker = ['14Oct2018_2147_150msWind'];
if (0)
PARAMS.savemarker = '14Oct2018_2147_150msWind';
PARAMS.savemarker = '14Oct2018_2147_75msWind';
end
movingwin = [0.1 0.01];

COHSTRUCT = lt_neural_LFP_GetCohStruct(LFPSTRUCT, PARAMS, SummaryStruct, movingwin);

%% ############ SAVE OVERALL PARAMS (inclyding time bins)
if exist('All_tbins', 'var')
    tbins = All_tbins{1};
    ffbins = All_ffbins{1};
    PARAMS.tbins = tbins;
    PARAMS.ffbins = ffbins;
    PARAMS.ffbinsedges = [10 25 32 80]; % edges, to plot timecourse in frequency bands
else
    assert(~isempty(COHSTRUCT.bird(1).experiment(1).setnum(1).motif));
    PARAMS.tbins = COHSTRUCT.bird(1).experiment(1).setnum(1).motif(1).t_relons;
    PARAMS.ffbins = COHSTRUCT.bird(1).experiment(1).setnum(1).motif(1).ffbins;
    PARAMS.ffbinsedges = [10 25 32 80]; % edges, to plot timecourse in frequency bands
end

%% ################# ANALYSES START FROM HERE: LOAD DATA
% ==== SAVE
save('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_14Oct2018_2147.mat', ...
    'MOTIFSTATS_Compiled', '-v7.3');

% ==== LOAD 
% 1) all data (100ms window)
clear all; close all;
load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_14Oct2018_2147.mat');
load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_14Oct2018_2147.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_14Oct2018_2147.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/14Oct2018_2147/COHSTRUCT.mat');

% 1) all data  [LATEST, 3 birds, up to wh72, RALMANLearn3]
clear all; close all;
load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_15Dec2018_1321.mat');
load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_15Dec2018_1321.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_15Dec2018_1321.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/15Dec2018_1321/COHSTRUCT.mat');

% 1) all data  [LATEST, 3 birds, up to wh72, RALMANLearn5]
clear all; close all;
load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_03Jan2019_2009.mat');
load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_03Jan2019_2009.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_03Jan2019_2009.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/03Jan2019_2009/COHSTRUCT.mat');


% 2) all data (150ms window)
clear all; close all;
load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_14Oct2018_2147.mat');
load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_14Oct2018_2147.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_14Oct2018_2147.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/14Oct2018_2147_150msWind/COHSTRUCT.mat');

% 3) just wh44
clear all; close all;
load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_05Dec2018_1758.mat');
load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_05Dec2018_1758.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_05Dec2018_1758.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/05Dec2018_1758/COHSTRUCT.mat');

% 4) just wh72 (done 12/12)
clear all; close all;
load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_12Dec2018_0036.mat');
load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_12Dec2018_0036.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_12Dec2018_0036.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/12Dec2018_0036/COHSTRUCT.mat');


% 4) just wh72 (done 12/20) - includes up to RALMANLearn5 - also includes
% all syls
clear all; close all;
load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_20Dec2018_0951.mat');
load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_20Dec2018_0951.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_20Dec2018_0951.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/20Dec2018_0951/COHSTRUCT.mat');


% 4) just wh72 RA LMAN LEran4
clear all; close all;
load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_18Dec2018_1837.mat');
load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_18Dec2018_1837.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_18Dec2018_1837.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/18Dec2018_1837/COHSTRUCT.mat');

% 4) just wh72 RALMANLearn4 and 5 - each with one fake RA and one real LMAN
clear all; close all;
load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_20Dec2018_1723.mat');
load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_20Dec2018_1723.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_20Dec2018_1723.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/20Dec2018_1723/COHSTRUCT.mat');

% 4) just wh72 RALMANLearn4 and 5 - each with one fake LMAN and one fake RA
clear all; close all;
load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_21Dec2018_0028.mat');
load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_21Dec2018_0028.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_21Dec2018_0028.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/21Dec2018_0028/COHSTRUCT.mat');


% 4) all fake data (3 birds, curated)
clear all; close all;
load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_03Jan2019_0153.mat');
load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_03Jan2019_0153.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_03Jan2019_0153.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/03Jan2019_0153/COHSTRUCT.mat');

% ====== update params to match COHSTRUCT
settmp = 1;
if isempty(COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif)
    settmp=2;
    assert(~isempty(COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif));
end
PARAMS.tbins = COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif(1).t_relons;
PARAMS.ffbins = COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif(1).ffbins;
PARAMS.ffbinsedges = [10 25 32 80]; % edges, to plot timecourse in frequency bands

%% #################### [SHORTCUT - PREPROCESING]
% combines all things would do in this script into its own script. 
% Just run this after loading new data.

lt_neural_Coher_PreProcess;

%% ################ NOTE DOWN BRAIN REGION PAIRS

COHSTRUCT = lt_neural_Coher_GetBrRegPairs(COHSTRUCT, MOTIFSTATS_pop, SummaryStruct);

    %% ======== FOR LEARNING, GET SWITCHES
SwitchStruct = lt_neural_LEARN_getswitch(SummaryStruct);

    
%% ======= PLOT RAW LFP, PHI, SPECTROGRAM, COHEROGRAM

close all;
birdtoplot = 'wh72pk12';
expttoplot = 'RALMANLearn3';
swnum = 7;
% birdtoplot = 'wh44wh39';
% expttoplot = 'RALMANlearn2';
% swnum = 1;
motiftoplot = 'jr(b)'; % ------ WILL plot target syl, unless say otherwise.
chanpairtoplot = [];
bregiontoget = 'LMAN-RA'; % WILL take random channel (for bregion pair) if no chanpairtoplot specified.
fs = 1500;
Nplot = 5; % n trials base and dur WN.
plotonlyLFP = 1; % then doesn't plot cohgram, spectra etc.
lt_neural_Coher_LFP_plottrials(SwitchCohStruct, [], SwitchStruct, ...
    PARAMS, birdtoplot, expttoplot, swnum, motiftoplot, chanpairtoplot, fs, ...
    Nplot, bregiontoget, plotonlyLFP);



%% ====== GET DATA FOR LEARNING (RELATIVE TO LEARNING SWITCHES)

if (0)
pairtoget = 'LMAN-RA';
SwitchCohStruct = lt_neural_Coher_LearnExtr(COHSTRUCT, MOTIFSTATS_pop, ...
    SwitchStruct, pairtoget, LFPSTRUCT, PARAMS);
else
    % --- doesn't extract cohmat, just saves path to cohmat.
    baseuseallinds =0 ;
pairtoget = 'LMAN-RA';
% pairtoget = 'LMANoutside-RA';
% pairtoget = {'LMANoutside-RA', 'LMAN-RAoutside', 'LMANoutside-RAoutside'};
SwitchCohStruct = lt_neural_Coher_LearnExtr2(COHSTRUCT, MOTIFSTATS_pop, ...
    SwitchStruct, pairtoget, LFPSTRUCT, PARAMS, baseuseallinds);    
end

PARAMS.bregionpair_toget = pairtoget;


%% ######################## [WN HITS] CHECK 1) fraction and 2)timing 

close all;
lt_neural_LFP_WNCheck(SwitchCohStruct, MOTIFSTATS_pop, SwitchStruct, ...
    PARAMS);

%% ########################## COH CORRELATE WITH PITCH?
%% ======== EXTRACT SCALARS

twind = [-0.08 -0.03];
fwind = [25 35];
twind = [-0.08 -0.03];
fwind = [22 32];
% twind = [-0.05 -0.0];
% fwind = [25 40];
% twind = [-0.09 -0.02];
% fwind = [15 40];
twind = [-0.07 -0.03]; % all combined
fwind = [22 32];

[SwitchCohStruct, PARAMS] = lt_neural_LFP_PitchCorr(COHSTRUCT, SwitchCohStruct,...
    PARAMS, twind, fwind);

%% ======= PLOT ALL DATA [COH SCALAR vs FF]
close all;
switchestoplot = [1];
birdstoplot = {};

numbirds = length(SwitchCohStruct.bird);
for i=1:numbirds
    bname = SummaryStruct.birds(i).birdname;
    numexpt = length(SwitchCohStruct.bird(i).exptnum);
    for ii=1:numexpt
        exptname = SummaryStruct.birds(i).exptnum_pop(ii).exptname;
        numswitch = length(SwitchCohStruct.bird(i).exptnum(ii).switchlist);
        for ss=1:numswitch
           
            if ~isempty(switchestoplot)
               if ~any(ismember(switchestoplot, ss))
                   continue
               end
            end
            
            figcount=1;
            subplotrows=6;
            subplotcols=3;
            fignums_alreadyused=[];
            hfigs=[];
            hsplots = [];

            nummotif = length(SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum);
           for mm=1:nummotif
               disp([num2str(i) '-' num2str(ii) '-' num2str(ss)]);
               datthis = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm);
               if isempty(datthis.chanpair)
                   continue
               end
               
               motifname = datthis.motifname;
               tmp = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningDirs(1:2:end);
               if any(strcmp(tmp, motifname))
                   istarg=1;
                   learndir = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningDirs{find(strcmp(tmp, motifname))+1};
               else
                   istarg=0;                   
               learndir = [];
               end
               
               % ==================== extract coherene scalars for this dat
               % ------------ COLLECT SCALAR
               t = datthis.tvals;
               ff = datthis.ffvals;
               cohscal = datthis.cohscalar;
               assert(size(t,2) == size(cohscal,1));
               
               % ------------ WN switch times
               swtime = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum;
               
               % ================ ONE PLOT FOR EACH CHANNEL PAIR
                nchanpairs = size(cohscal,2);
                hsplots = [];
                for cc = 1:nchanpairs
                    chpairthis = datthis.chanpair(cc,:);
                   
                    X = {};
                    Y = {};

                   % -- baseline
                   [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                   hsplots = [hsplots hsplot];
                   title({[bname '-' exptname '-s' num2str(ss)], ...
                       [motifname '-ch' num2str(chpairthis) '[BASE]']});
                   indtrial = t<swtime;
                   pcol = 'k';
                   
                   xlabel('coh scalar');
                   ylabel('ff');
                   x = cohscal(indtrial,cc);
                   y = ff(indtrial);
                   plot(x, y, 'x', 'Color', pcol);
                   lt_regress(y, x, 0, 0, 1, 1, 'm', 1);
                   
                   X{1} = x;
                   Y{1} = y;
                   
                   % -- WN
                   [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                   hsplots = [hsplots hsplot];
                   title({[bname '-' exptname '-s' num2str(ss)], ...
                       [motifname '-ch' num2str(chpairthis) '[WN]']});
                   indtrial = t>swtime;
                   pcol = 'm';
                   
                   xlabel('coh scalar');
                   ylabel('ff');
                   x = cohscal(indtrial,cc);
                   y = ff(indtrial);
                   plot(x, y, 'x', 'Color', pcol);
                   lt_regress(y, x, 0, 0, 1, 1, 'm', 1);
                   
                   X{2} = x;
                   Y{2} = y;

                   % ---- if target, note down learn dir
                   if istarg==1
                       lt_plot_text(min(x), max(y), ['targ, dir' num2str(learndir)], 'b');
                   end
                   
                   % -------------- OVERLAY DISTRIBUTIONS FOR BASE AND WN
                   [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                   title({[bname '-' exptname '-s' num2str(ss)], ...
                       [motifname '-ch' num2str(chpairthis)]});
                   xlabel('coh (mean, sem)');
                   ylabel('ff (mean, sem]');
                   xmean = cellfun(@mean, X);
                   ymean = cellfun(@mean, Y);
                   xerr = cellfun(@lt_sem, X);
                   yerr = cellfun(@lt_sem, Y);
                   lt_plot(xmean(1), ymean(1), {'Errors', yerr(1), 'Xerrors', xerr(1), 'Color', 'k'});
                   lt_plot(xmean(2), ymean(2), {'Errors', yerr(2), 'Xerrors', xerr(2), 'Color', 'r'});
                   xlim([0.2 0.8]);
                end
                linkaxes(hsplots, 'xy');
                              
           end
           pause;
           close all;
        end
    end
end



%% ############################################

%% ################################################# LEARNING CHANGES
%% ====== EXTRACT LFP ACROSS LEARNING
close all;
collectAllProcess =1; % then colelcts not just Cohmat, but all phi and spectra.
% if 1, then will not collect each trial but just means (as takes too much
% meory). (mean pre and post and also diff)
plotON = 0;
averagechanpairs= 0; % for each motif, average over all chan pairs [NOTE: this is not up to date]
onlyfirstswitch = 0;
removeBadSyls = 1; % LEAVE AT 1.
zscoreLFP = 3; % default 1, z-scores each t,ff bin separately.
% if 2, then doesn't zscore, instead normalizes as power proprotion
% (separately for each time slice and trial).
% if 3, then each ff window z-scored separately (good to see moudlation).
% if 4, then first 1) zscores within f, and 2) normalizes to all f (i.e.
% proportion
collectDiffMats = 0; % if 1, then collects differences (WN minus base). redundant, so leave at 0.

if (0)
[OUTSTRUCT, OUTSTRUCT_CohMatOnly] = lt_neural_LFP_Learn_Extr(SwitchStruct, SwitchCohStruct, ...
    plotON, averagechanpairs, PARAMS, onlyfirstswitch, removeBadSyls, ...
    collectAllProcess, zscoreLFP, collectDiffMats);
end


[OUTSTRUCT] = lt_neural_LFP_Learn_Extr(SwitchStruct, SwitchCohStruct, ...
    plotON, averagechanpairs, PARAMS, onlyfirstswitch, removeBadSyls, ...
    collectAllProcess, zscoreLFP, collectDiffMats);

%% for each case get wn minus base, mean coh scalar
cohscal_diff = [];
for i=1:length(OUTSTRUCT.bnum)
    
    indsbase = OUTSTRUCT.indsbase_epoch{i};
    indswn = OUTSTRUCT.indsWN_epoch{i};
    
    cohscal = OUTSTRUCT.cohscal{i};
    
    cohdiff = mean(cohscal(indswn)) - mean(cohscal(indsbase));
    
    cohscal_diff = [cohscal_diff; cohdiff];
end

OUTSTRUCT.cohscal_diff = cohscal_diff;
%% ====== don't want raw data?
clear SwitchCohStruct;
clear OUTSTRUCT_CohMatOnly



%% ############ [FILTER BANK - EXTRACT]
% Extracts for each motif filter bank data in the range dsried. 
% First step before any anaylsis using filter banks
% analogous to LFPSTRUCT extraction above.

% ====================== PARAMS
BirdsToPlot = {'pu69wh78', 'wh72pk12', 'wh44wh39'};
extrapad = 0.05; % collect extra (sec) on edges.
saveON = 1;
skipifdone = 1; % 0 to overwrite all data.

lt_neural_LFP_FiltBankExtr(SummaryStruct, SwitchCohStruct, MOTIFSTATS_pop, ...
    PARAMS, BirdsToPlot, extrapad, saveON, skipifdone);


%% ============== [FILTER BANK]
% ===== recalculate coherence using filtered data.




%% ============== [FILTER BANK] Sanity check - compare a random trial to raw data and extracted lfp.




%% ====== [RECALCULATE COHERENCE USING LFP DATA] [OR XCORR OF LFP]
% -- for each experiement, recalculate using currently extracted windows.
close all;
usecorrcoeff = 1; % for lfp xcorr/
thingstodo = {'cohere'};
% thingstodo = {'lfpxcorr'};
% thingstodo = {'waveletcoh'};
[OUTSTRUCT, PARAMS, LFPXCORR_Base, LFPXCORR_WN, LFPXCORR_freqsall] =  ...
    lt_neural_LFP_RecalcCoh(OUTSTRUCT, SwitchCohStruct, LFPSTRUCT, SwitchStruct, ...
    PARAMS, usecorrcoeff, thingstodo);





%% ====== COLLECT LEARN DIR AT TARGET

OUTSTRUCT = lt_neural_LFP_GetLearnDir(OUTSTRUCT, SwitchStruct);
    

%% ====== extract specific switch types
% [HERE JUST DISPLAYS...]
skipifnodata =1; % if 1, then olnly shows a switch if ther eis data extracted...
for j=1:length(SwitchStruct.bird)
    
    for jj=1:length(SwitchStruct.bird(j).exptnum)
        disp( ' ===========' );
        
        for ss=1:length(SwitchStruct.bird(j).exptnum(jj).switchlist)
            
            if skipifnodata==1
            try SwitchCohStruct.bird(j).exptnum(jj).switchlist(ss).motifnum(1);
            catch
                continue
            end
            end
            learningContingencies = SwitchStruct.bird(j).exptnum(jj).switchlist(ss).learningContingencies;
            
            strtoplot = '';
            for k=1:length(learningContingencies)/2
               strtoplot = [strtoplot ' -- ' [learningContingencies{2*k-1} ' [' num2str(learningContingencies{2*k}) ']']];
            end
            bname = SwitchStruct.bird(j).birdname;
            ename = SwitchStruct.bird(j).exptnum(jj).exptname;
            disp(' ');
            disp([bname '-' ename '-sw' num2str(ss)]);
            
            % ============ DISPLAY WHAT THE SAME TYPE SYLS ARE
            indstmp = OUTSTRUCT.bnum==j & OUTSTRUCT.enum==jj ...
                & OUTSTRUCT.switch==ss & OUTSTRUCT.issame==1 & ...
                OUTSTRUCT.istarg==0;
            
            samesyls = unique(OUTSTRUCT.motifname(indstmp));
            strtoplot = [strtoplot ' |||||||||| SAMETYPES: '];
            for k=1:length(samesyls)
                strtoplot = [strtoplot ' ' samesyls{k} ','];
            end
            disp(strtoplot);
            
            
        end
    end
end

%% ============== FIND DATA THAT MATCH CRITERIA OF SWITCH TYPE
% == will consider it a match if EVERY target motif for a given switch
% pasess criterion.

swtoget = {[0 1], [0 -1]}; % passes if matches ANY of these
firstswitchfortarget_withinday = 0; % if 1, then onlky keeps if all targets 
% for a given switch did not have a previous switch on the same day
indtoget_b_e_s = lt_neural_LEARN_FilterSwitches(SwitchStruct, swtoget, ...
    firstswitchfortarget_withinday);

% indsgrp = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch});
% for i=1:length(OUTSTRUCT)
%     bnum = OUTSTRUCT.bnum(i);
%     enum = OUTSTRUCT.enum(i);
%     swnum = OUTSTRUCT.switch(i);
%     
%     targmotifs = = OUTSTRUCT.motifname{i};
%     
%     lc = SwitchStruct.bird(bnum).exptnum(enum).switchlist(swnum).learningContingencies;
%     strcmp(lc, motifname)
% end

%% ====== SUMMARY PLOT OF COHERENCE LARNING

% SUMMARY PLOTS (compare diff syl types ...)
close all;
sumplottype = 'switches'; % i.e. what is datapoint?
% switches
% chanpairs
plotAllSwitchRaw = 0;
clim = [-0.1 0.1];
fieldtoplot = 'Spec1Mean_WNminusBase';
% Spec1Mean_WNminusBase
% Spec2Mean_WNminusBase
% CohMean_WNminusBase
% === v1 - % NOTE: this is obsolete. see below. TO RUN THIS with coher,
% make sure in extraction of OUTSTRUCT included collectDiffMats.
lt_neural_Coher_Learn_PlotSum(OUTSTRUCT, PARAMS, SwitchStruct, sumplottype, ...
    plotAllSwitchRaw, clim, fieldtoplot);


% ============== THIS IS BETTER - subsumes the above, more compacta ndf
% flexible.
fieldtoplot = 'coher';
% 'coher'
% 'spec'
birdstoplot = [3];
expttoplot = [];
% swtoplot = [1 7 9 11];
swtoplot = [];

% --- to get specific switch types. ... [is done in addition to above
% fitlers]
% swtoget = {}; % passes if matches ANY of these
swtoget = {[0 1], [0 -1]}; % passes if matches ANY of these
firstswitchfortarget_withinday = 1; % if 1, then onlky keeps if all targets 
% for a given switch did not have a previous switch on the same day
% swtoget = {[1 0], [-1 0]}; % passes if matches ANY of these
% firstswitchfortarget_withinday = 1; % if 1, then onlky keeps if all targets 
% for a given switch did not have a previous switch on the same day
firstswitchofday=1;
indtoget_b_e_s = lt_neural_LEARN_FilterSwitches(SwitchStruct, swtoget, ...
    firstswitchfortarget_withinday, firstswitchofday);

% NOTE: to get raw cohgram (before subtract) cd /blucurrently need to do
% breakpoint using spec and evaluate cohgram version instead. Should
% modify to plot cohgram.
timewindowtoplot = [-0.08 0]; % for spectra.
ffbinsedges = [20 35 80 130]; % edges, to plot timecourse in frequency bands
lt_neural_Coher_Learn_PlotSum2(OUTSTRUCT, PARAMS, SwitchStruct, sumplottype, ...
    plotAllSwitchRaw, clim, fieldtoplot, birdstoplot, timewindowtoplot, ...
    zscoreLFP, expttoplot, swtoplot, ffbinsedges, indtoget_b_e_s);


% ================ ONE PLOT FOR EACH SWITCH
close all;
fieldtoplot = 'coher';
% 'coher'
% 'spec'
% NOTE: to get raw cohgram (before subtract) currently need to do
% breakpoint using spec and evaluate cohgram version instead. Should
% modify to plot cohgram.
timewindowtoplot = [-0.08 0]; % for spectra.
birdstoplotTMP = [];
expttoplotTMP = [];
swtoplotTMP = [];
for i=1:max(OUTSTRUCT.bnum)
    if ~isempty(birdstoplotTMP)
        if ~any(birdstoplotTMP==i)
            continue
        end
    end
    for ii=1:max(OUTSTRUCT.enum)
        if ~isempty(expttoplotTMP)
            if ~any(expttoplotTMP==ii)
                continue
            end
        end
        for ss=1:max(OUTSTRUCT.switch)
            if ~isempty(swtoplotTMP)
                if ~any(swtoplotTMP==ss)
                    continue
                end
            end
            
            if ~any(OUTSTRUCT.bnum==i & OUTSTRUCT.enum==ii & OUTSTRUCT.switch==ss)
                continue
            end
            
            % === plot for this case
            birdstoplot = [i];
            expttoplot = [ii];
            swtoplot = [ss];
            lt_neural_Coher_Learn_PlotSum2(OUTSTRUCT, PARAMS, SwitchStruct, sumplottype, ...
                plotAllSwitchRaw, clim, fieldtoplot, birdstoplot, timewindowtoplot, ...
                zscoreLFP, expttoplot, swtoplot, ffbinsedges);
            
            birdname = SwitchStruct.bird(i).birdname;
            exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
            lt_subtitle([birdname '-' exptname '-s' num2str(ss)]);
            
        end
    end
end




%% ####################### [MOTIFPLOT] SUMMARIZE SCALAR CHANGE IN LEARNING
% I.E. EACH syl in order...
% ======= SUMMARIZE SCALAR RESULT, OVER ALL SWITCHES, MOTIFS, AND CHANNELS

close all;

swtoget = {[0 1], [0 -1]}; % passes if matches ANY of these
% swtoget = {[1 0], [-1 0]}; % passes if matches ANY of these
firstswitchfortarget_withinday = 1; % if 1, then onlky keeps if all targets 
% for a given switch did not have a previous switch on the same day
firstswitchofday =1;
indtoget_b_e_s = lt_neural_LEARN_FilterSwitches(SwitchStruct, swtoget, ...
    firstswitchfortarget_withinday, firstswitchofday);

lt_neural_Coher_SumPlotMotifs(OUTSTRUCT, SwitchStruct, MOTIFSTATS_Compiled, PARAMS, ...
    indtoget_b_e_s);


%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% +++++++++++++++++++++ [COMBINIGN MULTIPLE MOTIFPLOT STRUCTURES];

%% ==================== 1) [SAVE SUMMARY, FOR COMPARISON ACROSS EXTRACTED DATASETS]
% === Saves once for each bird/expt/switch - 
% === Therefore b/e/s must be matched across strtuctures that are comparing
% === saves: 1) mean coherogram (pre and post) and 2) scalar (pre
% and post), for each syllable (once for each channel pair). 3) is target?
% 4) same type? 5) motif name

% NOTE: overwrites any old data extracted using this dataset
% NOTE: currently only coded to take mean across channels
save_xchan_means = 1; % first collapses across channel pairs.

lt_neural_Coher_SaveMotifDat(OUTSTRUCT, SwitchStruct, MOTIFSTATS_Compiled, PARAMS, ...
    save_xchan_means);

%% #################### [PLOTS AND ANALYSES]
close all;
lt_neural_Coher_XDAT_Master;

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ######################### [MOTIFPLOT] SAME, DIFF MOTIFS

close all;

swtoget = {[0 1], [0 -1]}; % passes if matches ANY of these
% swtoget = {[1 0], [-1 0]}; % passes if matches ANY of these
firstswitchfortarget_withinday = 1; % if 1, then onlky keeps if all targets 
% for a given switch did not have a previous switch on the same day
firstswitch=1;
indtoget_b_e_s = lt_neural_LEARN_FilterSwitches(SwitchStruct, swtoget, ...
    firstswitchfortarget_withinday, firstswitch);

lt_neural_Coher_PlotSameMotif(OUTSTRUCT, SwitchStruct, MOTIFSTATS_Compiled, ...
    PARAMS, indtoget_b_e_s);

%% ############### [MOTIFPLOT] FOR EACH MOTIF, PLOT COH CHANGE DEPEDNING ON STATUS
close all;

% ########## FOR A GIVEN SYL, PLOT ITS CHANGE IN COHERNECE DEPENDING ON WHETHER IT IS TARG, SAME OR DIFF ACROSS EXPERIMENTS
averageOverChanPairs = 1; % default, since chan pairs are so similar.
statmethod = 'diff';
meanOverExpt = 1; % then for each motif gets one value for each status type (targ, same ..);

birdstoplot = [];
expttoplot = [];
swtoplot = [1 7 9 11];
swtoplot = [];


% === to get specific switch types. ... [is done in addition to above
% fitlers]
% swtoget = {}; % passes if matches ANY of these
swtoget = {[0 1], [0 -1]}; % passes if matches ANY of these
% swtoget = {[1 0], [-1 0]}; % passes if matches ANY of these
firstswitchfortarget_withinday = 1; % if 1, then onlky keeps if all targets 
% for a given switch did not have a previous switch on the same day
firstswitch=1;
indtoget_b_e_s = lt_neural_LEARN_FilterSwitches(SwitchStruct, swtoget, ...
    firstswitchfortarget_withinday, firstswitch);

lt_neural_Coher_SumPlotMotifs2(OUTSTRUCT, SwitchStruct, averageOverChanPairs, ...
    statmethod, indtoget_b_e_s, meanOverExpt, birdstoplot, expttoplot, swtoplot);



%% ########### [MOTIFPLOT] Change in coherence depending on position in motif


%% ################################# LFP CROSS-CORRELATIONS
% plot every single motif x channel pair. will pause after each experiemnt
% and then close on instruction
close all;
assert(exist('LFPXCORR_Base', 'var')==1, 'neeed to first run code to extract this...');
freqtoplot = 30;
lt_neural_LFP_LFPxcorr_plotall;

%% ############ [LFP XCORR] - CONSISTENT LAG?
% ==== quick, plots average lag across all datapoints.

tmp = cellfun(@(x)mean(x(:,3,2)), OUTSTRUCT.LFPXCORR_Base);
lt_figure; hold on; lt_plot_histogram(tmp);
xlabel('LMAN leads < --- > RA leads [sec]');
lt_plot_zeroline_vert;

%% ############## [LFP XCORR] - SUMMARY PLOTS

birdstoplot = [];
expttoplot = [];
swtoplot = [1];
freqtoplot = 30;

% ==== 1) get WN minus baseline differences
LFPXCORR_FracChange = nan(size(LFPXCORR_WN,1), size(LFPXCORR_WN{1},2)); % cases x freqs
for j=1:length(LFPXCORR_WN)
    
    % ---- get meadian within WN or base, 1 for each freq
    b = median(LFPXCORR_WN{j}(:,:,1), 1);
    a = median(LFPXCORR_Base{j}(:,:,1), 1);
%     b = mean(LFPXCORR_WN{j}(:,:,1), 1);
%     a = mean(LFPXCORR_Base{j}(:,:,1), 1);
    
    xcorr_fracchange = (b-a)./a;
    LFPXCORR_FracChange(j,:) = xcorr_fracchange;
end


% ==== put things into OUTSTRUCT.
OUTSTRUCT.LFPXCORR_Base = LFPXCORR_Base;
OUTSTRUCT.LFPXCORR_WN = LFPXCORR_WN;
OUTSTRUCT.LFPXCORR_FracChange = LFPXCORR_FracChange;
OUTSTRUCT.LFPXCORR_freqsall = LFPXCORR_freqsall;

lt_neural_LFP_LFPxcorr_Summary(OUTSTRUCT, PARAMS, SwitchStruct, ...
    birdstoplot, expttoplot, swtoplot, freqtoplot);

%% ###################################################################
%% ################ COHERENCE SCALAR PLOTS
% - Can predict change in coh based on ffvs coh correlation?
% - Is there signifincat correlation between coh and FF?
% - Shuffle stats for change in coherence across experiments.
% NOTE: USES SINGLE SCALAR VALUE EXTRACTED USING lt_neural_LFP_PitchCorr

%% ====== EXTRACT BASELINE CORRELATIONS.
useonlybaseepoch = 0; % if 1, then just limited epoch. if 0, then entie baseline data
% default: 0;
corrtype = 'spearman'; % spearman. or pearson or Kendall
cohdiff_usedprime = 0; % if 0, then just mean diff. if 1, then dprime
nboot = 10; % to get bootstrap SE

OUTSTRUCT = lt_neural_Coher_CohScalExtract(OUTSTRUCT, useonlybaseepoch, ...
    corrtype, cohdiff_usedprime, nboot);


%% ===== [SUMMARY PLOT] LEARNING CHANGE IN COHSCALAR AS FUNCTION OF BASE CORR? 
% i.e. can I predict the change in coherence in a given expt (switch) based
% on relationship between baseline coh/ff and direction of training
% NOTE: can test variation within and between experiments (centerdata)
% NOTE: can test both normative learn dir and actual ff change (useFF ...);
close all;
cohdiff_norm_to_global = 1; % if 1, then finds (for a given channel)
rho_norm_to_global = 0;
centerdata = 0; % then for each switch, centers data (i.e. across channels)
% NOTE: this is useful if want to ask about, within a tgiven expt, is there
% correlation. 
onlyfirstswitch = 0; % if 0, then all siwtches.
useFFchangeInsteadOfLearnDir=0; % then sign(ff(WN)-ff(base)); if 0, then learndir

lt_neural_Coher_CohScalPlot(OUTSTRUCT, SwitchStruct, cohdiff_norm_to_global, ...
    rho_norm_to_global, centerdata, onlyfirstswitch, useFFchangeInsteadOfLearnDir);


%% ###################################################################
%% ################ COHERENCE SCALAR PLOTS [SYSTEMATIC]
corrtype = 'spearman'; % spearman. or pearson or Kendall
Nshuff = 200;

% GOES THRU ALL T, FF BINS
% REQUIRES OUTSTRUCT_CohMatOnly
tbins = PARAMS.tbins;
fbins = PARAMS.ffbins;

FFcorrCoh = nan(length(tbins), length(fbins), length(OUTSTRUCT.bnum));
FFcorrCoh_pctileVsShuff = nan(length(tbins), length(fbins), length(OUTSTRUCT.bnum));
FFcorrCoh_zscoreVsShuff = nan(length(tbins), length(fbins), length(OUTSTRUCT.bnum));
FFcorrCoh_shuffCI = nan(length(tbins), length(fbins), 2, length(OUTSTRUCT.bnum));

savedir = '/bluejay5/lucas/analyses/neural/COHERENCE/SCALAR';
savedir = [savedir '/' PARAMS.savemarker];
if ~exist(savedir)
mkdir(savedir);
end


for i=1:length(OUTSTRUCT.bnum)
    
    inds_base = OUTSTRUCT.indsbase{i};
    
    cohmatall = OUTSTRUCT_CohMatOnly{i};
    ff = OUTSTRUCT.ffvals{i};
    
    % ================ GO THRU EACH T, FF BIN
    for j=1:length(tbins)
        for jj=1:length(fbins)
            
            cohthis = squeeze(cohmatall(j, jj, inds_base));
            ffthis = ff(inds_base);
            
            % ---------------- GET CORR
            rho_dat = corr(ffthis', cohthis, 'type', corrtype);
            
            % ---------------- GET SHUFFLE DISTRIBUTION OF CORR
            rho_shuff_all = [];
            for nn=1:Nshuff
                disp(['case ' num2str(i) ', bin1=' num2str(j) ', bin2=' num2str(jj) ', shuff' num2str(nn)]);
                
                % --- shuffle trials
                ffthis_shuff = ffthis(randperm(length(ffthis)));
                rho_shuff = corr(ffthis_shuff', cohthis, 'type', corrtype);
                rho_shuff_all = [rho_shuff_all; rho_shuff];
                
            end
            
            % ================= OUTPUT STATS
            FFcorrCoh(j, jj, i) = rho_dat;
            
            p = (sum(abs(rho_shuff_all)>=abs(rho_dat))+1)./(length(rho_shuff_all)+1);
            FFcorrCoh_pctileVsShuff(j, jj, i) = p;
            
            shuffmean = mean(rho_shuff_all);
            shuffstd = std(rho_shuff_all);
            rho_z = (rho_dat - shuffmean)/shuffstd;
            FFcorrCoh_zscoreVsShuff(j, jj, i) = rho_z;
            FFcorrCoh_shuffCI(j, jj, :, i) = prctile(rho_shuff_all, [2.75 97.5]);
            
        end
    end
save([savedir '/FFcorrCoh'], 'FFcorrCoh');
save([savedir '/FFcorrCoh_pctileVsShuff'], 'FFcorrCoh_pctileVsShuff');
save([savedir '/FFcorrCoh_zscoreVsShuff'], 'FFcorrCoh_zscoreVsShuff');
save([savedir '/FFcorrCoh_shuffCI'], 'FFcorrCoh_shuffCI');

end


% clear OUTSTRUCT_CohMatOnly;

%% ===================== 

tmp = FFcorrCoh_pctileVsShuff<0.05;

%% ====================== PLOT SHUFFLE SUMMARY [COH SCALAR]
indsgood = ~isnan(squeeze(FFcorrCoh(1,1,:))) & OUTSTRUCT.switch==1;

figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


% ============ 1) heat map of mean correlation
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('mean cor coeff, (not abs)');
cormean = nanmean(FFcorrCoh(:,:,indsgood), 3);
imagesc(tbins, fbins, cormean');
colorbar('East');
axis tight;


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('mean cor coeff, (z,score not abs)');
cormean = nanmean(FFcorrCoh_zscoreVsShuff(:,:,indsgood), 3);
imagesc(tbins, fbins, cormean');
colorbar('East');
axis tight;


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('mean cor coeff, (mean of abs val)');
cormean = nanmean(abs(FFcorrCoh(:,:,indsgood)), 3);
imagesc(tbins, fbins, cormean');
colorbar('East');
axis tight;


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('mean cor coeff, (mean of abs val of zscore)');
cormean = nanmean(abs(FFcorrCoh_zscoreVsShuff(:,:,indsgood)), 3);
imagesc(tbins, fbins, cormean');
colorbar('East');
axis tight;


%% ###################################### PHI
%% ==== MEAN PHI ACROSS ALL SYLS

PhiMeanAll = nan(length(PARAMS.tbins), length(PARAMS.ffbins), length(OUTSTRUCT.bnum));

for i=1:length(OUTSTRUCT.bnum)
   disp(i);
    % ======= collect mean baseline phi
    inds_base = OUTSTRUCT.indsbase_epoch{i};
    phiall = OUTSTRUCT.PhiMat{i};
    
    phimat = phiall(:,:, inds_base);
    
    phimean = nan(size(phimat,1), size(phimat,2));
    for j=1:size(phimat,1)
        for jj=1:size(phimat,2)
       
            phimean(j, jj) = circ_mean(squeeze(phimat(j, jj, :)));
        end
    end
    
    PhiMeanAll(:,:, i) = phimean;
    
    if (0)
       figure; hold on;
       imagesc(phimean');
    end
    
end


%% ========== PLOT MEAN PHI
figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

t=PARAMS.tbins;
f = PARAMS.ffbins;
clim = [-pi pi];

% =========== GLOBAL MEAN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('global mean');

phimean = nan(size(PhiMeanAll,1), size(PhiMeanAll,2));
for j=1:size(PhiMeanAll,1)
    for jj=1:size(PhiMeanAll,2)
        phimean(j, jj) = circ_mean(squeeze(PhiMeanAll(j, jj, :)));
    end
end

imagesc(t, f, phimean', clim);
colorbar;
axis tight;


% =========== GLOBAL MEAN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('global mean [smaller clim range]');

phimean = nan(size(PhiMeanAll,1), size(PhiMeanAll,2));
for j=1:size(PhiMeanAll,1)
    for jj=1:size(PhiMeanAll,2)
        phimean(j, jj) = circ_mean(squeeze(PhiMeanAll(j, jj, :)));
    end
end

imagesc(t, f, phimean', [-1.5 1.5]);
colorbar;
axis tight;


% =========== GLOBAL MEAN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('global mean [smaller clim range]');

phimean = nan(size(PhiMeanAll,1), size(PhiMeanAll,2));
for j=1:size(PhiMeanAll,1)
    for jj=1:size(PhiMeanAll,2)
        phimean(j, jj) = circ_mean(squeeze(PhiMeanAll(j, jj, :)));
    end
end

imagesc(t, f, phimean', [-0.5 0.5]);
colorbar;
axis tight;


% ============ PLOT ANGLE SPECTRUM AT DIFFERENT TIMEPOINTS
timebins = PARAMS.tbins(1):0.02:PARAMS.tbins(end);
tbininds = discretize(t, timebins);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('angle spectra at diff time bins');
plotcols = lt_make_plot_colors(length(timebins)-1, 1, [1 0 0]);

for i=1:length(timebins)-1
    indsthis = tbininds==i;
    
    phivec = circ_mean(phimean(indsthis, :));
    
    plot(f, phivec, '-', 'Color', plotcols{i});
    
    lt_plot_text(f(end), phivec(end), ['t=' num2str(timebins(i)) ' to ' num2str(timebins(i+1))], plotcols{i});
    
end
lt_plot_zeroline;
    ylim([-pi pi]);
    
if (0) % OBSOLETE - above plots all on one plot. here plots them on separate plots.
% ============ PLOT ANGLE SPECTRUM AT DIFFERENT TIMEPOINTS
timebins = PARAMS.tbins(1):0.02:PARAMS.tbins(end);
tbininds = discretize(t, timebins);
for i=1:length(timebins)-1
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    indsthis = tbininds==i;
    title(['t=' num2str(timebins(i)) ' to ' num2str(timebins(i+1))]);
    
    phivec = circ_mean(phimean(indsthis, :));
    
    plot(f, phivec, '-k');
    lt_plot_zeroline;
    ylim([-pi pi]);
    
end
end

%% ===== [PHI] CALCULATE THINGS

% tbintoget = [11 12];
% fbintoget = [3 4];
tbintoget = [11]; % NOTE: currently only support single index. 
fbintoget = [4];

% ============= ORDER MATTERS FOR PHI;
assert(all(strcmp(OUTSTRUCT.bregionpairs_originalorder, 'LMAN-RA')), 'assumes that is LMAN-->RA');


% ========== 1) PLV AND MEAN PHI
OUTSTRUCT.Phi_mean_BaseWN = nan(length(OUTSTRUCT.PhiMat), 2);
OUTSTRUCT.Phi_PLV_BaseWN = nan(length(OUTSTRUCT.PhiMat), 2);
for i=1:length(OUTSTRUCT.PhiMat)
   disp(i);
    phimat = OUTSTRUCT.PhiMat{i};
    inds_base = OUTSTRUCT.indsbase_epoch{i};
    inds_WN = OUTSTRUCT.indsWN_epoch{i};
    
    % ========== 1) PLV
    % - BASE
    indsthis = inds_base;
    
    phivec= squeeze(phimat(tbintoget, fbintoget, indsthis));
    [plv, mean_cart, mean_pol, bootstats] = lt_neural_QUICK_PhaseLockVal(phivec, 1);
    
    OUTSTRUCT.Phi_mean_BaseWN(i,1) = mean_pol;
    OUTSTRUCT.Phi_PLV_BaseWN(i,1) = plv;

    % - BASE
    indsthis = inds_WN;
    
    phivec= squeeze(phimat(tbintoget, fbintoget, indsthis));
    [plv, mean_cart, mean_pol, bootstats] = lt_neural_QUICK_PhaseLockVal(phivec, 1);
    
    OUTSTRUCT.Phi_mean_BaseWN(i,2) = mean_pol;
    OUTSTRUCT.Phi_PLV_BaseWN(i,2) = plv;
    
end


%% ==== [QUICK] plots PLV across target syls. [dirty, since not controlled for sampel size]
indsthis = OUTSTRUCT.istarg==1;
figure; hold on; 
lt_plot_histogram(OUTSTRUCT.Phi_PLV_BaseWN(indsthis,1), [], 1, 1, [], 1, 'k');
lt_plot_histogram(OUTSTRUCT.Phi_PLV_BaseWN(indsthis,2), [], 1, 1, [], 1, 'r');
figure; hold on;
plot([1 2], [OUTSTRUCT.Phi_PLV_BaseWN(indsthis,1) OUTSTRUCT.Phi_PLV_BaseWN(indsthis,2)], '-ok');
xlim([0 3]);


%% === [QUICK] Plots random data, single trial phi matrix.

lt_figure; hold on;
colormap('parula')
for j=1:8
    
i = randi(200);
ii = randi(750);
try
phimat = OUTSTRUCT.PhiMat{ii}(:,:,i);
catch err
i = randi(200);
ii = randi(750);
end

lt_subplot(4,2,j); hold on;

imagesc(PARAMS.tbins, PARAMS.ffbins, phimat', [-pi pi]);
colorbar
axis tight;
end


%% ====== [SUMMARY PLOT] PHI DURING LEARNING
% 1) Changes in Phi consistency across trials
% 2) Changes in mean Phi (Wn vs. base)
% NOTE: Assume that order of data is LMAN-->RA (i.e. negative phi means LMAN leads). Code
% will make assertions so that this must be true.

close all;
cohdiff_norm_to_global = 1; % if 1, then finds (for a given channel)
rho_norm_to_global = 0;
centerdata = 0; % then for each switch, centers data (i.e. across channels)
% NOTE: this is useful if want to ask about, within a tgiven expt, is there
% correlation. 
onlyfirstswitch = 1; % if 0, then all siwtches.
useFFchangeInsteadOfLearnDir=0; % then sign(ff(WN)-ff(base)); if 0, then learndir

lt_neural_Coher_PhiSumPlot(OUTSTRUCT, SwitchStruct, cohdiff_norm_to_global, ...
    rho_norm_to_global, centerdata, onlyfirstswitch, useFFchangeInsteadOfLearnDir);



%% ====== [SUMMARY PLOT] SHUFFLE test of coherence change at target

% === 1) plot distribution of coh change and overlay tareget [each expt]



%% ===== PLOT PHASE DISTRIBUTIONS
close all;
figcount=1;
subplotrows=4;
subplotcols=5;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


assert(all(strcmp(OUTSTRUCT.bregionpair, 'LMAN-RA')), 'assumes chan1 is LAMN, chang 2 is RA...');
tbins = PARAMS.tbins;
ffbins = PARAMS.ffbins;

% ---- params
ffbintoplot = 4;
tbintoplot = 11;
nrand = 10; % random trials to plot
nanglebins = 16;
ntoplot = 10; % random channels to plot
equalizeN = 1; % equalizes sample size by random sample from that with more samples.
plotraw = 1; % if 1, then takes random subset ot plot. if 0, then takes all.

% === targets
indslist = find(OUTSTRUCT.istarg==1)';
if plotraw==1
indslist = indslist(randperm(length(indslist), ntoplot));
else
end

% ===== to collect
PLV_all = nan(length(indslist), 2); % [base, wn]

for j=1:length(indslist)
    indthis = indslist(j);
    
    trialsbase = OUTSTRUCT.indsbase_epoch{indthis};
    trialsWN = OUTSTRUCT.indsWN_epoch{indthis};
    phimat = OUTSTRUCT.PhiMat{indthis};    
    maxtrials = min([length(trialsbase) length(trialsWN)]);
    
    
    % ================== BASELINE
    phithis = squeeze(phimat(:, ffbintoplot, trialsbase));
    if equalizeN==1
        phithis = phithis(:, randperm(size(phithis,2), maxtrials));
    end
    basewn_label = 'BASELINE';
    
    if plotraw==1
    % -- 1) plot a few trials
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('random trials');
    ylabel(basewn_label);
    indrand = randperm(size(phithis,2), nrand);
    plot(tbins, phithis(:, indrand), '-x');
    axis tight;
    
    % -- 2) plot all trials
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['all trials, N=' num2str(size(phithis,2))]);
    ylabel('mean, std');
    plot(tbins, phithis, 'kx');
    plot(tbins, phithis+2*pi, 'kx');
%     phimean = mean(phithis,2);
%     phistd = std(phithis, [], 2);
%     lt_plot(tbins, phimean, {'Errors', phistd, 'Color', 'r'});
    axis tight
    ylim([-pi 2*pi])
    
    % -- 3) plot circular distribution of angles
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['all trials, N=' num2str(size(phithis,2))]);
    ylabel(['at tbin:' num2str(tbintoplot)]);
    theta = phithis(tbintoplot, :);
%     theta(theta<0) = theta(theta<0)+2*pi; % so that all are positive
    rose(theta);
   
    % ---- 4) show vector sum
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['all trials, N=' num2str(size(phithis,2))]);
    
    [x,y] = pol2cart(theta, ones(size(theta)));
    line([0 0], [mean(x) mean(y)], 'LineWidth', 3);
    xlim([-1 1]), ylim([-1 1]);
    grid on;
    
    
    % ---- get shuffle distribution of phase lock values.
    % DAT
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('plv (shuff, dat)');
    plv_dat = lt_neural_QUICK_PhaseLockVal(theta);
    % SHUFF
    plv_shuff_all = [];
    nshuff = 100;
    ntrials = length(theta);
    for k=1:nshuff
        % -- sample uniformly from circle.
        theta_shuff = 2*pi*rand(size(theta))-pi;
        plv_shuff = lt_neural_QUICK_PhaseLockVal(theta_shuff);
        plv_shuff_all = [plv_shuff_all; plv_shuff];
    end
    lt_plot_histogram(plv_shuff_all);
    line([plv_dat plv_dat], ylim, 'Color', 'r');
    xlim([-0.05 0.7]);
    end
    
    % ========
    theta = phithis(tbintoplot, :);
    plv_dat = lt_neural_QUICK_PhaseLockVal(theta);
    PLV_all(j, 1) = plv_dat;
    
    
    % ================== WN
    phithis = squeeze(phimat(:, ffbintoplot, trialsWN));
    if equalizeN==1
        phithis = phithis(:, randperm(size(phithis,2), maxtrials));
    end
    basewn_label = 'WN';
    
    if plotraw==1
    % -- 1) plot a few trials
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('random trials');
    ylabel(basewn_label);
    indrand = randperm(size(phithis,2), nrand);
    plot(tbins, phithis(:, indrand), '-x');
    axis tight;
    
    % -- 2) plot all trials
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['all trials, N=' num2str(size(phitmp,2))]);
    ylabel('mean, std');
    plot(tbins, phithis, 'kx');
    plot(tbins, phithis+2*pi, 'kx');
%     phimean = mean(phithis,2);
%     phistd = std(phithis, [], 2);
%     lt_plot(tbins, phimean, {'Errors', phistd, 'Color', 'r'});
    axis tight
    ylim([-pi 2*pi])
    
    % -- 3) plot circular distribution of angles
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['all trials, N=' num2str(size(phithis,2))]);
    ylabel(['at tbin:' num2str(tbintoplot)]);
    theta = phithis(tbintoplot, :);
%     theta(theta<0) = theta(theta<0)+2*pi; % so that all are positive
    rose(theta);
   
    % ---- 4) show vector sum
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['all trials, N=' num2str(size(phithis,2))]);
    
    [x,y] = pol2cart(theta, ones(size(theta)));
    line([0 0], [mean(x) mean(y)], 'LineWidth', 3);
    xlim([-1 1]), ylim([-1 1]);
    grid on;
    
    
    % ---- get shuffle distribution of phase lock values.
    % DAT
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('plv (shuff, dat)');
    plv_dat = lt_neural_QUICK_PhaseLockVal(theta);
    % SHUFF
    plv_shuff_all = [];
    nshuff = 100;
    ntrials = length(theta);
    for k=1:nshuff
        % -- sample uniformly from circle.
        theta_shuff = 2*pi*rand(size(theta))-pi;
        plv_shuff = lt_neural_QUICK_PhaseLockVal(theta_shuff);
        plv_shuff_all = [plv_shuff_all; plv_shuff];
    end
    lt_plot_histogram(plv_shuff_all);
    line([plv_dat plv_dat], ylim, 'Color', 'r');
    xlim([-0.05 0.7]);
    end
    
    theta = phithis(tbintoplot, :);
    plv_dat = lt_neural_QUICK_PhaseLockVal(theta);
    PLV_all(j, 2) = plv_dat;

end

lt_figure; hold on;
x = [1 2];
plot(x, PLV_all', '-ok');
lt_plot(x+0.2, mean(PLV_all,1), {'Errors', lt_sem(PLV_all), 'Color', 'r'});
xlim([0 3]);
xlabel('base - WN');
ylabel('PLV (equal N)');
signrank(PLV_all(:,1), PLV_all(:,2))

%% ===== PLOT LFP FOR MULTIPLE TRIALS OF ALL MOTIFS

close all;



%% ===== COMPUTE SPECTROGRAMS USING EXTRACTED

