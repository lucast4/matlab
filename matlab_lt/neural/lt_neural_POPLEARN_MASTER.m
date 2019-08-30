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




%% ==== GET PARAMS
clear PARAMS;

PARAMS.motif_predur = MOTIFSTATS_Compiled.birds(1).MOTIFSTATS.params.motif_predur;
PARAMS.motif_postdur = MOTIFSTATS_Compiled.birds(1).MOTIFSTATS.params.motif_postdur;
PARAMS.alignbyonset = MOTIFSTATS_Compiled.birds(1).MOTIFSTATS.neurons(1).motif(1).Params.REGEXP.alignByOnset;
assert(~isempty(PARAMS.alignbyonset));
% PARAMS.savemarker = '14Oct2018_2147';
PARAMS.savemarker = input('what is save marker? (e.g. 14Oct2018_2147)? ', 's');



%% ==== get times of WN hits

timewindhit = [-0.15 0.1]; % window, rel syl onset, during which any WN onset or offset will be considered a hit.
% leave expty [] to use pre and postdur from data extraction.

MOTIFSTATS_Compiled = lt_neural_POPLEARN_ExtrWNtime(timewindhit, SummaryStruct,...
    MOTIFSTATS_Compiled, PARAMS);

%% =============== [WN hits] - resave structure

% ==== SAVE
save(['/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_' ...
    PARAMS.savemarker '.mat'], ...
    'MOTIFSTATS_Compiled', '-v7.3');


%% ==== REMOVE DIR SONG
MOTIFSTATS_Compiled = lt_neural_QUICK_MotCom_RemoveDIR(MOTIFSTATS_Compiled);


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
BirdsToPlot = {'pu69wh78', 'wh44wh39', 'wh72pk12', 'gr48bu5'};
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


%% #################  [GOOD] [PLOT EXAMPLE RAW TRIALS]

birdplot = 'pu69wh78';
exptplot = 'RALMANlearn1';
swplot = 1;
motifplot = []; % [string] leave blank for target


% #####################################################
i = find(strcmp({SummaryStruct.birds.birdname}, birdplot));
ii = find(strcmp({SummaryStruct.birds(i).exptnum_pop.exptname}, exptplot));
ss = swplot;

indsthis = find(OUTSTRUCT.bnum==i & OUTSTRUCT.enum==ii & OUTSTRUCT.switch==ss & ...
    OUTSTRUCT.istarg==1);

% --- target syl motif number?
mm = unique(OUTSTRUCT.motifnum(indsthis));
motifplot = OUTSTRUCT.motifname{indsthis(1)};
if length(mm)>1
    disp('MULTIPLE TARGS, taking first one');
    mm = mm(1);
    indsthis = find(OUTSTRUCT.bnum==i & OUTSTRUCT.enum==ii & OUTSTRUCT.switch==ss & ...
        OUTSTRUCT.istarg==1 & OUTSTRUCT.motifnum==mm);
    motifplot = OUTSTRUCT.motifname{indsthis(1)};
end
% assert(length(mm)==1, 'multipel targets?');


% ======== 1) EXTRACT RAW DATA (NOT JUST LFP)
%     close all;

%     birdplot = 'pu69wh78';
%     exptplot = 'RALMANlearn1';
%     swplot = 1;
%     motifplot = []; % [string] leave blank for target
extrapad = 0.05; % seconds, pre and post...
[DatAll, t_onoff, fs, bregionlist, chanlist_toget, i, ii, mm] = ...
    lt_neural_LFP_PlotEgRaw_Extract(PARAMS, SwitchStruct, MOTIFSTATS_pop, ...
    SummaryStruct, SwitchCohStruct, birdplot, exptplot, swplot, motifplot, ...
    extrapad);


% ######################3 EXTRACT DATA AND PLOT
% chanpairs = OUTSTRUCT.chanpair(indsthis, :);
% 
% tvals = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).tvals;
% ffvals = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).ffvals;
% cohscal_allpairs = OUTSTRUCT.cohscal(indsthis);
% indsbase = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).indsbase_epoch;
% indsWN = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).indsWN_epoch;
% indsbase_all = OUTSTRUCT.indsbase{indsthis(1)};
% indsWN_all = OUTSTRUCT.indsWN{indsthis(1)};



% ================== FIND TRIAL TO PLOT
% trialmode = 'random';
trialmode = 'highcoh';
% trialmode = 'lowcoh';
N = 5;

ntrials = size(DatAll{1},2);

% ---- list of all possible trials
indsthis = find(OUTSTRUCT.bnum==i & OUTSTRUCT.enum==ii & OUTSTRUCT.switch==ss & ...
    OUTSTRUCT.istarg==1);
triallist = [find(OUTSTRUCT.indsbase{indsthis(1)}) find(OUTSTRUCT.indsWN{indsthis(1)})];

% ---- what is vector of coherence
cohscal = OUTSTRUCT.cohscal(indsthis);
cohscal = mean(cell2mat(cohscal),1); % take average over all chan pairs

[~, indsort] = sort(cohscal);
indsort = indsort(ismember(indsort, triallist));
indslow = indsort(1:N);
indshi = indsort(end-N+1:end);
    
    
if strcmp(trialmode, 'random')
    trialstoget = triallist(randperm(length(triallist), N));
elseif strcmp(trialmode, 'highcoh')
    trialstoget = indshi;
elseif strcmp(trialmode, 'lowcoh')
    trialstoget = indslow;
end

% ======================================= PLOT ALL THINGS COMBINED
% _------------ MAKE SURE TO FIRST HAVE COHEROGRAMS...
close all;
filt_low = 20;
filt_hi = 35;
% filt_low = 3200;
% filt_hi = 3300;
% trialtoplot = 100;
plotUnfiltered = 0; % defualt 0, plots high and low pass. if 1, then plots unfiltered instead of high pass.
lt_neural_LFP_PlotEgRaw_v2(DatAll, t_onoff, fs, bregionlist, chanlist_toget, ...
    i, ii, swplot, mm, trialstoget, filt_low, filt_hi, SwitchCohStruct, ...
    SwitchStruct, OUTSTRUCT, MOTIFSTATS_pop, PARAMS, ...
    SummaryStruct, OUTSTRUCT_CohMatOnly, plotUnfiltered);


% ================================= PLOT HEAT MAP OF ALL TRIALS
% _------------ MAKE SURE TO FIRST HAVE COHEROGRAMS...
close all;
filt_low = 20;
filt_hi = 35;
noiseband = [3100 3400];
plotUnfiltered = 0; % defualt 0, plots high and low pass. if 1, then plots unfiltered instead of high pass.
trialOrder = 'coh';
% trialOrder = 'trials'; % in order of actual trials
cohwind_f = PARAMS.cohscal_fwind;

lt_neural_LFP_PlotRawHeatmap(DatAll, t_onoff, fs, bregionlist, chanlist_toget, ...
    i, ii, swplot, mm, filt_low, filt_hi, SwitchCohStruct, ...
    SwitchStruct, OUTSTRUCT, MOTIFSTATS_pop, PARAMS, ...
    SummaryStruct, OUTSTRUCT_CohMatOnly, plotUnfiltered, trialOrder, ...
    noiseband, cohwind_f);


%% =========== [do above] but iterate over all expts (i.e. all targ syls)

lt_neural_LFP_PlotEgRaw_v2_allbirds;

%% ################## [GOOD, RAW PLOTS]
lt_neural_Coher_AllRaw;


%% #################### [DIAGNOSTIC] - display times of songs for each experiment

for i=1:length(SwitchCohStruct.bird)
    bname = SwitchStruct.bird(i).birdname;
   for ii=1:length(SwitchCohStruct.bird(i).exptnum)
       ename = SwitchStruct.bird(i).exptnum(ii).exptname;
       for iii=1:length(SwitchCohStruct.bird(i).exptnum(ii).switchlist)
           if isempty(SwitchCohStruct.bird(i).exptnum(ii).switchlist(iii).motifnum)
               continue
           end
          indsbase = SwitchCohStruct.bird(i).exptnum(ii).switchlist(iii).motifnum(1).indsbase;
          indswn = SwitchCohStruct.bird(i).exptnum(ii).switchlist(iii).motifnum(1).indsWN;
           
          disp(' ');
          disp('===================================== ');
          disp([bname ' - ' ename ' - sw' num2str(iii)]);
          t = SwitchCohStruct.bird(i).exptnum(ii).switchlist(iii).motifnum(1).tvals;
          t = [t(indsbase) t(indswn)];
          disp([datestr(t(1), 'ddmmmyyyy-HHMM') ' to ' datestr(max(t), 'ddmmmyyyy-HHMM')]);
       end
   end
end


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
% 1) all data, 2/5/19, up to gr48, RALMANLearn6
clear all; close all;
load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_05Feb2019_2142.mat');
load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_05Feb2019_2142.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_05Feb2019_2142.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/05Feb2019_2142/COHSTRUCT.mat');


% 1) gr48 - RALMANLearn7 ONLY
clear all; close all;
load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_26Mar2019_2352.mat');
load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_26Mar2019_2352.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_26Mar2019_2352.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/26Mar2019_2352/COHSTRUCT.mat');

% 1) gr48 - RALMANLearn9 ONLY
clear all; close all;
load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_30Mar2019_2048.mat');
% load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_26Mar2019_2352.mat');
% load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_26Mar2019_2352.mat');
% load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/26Mar2019_2352/COHSTRUCT.mat');


% 1) MULTI-DAY LEARNING ONLY
clear all; close all;
load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_28Mar2019_1106.mat');
load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_28Mar2019_1106.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_28Mar2019_1106.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/28Mar2019_1106/COHSTRUCT.mat');

% ===========================

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


% 1) all data  [LATEST, 4 birds, up to gr48, RALMANLearn3]
clear all; close all;
load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_20Jan2019_2134.mat');
load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_20Jan2019_2134.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_20Jan2019_2134.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/20Jan2019_2134/COHSTRUCT.mat');


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

% 4) all fake data (3 birds, curated) JAN 6 - JUST WH72 RALMANLEARN6
clear all; close all;
load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_03Jan2019_0153.mat');
% load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_03Jan2019_0153.mat');
% load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_03Jan2019_0153.mat');
% load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/03Jan2019_0153/COHSTRUCT.mat');


% 4) all fake data (for up to gr48 - RALMANLearn6) - LATEST BIRDS
clear all; close all;
load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_16Feb2019_0037.mat');
load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_16Feb2019_0037.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_16Feb2019_0037.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/16Feb2019_0037/COHSTRUCT.mat');

% 4) gr48bu5 - RALMANLearn2
clear all; close all;
load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_17Jan2019_1814.mat');
load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_17Jan2019_1814.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_17Jan2019_1814.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/17Jan2019_1814/COHSTRUCT.mat');


% 4) gr48bu5 - RALMANLearn2, RALMANLearn3 [18 Jan, around 1am]
clear all; close all;
load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_18Jan2019_0115.mat');
load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_18Jan2019_0115.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_18Jan2019_0115.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/18Jan2019_0115/COHSTRUCT.mat');

% 4) gr48bu5 - RALMANLearn4
clear all; close all;
load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_23Jan2019_1515.mat');
load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_23Jan2019_1515.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_23Jan2019_1515.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/23Jan2019_1515/COHSTRUCT.mat');


% 4) gr48bu5 - RALMANLearn3 [18 Jan, around 1am]
clear all; close all;
load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_23Jan2019_1624.mat');
load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_23Jan2019_1624.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_23Jan2019_1624.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/23Jan2019_1624/COHSTRUCT.mat');


% 4) gr48bu5 - RALMANLearn6 
clear all; close all;
load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_01Feb2019_1705.mat');
load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_01Feb2019_1705.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_01Feb2019_1705.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/01Feb2019_1705/COHSTRUCT.mat');

% 4) gr48bu5 - RALMANLearn5
clear all; close all;
load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_01Feb2019_2024.mat');
load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_01Feb2019_2024.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_01Feb2019_2024.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/01Feb2019_2024/COHSTRUCT.mat');

% ====== update params to match COHSTRUCT
settmp = 1;
if isempty(COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif)
    settmp=2;
    assert(~isempty(COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif));
end
PARAMS.tbins = COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif(1).t_relons;
PARAMS.ffbins = COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif(1).ffbins;
PARAMS.ffbinsedges = [20 35 80 130]; % edges, to plot timecourse in frequency bands

%% #################### [SHORTCUT - PREPROCESING]
% combines all things would do in this script into its own script. 
% Just run this after loading new data.

disp('CHECK SCRIPT - some things not auto correct');
pause;
lt_neural_Coher_PreProcess;

lt_neural_Coher_PreProcess_multiday; % for multiday learning 
% only difference is loads a different switchstruct.

%% ##################### [SCRIPT WITH LIST OF STUFF TO DO IN ORDER]

lt_neural_Coher_ScriptInOrder;


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
removeBadTrials = 1; % will only affect epoch inds.
SwitchCohStruct = lt_neural_Coher_LearnExtr2(COHSTRUCT, MOTIFSTATS_pop, ...
    SwitchStruct, pairtoget, LFPSTRUCT, PARAMS, baseuseallinds, removeBadTrials);    
end

PARAMS.bregionpair_toget = pairtoget;


%% ######################## [WN HITS] CHECK 1) fraction and 2)timing 

close all;
lt_neural_LFP_WNCheck(SwitchCohStruct, MOTIFSTATS_pop, SwitchStruct, ...
    PARAMS);

%% ######################## [WN LATENCY] Coherence change temporally aligned to WN
% - comapred to aligned to syl onset...
% - need to have extracted WN first before runnign this.
close all;

prctiletouse = 2.5;
dozscore=0;

prewind_relWN = [-0.1 -0.05]; % rel the percentiel you want to use.
onlyusegoodtargsyls = 1; % default = 1, i.e. gets timing only from that that will actually analyze.

% =================== 2) RUN AANLYSIS./
lt_neural_LFP_AlignToWN(OUTSTRUCT, SwitchCohStruct, MOTIFSTATS_pop, ...
    SwitchStruct, PARAMS, prctiletouse, dozscore, prewind_relWN, ...
    onlyusegoodtargsyls);


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
% twind = [-0.02 0]; % all combined
fwind = [22 32];

% ============= EXTRACT RELATIVE TO WN ONSET?
useWNtiming=0;
WNprctile = 2.5;
prewind_relWN = [-0.1 -0.05]; % rel the percentiel you want to use.

% === linearly interpolate coherence?
interpol =0; % across time  [NOT YET DONE].

onlyusegoodtargsyls = 1; % default = 1, i.e. gets timing only from that that will actually analyze.

[SwitchCohStruct, PARAMS] = lt_neural_LFP_PitchCorr(COHSTRUCT, SwitchCohStruct,...
    SwitchStruct, PARAMS, twind, fwind, useWNtiming, WNprctile, prewind_relWN, ...
    interpol, onlyusegoodtargsyls);

%% ======= PLOT ALL DATA [COH SCALAR vs FF]
close all;
switchestoplot = [1];
birdstoplot = {'pu69wh78'};

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
zscoreLFP = 3; % default 1, z-scores each t,ff bin separately.
% if 2, then doesn't zscore, instead normalizes as power proprotion
% (separately for each time slice and trial).
% if 3, then each ff window z-scored separately (good to see moudlation).
% if 4, then first 1) zscores within f, and 2) normalizes to all f (i.e.
% proportion
collectDiffMats = 0; % if 1, then collects differences (WN minus base). redundant, so leave at 0.
removeBadChans = 1; % default: 1
removeBadSyls = 1; % LEAVE AT 1.
typesToRemove = {'wn'}; % only remove syls that are bad becuase preceded by WN
% typesToRemove = {'wn', 'noise'};
saveSpec = 0; % saves spectrograms ...

if (0)
[OUTSTRUCT, OUTSTRUCT_CohMatOnly] = lt_neural_LFP_Learn_Extr(SwitchStruct, SwitchCohStruct, ...
    plotON, averagechanpairs, PARAMS, onlyfirstswitch, removeBadSyls, ...
    collectAllProcess, zscoreLFP, collectDiffMats, removeBadChans, typesToRemove);
end

[OUTSTRUCT] = lt_neural_LFP_Learn_Extr(SwitchStruct, SwitchCohStruct, ...
    plotON, averagechanpairs, PARAMS, onlyfirstswitch, removeBadSyls, ...
    collectAllProcess, zscoreLFP, collectDiffMats, removeBadChans, typesToRemove, ...
    saveSpec);


%% ===== get wn timing.
onlyusegoodtargsyls = 1;
prctiletouse = 2.5;
useallwn = 1;
[OUTSTRUCT] = lt_neural_Coher_GetWNtiming(OUTSTRUCT, onlyusegoodtargsyls, ...
    useallwn, prctiletouse, SwitchStruct, SwitchCohStruct, PARAMS);

%% ====== don't want raw data?
clear SwitchCohStruct;
clear OUTSTRUCT_CohMatOnly



%% ############ [FILTER BANK - EXTRACT]
% Extracts for each motif filter bank data in the range dsried. 
% First step before any anaylsis using filter banks
% analogous to LFPSTRUCT extraction above.

disp('NOTE: must first run lt_neural_tools_LFPfiltbank in MasterScript_v2');
pause; 

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

% == ssave to disk?
saveON =1;
savename = 'cohere_150ms'; % will save specifically for this PARAMS.savemarker

% === overwrite OUTSTRUCT?
overwriteOUT = 0; 

% === norm by subtracting shift predictor?
normtype = 'shifted';
% normtype = '';

onlygoodexpt = 1;

[OUTSTRUCT, PARAMS, LFPXCORR_Base, LFPXCORR_WN, LFPXCORR_freqsall] =  ...
    lt_neural_LFP_RecalcCoh(OUTSTRUCT, SwitchCohStruct, LFPSTRUCT, SwitchStruct, ...
    PARAMS, usecorrcoeff, thingstodo, saveON, savename, overwriteOUT, normtype, ...
    onlygoodexpt);


%% ======== [RECALCULATE COHERENCE, LOAD PREVIOUSLY CALCUALTED VERSIONS]

cd(['/bluejay0/bluejay2/lucas/analyses/neural/COHERENCE/RecalcCoh/' PARAMS.savemarker])

tmp = load('savestruct_cohere_PlusShuff.mat');
% tmp = load('savestruct_cohere_150ms.mat');

OUTSTRUCT.CohMean_Base = tmp.savestruct.dat.CohMean_Base;
OUTSTRUCT.CohMean_WN = tmp.savestruct.dat.CohMean_WN;
OUTSTRUCT_CohMatOnly = tmp.savestruct.dat.CohAlltrials;
OUTSTRUCT_CohMatOnly_shift = tmp.savestruct.dat.CohAlltrials_shuff;

PARAMS.tbins_old = PARAMS.tbins;
PARAMS.ffbins_old = PARAMS.ffbins;
PARAMS.tbins = tmp.savestruct.dat.t;
PARAMS.ffbins = tmp.savestruct.dat.f;

clear tmp;

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


%% ========= [COHERENCE MATRIX] REALIGN BY WN ONSET
% Realign all cohernece matrices by WN time of target.
useallwn = 1; % default is 0 (just epoch) but for some experiments not enough data in those epoch?
[OUTSTRUCT, PARAMS] = lt_neural_Coher_RealignbyWN(OUTSTRUCT, SwitchCohStruct, ...
    SwitchStruct, PARAMS, useallwn);


%% ==================== [RECALCULATE COHERENCE MATRIX and SCLALAR
% ===== ALSO RECALCULATE LFP
% BASED ON NEW SET OF TRIALS]

% ===== params for recalc of sclar
% twind = [-0.07 -0.03]; % all combined
% fwind = [22 32];

twind = [-0.07 -0.03]; % all combined
fwind = [22 36];

% twind = [-0.05 -0.03]; % all combined
% fwind = [22 32];
twind = [-0.07 -0.03]; % all combined
fwind = [22 60];

extractLFP = 0;
% need to have previously run Learn_Extr with extraction of LFP
lfpUseMedian = 1;
% specScaleType = 3; % see zscore above.
specScaleType = 1; % see zscore above.

% - EXTRACT RELATIVE TO WN ONSET?
useWNtiming=1;
WNprctile = 2.5;
prewind_relWN = [-0.1 -0.05]; % rel the percentiel you want to use.

% wntouse = 'half';
wntouse = 'quarter';
% wntouse = 'third';
% wntouse = 'firsthalf';
% wntouse = 'half';

% RemoveIfTooFewTrials =1; % default
RemoveIfTooFewTrials = 1;

% ========= remove bad syl
% removebadsyl=1; % default
removebadsyl=1;
removebadtrialtype = 'lfp'; % remove cases which by eye looked like large shared noise fluctuations...
% removebadtrialtype = 'all';
% removebadtrialtype = 'spikes';
% removebadtrialtype = '';

% --- norm to shift?
normtoshuff = 0;
normtype = 'minus';
% normtype = 'minusglob'; % then uses mean across base and train...
% normtype = 'zscore';
% normtype = 'zscorerelbase'; % if this, then must set normtoshuff to 0 [also, will only do this for scalar]
% normtype = ''; % if this, then must set normtoshuff to 0 [also, will only do this for scalar]

% == coherence use median? (across trials)
cohUseMedian =0; % NOT READY - doesn't do anything.

if normtoshuff==1
[OUTSTRUCT, PARAMS] = lt_neural_Coher_RecalcMat(SwitchStruct, OUTSTRUCT, OUTSTRUCT_CohMatOnly, ...
    SwitchCohStruct, PARAMS, twind, fwind, wntouse, useWNtiming, ...
    prewind_relWN, COHSTRUCT, RemoveIfTooFewTrials, removebadtrialtype, extractLFP, lfpUseMedian, specScaleType, ...
    cohUseMedian, removebadsyl, normtoshuff, normtype, OUTSTRUCT_CohMatOnly_shift);
else    
[OUTSTRUCT, PARAMS] = lt_neural_Coher_RecalcMat(SwitchStruct, OUTSTRUCT, OUTSTRUCT_CohMatOnly, ...
    SwitchCohStruct, PARAMS, twind, fwind, wntouse, useWNtiming, ...
    prewind_relWN, COHSTRUCT, RemoveIfTooFewTrials, removebadtrialtype, extractLFP, lfpUseMedian, ...
    specScaleType, cohUseMedian, removebadsyl, normtoshuff, normtype);
end


% ===== USE THIS VERSION IF ONLY GETTING LFP
if (0)
[OUTSTRUCT, PARAMS] = lt_neural_Coher_RecalcMat(SwitchStruct, OUTSTRUCT, '', ...
    SwitchCohStruct, PARAMS, twind, fwind, wntouse, useWNtiming, ...
    prewind_relWN, COHSTRUCT, RemoveIfTooFewTrials, removebadtrialtype, extractLFP, lfpUseMedian, specScaleType, ...
    cohUseMedian, removebadsyl, normtoshuff, normtype);
end

%% ====== EXTRACT SCALARS FOR LFP POWER CHANGE




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
birdstoplot = [];
expttoplot = [];
% swtoplot = [1 7 9 11];
swtoplot = [];
useAbsVal = 0; 

% --- to get specific switch types. ... [is done in addition to above
% fitlers]
% swtoget = {}; % passes if matches ANY of these
swtoget = {[0 1], [0 -1]}; % passes if matches ANY of these
firstswitchfortarget_withinday = 1; % if 1, then onlky keeps if all targets 
% for a given switch did not have a previous switch on the same day
% swtoget = {[1 0], [-1 0]}; % passes if matches ANY of these
% swtoget = {[1 0], [-1 0], [1 -1], [-1 1]}; % passes if matches ANY of these
% firstswitchfortarget_withinday = 1; % if 1, then onlky keeps if all targets 
% for a given switch did not have a previous switch on the same day
firstswitchofday=1;
handRemoveLFP = 1;
indtoget_b_e_s = lt_neural_LEARN_FilterSwitches(SwitchStruct, swtoget, ...
    firstswitchfortarget_withinday, firstswitchofday, [], handRemoveLFP);

% NOTE: to get raw cohgram (before subtract) cd /blucurrently need to do
% breakpoint using spec and evaluate cohgram version instead. Should
% modify to plot cohgram.
timewindowtoplot = [-0.08 0]; % for spectra.
ffbinsedges = [20 35 80 130]; % edges, to plot timecourse in frequency bands
lt_neural_Coher_Learn_PlotSum2(OUTSTRUCT, PARAMS, SwitchStruct, sumplottype, ...
    plotAllSwitchRaw, clim, fieldtoplot, birdstoplot, timewindowtoplot, ...
    zscoreLFP, expttoplot, swtoplot, ffbinsedges, indtoget_b_e_s, useAbsVal);


% ================ ONE PLOT FOR EACH SWITCH
close all;
fieldtoplot = 'coher';
% 'coher'
% 'spec'
% NOTE: to get raw cohgram (before subtract) currently need to do
% breakpoint using spec and evaluate cohgram version instead. Should
% modify to plot cohgram.
timewindowtoplot = [-0.08 0]; % for spectra.
birdstoplotTMP = [4];
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
            sumplottypeTMP = 'chanpairs'; % i.e. what is datapoint?
            lt_neural_Coher_Learn_PlotSum2(OUTSTRUCT, PARAMS, SwitchStruct, sumplottypeTMP, ...
                plotAllSwitchRaw, clim, fieldtoplot, birdstoplot, timewindowtoplot, ...
                zscoreLFP, expttoplot, swtoplot, ffbinsedges);
            
            birdname = SwitchStruct.bird(i).birdname;
            exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
            lt_subtitle([birdname '-' exptname '-s' num2str(ss)]);
            
        end
    end
end




%% ################### SUMMARY PLOT OF PITCH LAERNING
close all
onlygoodexpt = 1;
lt_neural_POPLEARN_Coher_PitchLearn(OUTSTRUCT, SwitchStruct, PARAMS, ...
    onlygoodexpt);

%% ################# [MOTIFPLOT - MAIN] SUMMARIZE SCALAR CHANGE IN LEARNING
% I.E. EACH syl in order...
% ======= SUMMARIZE SCALAR RESULT, OVER ALL SWITCHES, MOTIFS, AND CHANNELS

close all;

swtoget = {[0 1], [0 -1]}; % passes if matches ANY of these
% swtoget = {[1 0], [-1 0]}; % passes if matches ANY of these
firstswitchfortarget_withinday =1; % if 1, then onlky keeps if all targets 
% for a given switch did not have a previous switch on the same day
firstswitchofday =1;
indtoget_b_e_s = lt_neural_LEARN_FilterSwitches(SwitchStruct, swtoget, ...
    firstswitchfortarget_withinday, firstswitchofday);
% indtoget_b_e_s = [];

lt_neural_Coher_SumPlotMotifs(OUTSTRUCT, SwitchStruct, MOTIFSTATS_Compiled, PARAMS, ...
    indtoget_b_e_s);


%% ################# [MOTIFPLOT - SPECTRAL POWER] 
% SUMMARIZE SCALAR CHANGE IN LEARNING
% I.E. EACH syl in order...
% ======= SUMMARIZE SCALAR RESULT, OVER ALL SWITCHES, MOTIFS, AND CHANNELS

close all;

swtoget = {[0 1], [0 -1]}; % passes if matches ANY of these
% swtoget = {[1 0], [-1 0]}; % passes if matches ANY of these
firstswitchfortarget_withinday =1; % if 1, then onlky keeps if all targets 
% for a given switch did not have a previous switch on the same day
firstswitchofday =1;
handRemoveLFP = 1;
indtoget_b_e_s = lt_neural_LEARN_FilterSwitches(SwitchStruct, swtoget, ...
    firstswitchfortarget_withinday, firstswitchofday, [], handRemoveLFP);
% indtoget_b_e_s = [];

fieldtoplot = 'specdiff_chan1_all';
lt_neural_Coher_SumPlotMotifs_Spec(OUTSTRUCT, SwitchStruct, MOTIFSTATS_Compiled, PARAMS, ...
    indtoget_b_e_s, fieldtoplot);


%% ===================== [COMPARE RA AND LMAN FOR CHANGE IN LFP POWER]




%% ################## [GOOD] [SCALAR & PITCH LERANING- COMPARE TO LEARNING RATE]
close all;
nsegs = 4; % quartiles, then 4... (for both coh and learning)
learnConvertToRate = 1; % then learning (z) is z dividided by mean time of bin (hours).
lt_neural_Coher_PitchLearnCoh(OUTSTRUCT, PARAMS, SwitchCohStruct, ...
    SwitchStruct, nsegs, learnConvertToRate);



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
% statmethod = 'diff'; % diff minus mean across all motifs in each expt.
statmethod = 'minDiffType'; % minsu diff type.
% statmethod = 'minusbase'; % simply minus base.
meanOverExpt = 1; % then for each motif gets one value for each status type (targ, same ..);

birdstoplot = [];
expttoplot = [];
% swtoplot = [1 7 9 11];
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

%% ====== [COHSCALAR CORRELATIONS] EXTRACT BASELINE [and WN] CORRELATIONS.
useonlybaseepoch = 0; % if 1, then just limited epoch. if 0, then entie baseline data
% default: 0;
corrtype = 'spearman'; % spearman. or pearson or Kendall
cohdiff_usedprime = 0; % if 0, then just mean diff. if 1, then dprime
nboot = 15; % to get bootstrap SE

OUTSTRUCT = lt_neural_Coher_CohScalExtract(OUTSTRUCT, useonlybaseepoch, ...
    corrtype, cohdiff_usedprime, nboot);


%% ======= [COHSCALAR CORR] OVERVIEW PLOT - 
% === each switch, each motif, each chan pair, pre and post WN
close all;
onlyplotgood=1;
lt_neural_Coher_CohScalMotifPlot(OUTSTRUCT, SwitchStruct, MOTIFSTATS_Compiled, PARAMS, ...
    onlyplotgood);


%% ===== [COHSCALAR CORRELATIONS - SUMMARY PLOT] LEARNING CHANGE IN COHSCALAR AS FUNCTION OF BASE CORR? 
% i.e. can I predict the change in coherence in a given expt (switch) based
% on relationship between baseline coh/ff and direction of training
% NOTE: can test variation within and between experiments (centerdata)
% NOTE: can test both normative learn dir and actual ff change (useFF ...);
close all;
cohdiff_norm_to_global = 1; % if 1, then finds (for a given channel pair) all motifs.
rho_norm_to_global = 0; % same as above, but for rho.
centerdata = 0; % then for each switch, centers data (i.e. across channels)
% NOTE: this is useful if want to ask about, within a tgiven expt, is there
% correlation 
onlyfirstswitch = 0; % if 0, then all siwtches.
useFFchangeInsteadOfLearnDir=0; % then sign(ff(WN)-ff(base)); if 0, then learndir

% === only initianl elarning
onlygoodexpt=1;

lt_neural_Coher_CohScalPlot(OUTSTRUCT, SwitchStruct, cohdiff_norm_to_global, ...
    rho_norm_to_global, centerdata, onlyfirstswitch, useFFchangeInsteadOfLearnDir, ...
    onlygoodexpt);

%% ===== [COHSCALAR] 


%% ###################################################################
%% ################ COHERENCE SCALAR PLOTS [SYSTEMATIC]
corrtype = 'spearman'; % spearman. or pearson or Kendall
Nshuff = 500;
saveON = 0;

% GOES THRU ALL T, FF BINS
% REQUIRES OUTSTRUCT_CohMatOnly
tbins = PARAMS.tbins;
fbins = PARAMS.ffbins;

tbinstokeep = [-0.15 0.06];
fbinstokeep = [15 125];

FFcorrCoh = nan(length(tbins), length(fbins), length(OUTSTRUCT.bnum));
FFcorrCoh_pval = nan(length(tbins), length(fbins), length(OUTSTRUCT.bnum));
FFcorrCoh_pctileVsShuff = nan(length(tbins), length(fbins), length(OUTSTRUCT.bnum));
FFcorrCoh_zscoreVsShuff = nan(length(tbins), length(fbins), length(OUTSTRUCT.bnum));
FFcorrCoh_shuffCI = nan(length(tbins), length(fbins), 2, length(OUTSTRUCT.bnum));

savedir = '/bluejay5/lucas/analyses/neural/COHERENCE/SCALAR';
savedir = [savedir '/' PARAMS.savemarker];
if ~exist(savedir)
    mkdir(savedir);
end

tic
for i=1:length(OUTSTRUCT.bnum)
    
    inds_base = OUTSTRUCT.indsbase{i};
    cohmatall = OUTSTRUCT_CohMatOnly{i};
    ff = OUTSTRUCT.ffvals{i};
    
    if all(isnan(ff))
        disp('SKIP! no defined FF');
        continue
    end
    
    % ================ GO THRU EACH T, FF BIN
    for j=1:length(tbins)
        if tbins(j)<tbinstokeep(1) | tbins(j)>tbinstokeep(2)
            continue
        end
            
        for jj=1:length(fbins)
        if fbins(jj)<fbinstokeep(1) | fbins(jj)>fbinstokeep(2)
            continue
        end
        
        disp(['case ' num2str(i) ', bin1=' num2str(j) ', bin2=' num2str(jj)]);
            
            cohthis = squeeze(cohmatall(j, jj, inds_base));
            ffthis = ff(inds_base);
            
            % ---------------- GET CORR
            [rho_dat, pval] = corr(ffthis', cohthis, 'type', corrtype);
            
            % ---------------- GET SHUFFLE DISTRIBUTION OF CORR
            rho_shuff_all = nan(1, Nshuff);
            nrends = length(ffthis);
            for nn=1:Nshuff
                %                 disp(['case ' num2str(i) ', bin1=' num2str(j) ', bin2=' num2str(jj) ', shuff' num2str(nn)]);
                
                % --- shuffle trials
                ffthis_shuff = ffthis(randperm(nrends));
                rho_shuff = corr(ffthis_shuff', cohthis, 'type', corrtype);
                rho_shuff_all(nn) = rho_shuff;
            end
            %             rho_shuff_all = [];
            %             for nn=1:Nshuff
            % %                 disp(['case ' num2str(i) ', bin1=' num2str(j) ', bin2=' num2str(jj) ', shuff' num2str(nn)]);
            %
            %                 % --- shuffle trials
            %                 ffthis_shuff = ffthis(randperm(length(ffthis)));
            %                 rho_shuff = corr(ffthis_shuff', cohthis, 'type', corrtype);
            %                 rho_shuff_all = [rho_shuff_all; rho_shuff];
            %             end
            
            % ================= OUTPUT STATS
            FFcorrCoh(j, jj, i) = rho_dat;
            FFcorrCoh_pval(j, jj, i) = pval;
%             keyboard
            p = (sum(abs(rho_shuff_all)>=abs(rho_dat))+1)./(length(rho_shuff_all)+1);
            FFcorrCoh_pctileVsShuff(j, jj, i) = p;
            
            FFcorrCoh_zscoreVsShuff(j, jj, i) = (rho_dat - mean(rho_shuff_all))/std(rho_shuff_all);
            FFcorrCoh_shuffCI(j, jj, :, i) = prctile(rho_shuff_all, [2.75 97.5]);
            
%             if p<0.05 & rho_dat>FFcorrCoh_shuffCI(j, jj, 1, i) & rho_dat<FFcorrCoh_shuffCI(j, jj, 2, i)
%                 keyboard
%             end
%             
        end
    end
    if saveON ==1
        if mod(i, 10)==0
            save([savedir '/FFcorrCoh'], 'FFcorrCoh');
            save([savedir '/FFcorrCoh_pval'], 'FFcorrCoh_pval');
            save([savedir '/FFcorrCoh_pctileVsShuff'], 'FFcorrCoh_pctileVsShuff');
            save([savedir '/FFcorrCoh_zscoreVsShuff'], 'FFcorrCoh_zscoreVsShuff');
            save([savedir '/FFcorrCoh_shuffCI'], 'FFcorrCoh_shuffCI');
        end
    end
end
toc

% clear OUTSTRUCT_CohMatOnly;

%% ===================== 

tmp = FFcorrCoh_pctileVsShuff<0.05;


%% ====================== [COHSCALAR] - OVERVIEW PLOT
close all;
onlygood = 1;
plotEachSylChan = 0;
useshuffpval =1; % then asks if is outside shuffle [2.5 97.5]
lt_neural_Coher_COHSCALAR_Overview(OUTSTRUCT, FFcorrCoh, FFcorrCoh_pval, ...
    FFcorrCoh_pctileVsShuff, FFcorrCoh_zscoreVsShuff, FFcorrCoh_shuffCI, ...
    SwitchStruct, onlygood, PARAMS, plotEachSylChan, useshuffpval);

%% ====================== PLOT SHUFFLE SUMMARY [COH SCALAR]
indsgood = ~isnan(squeeze(FFcorrCoh(10,10,:))) & OUTSTRUCT.switch==1;

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

% ==== plot specific bins
lt_figure; hold on;
tmp = squeeze(FFcorrCoh(11,4,:));
tmp = tmp(~isnan(tmp));
[~, xcenters] = lt_plot_histogram(tmp, '', 1, 0, [], 1, 'k');

tmp = squeeze(FFcorrCoh(8,13,:));
tmp = tmp(~isnan(tmp));
lt_plot_histogram(tmp, xcenters, 1, 0, [], 1, 'r');

lt_figure; hold on;
tmp = squeeze(FFcorrCoh(11,4,:));
tmp1 = tmp(~isnan(tmp));

tmp = squeeze(FFcorrCoh(8,13,:));
tmp2 = tmp(~isnan(tmp));
plot(tmp1, tmp2, 'kx');
lt_regress(tmp1, tmp1, 1, 0, 1, 1, 'r');
lt_plot_makesquare_plot45line(gca, 'b');

%% ###################################### PHI
%% ==== MEAN PHI ACROSS ALL SYLS

PhiMeanAll = nan(length(PARAMS.tbins), length(PARAMS.ffbins), length(OUTSTRUCT.bnum));
PhiMeanAll_WN = nan(length(PARAMS.tbins), length(PARAMS.ffbins), length(OUTSTRUCT.bnum));

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
    
    
    
    % ======= collect mean WN
    indsWN = OUTSTRUCT.indsWN_epoch{i};
    phiall = OUTSTRUCT.PhiMat{i};
    
    phimean = nan(size(phimat,1), size(phimat,2));
    for j=1:size(phiall,1)
        for jj=1:size(phiall,2)
            phimean(j, jj) = circ_mean(squeeze(phiall(j, jj, indsWN)));
        end
    end
    PhiMeanAll_WN(:,:, i) = phimean;
end



%% ============ [PHI] COMPUTE SCALAR BY AVERAGING WITHIN A BIN
tbins = [-0.07 -0.03];
fbins = [20 35];

inds_t = PARAMS.tbins>tbins(1) & PARAMS.tbins<tbins(2);
inds_f = PARAMS.ffbins>fbins(1) & PARAMS.ffbins<fbins(2);

% ============= BASELINE
% tmp = squeeze(mean(mean(PhiMeanAll(inds_t, inds_f, :),2),1));
tmp = squeeze(circ_mean(circ_mean(PhiMeanAll(inds_t, inds_f, :),[],2),[],1));
OUTSTRUCT.PhiScalar = tmp;


% ============= WN
% tmp = squeeze(mean(mean(PhiMeanAll(inds_t, inds_f, :),2),1));
tmp = squeeze(circ_mean(circ_mean(PhiMeanAll_WN(inds_t, inds_f, :),[],2),[],1));
OUTSTRUCT.PhiScalar_WN = tmp;

%% =========== [PHI - MOTIF PLOT]

% I.E. EACH syl in order...
% ======= SUMMARIZE SCALAR RESULT, OVER ALL SWITCHES, MOTIFS, AND CHANNELS

close all;

swtoget = {[0 1], [0 -1]}; % passes if matches ANY of these
% swtoget = {[1 0], [-1 0]}; % passes if matches ANY of these
firstswitchfortarget_withinday =1; % if 1, then onlky keeps if all targets 
% for a given switch did not have a previous switch on the same day
firstswitchofday =1;
indtoget_b_e_s = lt_neural_LEARN_FilterSwitches(SwitchStruct, swtoget, ...
    firstswitchfortarget_withinday, firstswitchofday);
% indtoget_b_e_s = [];

lt_neural_Coher_PHI_SumPlotMotifs(OUTSTRUCT, SwitchStruct, MOTIFSTATS_Compiled, PARAMS, ...
    indtoget_b_e_s);

%% ========== [PHI] PLOT MEAN PHI
close all;
birdtoplot = [];
datlevel = 'motif';
lt_neural_Coher_PHI_PlotSum(OUTSTRUCT, PARAMS, PhiMeanAll, birdtoplot);

%% ===== [PHI] CALCULATE THINGS

% tbintoget = [11 12];
% fbintoget = [3 4];
tbintoget = [11]; % NOTE: currently only support single index.
fbintoget = [4];

useSingleBinVersion = 0; % if 1, then uses multi bin mean extgracted above. otehrwise uses single bin

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
    
    if useSingleBinVersion==1
        OUTSTRUCT.Phi_mean_BaseWN(i,1) = mean_pol;
    else
        OUTSTRUCT.Phi_mean_BaseWN(i,1) = OUTSTRUCT.PhiScalar(i);
    end
    OUTSTRUCT.Phi_PLV_BaseWN(i,1) = plv;
    
    % - WN
    indsthis = inds_WN;
    
    phivec= squeeze(phimat(tbintoget, fbintoget, indsthis));
    [plv, mean_cart, mean_pol, bootstats] = lt_neural_QUICK_PhaseLockVal(phivec, 1);
    
    if useSingleBinVersion==1
        OUTSTRUCT.Phi_mean_BaseWN(i,2) = mean_pol;
    else
        OUTSTRUCT.Phi_mean_BaseWN(i,2) = OUTSTRUCT.PhiScalar_WN(i);
    end
    OUTSTRUCT.Phi_PLV_BaseWN(i,2) = plv;
    
end


%% ==== [QUICK] plots PLV across target syls. [dirty, since not controlled for sampel size]
indsthis = OUTSTRUCT.istarg==1;
figure; hold on; 
lt_plot_histogram(OUTSTRUCT.Phi_PLV_BaseWN(indsthis,1), [], 1, 1, [], 1, 'k');
lt_plot_histogram(OUTSTRUCT.Phi_PLV_BaseWN(indsthis,2), [], 1, 1, [], 1, 'r');

figure; hold on;
Y = [OUTSTRUCT.Phi_PLV_BaseWN(indsthis,1) OUTSTRUCT.Phi_PLV_BaseWN(indsthis,2)];
plot([1 2], [OUTSTRUCT.Phi_PLV_BaseWN(indsthis,1) OUTSTRUCT.Phi_PLV_BaseWN(indsthis,2)], '-ok');
xlim([0 3]);
p = signrank(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'srank', 1);

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
% cohdiff_norm_to_global = 1; % if 1, then finds (for a given channel)
% rho_norm_to_global = 0;
% centerdata = 0; % then for each switch, centers data (i.e. across channels)
% NOTE: this is useful if want to ask about, within a tgiven expt, is there
% correlation. 
% onlyfirstswitch = 1; % if 0, then all siwtches.
useFFchangeInsteadOfLearnDir=0; % then sign(ff(WN)-ff(base)); if 0, then learndir
onlygood = 1;

lt_neural_Coher_PhiSumPlot(OUTSTRUCT, SwitchStruct, ...
    useFFchangeInsteadOfLearnDir, onlygood);



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


%% #################################################################


% ======= GET XCOV OF DIFFERNET TYPES:
% 1) SPIKE XCOV (SHUFFLE SUBTRACTED)
% 2) XCOV OF SMOOTHED FR
% 3) SPIKE XCOV (NOT SHUFFLE SUBTRACTED)
% 4) SPIKE XCOV (USING FREQUENCY DOMAIN)

BirdsToSkip = {'wh72pk12'};

% ====== binning spikes:
binsize_spk = 0.0025; % default, 5ms bins for cross corr
% xcov_dattotake = [-0.12 0.02]; % rel syl onset.
% xcov_dattotake = [-0.05 0.02]; % rel syl onset.
xcov_dattotake = [-0.06 0]; % rel syl onset. [BEST WINDOW]
% xcov_dattotake = [-0.08 0.02]; % rel syl onset.
% xcovwindmax = 0.06; % seconds
xcovwindmax = 0.03; % seconds
normmethod = 'unbiased';
% normmethod = 'coeff';

% ======= norm bined spikes to prob in bin? (i.e. sum to 1)
normspiketoprob = 0; % if 0, then uses spike count

% ======= bregionpairtoget
bregionpairtoget = 'LMAN-RA';
% bregionpairtoget = 'LMANoutside-RAoutside';

% ======== to remove units that were extracted only for LFP
removeIfLFPOnly = 1;

% ====== FOR GETTING RUNNING XCOV
getXgram=1;
windsize = 0.06;
windshift = 0.005;
% windsize = 0.08;
% windshift = 0 = {};.005;
% ---

goodexptonly = 1;

% ============== CONTROL - TEMPORAL JITTER
% SEE Smith & Kohn, JNeurosci 2008
% a jitter window, shuffles spikes across trials. breaks synchrony, but
% keeps other stats.
dojitter = 1;
jitterwindSize = 0.05;

% ================= RUN
[SwitchXCovStruct, PARAMS] = lt_neural_POPLEARN_XcovCalc(BirdsToSkip, ...
    binsize_spk, xcov_dattotake, xcovwindmax, normmethod, normspiketoprob, ...
    bregionpairtoget, removeIfLFPOnly, getXgram, windsize, windshift, PARAMS, ...
    SwitchCohStruct, SwitchStruct, MOTIFSTATS_pop, SummaryStruct, goodexptonly, ...
    dojitter, jitterwindSize);


%% ====================== [XCOV] EXTRACTION

if isfield(PARAMS, 'Xcov_ccLags_beforesmooth')
   % --- reset lag values to before smoothing... (so that can redo
   % extraction here)
    PARAMS.Xcov_ccLags=PARAMS.Xcov_ccLags_beforesmooth;
end

% ====== usealltrial
% --- option 1: [this trumps others]
usealltrials = 0; % otehrwise uses epchs.
% --- option 2: [only does this if usealltrials=0]
useallbase = 0; % 1: all; 0=last half.
% wntouse = 'all';
% wntouse = 'half';
% wntouse = 'third';
wntouse = 'quarter';

% ============== smoothe in lag space? 
dosmooth = 1; % interpolate, then smooth
dosmooth_sigma = 2*PARAMS.Xcov.binsize_spk; % sec

% ===== what version of xcov?
xcovver = 'zscore'; % i.e. each lag bin, zscore relative to shuffle.
% xcovver = 'scale'; % % NOT DONE YET! scales by relative magnitude of shuffle functions.
% xcovver = ''; % empty means difference of means.
% xcovver= 'coherency';

% ===== get xcovgram?
getxgram = 1;
getxgram_epochbins = 4; % if empty, then ignore. otherwise divides training into this many even bins

% ====== removebad syl?
removebadsyl = 1;
removebadtrials =1;
removebadchans = 1;
plotraw = 0; % will keyboard on targets, and can run to see how compute xcov [i.e. different methods]

% ======= AUTO PARAMS
if strcmp(xcovver, 'coherency')
    dosmooth = 0; % have not coded this yet, should take smoothing AFTER get coherency.
end
    
% ===== split [base WN] and [epochs] into hi and low ff trials.
getHiLoFFSplit = 1; 

% ====
hilosplit_shuffver=1;

% =============== RUN
[OUTSTRUCT_XCOV, PARAMS, NanCountAll] = lt_neural_POPLEARN_XcovExtr(SwitchXCovStruct, ...
    SwitchStruct, PARAMS, SwitchCohStruct, OUTSTRUCT, usealltrials, ...
    useallbase, dosmooth, dosmooth_sigma, removebadsyl, ...
    PARAMS.Xcov.Xgram.windlist, plotraw, xcovver, wntouse, removebadtrials, ...
    getxgram, removebadchans, getxgram_epochbins, getHiLoFFSplit, hilosplit_shuffver);

assert(all(PARAMS.xcenters_gram == mean(PARAMS.Xcov.Xgram.windlist,2)))

% ================ plot nans
lt_figure; hold on;
title('trials with nan for xcorr [real dat]');
xlabel('N total');
ylabel('N that is nan');
plot(NanCountAll(:,2), NanCountAll(:,1), 'ok');


% ========== POST
OUTSTRUCT_XCOV = lt_neural_POPLEARN_XCovExtr_post(OUTSTRUCT_XCOV, ...
    SwitchCohStruct, SummaryStruct, SwitchStruct);

%% =================== [EXTRACT] SPLIT BY FF (HIGH VS. LOW)
if (0) % NOW INCORPORATED INTO THE ABOUT OUTSTRUCT CODE
% ==== baseline and WN splits intow high and low FF
% ==== also does for each bin of epoch, if epochs exist.
% ===== MUST RUN THIS AFTER GETTING OUTSTRUCT_XCOV, since it determines
% whether use epoch bins.
normspiketoprob = 0; % must match previuos... above
removeIfLFPOnly = 1; % must match previuos... above

[SwitchXCovStruct, OUTSTRUCT_XCOV, PARAMS] = lt_neural_POPLEARN_XcovCalc_FFsplit(...
    SwitchXCovStruct, normspiketoprob, PARAMS, SwitchStruct, MOTIFSTATS_pop, ...
    SummaryStruct, SwitchCohStruct, removeIfLFPOnly);
end

%% ================= SMOOTH ALL XCOV
% NOTE: no need to run unless doing coherence xcovver above
% IF do run, then need to include epoch analyses here
close all;
dosmooth = 1; % interpolate, then smooth
dosmooth_sigma = 2*PARAMS.Xcov.binsize_spk; % sec

[OUTSTRUCT_XCOV, PARAMS] = lt_neural_POPLEARN_XcovSmth(...
    OUTSTRUCT_XCOV, PARAMS, dosmooth_sigma);



    
%% [SAVE] --- especially useful if computed covariance heat maps

savedir = ['/bluejay0/bluejay2/lucas/analyses/neural/POPLEARN/OUTSTRUCT_XCOV/' PARAMS.savemarker];
if ~exist(savedir, 'dir')
    mkdir(savedir);
end
save([savedir '/OUTSTRUCT_XCOV.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
save([savedir '/SwitchXCovStruct.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS.mat'], 'PARAMS', '-v7.3');

% === save variant
savemarker = 'zscore';
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
save([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3');

savemarker = '80mswind_zscore_smooth';
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
save([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3');

savemarker = '40mswind_zscore_smooth';
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
save([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3');

savemarker = 'zscore_25msbin';
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
save([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3');

savemarker = 'zscore_25msbin_smoothed';
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
save([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3');

savemarker = 'frprob';
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
save([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3');


% ============= BELOW: ALL USING 
% ==============
savemarker = '80mswind_zscore_smooth_021319'; % latest (clean chans, syls, trials).
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
save([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3');


% ========== 80ms, zscore, smoothed, using all base and 1/3 of WN trials.
savemarker = '80mswind_zscore_smooth_021319_2'; % latest (clean chans, syls, trials).
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
save([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3');

% ========== 100ms, zscore, smoothed, using 1/2 base and 1/3 of WN trials.
savemarker = '100mswind_zscore_smooth_021419'; % latest (clean chans, syls, trials).
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
save([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3');

% ========== 100ms, zscore, smoothed, using 1/2 base and 1/4 of WN trials
% [GOOD]
savemarker = '100mswind_zscore_smooth_021419_2'; % latest (clean chans, syls, trials).
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
save([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3');


% ========== 100ms, zscore, smoothed, using 1/2 base and 1/4 of WN trials
% [GOOD] [extending xcenter bins into syllable onset] [all up to gr48-
% Learn6]
savemarker = '100mswind_zscore_smooth_021619'; % latest (clean chans, syls, trials).
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
save([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3');


% ========== 100ms, zscore, smoothed, using 1/2 base and 1/4 of WN trials
% [GOOD] [extending xcenter bins into syllable onset] [all up to gr48-
% Learn6] [ALSO GETTING AUTOCOVARIANCE] [OUTSTRUCT is COHERENCY]
savemarker = '100mswind_zscore_smooth_021919'; % latest (clean chans, syls, trials).
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
save([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3');


savemarker = '100mswind_zscore_smooth_030519'; % latest, and 4 epochs for Xcovgram
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
% save([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3');


savemarker = '100mswind_zscore_smooth_030519_3bins'; % latest, and 3 epochs for Xcovgram
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
% save([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3');

savemarker = '100mswind_zscore_smooth_030519_5bins'; % latest, and 5 epochs for Xcovgram
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
% save([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3');

savemarker = '100mswind_zscore_smooth_030519_6bins'; % latest, and 6 epochs for Xcovgram
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
% save([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3');

savemarker = '100mswind_zscore_smooth_030519_2bins'; % latest, and 6 epochs for Xcovgram
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
% save([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3'

% =============== [60 m s bins] [3 epochs]
savemarker = '60mswind_zscore_smooth_031019_3bins'; % latest, and 6 epochs for Xcovgram
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
save([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3')

savemarker = '60mswind_zscore_smooth_031019_5bins'; % latest, and 6 epochs for Xcovgram
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
% save([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3')


% =============== [60 m s bins] [3 epochs] [hi lo ff split]
savemarker = '60mswind_zscore_smooth_031019_3bins_hiloFF';  % latest, and 6 epochs for Xcovgram
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
% save([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3')


% ############################### [3/19/19 - LATEST, 60MS, after check
% rasteres and raw neural and all clearned). [3 epochs; hi lo split]
savemarker = '60mswind_031919_3bins_hiloFF';  % latest, 
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
save([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3')

% ====================== 
savemarker = '60mswind_031919_3binsAll';  % 3 bins for both xcovgram and epochs
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
save([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3')

% ====================== 
savemarker = '60mswind_031919_4binsAll';  % 4 bins for both xcovgram and epochs
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3')


savemarker = '60mswind_031919_3bins_hiloFF_shuffver2';  % using shufflever = 2
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3')



% ############################### [3/19/19 - LATEST, 60MS, after check
% rasteres and raw neural and all clearned). [3 epochs; hi lo split]
% [JITTER: 20MS]
savemarker = '60mswind_031919_3bins_hiloFF_Jitter20ms';  % latest, and 6 epochs for Xcovgram
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
save([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3')

savemarker = '60mswind_031919_3bins_hiloFF_Jitter50ms';  % latest, and 6 epochs for Xcovgram
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3')


% ############################### [3/25/19 - 60ms bins - gr48bu5 -
% RALMANLearn7]
savemarker = 'gr48bu5_RALMANLearn7';  % latest, and 6 epochs for Xcovgram
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
save([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3')

% ############################### [3/29/19 - 60ms bins - gr48bu5 -
% RALMANLearn9]
savemarker = 'gr48bu5_RALMANLearn9';  % latest, and 6 epochs for Xcovgram
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
save([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3')

% ############################### [3/28/19 - 60ms bins - multi-day learning
savemarker = 'multiday_032819';  % latest, and 6 epochs for Xcovgram
save([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV', '-v7.3');
save([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct', '-v7.3');
save([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS', '-v7.3')


%% [LOAD] --- 
% NOTE: the same switchxcov struct could apply for adifferent outstructs.
if (0)
savedir = ['/bluejay0/bluejay2/lucas/analyses/neural/POPLEARN/OUTSTRUCT_XCOV/' PARAMS.savemarker];
load([savedir '/OUTSTRUCT_XCOV.mat']);
load([savedir '/SwitchXCovStruct.mat']);
load([savedir '/PARAMS.mat']);
end


if (0)
    
savemarker = '40mswind_zscore_smooth';
load([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat']);
load([savedir '/SwitchXCovStruct_' savemarker '.mat']);
load([savedir '/PARAMS_' savemarker '.mat']);

savemarker = '60mswind_zscore_smooth_031019_3bins';
load([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat']);
load([savedir '/SwitchXCovStruct_' savemarker '.mat']);
load([savedir '/PARAMS_' savemarker '.mat']);

% =============== [60 m s bins] [3 epochs] [hi lo ff split]
savemarker = '60mswind_zscore_smooth_031019_3bins_hiloFF';
load([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat']);
% load([savedir '/SwitchXCovStruct_' savemarker '.mat']);
load([savedir '/PARAMS_' savemarker '.mat']);

% ================ [60ms bins - after recurating data][3 epochs] [hi lo ff split]
% "quarter" for xcov gram, 4 bins for epoch
% ################### GOOD
savemarker = '60mswind_031919_3bins_hiloFF';  % latest, and 6 epochs for Xcovgram
load([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat']);
load([savedir '/SwitchXCovStruct_' savemarker '.mat']);
load([savedir '/PARAMS_' savemarker '.mat']);

% ====================== [SAME, but 3 bins for both xcovgram and epochs]
savemarker = '60mswind_031919_3binsAll';  % 
load([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV');
load([savedir '/SwitchXCovStruct_' savemarker '.mat'], 'SwitchXCovStruct');
load([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS')

% ====================== [SAME, but 4 bins for both xcovgram and epochs]
savemarker = '60mswind_031919_4binsAll';  % 4 bins for both xcovgram and epochs
load([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat'], 'OUTSTRUCT_XCOV');
load([savedir '/PARAMS_' savemarker '.mat'], 'PARAMS')


% ############################### [3/19/19 - LATEST, 60MS, after check
% rasteres and raw neural and all clearned). [3 epochs; hi lo split]
% [JITTER: 20MS]
savemarker = '60mswind_031919_3bins_hiloFF_Jitter20ms';  % latest, and 6 epochs for Xcovgram
load([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat']);
load([savedir '/SwitchXCovStruct_' savemarker '.mat']);
load([savedir '/PARAMS_' savemarker '.mat']);


% ================ [60ms bins - gr48bu5 - RALMANLearn7]
savemarker = 'gr48bu5_RALMANLearn7';  % latest, and 6 epochs for Xcovgram
load([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat']);
load([savedir '/SwitchXCovStruct_' savemarker '.mat']);
load([savedir '/PARAMS_' savemarker '.mat']);


end

if (0)
savemarker = '80mswind_zscore_smooth';
load([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat']);
load([savedir '/SwitchXCovStruct_' savemarker '.mat']);
load([savedir '/PARAMS_' savemarker '.mat']);

savemarker = '100mswind_zscore_smooth_021619';
load([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat']);
load([savedir '/SwitchXCovStruct_' savemarker '.mat']);
load([savedir '/PARAMS_' savemarker '.mat']);

savemarker = '100mswind_zscore_smooth_030519'; % latest, and 4 epochs for Xcovgram
load([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat']);
load([savedir '/SwitchXCovStruct_100mswind_zscore_smooth_021619.mat']);
load([savedir '/PARAMS_' savemarker '.mat']);

savemarker = '100mswind_zscore_smooth_030519_3bins'; % latest, and 3 epochs for Xcovgram
load([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat']);
load([savedir '/SwitchXCovStruct_100mswind_zscore_smooth_021619.mat']);
load([savedir '/PARAMS_' savemarker '.mat']);

savemarker = '100mswind_zscore_smooth_030519_5bins'; % latest, and 5 epochs for Xcovgram
load([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat']);
load([savedir '/SwitchXCovStruct_100mswind_zscore_smooth_021619.mat']);
load([savedir '/PARAMS_' savemarker '.mat']);

savemarker = '100mswind_zscore_smooth_030519_6bins'; % latest, and 5 epochs for Xcovgram
load([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat']);
load([savedir '/SwitchXCovStruct_100mswind_zscore_smooth_021619.mat']);
load([savedir '/PARAMS_' savemarker '.mat']);

savemarker = '100mswind_zscore_smooth_030519_2bins'; % latest, and 5 epochs for Xcovgram
load([savedir '/OUTSTRUCT_XCOV_' savemarker '.mat']);
load([savedir '/SwitchXCovStruct_100mswind_zscore_smooth_021619.mat']);
load([savedir '/PARAMS_' savemarker '.mat']);

end


%% ========== [REALIGN XCOVGRAMS TO WN ONSET]





%% =========== [SEE IF ALIGNT TO WN?]
% 2) Plot separating expereiments by timing.
close all;
prewind_relWN = [-0.1 -0.05]; % rel the percentiel you want to use.
prewind_relSyl = [-0.07 -0.03];
lt_neural_LFP_AlignToWN_Xcov(OUTSTRUCT_XCOV, OUTSTRUCT, SwitchStruct, ...
    PARAMS, prewind_relWN, prewind_relSyl);

%% ================== AD HOC COMBINING OF EPOCHS
% == PU69 RL2, since there is gap in trials, need to combine.
numepochs = size(OUTSTRUCT_XCOV.XcovgramWN_FFsplits_Epochs{1}, 3);
disp(['nume pochs = ' num2str(numepochs) '?']);
pause;
assert(size(OUTSTRUCT_XCOV.XcovgramWN_FFsplits_Epochs,2)==2, 'then need to modify code for OUTSTRUCT_XCOV.XcovgramWN_FFsplits_Epochs below');



btoget = find(strcmp({SummaryStruct.birds.birdname}, 'pu69wh78'));
etoget = find(strcmp({SwitchStruct.bird(btoget).exptnum.exptname}, 'RALMANlearn2'));

indsthis = find(OUTSTRUCT_XCOV.bnum==btoget & OUTSTRUCT_XCOV.enum==etoget);
 
for j=indsthis'

    % ====================
    tmp = OUTSTRUCT_XCOV.XcovgramWN_epochs{j};    
    tmp(:,:, numepochs) = nanmean(tmp, 3);
    tmp(:,:, 1:numepochs-1) = nan;
    OUTSTRUCT_XCOV.XcovgramWN_epochs{j} = tmp;
    
    % ====================
    tmp = OUTSTRUCT_XCOV.XcovgramWN_FFsplits_Epochs{j, 1};    
    tmp(:,:, numepochs) = nanmean(tmp, 3);
    tmp(:,:, 1:numepochs-1) = nan;
    OUTSTRUCT_XCOV.XcovgramWN_FFsplits_Epochs{j,1} = tmp;

    tmp = OUTSTRUCT_XCOV.XcovgramWN_FFsplits_Epochs{j,2};    
    tmp(:,:, numepochs) = nanmean(tmp, 3);
    tmp(:,:, 1:numepochs-1) = nan;
    OUTSTRUCT_XCOV.XcovgramWN_FFsplits_Epochs{j,2} = tmp;

    % ====================
    tmp = OUTSTRUCT_XCOV.trialedges_epoch{j};    
    tmp([numepochs numepochs+1]) = [tmp(1) tmp(end)];
    tmp(1:numepochs-1) = nan;
    OUTSTRUCT_XCOV.trialedges_epoch{j} = tmp;

end


%% ========= 2) EXTRACT SCALAR FOR TWO PEAKS FOR EACH EXPERIMENT
twindows = {[-0.06 -0.03], [-0.95 -0.035], [-0.06 -0.03]}; % one array for each window [x centers...] [inclusive]
lagwindows = {[-0.01 -0.005], [0.02 0.05], [0.005 0.01]};

% 5ms bins
twindows = {[-0.08 -0.03], [-0.08 -0.03]}; % one array for each window [x centers...] [inclusive]
lagwindows = {[-0.01 0], [0.03 0.045]};

% 2.5ms bins, 80ms, default.
twindows = {[-0.07 -0.03], [-0.08 -0.03]}; % one array for each window [x centers...] [inclusive]
lagwindows = {[-0.008 0.001], [0.028 0.047]};


% 2.5ms bin, smoothed, 80ms window, good?
twindows = {[-0.06 -0.03], [-0.08 -0.03]}; % one array for each window [x centers...] [inclusive]
lagwindows = {[-0.012 0.002], [0.028 0.047]};

% 2.5ms bins, 100ms, default.
twindows = {[-0.07 -0.03], [-0.07 -0.03]}; % one array for each window [x centers...] [inclusive]
lagwindows = {[-0.011 0.005], [0.03 0.045]};

lagwindows = {[-0.014 0.006], [0.03 0.045]}; % OK
% lagwindows = {[-0.008 0.001], [0.03 0.045]};

% lagwindows = {[-0.015 0.005], [0.026 0.048]};
% twindows = {[-0.06 -0.03], [-0.07 -0.03]}; % one array for each window [x centers...] [inclusive]

% ========= GOOD: [100ms window]
alignto = 'sylonset'; 
twindows = {[-0.07 -0.03], [-0.07 -0.03]}; % one array for each window [x centers...] [inclusive]
lagwindows = {[-0.014 0.006], [0.03 0.045]}; % OK

% ========= GOOD: [40ms window]
alignto = 'sylonset'; 
twindows = {[-0.08 -0.02], [-0.1 -0.01]}; % one array for each window [x centers...] [inclusive]
lagwindows = {[-0.01 0.005], [-0.014 0.006]}; % OK

% ========= GOOD: [60ms window]
% alignto = 'sylonset'; 
% twindows = {[-0.07 -0.03], [-0.1 -0.01]}; % one array for each window [x centers...] [inclusive]
% lagwindows = {[-0.01 0.005], [-0.014 0.006]}; % OK
alignto = 'sylonset'; 
twindows = {[-0.07 -0.03]}; % one array for each window [x centers...] [inclusive]
lagwindows = {[-0.01 0.005]}; % OK

% ========= GOOD: [60ms window]
% alignto = 'sylonset'; 
% twindows = {[-0.07 -0.03], [-0.1 -0.01]}; % one array for each window [x centers...] [inclusive]
% lagwindows = {[-0.01 0.005], [-0.014 0.006]}; % OK
alignto = 'sylonset'; 
twindows = {[-0.06 -0.03]}; % one array for each window [x centers...] [inclusive]
lagwindows = {[-0.009 0.006]}; % OK

% % ========= GOOD: [60ms window]
% alignto = 'sylonset'; 
% twindows = {[-0.09 -0.02], [-0.1 -0.01]}; % one array for each window [x centers...] [inclusive]
% lagwindows = {[-0.014 0.006], [-0.014 0.006]}; % OK
% 
% 
% alignto = 'wnonset'; % not good.
% twindows = {[-0.09 -0.06], [-0.09 -0.06]}; % one array for  0 each window [x centers...] [inclusive]
% lagwindows = {[-0.015 0.005], [0.03 0.05]}; % OK
% 
% alignto = 'wnonset'; % using median of WN
% twindows = {[-0.11 -0.06], [-0.11 -0.065]}; % one array for each window [x centers...] [inclusive]
% lagwindows = {[-0.015 0.005], [0.03 0.05]}; % OK

% NOTE: 
% pu69, overnight, max preceding: -0.08 (center)
% alignto = 'sylonset'; 
% alignto = 'wnonset';
OUTSTRUCT_XCOV = lt_neural_POPLEARN_XCov_ExtrScal(OUTSTRUCT_XCOV, OUTSTRUCT, ...
    PARAMS, twindows, lagwindows, alignto);


%% ========= 2) EXTRACT TIME SLICE (OVERWRITE ORIGINAL TIME SLICE...)
% NOTE: THIS IS THE ONLY THING CURRENTLY MODIFIED TO ALSO DO FOR EPOCH
% ANALYSES.

twindow = [-0.085 -0.04]; % DEFAULT, 80ms, smoothed... one array for each window [x centers...] [inclusive]
% twindow = [-0.08 -0.04]; % one array for each window [x centers...] [inclusive]
twindow = [-0.07 -0.03]; % one array for each window [x centers...] [inclusive]

% NOTE: 
% pu69, overnight, max preceding: -0.08 (center)
% alignto = 'sylonset'; 
alignto = 'wnonset';

% ======= GET ALL DATA (FOR BASELINE ANALYSIS ONLY)
twindow = [-0.07 0.005]; % one array for each window [x centers...] [inclusive]
alignto = 'sylonset'; 


% ======= GOOD [100ms]
twindow = [-0.07 -0.03]; % one array for each window [x centers...] [inclusive]
alignto = 'sylonset'; 

% ======= GOOD [40ms window]
twindow = [-0.1 -0.01]; % one array for each window [x centers...] [inclusive]
alignto = 'sylonset'; 

% ======= GOOD [60ms window]
twindow = [-0.06 -0.03]; % one array for each window [x centers...] [inclusive]
alignto = 'sylonset'; 
% twindow = [-0.09 -0.02]; % one array for each window [x centers...] [inclusive]
% alignto = 'sylonset'; 


% twindow = [-0.1 -0.05]; % one array for each window [x centers...] [inclusive]
% % twindow = [-0.09 -0.07]; % one array for each window [x centers...] [inclusive]
% alignto = 'wnonset';

OUTSTRUCT_XCOV = lt_neural_POPLEARN_XCov_ExtrSlice(OUTSTRUCT_XCOV, ...
    OUTSTRUCT, PARAMS, twindow, alignto);

%% ============ [REALIGN TO WN ONSET]

[OUTSTRUCT_XCOV, PARAMS] = lt_neural_POPLEARN_XCov_ReAlign(OUTSTRUCT_XCOV, ...
    OUTSTRUCT, PARAMS);


%% =========== [XCOV] PLOT EACH INDIVIDUAL EXPERIMENT

close all;
plotNotminShuff=0; % 1, then not shift subtracted
onlygoodexpt = 0;
clim = [-0.1 0.1];
lt_neural_POPLEARN_XCov_PlotAll(OUTSTRUCT_XCOV, SwitchStruct, PARAMS, ...
    plotNotminShuff, onlygoodexpt, clim);


%% ============ [XCOV] PLOT EACH SYL/CHAN IN EACH EXPT

close all;
plotNotminShuff=0; % 1, then not shift subtracted
onlygoodexpt = 0;
clim = [-0.1 0.1];
btoplot = 2;
etoplot = 2;
stoplot = 4;
lt_neural_POPLEARN_XCov_PlotEach(OUTSTRUCT_XCOV, SwitchStruct, PARAMS, ...
    plotNotminShuff, onlygoodexpt, clim, btoplot, etoplot, stoplot);



%% ============ [XCOV] PLOT EACH SYL/CHAN IN EACH EXPT [FF SPLIT VERSION]

close all;
plotNotminShuff=0; % 1, then not shift subtracted
onlygoodexpt = 0;
clim = [-0.1 0.1];
btoplot = 2;
etoplot = 1;
stoplot = 1;
epochstoplot = [1:3]; % will average over these.

lt_neural_POPLEARN_XCov_PlotEachSplit(OUTSTRUCT_XCOV, SwitchStruct, PARAMS, ...
    plotNotminShuff, onlygoodexpt, clim, btoplot, etoplot, stoplot, ...
    epochstoplot);


%% ========== [XCOV] SUMMARIZE FFSPLIT XCOV SLICE DIFFERENCES
close all;
epochstoplot = 3:4;
lt_neural_POPLEARN_XCov_FFsplitSlice(OUTSTRUCT_XCOV, PARAMS, epochstoplot)


%% =========== [XCOV] sUMMARIZE LEARNING

close all;
onlygoodexpt = 0;
lt_neural_POPLEARN_Xcov_PitchLearn(OUTSTRUCT, OUTSTRUCT_XCOV, SwitchStruct, ...
    PARAMS, onlygoodexpt);



%% ############################ [EPOCHS ANALYSIS - I.E. TRAINIGN BINS]
%% =========== [EXTRACTION] GET LEARNING INTO BINS
% === ALSO PLOTS SAMPLE SIZES...
onlygoodexpt = 1;
% windowprem = [-0.1 0.01]; % for counting spikes
windowprem = [-0.1 0]; % for counting spikes
OUTSTRUCT_XCOV = lt_neural_POPLEARN_Xcov_ExtrLearn(OUTSTRUCT, OUTSTRUCT_XCOV, ... 
    SwitchStruct, PARAMS, SwitchCohStruct, MOTIFSTATS_pop, windowprem);


%% =========== [EPOCH SPLIT STATS - PLOT DISTRIBUTIONS] [FRATE EFFECT ALSO]
    % 1) across epoch
    % 2) across song bout
    % 3) mean frate difference
    
close all;
onlygoodexpt = 1;
plotRaw=0;
epochWN = [3:4];
lt_neural_POPLEARN_Xcov_Epochs_Distr(OUTSTRUCT_XCOV, PARAMS, ...
    onlygoodexpt, SwitchStruct, plotRaw, epochWN)

%% =========== [XCOV SLICE EPOCHS] - change over epohcs
% must have exctracted epoch data above.
close all;
dattype = 'chan';
lt_neural_POPLEARN_Xcov_Epochs(OUTSTRUCT, OUTSTRUCT_XCOV, PARAMS, ...
    onlygoodexpt, SwitchStruct, dattype)

%% ============ [XCOV EPOCHS] - SCALAR, COMAPARE TO LEARNING (Z)
close all;
dattype = 'chan';
% dattype = 'switch';
scalwind = 1; % t/f window for scalar. (i.e which peak)
syltype = 1; %, 1,2,3 = targ, same, diff.
doFFsplit=0;

% casestokeep = 'goodlearn';
% casestokeep = 'goodneural';
% casestokeep = 'goodlearnneural';
% casestokeep = 'badlearn';
casestokeep = 'all';


% ========== TO LOOK AT OVERALL LEARNING/LEARNING RATE
% mintraindur = 0; % hours. (have to have data)
% mintotaltrain = 0; %q hours (end of data minus base) (doesn't have to actually have data throughotu)


% ========== TO LOOK AT WITHIN EXPT CHANGES
% mintraindur = 2.5; % hours. (have to have data)
% mintotaltrain = 2; % hours (end of data minus base) (doesn't have to actually have data throughotu)
mintraindur = 1; % hours. (have to have data)
mintotaltrain = 1; % hours (end of data minus base) (doesn't have to actually have data throughotu)
% mintraindur = 0; % hours. (have to have data)
% mintotaltrain = 0; % hours (end of data minus base) (doesn't have to actually have data throughotu)

% === to look at within changes, scaled bin-to-bin increments
% mintraindur = 2.5; % hours. (have to have data)
% mintotaltrain = 2; % hours (end of data minus base) (doesn't have to actually have data throughotu)
% casestokeep = 'goodlearnneural';

% doshift =0; % if 1, then asks whetehr neural in bin n predicts learn in bin n+1;
adhoc_replacelearnwithWind2 = 0; % NOTE: have not modified FFsplit code to run this.

% === NOTES ABOUT FFSPLIT
% - have not coded for dattype=switch
% - have not used if, so if do not have FF splits extracted, this entier code wil not run
% - all plots are AFTER removing too short experiments
% - to plot all the varaints (ie.. towards baseline bias/away bias) have to
% go in by hand and uncomment various thigns int he code (in the
% subscript...)
% - the "casestokeep" does not apply for ffsplit analyses.

% ---- FF SPLIT? USE THESE PARAMS
% FFsplit_pu69learn2_combine = 1; % combines all bins into one bin (since only looked at trials at end of learning, since noisy neural during learning)
% mintraindur = 0; % hours. (have to have data)
% mintotaltrain = 3; % hours (end of data minus base) (doesn't have to actually have data throughotu)
% doFFsplit=1;
% useAd_Nonad_Average_forBaseline = 0; % if 0, then (adaptive minus adaptive), (nonad - nonad)
% % if 1, then for each case subtract the baseline mean over non and ad)
% % DEFAULT: 0

lt_neural_POPLEARN_Xcov_EpochScal(OUTSTRUCT_XCOV, PARAMS, ...
    onlygoodexpt, SwitchStruct, dattype, scalwind, syltype, casestokeep, ...
    mintraindur, mintotaltrain, adhoc_replacelearnwithWind2, ...
    FFsplit_pu69learn2_combine, doFFsplit, useAd_Nonad_Average_forBaseline)






%% ================= [NEURAL BIAS VS DIRECTED SONG BIAS]
% ======== DOES DIRECTION OF BASELINE BIAS CORRELATE WITH BASELINE AFP
% BIAS?
close all;
dirstruct = lt_DirSong_MotifID;
useff_zscore = 1; % if 0, then mean FF diff, if 1, then dir zscore from undir. (each day, then take mean)
%  if 2, then uses percent change.
targsylonly = 0; % default 0.
cvdiffmethod = 'diff'; % cv undir minus cv dir
cvdiffmethod = 'percent'; % cv undir minus cv dir
onlygoodexpt = 1;
lt_neural_POPLEARN_Xcov_DirBias(OUTSTRUCT_XCOV, PARAMS, ...
    onlygoodexpt, SwitchStruct, dirstruct, useff_zscore, targsylonly,  ...
    cvdiffmethod);



%% #####################################################################
%% ============ [XCOV] PLOT SUMMARY ACROSS EXPERIMENTS
close all;
% datlevel = 'switch';
datlevel = 'neurpair';
fracchange = 0; % if 1, then frac change in xcov. if 0, then difference
plotNotminShuff=0;

% --- to get specific switch types. ... [is done in addition to above
% fitlers]
birdtoget = [1 2 4];
swtoget = {[0 1], [0 -1]}; % passes if matches ANY of these
firstswitchfortarget_withinday = 1; % if 1, then onlky keeps if all targets
% for a given switch did not have a previous switch on the same day
% swtoget = {[1 0], [-1 0]}; % passes if matches ANY of these
% firstswitchfortarget_withinday = 1; % if 1, then onlky keeps if all targets
% for a given switch did not have a previous switch on the same day
firstswitchofday=1;
indtoget_b_e_s = lt_neural_LEARN_FilterSwitches(SwitchStruct, swtoget, ...
    firstswitchfortarget_withinday, firstswitchofday, birdtoget);
% indtoget_b_e_s = [];

% clim = [-1.5e-4 1.5e-4];
clim = [-0.15 0.15];
clim = [-0.05 0.05];

% ======== align to WN or syl?
alignto = 'syl';
% alignto = 'wn';

lt_neural_POPLEARN_XCov_PlotSum(OUTSTRUCT_XCOV, SwitchStruct, PARAMS, ...
    indtoget_b_e_s, datlevel, fracchange, plotNotminShuff, clim, ...
    alignto);


lt_subtitle(['[' normmethod '], alltrials=' num2str(usealltrials)]);

%% =========== [XCOV] TIMECOURSES

% lagwindows = {[-0.01 0], [0.03 0.055]};


close all;
onlygoodexpt = 1;
% lagwindows = [-0.011 -0.001 0.029 0.056]; % exclusive... , will take all intervals.
% lagwindows = [-0.011 -0.001 0.029 0.046]; % exclusive... , will take all intervals.
% lagwindows = [-0.0126 0.001 0.029 0.046]; % 2.5ms windows
% lagwindows = [-0.011 0.001 0.027 0.045]; % 2.5ms, smoothed DEFAULT
lagwindows = [-0.015 0.007 0.029 0.046]; % OK % 100ms window
ffbinsedges_indstoplot = [1 3]; % which ones of lagwindow intervals to p[lot?

lagwindows = [-0.016 0.006]; % OK % 60ms window
ffbinsedges_indstoplot = [1]; % which ones of lagwindow intervals to p[lot?
dattype = 'chan';
% dattype = 'switch';
% clim = [-0.015 0.015];
plotindivswitch=0;
% clim = [-1.5e-4 1.5e-4];
clim = [-0.1 0.1];
XLIM = [-0.1 -0.04];
YLIMGRAM = [-0.05 0.05];
clim = [-0.3 0.3];
XLIM = [-0.1 -0.01];
YLIMGRAM = [-0.05 0.05];
plotlines = 1; % then plots each channel on the summary figure.

% === do ffsplit?
ffsplitparams = struct;
ffsplitparams.dosplit = 1; % on or off
ffsplitparams.epochtoplot = 3:4;
ffsplitparams.adaptivebin = 2; % usually 1=non-adaptive; 2=adaptive, 3= take difference

lt_neural_POPLEARN_Xcov_PlotTcourse(OUTSTRUCT_XCOV, SwitchStruct, ...
    onlygoodexpt, PARAMS, dattype, lagwindows, clim, ffbinsedges_indstoplot, ...
    plotindivswitch, XLIM, YLIMGRAM, plotlines, ffsplitparams);


%% =========== [XCOV] PLOT INVIDIDUAL EXPEIRMNTS, BUT XCOV GRAM
close all;
clim = [-0.15 0.15];
birdtoplot = [1];
lt_neural_POPLEARN_XCov_PlotAllGram(OUTSTRUCT_XCOV, SwitchStruct, PARAMS, ...
    clim, birdtoplot)


%% =========== [XCOV] PLOT EACH MOTIF/CHAN, FOR A GIVEN EXPERIMENT
% HEAT MAP

close all;

bb = 1;
ee = 1;
sw = 1;
clim = [-0.15 0.15];

bname = SwitchStruct.bird(bb).birdname;
ename = SwitchStruct.bird(bb).exptnum(ee).exptname;

disp('==============');
disp([bname ' --- ' ename ' --- sw' num2str(sw)]);

% === unique channels pairs
[indsgrp, indsgrpU] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, OUTSTRUCT_XCOV.enum, ...
    OUTSTRUCT_XCOV.switch, OUTSTRUCT_XCOV.chanpair, OUTSTRUCT_XCOV.motifnum});

figcount=1;
subplotrows=6;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for i=1:length(indsgrpU)
    indsthis = OUTSTRUCT_XCOV.bnum==bb & OUTSTRUCT_XCOV.enum==ee & ...
        OUTSTRUCT_XCOV.switch==sw & indsgrp==indsgrpU(i);
    if ~any(indsthis)
        continue
    end
    
    % ======= go thru each motif and plot
    assert(sum(indsthis)==1);
    xgram_base = OUTSTRUCT_XCOV.XcovgramBase{indsthis};
    xgram_wn = OUTSTRUCT_XCOV.XcovgramWN{indsthis};
    motif = OUTSTRUCT_XCOV.motifname{indsthis};
    issame = OUTSTRUCT_XCOV.issame(indsthis);
    istarg = OUTSTRUCT_XCOV.istarg(indsthis);
    if istarg==1
        ptit = 'TARG';
    elseif issame==1
        ptit = 'SAME';
    else
        ptit = 'DIFF';
    end
        
    % ============== 1) BASELINE
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('BASE');
    ylabel([motif '[' ptit ']']);
    lt_neural_Coher_Plot(xgram_base, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags,...
        1, '', clim);
    lt_plot_zeroline;
    line([-0.05 -0.05], ylim, 'Color', 'm');
    
    % ============== 2) WN
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('WN');
%     ylabel([motif '[' ptit ']']);
    lt_neural_Coher_Plot(xgram_wn, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags,...
        1, '', clim);
    lt_plot_zeroline;
    line([-0.05 -0.05], ylim, 'Color', 'm');
    
    % ============== 3) WN MINUS BASE
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('WN-base');
%     ylabel([motif '[' ptit ']']);
    lt_neural_Coher_Plot(xgram_wn-xgram_base, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags,...
        1, '', clim);
    lt_plot_zeroline;
    line([-0.05 -0.05], ylim, 'Color', 'm');
end

%% =========== [XCOV SCALAR] - plot each experiment [MOTIFS]
close all;
windthis =1;
Yscalar = cellfun(@(x)(x(windthis,2)-x(windthis,1)), OUTSTRUCT_XCOV.Xcovscal_window_BaseWN);
% Yscalar = cellfun(@(x)mean(x(windthis,:),2), OUTSTRUCT_XCOV.Xcovscal_window_BaseWN);
% clim = [-0.2 0.2];
clim = [-0.2 0.2];
% clim = [-1.5e-4 1.5e-4];
onlygoodexpt = 1;
expttype = 'xcov_spikes';
% plotlevel = 'switch';
plotlevel = 'chanpair';

useoldInds = 0; % default 0. only use 1 if you haven't extracted inds_base_epoch to OUTSTRUCT_XCOV yet...
% there will be error if you need to switch to 1.

lt_neural_POPLEARN_XCov_PlotScal(Yscalar, OUTSTRUCT_XCOV, SwitchCohStruct, SwitchStruct, ...
    MOTIFSTATS_Compiled, PARAMS, clim, onlygoodexpt, expttype, plotlevel, ...
    OUTSTRUCT, useoldInds);

%% ========= [GOOD - PUBLICATION] Linear mixed effects of change in xcov
close all;

windthis =1;
Yscalar = cellfun(@(x)(x(windthis,2)-x(windthis,1)), OUTSTRUCT_XCOV.Xcovscal_window_BaseWN);

clim = [-0.2 0.2];
onlygoodexpt = 1;
expttype = 'xcov_spikes';

onlyIfSameType=0; % if 1, then only expt that has at least one asme type.
% will also include same-type in model.

lt_neural_POPLEARN_XCov_LME(Yscalar, OUTSTRUCT_XCOV, ...
    SwitchStruct, MOTIFSTATS_Compiled, PARAMS, clim, ...
    onlygoodexpt, expttype, onlyIfSameType)

%% ========= [FR SCALAR, CHANGE FROM BASELINE]

epochtoplot = 3;
normtype = 'usefrac';
% useminus
j=1;

fr1_base = {};
fr2_base = {};
fr1_wn = {};
fr2_wn = {};

for j=1:length(OUTSTRUCT_XCOV.epochSplitStatsAll)
    fr1_base = [fr1_base; OUTSTRUCT_XCOV.epochSplitStatsAll{j}(1).nspksByNeuron(1)];
    fr2_base = [fr2_base; OUTSTRUCT_XCOV.epochSplitStatsAll{j}(1).nspksByNeuron(2)];
    fr1_wn = [fr1_wn; OUTSTRUCT_XCOV.epochSplitStatsAll{j}(epochtoplot+1).nspksByNeuron(1)];
    fr2_wn = [fr2_wn; OUTSTRUCT_XCOV.epochSplitStatsAll{j}(epochtoplot+1).nspksByNeuron(2)];
end

if strcmp(normtype, 'useminus')
frchange_neur1 = cellfun(@(x)mean(x), fr1_wn) - cellfun(@(x)mean(x), fr1_base);
frchange_neur2 = cellfun(@(x)mean(x), fr2_wn) - cellfun(@(x)mean(x), fr2_base);
elseif strcmp(normtype, 'usefrac')
frchange_neur1 = cellfun(@(x)mean(x), fr1_wn)./cellfun(@(x)mean(x), fr1_base);
frchange_neur2 = cellfun(@(x)mean(x), fr2_wn)./cellfun(@(x)mean(x), fr2_base);    
end
% ============== PLOT SCALAR
close all;
clim = [min([frchange_neur1; frchange_neur2]) ...
    max([frchange_neur1; frchange_neur2])];

onlygoodexpt = 0;
expttype = 'xcov_spikes';
plotlevel = 'chanpair';

useoldInds = 0; % default 0. only use 1 if you haven't extracted inds_base_epoch to OUTSTRUCT_XCOV yet...
% there will be error if you need to switch to 1.

assert(all(strcmp(OUTSTRUCT_XCOV.bregionpair, 'LMAN-RA')), 'not true that neurons 1 and 2 and LMAN and RA');
   
% ============= 1) CHANGE IN FR FOR NERUON 1 (is always LMAN)
lt_neural_POPLEARN_XCov_PlotScal(frchange_neur1, OUTSTRUCT_XCOV, SwitchCohStruct, SwitchStruct, ...
    MOTIFSTATS_Compiled, PARAMS, clim, onlygoodexpt, expttype, plotlevel, ...
    OUTSTRUCT, useoldInds);

% ============= 1) CHANGE IN FR FOR NERUON 2 (is always RA)
close all;
lt_neural_POPLEARN_XCov_PlotScal(frchange_neur2, OUTSTRUCT_XCOV, SwitchCohStruct, SwitchStruct, ...
    MOTIFSTATS_Compiled, PARAMS, clim, onlygoodexpt, expttype, plotlevel, ...
    OUTSTRUCT, useoldInds);


%% ========== [XCOV SCALAR] - MOTIF POSITION...
% NOTE: is currently plotting at level of channel pair. to change is easy -
% eitehr replace which indsgrp is used in the code by hand or add as
% variable.
close all;
windthis =1;
Yscalar = cellfun(@(x)(x(windthis,2)-x(windthis,1)), OUTSTRUCT_XCOV.Xcovscal_window_BaseWN);

onlygoodexpt = 1;
expttype = 'xcov_spikes';

lt_neural_POPLEARN_Xcov_MotifPosition(Yscalar, OUTSTRUCT_XCOV, SwitchCohStruct, ...
    SwitchStruct, MOTIFSTATS_Compiled, PARAMS, onlygoodexpt, expttype);



%% ========= ACROSS-EXPERIMENT CORRELATIONS




%% ========== [XCOV, XCOV PEAK OVER TIME]



%% ========= [XCOV AT BASELINE... - PLOT DISTRUBUTION OVER SYLS]

close all;
plotorigxcov = 0; % if 0, then recalcualed, otherwise 1 origianl. [default 0]
lt_neural_POPLEARN_Xcov_Sumcorr(OUTSTRUCT_XCOV, SummaryStruct, ...
    PARAMS, SwitchStruct, plotorigxcov);


%% ========= [PLOT SLIDING DOT PRODUCT USED TO COMPUTE XCOV]
close all;
% === option1:
birdtoplot = 'wh44wh39';
motiftoplot = '(j)n';
pairtoplot = [1 2];
% === option 2(overrides 1);
plotrandpair = 0;



recalcver = 0;
% 0 = zscore;
% 1 = first do minus, then do for all shuff, then get z rel shuff. [should
% be identical to 0]

% === i.e. shuffle. data, and zscore.
if plotrandpair==1
indthis = randi(length(OUTSTRUCT_XCOV.bnum));

bnum = OUTSTRUCT_XCOV.bnum(indthis);
motiftoplot = OUTSTRUCT_XCOV.motifname{indthis};
pairtoplot = OUTSTRUCT_XCOV.neurpair(indthis, :);
else

% ####################
bnum = find(strcmp({SummaryStruct.birds.birdname}, birdtoplot));
indthis = find(OUTSTRUCT_XCOV.bnum==bnum & ...
    strcmp(OUTSTRUCT_XCOV.motifname, motiftoplot) & ...
    ismember(OUTSTRUCT_XCOV.neurpair, pairtoplot, 'rows'));
assert(length(indthis)==1);
end

enum = OUTSTRUCT_XCOV.enum(indthis);

% #####################################
y = OUTSTRUCT_XCOV.Xcov_DotProd_trials{indthis}.base_dat;
yshuff = OUTSTRUCT_XCOV.Xcov_DotProd_trials{indthis}.base_shuff;
yxcov = OUTSTRUCT_XCOV.XcovBase_orig(indthis, :);
indsbase = OUTSTRUCT_XCOV.inds_base_epoch{indthis}

if recalcver==1
    % --- first convert both data and shuffle to xcov (i.e. subtraact the
    % man of shuffle).
    
        zmean = mean(yshuff,1);
%         zstd = std(yshuff, [], 1);

        % -- first convert all to xcov (i..e minus mean of shuff).
        yshuff = (yshuff - zmean);
        y = y - zmean;
    
end

x = PARAMS.Xcov_ccLags;

lt_figure; hold on;
hsplots = [];

% ==== 
hsplot = lt_subplot(3,2,1); hold on;
hsplots = [hsplots hsplot];
xlabel('lag');
ylabel('dot prod');
title('data trials');
plot(x, y', '-r');
lt_plot_zeroline;

hsplot = lt_subplot(3,2,2); hold on;
hsplots = [hsplots hsplot];
xlabel('lag');
ylabel('dot prod');
title('data trials');
plot(x, yshuff', '-k');
lt_plot_zeroline;

% === overlay
hsplot = lt_subplot(3,2,3); hold on;
hsplots = [hsplots hsplot];
title([num2str(bnum) '-e' num2str(enum) ','  motiftoplot '-neur:' num2str(pairtoplot)]);
plot(x, yshuff', '-', 'Color', [0.7 0.7 0.7]);
shadedErrorBar(x, mean(yshuff,1), std(yshuff), {'Color' , 'g'},1);
shadedErrorBar(x, mean(yshuff,1), lt_sem(yshuff), {'Color' , 'k'},1);
plot(x, mean(y,1), '-r');
lt_plot_zeroline;
xlim([-0.03 0.03]);

hsplot = lt_subplot(3,2,4); hold on;
% hsplots = [hsplots hsplot];
title('xcov');
ylabel('b=orig; m=recalc');
plot(x, yxcov, '-b');
lt_plot_zeroline;

% ======== recalculate xcov
% if recalcver==0
yxcov_recalc = (mean(y,1) - mean(yshuff,1))./std(yshuff, [], 1);
% elseif recalcver==1
%         
%         % --- then convert data to zscore rel shuff
%         yxcov_recalc = (mean(y2,1) - mean(yshuff2,1))./std(yshuff2,[], 1);
% end
plot(x, yxcov_recalc, 'm');

linkaxes(hsplots, 'xy');
xlim([-0.03 0.03]);

%% ========= [RAW, RASTERS] FOR A SINGLE CASE, PLOT RANDOM SUBSET OF TRIALS

close all;
% ======================== PLOT PAIRED RASTERS, ALIGNED TO SYL
BirdExptPairsToPlot = {'gr48bu5', 'RALMANLearn9'}; % ordered p[aors {bird, expt, bird, expt } ....
% BirdExptPairsToPlot = {'wh44wh39', 'RALMANlearn3'}; % ordered p[aors {bird, expt, bird, expt } ....
% BirdExptPairsToPlot = {'gr48bu5', 'RALMANLearn2'}; % ordered p[aors {bird, expt, bird, expt } ....
SwitchToPlot = [1]; % array of swiches.
neurpair_globID = [1 2]; % only one pair allowed: e.g. [9 17], actual global ID.
TypeOfPairToPlot = {'LMAN-RA'}; % e.g. 'LMAN-RA' (in alphabetical order)
motiftoplot = {'a(b)'}; % if empty, then plots target(s).
xgramxlim = [-0.085 -0.03];
xgramylim = [-0.045 0.045];
clim = [-0.25 0.25]; % for cov gram;
Nrand = 50; % number of trials from base and WN to extract
plotsimple = 1;
plotrawneural = 1; % will plot a few raw neural examples
lt_neural_POPLEARN_PairRast(BirdExptPairsToPlot, SwitchToPlot, ...
    neurpair_globID, motiftoplot, OUTSTRUCT, OUTSTRUCT_XCOV, SwitchStruct, ...
    SwitchCohStruct, MOTIFSTATS_Compiled, MOTIFSTATS_pop, SummaryStruct, ...
    PARAMS, xgramxlim, xgramylim, clim, Nrand, plotsimple, plotrawneural);


%% ===================== [RAW RASTERS] PLOT ALL CASES, SAVING
close all;
for i=1:length(OUTSTRUCT_XCOV.bnum)
    
    bnum = OUTSTRUCT_XCOV.bnum(i);
    enum = OUTSTRUCT_XCOV.enum(i);
    sw = OUTSTRUCT_XCOV.switch(i);
    neurpair = OUTSTRUCT_XCOV.neurpair(i,:);
    motif = OUTSTRUCT_XCOV.motifname{i};
    
    bname = SwitchStruct.bird(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    
    
    BirdExptPairsToPlot = {bname, ename}; % ordered p[aors {bird, expt, bird, expt } ....
    SwitchToPlot = [sw]; % array of swiches.
    neurpair_globID = neurpair; % only one pair allowed: e.g. [9 17], actual global ID.
    TypeOfPairToPlot = {'LMAN-RA'}; % e.g. 'LMAN-RA' (in alphabetical order)
    motiftoplot = {motif}; % if empty, then plots target(s).
    xgramxlim = [-0.085 -0.03];
    xgramylim = [-0.045 0.045];
    clim = [-0.25 0.25]; % for cov gram;
    Nrand = 40; % number of trials from base and WN to extract
    plotsimple = 1;
    plotrawneural = 1; % will plot a few raw neural examples
    lt_neural_POPLEARN_PairRast(BirdExptPairsToPlot, SwitchToPlot, ...
        neurpair_globID, motiftoplot, OUTSTRUCT, OUTSTRUCT_XCOV, SwitchStruct, ...
        SwitchCohStruct, MOTIFSTATS_Compiled, MOTIFSTATS_pop, SummaryStruct, ...
        PARAMS, xgramxlim, xgramylim, clim, Nrand, plotsimple, plotrawneural);
    
    dirthis = [bname '-' ename '-sw' num2str(sw) '-' motif '-n' num2str(neurpair(1)) '_' num2str(neurpair(2))];
    cd('/bluejay0/bluejay2/lucas/analyses/neural/LEARN/LMAN-RA');
    mkdir(dirthis);
    cd(dirthis);
%     savefig(gcf, [bname '-' ename '-sw' num2str(sw) '-' motif '-n' num2str(neurpair(1)) '_' num2str(neurpair(2)) '.fig'], 'compact');
    lt_save_all_figs;
    close all;
    cd ..
end


%% ======================= [RAW RASTERS]
% ITERATE AND PLOT ALL CASES
cd('/bluejay0/bluejay2/lucas/analyses/neural/LEARN/LMAN-RA');
% cd('/bluejay0/bluejay2/lucas/analyses/neural/LEARN/LMAN-RA_beforeSaveRawNeural/LMAN-RA');
close all;

birdtoplot = 'gr48bu5';
expttoplot = 'RALMANLearn5';
swtoplot = 1;
motifplot = 'j(b)';

for i=1:length(SwitchStruct.bird)
    bname = SwitchStruct.bird(i).birdname;
    if ~isempty(birdtoplot)
        if ~strcmp(birdtoplot, bname)
            continue
        end
    end
    
    for ii=1:length(SwitchStruct.bird(i).exptnum)
        ename = SwitchStruct.bird(i).exptnum(ii).exptname;
        if ~isempty(expttoplot)
            if ~strcmp(expttoplot, ename)
                continue
            end
        end
        
        for ss =1:length(SwitchStruct.bird(i).exptnum(ii).switchlist)
            if ~isempty(swtoplot)
                if swtoplot ~= ss
                    continue
                end
            end
            
            % ====== get list of motifs
            indsthis = OUTSTRUCT_XCOV.bnum==i & OUTSTRUCT_XCOV.enum==ii & ...
                OUTSTRUCT_XCOV.switch==ss;
            
            if ~any(indsthis)
                continue
            end
            
            motiflist = unique(OUTSTRUCT_XCOV.motifname(indsthis));
            
            for mm=1:length(motiflist)
                motifthis = motiflist{mm};
                
                if ~isempty(motifplot)
                    if ~strcmp(motifplot, motifthis)
                        continue
                    end
                end
                
                if (0) % version 1 - load saved files
                    fname_prefix = [bname '-' ename '-sw' num2str(ss) '-' motifthis '*.fig'];
                    
                    % -- find all files with this prefix
                    flist = dir(fname_prefix);
                    for j=1:length(flist)
                        openfig(flist(j).name);
                    end
                    pause;
                    close all;
                else % version 1 - got to directory
                    fname_prefix = [bname '-' ename '-sw' num2str(ss) '-' motifthis '*'];
                    tmp = dir(fname_prefix);
                    % --- go to all folders and load
                    for j=1:length(tmp)
                       cd(tmp(j).name);
                       openfig('figsall.fig');
                       cd ..
                    end
                    pause; 
                    close all;
                end
            end
        end
    end
end


%% ##################################### [MORE CONSISTENT SYLLABLE LOCKED ACTIVITY?]
%% ============= [CONSISTENT ACTIVITY] EXTRACTION
% onlykeepgoodcase = 1;
% onlyadjacentpairs = 1; % if 0, then all trial pairs...
% corrwind = [-0.1 0.01]; % rel syl onset, to get pairwise trial by trial correlations.
% corrwind = [-0.12 0]; % rel syl onset, to get pairwise trial by trial correlations.
% baseallinds = 1; % 0 is epocjh, 1 is all % default:1
% wnallinds = 0;
% if strcmp(wntouse, 'quarter')
% usehalfwntmp = 1; % entered above...
% end
% % usehalfwntmp = 0; % entered above...
onlykeepgoodcase = 1;
onlyadjacentpairs = 1; % if 0, then all trial pairs...
corrwind = [-0.1 0.0]; % rel syl onset, to get pairwise trial by trial correlations.
baseallinds = 0; % 0 is epocjh, 1 is all % default:1
wnallinds = 0;
if strcmp(wntouse, 'quarter')
usehalfwntmp = 1; % entered above...
end
% usehalfwntmp = 0; % entered above...
usesameindsasXcov = 1; % this trumps everything else --> will use same inds as xcov analyses.
OUTSTRUCT_units = lt_neural_POPLEARN_Cons_Extr(OUTSTRUCT, OUTSTRUCT_XCOV, SwitchStruct, ...
    SwitchCohStruct, MOTIFSTATS_Compiled, MOTIFSTATS_pop, SummaryStruct, ...
    PARAMS, onlykeepgoodcase, onlyadjacentpairs, corrwind, baseallinds, ...
    wnallinds, usehalfwntmp, usesameindsasXcov);


%% ============ [CONSISTENT ACTIVITY - LFP ONLY]

onlykeepgoodcase = 1;
onlyadjacentpairs = 1; % if 0, then all trial pairs...
corrwind = [-0.1 0.01]; % rel syl onset, to get pairwise trial by trial correlations.
corrwind = [-0.1 0]; % rel syl onset, to get pairwise trial by trial correlations.
baseallinds = 0; % 0 is epocjh, 1 is all % default:1
wnallinds = 0;
if strcmp(wntouse, 'quarter')
usehalfwntmp = 1; % entered above...
end
usehalfwntmp = 0; % entered above...
lfpfilt = [15 35];
OUTSTRUCT_lfp = lt_neural_POPLEARN_Cons_Extr_LFP(OUTSTRUCT, OUTSTRUCT_XCOV, SwitchStruct, ...
    SwitchCohStruct, MOTIFSTATS_Compiled, MOTIFSTATS_pop, SummaryStruct, ...
    PARAMS, onlykeepgoodcase, onlyadjacentpairs, corrwind, baseallinds, ...
    wnallinds, usehalfwntmp, lfpfilt);

OUTSTRUCT_lfp.neurID =  OUTSTRUCT_lfp.chanthis;

%% ============= [PLOTS] OVERVIEW
% === FOR EACH EXPERIMENT, PLOT ALL DATA
close all;
birdtoplot = [];
onlygoodexpt = 1;
lt_neural_POPLEARN_Cons_PlotOver(OUTSTRUCT_units, SwitchStruct, SummaryStruct, ...
    birdtoplot, onlygoodexpt);


%% ============ [PLOTS] SUMMARY

% === FOR EACH EXPERIMENT, PLOT ALL DATA
close all;
birdtoplot = [];
expttoplot = [];
onlygoodexpt = 1;
plotlevel = 'switch';
% plotlevel = 'chanpair';
lt_neural_POPLEARN_Cons_PlotOverSum(OUTSTRUCT_units, OUTSTRUCT_XCOV, ...
    SwitchStruct, SummaryStruct, birdtoplot, onlygoodexpt, plotlevel, ...
    corrwind, expttoplot);

%% ============ [PLOTS-LFP] SUMMARY
close all;
birdtoplot = [1 2 4];
expttoplot = [];
onlygoodexpt = 1;
% plotlevel = 'switch';
plotlevel = 'chanpair';
lt_neural_POPLEARN_Cons_PlotOverSum_LFP(OUTSTRUCT_lfp, OUTSTRUCT, ...
    SwitchStruct, SummaryStruct, birdtoplot, onlygoodexpt, plotlevel, ...
    corrwind, expttoplot);


%% ##############################################################
%% ################## [COMPARE SYL-LOCKED ACTIVITY] 
% --- clear things that take up memory but not needed for following
% analyses.
clear SwitchXCovStruct
clear LFPSTRUCT

%% ============ [PLOTS] - trial-mean LFP and MU, and their relations
% ============== WILL ALWAYS USE TRIALS IN OUTSTRUCT_XCOV (THOSE ARE GOOD)

% =============== 1) FOR EACH TARGET SYL PLOT ALL LFP AND NEURONS (IN RA
% AND LMAN)
close all;
onlygoodexpt = 1;
xtoplot = [-0.12 0.02];
xtoplot = [-0.1 0.01];
% xtoplot = [-0.15 0.01];
plotraw = 0; % 1: plots each syl/expt, chan etc...
disp('NOTE: remember to include wh72 in this analysis');
plotAlsoWN = 1; % matters even if plotraw==0. assumes that normal plot is for baseline.
fpass = [12 150]; % for bandpass filtering LFP.
% fpass = [20 35]; % for bandpass filtering LFP.
lt_neural_POPLEARN_SylLocked_Over(OUTSTRUCT, OUTSTRUCT_XCOV, SwitchStruct, ...
    SwitchCohStruct, MOTIFSTATS_Compiled, MOTIFSTATS_pop, SummaryStruct, ...
    PARAMS, onlygoodexpt, xtoplot, plotraw, plotAlsoWN, fpass);


%% ============ [PLOT] individual trials, LFP, SPIKING

close all;
onlygoodexpt = 1;
% xtoplot = [-0.12 0.02];
% xtoplot = [-0.14 0.01];
xtoplot = [-0.1 0.0];
xtoplot = [-0.1 0.0];
plotraw = 0; % 1: plots each trial (examples...)
disp('NOTE: remember to include wh72 in this analysis');
zscore_lfp = 0; % if 1, then z-scores concatenating all trials (a given chan) within time segment (after filtering)
% fpass = [40 100]; % for bandpass filtering LFP.
% fpass = [20 35]; % for bandpass filtering LFP.
% fpass = [5 350]; % for bandpass filtering LFP.
fpass = [25 60]; % for bandpass filtering LFP.
% sta_wind = [-0.05 0.05]; % relative to spike, in sec % will only get spikes that are within data...
sta_wind = [-0.03 0.03]; % relative to spike, in sec % will only get spikes that are within data...
% kernelSD = 0.015; % empyt for default (for spike smoothgin)
kernelSD = 0.005; % empyt for default (for spike smoothgin)
removeBadLFP = 0; % then things that look like large fluctuations..
% if only care about smooth FR, then set to 0. if care about LFP and want
% to remove things in a conservative way, then set to 1.
[DATSTRUCT_SPK, DATSTRUCT_LFP, PARAMS] = lt_neural_POPLEARN_SylLocked_PlotTrials(OUTSTRUCT, OUTSTRUCT_XCOV, SwitchStruct, ...
    SwitchCohStruct, MOTIFSTATS_pop, SummaryStruct, ...
    PARAMS, onlygoodexpt, xtoplot, plotraw, zscore_lfp, fpass, sta_wind, ...
    kernelSD, removeBadLFP);
% ==== older params
lt_neural_POPLEARN_SylLocked_params;


% ============ [ANALYSIS] NEURAL CORRELATION PREDICTS STRONGER LFP POPWER?
close all;
plotRawRand = 0;
lt_neural_POPLEARN_SylLocked_Coord(DATSTRUCT_SPK, DATSTRUCT_LFP, PARAMS, ...
    plotRawRand)


% =========== [ANALYSIS] TRIAL-BY-TRIAL CROSS CORRELATIONS [LFP AND
% SPIKING]
close all;
plotRawRand = 0;
lt_neural_POPLEARN_SylLocked_Timing(DATSTRUCT_SPK, DATSTRUCT_LFP, PARAMS, ...
    plotRawRand)


% =========== [LFP AND SPIKING SPECTRA - AND COHERENCE]
close all;
xtoplot;
lt_neural_POPLEARN_SylLocked_Spectra(DATSTRUCT_SPK, DATSTRUCT_LFP, PARAMS, ...
    xtoplot, SwitchCohStruct, SwitchStruct);




%% ============ [ANALYSIS] FR CORR CHANGE DURING LEARNIGN?
close all;
DATSTRUCT_POP = lt_neural_POPLEARN_SylLocked_CoorLearn(DATSTRUCT_SPK, DATSTRUCT_LFP, ...
    PARAMS, OUTSTRUCT_XCOV, SwitchStruct);


%% ============ [ANALYSIS] LMAN POP CORRELATION 
% 1) RELATED TO LMAN-RA XCOV? FROM TRIAL TO TRIAL
% 2) RELATED TO BASELINE AFP BIAS?

dirstruct = lt_DirSong_MotifID;
dirstruct = lt_DirSong_ExtractAFPbias(dirstruct);

close all;
bregionthis = 'LMAN';
plotbase = 0; % then usese LMANcoord-pitch correaltions from baseline. if 0 then
% uses from WN. either case will use AFP bias onlyu from baseline.
DATSTRUCT_POP = lt_neural_POPLEARN_SylLocked_LMANcoord(DATSTRUCT_POP, DATSTRUCT_SPK, DATSTRUCT_LFP, ...
    PARAMS, OUTSTRUCT_XCOV, SwitchStruct, dirstruct, OUTSTRUCT, ...
    bregionthis, plotbase);



%% ======= [EXTRACTION] REPLIT LMAN-RA COORDINATION, BASED ON WITHIN-REGION COORD
% NOTE: This REXTRACTS XCOV AS BEFORE, BUT INSTEAD OF USING FF TO SPLIT, IT
% USES [LMAN OR RA] COORDINATINO TO SPLIT
% Should be able to use the output here in any function that works for
% OUTSTRUCT_XCOV.
% - Should do separately for LMAN and RA. 
% - Code follows, for using with this structure.

bregionToSplitBy = 'RA'; % i.e. if LMAN then splits by LMAN ensemble coordination

% ============== smoothe in lag space? 
dosmooth = 1; % interpolate, then smooth
dosmooth_sigma = 2*PARAMS.Xcov.binsize_spk; % sec

% ===== what version of xcov?
xcovver = 'zscore'; % i.e. each lag bin, zscore relative to shuffle.
% xcovver = 'scale'; % % NOT DONE YET! scales by relative magnitude of shuffle functions.
% xcovver = ''; % empty means difference of means.
% xcovver= 'coherency';

% ===== get xcovgram?
getxgram = 1;
getxgram_epochbins = 3; % if empty, then ignore. otherwise divides training into this many even bins

% ====== removebad syl?
removebadsyl = 1;
removebadtrials =1;
removebadchans = 1;
plotraw = 0; % will keyboard on targets, and can run to see how compute xcov [i.e. different methods]

% ======= AUTO PARAMS
if strcmp(xcovver, 'coherency')
    dosmooth = 0; % have not coded this yet, should take smoothing AFTER get coherency.
end
    
% ===== split [base WN] and [epochs] into hi and low ff trials.
getHiLoFFSplit = 1; 

hilosplit_shuffver=1;


[OUTSTRUCT_MeanRhoSplit, PARAMS_MeanRhoSplit] = lt_neural_POPLEARN_XcovSplitByLMANcoord(...
    SwitchXCovStruct, SwitchStruct, PARAMS, SwitchCohStruct, dosmooth, dosmooth_sigma, ...
    xcovver, getxgram, removebadchans, getHiLoFFSplit, hilosplit_shuffver, OUTSTRUCT_XCOV, ...
    DATSTRUCT_POP, bregionToSplitBy);


%% ====== [TO SAVE]
savedir = ['/bluejay0/bluejay2/lucas/analyses/neural/POPLEARN/OUTSTRUCT_XCOV/' PARAMS.savemarker];

savemarker = '60mswind_031919_3bins_hiloFF';  % latest, and 6 epochs for Xcovgram
save([savedir '/OUTSTRUCT_MeanRhoSplit_' bregionToSplitBy '_' savemarker '.mat'], 'OUTSTRUCT_MeanRhoSplit', '-v7.3');
save([savedir '/PARAMS_MeanRhoSplit_' bregionToSplitBy '_' savemarker '.mat'], 'PARAMS_MeanRhoSplit', '-v7.3')

%% ======= [TO LOAD]

savedir = ['/bluejay0/bluejay2/lucas/analyses/neural/POPLEARN/OUTSTRUCT_XCOV/' PARAMS.savemarker];

savemarker = 'LMAN_60mswind_031919_3bins_hiloFF';  % latest, and 6 epochs for Xcovgram
savemarker = 'RA_60mswind_031919_3bins_hiloFF';  % latest, and 6 epochs for Xcovgram
load([savedir '/OUTSTRUCT_MeanRhoSplit_'  savemarker '.mat']);
load([savedir '/PARAMS_MeanRhoSplit_' savemarker '.mat']);



%% ====== [POST EXTRACTINON]
OUTSTRUCT_MeanRhoSplit = lt_neural_POPLEARN_XCovExtr_post(OUTSTRUCT_MeanRhoSplit, ...
    SwitchCohStruct, SummaryStruct, SwitchStruct);

%% ========= 2) EXTRACT SCALAR FOR TWO PEAKS FOR EACH EXPERIMENT

% ========= GOOD: [60ms window]
alignto = 'sylonset'; 
twindows = {[-0.07 -0.03]}; % one array for each window [x centers...] [inclusive]
lagwindows = {[-0.01 0.005]}; % OK


OUTSTRUCT_MeanRhoSplit = lt_neural_POPLEARN_XCov_ExtrScal(OUTSTRUCT_MeanRhoSplit, OUTSTRUCT, ...
    PARAMS_MeanRhoSplit, twindows, lagwindows, alignto);


%% ========= 2) EXTRACT TIME SLICE (OVERWRITE ORIGINAL TIME SLICE...)
% ======= GOOD [60ms window]
twindow = [-0.07 -0.03]; % one array for each window [x centers...] [inclusive]
alignto = 'sylonset'; 

OUTSTRUCT_MeanRhoSplit = lt_neural_POPLEARN_XCov_ExtrSlice(OUTSTRUCT_MeanRhoSplit, ...
    OUTSTRUCT, PARAMS_MeanRhoSplit, twindow, alignto);

%% ====== [PLOT] XCOV SLICES (ADAPTIVE VS. NONADAPTIVE)

close all;
epochstoplot = 3;
ploteachsyl = 1; % if 0, then only plots target syl.
lt_neural_POPLEARN_XCov_FFsplitSlice(OUTSTRUCT_MeanRhoSplit, ...
    PARAMS_MeanRhoSplit, epochstoplot, ploteachsyl)


%% ================= [NEURAL BIAS VS DIRECTED SONG BIAS]
% ======== DOES DIRECTION OF BASELINE BIAS CORRELATE WITH BASELINE AFP
% BIAS?
close all;
dirstruct = lt_DirSong_MotifID;
useff_zscore = 1; % if 0, then mean FF diff, if 1, then dir zscore from undir. (each day, then take mean)
%  if 2, then uses percent change.
targsylonly = 0; % default 0.
cvdiffmethod = 'diff'; % cv undir minus cv dir
cvdiffmethod = 'percent'; % cv undir minus cv dir
onlygoodexpt = 1;
giveFakeAFPbias = 1; % default 0. if 1, then makes bias exactly 0 if no detectable FF.
% allows plotting on same plot.

lt_neural_POPLEARN_Xcov_DirBias(OUTSTRUCT_MeanRhoSplit, PARAMS_MeanRhoSplit, ...
    onlygoodexpt, SwitchStruct, dirstruct, useff_zscore, targsylonly,  ...
    cvdiffmethod, giveFakeAFPbias);





%% ###################################################################
%% =========== [PLOT ALL CROSS CORRELATIONS - USING BINNED SPIKES]
% datapoint = combination of syllablwe and neuron pair.
close all;
onlygoodexpt = 1;
lt_neural_POPLEARN_XCov_BaseCorr(OUTSTRUCT_XCOV, PARAMS, onlygoodexpt, SwitchStruct)

%% ============= RAW cross correlation (actually sliding dot product)

y = mean(OUTSTRUCT_XCOV.XcovBase_NoMinShuff,1);
ysem = lt_sem(OUTSTRUCT_XCOV.XcovBase_NoMinShuff);
x = PARAMS.Xcov_ccLags;
lt_figure;  hold on;
xlabel('LMAN -- RA');
ylabel('sliding dot product');
shadedErrorBar(x, y, ysem, {'Color', 'k'},1);
xlim([-0.04 0.04]);
lt_plot_zeroline_vert;



%% ######################################################################
%% ################################ [SIMULATION]
%% ========== 1) trial to trial mean fluctuation - could that account for xcov changes?

lt_neural_POPLEARN_MODEL_meanfluc;

