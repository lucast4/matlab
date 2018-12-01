%% ++++++++++++++++++++++++++++++++++++INTRO
% Method 1 and Method 2 below



%% +++++++++++++++++++++++++++++++++ METHOD 1 - autolabels using evtafv4 some - uses tameplate matching FF for context analysis
% === useful when have ginormous data, e.g. many days with rap[id context
% switching.

%% LT 4/26/15 - general script with directions for how to perform context analyses
% This requires using EvTAF_v4_LT_v1 or higher - i.e. has note group
% functionality.  It allows automatic switching of note groups (a note
% group is a set of notes that are imposed during a given epoch).


%% Single day data extraction
% extracts data using evtafsim, also extracts notegroup information. puts
% all into one structure and saves into a common dir.
lt_make_batch(1);

% -- 1) what batch and config file to use to run evtafsim.  if needed, make a batch.
Params.batch='batch.keep';
% Params.config= '/bluejay4/lucas/birds/rd23gr89/config_042115_fixedfreqbounds.evconfig2'; % for CtxtDepPitch (before 4/26)
% Params.config= '/bluejay4/lucas/birds/rd23gr89/042715_CtxtDepPitch2_TargB2_day2/config.evconfig2'; % for CtxtDepPitch2 (4/26 +)
Params.config= '/bluejay4/lucas/birds/rd23gr89/042815_CtxtDepPitch2_TargB2__day3/config.evconfig2'; % for CtxtDepPitch2 (4/28 +) (better hits)

Params.expt_phrase='Stimulation'; % folder where will save. leave blankt o get automatically from dir name
Params.dataID='All'; % e.g. id of data (e.g. 'All' for all data in data). if blank, uses batch name.

% -- Run evtafsim and save information
AllData=lt_context_ExtractDayData(Params);



%% DATA EXTRACTION - ACROSS DAYS
clear all; close all;
% metadata params
Params.metadata.experiment = 'LightDepPitch';
Params.metadata.condition='';
Params.metadata.notes='';
Params.metadata.date_range={'05May2015','05May2015'};
Params.metadata.only_labeled_dirs=0;

% analysis params
Params.analysis.batch='batch.keep';
% last version config (template used at very end of SeqDepPitch)
Params.analysis.config = '/bluejay4/lucas/birds/pu35wh17/061215_SeqDepPitch_durWN_day9/config_061115.evconfig2';
Params.analysis.dataID='All'; % e.g. id of data (e.g. 'All' for all data in data). if blank, uses batch name.
Params.analysis.Make_BatchKeep=1; % will make a new batch.keep for each day


Dirs_Missing_notegroups=lt_context_ExtractData_AcrossDays(Params);

%% BELOW, STUFF TO COLLECT DATA AND PLOT - NEED TO CONSOLIDATE INTO VARIOUS FUNCTIONS





%% COLLECT AllData Structures from multiple days
clear all; close all;

% Params
Params_alldays.CollectAllData=1; % if 1, then collects all data starting from first day
Params_alldays.firstday='05May2015';
Params_alldays.lastday='14May2015';

Params_alldays.NoteToPlot=2; % this is the note whose detects we will analyze (i.e. this note should get all renditions of the syl)
Params_alldays.RunBin=10;

Params_alldays.BoundaryTimes={'05May2014-1423', '08May2014-1423'}; % in format of e.g. 05May2014-1423, these are times of switching in experiment (e.g. turning WN off and on, changing pitch contingency, etc)

Params_alldays.Edge_Num_Rends = 20; % num rends to call "edges" (defualt: queries)
Params_alldays.Probe_CSplus=[1 2]; % [from to] (actual NG nums) (e.g. from no light --> light on(probe))
Params_alldays.Probe_CSminus=[1 3]; % [from to] (actual NG nums) (e.g. from no light --> no light (probe))

Params_alldays.PhaseToCompare1=4; % e.g. [light + WN up] phase
Params_alldays.PhaseToCompare2=5; % e.g. [light + WN dn] phase 

Params_alldays.throw_out_if_epoch_diff_days=1; % throws out any transitions that overlap with O/N (potentially 2 per O/N)

lt_context_CompileAndPlot(Params_alldays);




%% ++++++++++++++++++++++++++++++++ METHOD 2 - hand labeled - exctracts data, uses pitch contour

%% ======================= Extracting data across days
% --- day directories must be in format
% [date]_[experimentname]_[context_name].
% -- Run this in bird folder to extract all days data to subfolder

clear all; close all;

Params.syl_list={'b'}; % single syls
Params.freq_range_list={[4150 5100]};
Params.pc_dur_list=[0.11];
Params.pc_harms_list=[1];

Params.batch='batch.labeled.all';
Params.experiment = 'Association2';

date_range={'02Apr2017','02Apr2017'}; % e.g. {'20Apr2015','20May2015'}. leave blank ('') for all days

% ---- 2) Collect note group information? 
CollectNoteGroup = 1; % set to 1 if want to use online NoteGroups. Otherwise will do context stuff using the 
% [context_name] condition in the folder name - i.e. for bulk context, not
% rapid switching experiments. 

lt_extract_AllDaysPC(Params, date_range, CollectNoteGroup)


%% %% Compiling data - go to "extract" folder first
clear all; close all;
Params_global.CompilePC.PC_window_list={'b', [22 74]}; % syl, value pairs [single syls]
Params_global.CompilePC.FirstDay='';
Params_global.CompilePC.LastDay='';
plotON=1; % pitch contours, all days, all syls
saveON=1;

% Regular expressions - first calculates FF etc, then performs regular
% expressions
Params_global.CompilePC.regexp_list={'c(b)'}; % e.g. {'dcc(b)', 'ab+(g)'} : dcc(b) means match dccb, and give me ind of b in parantheses.  ab+(g) means match ab... (anly length repeat), then g. give me ind of g
Params_global.CompilePC.regexp_fieldnames={'cB'}; % whatever
% want to call structure field (if this cell array not defined, then will
% attempt to use the regexp names.
    
[ALLDATSTRUCT, Params_global]= lt_extract_CompilePC(plotON, Params_global, saveON);


%% Convert to context1 format
% --- TO BE ABLE TO RUN USING CONTEXT PLOT SAME AS FOR METHOD 1
% Params_global.ConvertToContext1.NoteNum_codes={'dcc_b_', 1, 'bcc_b_', 2}; % {notestring, notenum} pairs - notestring either single syl (e.g. 'a') or regexp, using underscores (e.g. 'dcc_b_')
Params_global.ConvertToContext1.NoteNum_codes={'cB', 1}; % {notestring, notenum} pairs - notestring either single syl (e.g. 'a') or regexp, using underscores (e.g. 'dcc_b_')

% syl='b';

Params_global.ConvertToContext1.UseAutoNG = 1; % if 1, this uses .ltrec2 saved note group information (this was extracted in lt_extract_AllDaysPC)
% if 0, then uses the names of directories - assumes that each directory
% (for a given day) was a different context
% ([date]_[experimentname]_[context_name]); if 0, then need to fill this
% out. If 1, then don't need this - can leave the field as whatever, will
% delete the entry.
Params_global.ConvertToContext1.NoteGroupNum_codes={'away', 1, 'home', 2}; % {NoteGroup_name, NoteGroupNum} pairs - name must match what is in "condition" field. NoteGroupNum can be anything (keep it from 1, 2, ...);



[AllSongsDataMatrix, Params_alldays]= lt_context2_ConvertToContext1(ALLDATSTRUCT, Params_global);



%% PLOT
% USING SAME CONTEXT PLOT CODE FROM METHOD 1
close all;
Params_alldays.NoteToPlot=1;
Params_alldays.RunBin=10;

Params_alldays.BoundaryTimes={'15Nov2015-0000'}; % in format of e.g. 05May2014-1423, these are times of switching in experiment (e.g. turning WN off and on, changing pitch contingency, etc)
Params_alldays.Edge_Num_Rends = 40; % num rends to call "edges"

Params_alldays.throw_out_if_epoch_diff_days=0; % throws out any transitions that overlap with O/N (potentially 2 per O/N)
one_switch_a_day=1; % manual switching experiemnts.

plotIndTrials = 1;

lt_context_PLOT(AllSongsDataMatrix, Params_alldays, one_switch_a_day);


%% ########################## FOR OPTO ASSOCIATION 
% ====== SEE or60_analysis_Association1 for details.
% ==== assumes is interleaved laser paired with WN, and switching over
% phases.

%% === 1) EXTRACT SONG DATA
justDoExtraction = 1;
[SORTED_DATA, PARAMS_GLOB] = lt_context_PLOT(AllSongsDataMatrix, Params_alldays, one_switch_a_day, plotIndTrials, ...
    justDoExtraction);


%% === 2) GET SUMMAR ACROSS PHASES, LASER EFFECT.
% datatype = 'phase'; % one value per phase
datatype = 'day'; % one value per day

% ======================================
% === for each phase and each notegroup get one mean FF across trials
numNG = length(SORTED_DATA.ByNoteGroup);
numPhase = length(PARAMS_GLOB.Phases_DayBounds);

% AllFFmeans_NGxPhase = nan(numNG, numPhase);
AllFFmeans_NGxPhase = cell(numNG, numPhase);
for ng = 1:numNG
    for pp=1:numPhase
        
        % --- get "epochs" within this phase
        indsepoch = SORTED_DATA.ByNoteGroup(ng).Stats_OneDataPtPerEpoch.PhaseNum_array == pp;
        
        % --- get all trials within this pohase
        ffvals = cell2mat(cellfun(@transpose, SORTED_DATA.ByNoteGroup(ng).Stats_OneDataPtPerEpoch.RAW.FFvals(indsepoch), 'UniformOutput', 0));
        tvals = cell2mat(cellfun(@transpose, SORTED_DATA.ByNoteGroup(ng).Stats_OneDataPtPerEpoch.RAW.Tvals(indsepoch), 'UniformOutput', 0));
        %         ffvals = SORTED_DATA.ByNoteGroup(ng).Stats_OneDataPtPerEpoch.meanFF(indsepoch);
        
        %         AllFFmeans_NGxPhase(ng, pp) = mean(ffvals);
        if strcmp(datatype, 'phase')
            AllFFmeans_NGxPhase{ng, pp} = mean(ffvals);
        elseif strcmp(datatype, 'day')
            ffmeans = grpstats(ffvals, floor(tvals), {'mean'}); % get day means
            AllFFmeans_NGxPhase{ng, pp} = ffmeans;
        end
    end
end
    
% ====== for each phase, specify direction of learning
% disp('this is learning in notegroup 0. notegroup 0 must always be laser');
% disp('i.e. assumes that laser is in notegroup 0 and that is paired with some dir of laerning that depends on phase of expt');
disp(' THIS ASSUMES THAT THE NG ASSOCAITED WITH LASER DOES NOT CHANGE DURING THE EXPT...');
pause;

PARAMS_GLOB.Assoc_PhaseLearndirMapping = nan(numPhase,1);
PARAMS_GLOB.Assoc_PhaseLearndirMapping(3) = 1;
PARAMS_GLOB.Assoc_PhaseLearndirMapping(4) = -1;
PARAMS_GLOB.Assoc_PhaseLearndirMapping(6) = -1;
PARAMS_GLOB.Assoc_PhaseLearndirMapping(7) = -1;
PARAMS_GLOB.Assoc_PhaseLearndirMapping(8) = 1;
PARAMS_GLOB.Assoc_PhaseLearndirMapping(9) = 1;
PARAMS_GLOB.Assoc_PhaseLearndirMapping(11) = 1;
PARAMS_GLOB.Assoc_LaserNGnum = 0;
PARAMS_GLOB.Assoc_LaserOffNGnum = 1;

% ============= go thru all phases, collect data
ng_las = PARAMS_GLOB.Assoc_LaserNGnum+1;
ng_nolas = PARAMS_GLOB.Assoc_LaserOffNGnum+1;

phasestoget = find(~isnan(PARAMS_GLOB.Assoc_PhaseLearndirMapping));

FFmeans_Las_NoLas = [];
LaserLearnDir = [];
Phasenum = [];
for pp=phasestoget'
   
    ffthis = AllFFmeans_NGxPhase([ng_las ng_nolas], pp);
    trdir = PARAMS_GLOB.Assoc_PhaseLearndirMapping(pp);
    
    % ======= save output
    FFmeans_Las_NoLas = [FFmeans_Las_NoLas; ffthis'];
    LaserLearnDir = [LaserLearnDir; trdir];
    Phasenum = [Phasenum; pp];

end

% =============== BREAK OUT INDIVIDUAL DAYS (IF MULTIPE IN PHASES)
FFmeans_Las_NoLas_TMP = [];
LaserLearnDir_TMP = [];
Phasenum_TMP = [];
nrows = size(FFmeans_Las_NoLas,1);
for i=1:nrows

    ffmeans = [FFmeans_Las_NoLas{i,1} FFmeans_Las_NoLas{i,2}];
    ldir = LaserLearnDir(i);
    ph = Phasenum(i);
    
    % ===== stick into new output
FFmeans_Las_NoLas_TMP = [FFmeans_Las_NoLas_TMP; ffmeans];
LaserLearnDir_TMP = [LaserLearnDir_TMP; ones(size(ffmeans,1),1)*ldir];
Phasenum_TMP = [Phasenum_TMP; ones(size(ffmeans,1),1)*ph];
end
FFmeans_Las_NoLas = FFmeans_Las_NoLas_TMP;
LaserLearnDir = LaserLearnDir_TMP;
Phasenum = Phasenum_TMP;


% ############### [PLOT]
figcount=1;
subplotrows=4;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% ====== [LASER + WN UP]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('laser -- nolaser');
ylabel('mean ff (each phase)');
title('laser + WN up');

ff = FFmeans_Las_NoLas(LaserLearnDir==1, :);
x = [1 2];
plot(x, ff', '-ok');
xlim([0 3]);

% ====== [LASER + WN DN]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('laser -- nolaser');
ylabel('mean ff (each phase)');
title('laser + WN dn');

ff = FFmeans_Las_NoLas(LaserLearnDir==-1, :);
x = [1 2];
plot(x, ff', '-ok');
xlim([0 3]);


% ======== [COMBINED - PLOT DIFFERENCE (LASER MINUS NO LASER)]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('laser+DN ------ Laser+UP');
ylabel('diff in mean ff (laser - nolaser)');
title('phase = dtpt');

ffdiff = FFmeans_Las_NoLas(:,1) - FFmeans_Las_NoLas(:,2);
x = LaserLearnDir;
plot(x, ffdiff, 'ok');
% -- means
[ymean, ysem] = grpstats(ffdiff, x, {'mean', 'sem'});
lt_plot(unique(x)+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
[h, p] = ttest2(ffdiff(x==1), ffdiff(x==-1));
lt_plot_pvalue(p, 'ttest2', 1);
xlim([-2 2]);
lt_plot_zeroline;
axis tight


