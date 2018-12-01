
function lt_Opto_Stim_analy_SUMMARY_MultDayAnaly_v3(Params_metadata, Params_glob, KeepOutliers, SepBySyl)
%% LT 9/28/18 - does separate analyses for different desired syls.
% Running this function does analysis for one syl. is modified from
% original version do make this analysis in subdir for specific syl...

% SepBySyl = 1; % then saves this in its own subdir (with name given by syl)
if ~exist('SepBySyl', 'var')
    SepBySyl = 0;
end

if SepBySyl==1
    Params_glob.Delete_old_analysis_folder=0;
    % can't be 1, since it if SepBySyl=1 then this will look for and delete
    % delete higher level folder containing all the individual syls. 
end

if isempty(Params_glob.StimDur)
    Params_glob.StimDur=10;
end
%% LT version 3 - 7/9/15 - MODIFIED TO AUTOMATICALLY FIND DIR NAMES (uses lt_metadata analysis)
% IMPORTANT: dirs must be structured as (for e.g.)
% /bluejay4/lucas/birds/pk32/070115_Reversion1_preWN_STIMoff
% (date_experiment_somethingelse_STIMoff[or STIMon]).

% ALSO: Stim epochs must be catch trials (because will use StimCatch and
% StimNotCatch fields to filter data)

% WILL TAKE ALL LABELED SONGS

% RUN THIS IN BIRD FOLDER


% %% PARAMS ENTRY
% clear all; close all;
% 
% % =======, to find directories
% experiment='Reversion1'; % 1st underscore ...
% condition='';
% notes='';
% date_range={'08Jul2015', '08Jul2015'};
% only_labeled_dirs=1;
% prog_vers=2; % version of program to run - DONT CHANGE
% 
% 
% % ===== For opto analysis
% Params_glob.DataInfo='UnDir'; % To note information about data, 'Dir' and 'UnDir' means this is directed or undirected song.
% Params_glob.Fs=32000;
% 
% %         WHERE TO ALIGN SONGS
% Params_glob.SylTarg='c'; % aligns to onset of this
% Params_glob.PreNote='';
% 
% Params_glob.PreDur=0.10; % sec, how much data to take before targ
% Params_glob.PostSylOnsetDur=0.4; % sec, data post to take
% Params_glob.TargNoteNum=1; % template note directed to target
% Params_glob.TargTrigCatchRate=0; % catch trial fraction at target
% 
% %         STIM
% Params_glob.notenum_stim=0; % of Stim - CONFIRMED THIS WORKS
% Params_glob.StimDelay=55; % in ms, delay from trigger to stim on
% Params_glob.StimDur=100; % duration of Stim, in ms
% Params_glob.StimType='const'; % either 'constant' or 'pulse'
% 
% %     	for pitch contour
% Params_glob.PC_freqrange=[2400 3650]; % for both pitch contour and find note
% Params_glob.pc_harms=1;
% 
% %     	To sort into stim vs. notstim trials
% Params_glob.StimLatencyWindow=[0, 100]; % [a b] a and b are times (in negative ms) preceding syl onset. If stim is within this window, then will call a "stim trial."
% 
% % Time Window List
% % DURING WN (on window 3)
% Params_glob.TimeWindowList{1}=[115 155];
% Params_glob.TimeWindowList{2}=[197 240];
% Params_glob.TimeWindowList{3}=[264 272];
% 
% % Plotting over time
% Params_glob.SmthBin=20; % smooth # rends
% 

%% RUN - GET DIRECTORIES

BaseDir=pwd;

MetadataStruct=lt_metadata_collect;

% SECOND, extracts dirs that meet criteria given above.
prog_vers=2;
[MetadataStruct_filtered]=lt_metadata_find_dirs(MetadataStruct, Params_metadata.experiment, Params_metadata.condition, Params_metadata.notes, ...
    Params_metadata.date_range, Params_metadata.only_labeled_dirs, prog_vers);


% ========== OTHER PARAMS
NumDirs=length(MetadataStruct_filtered);


%% RUN - go through all dirs and run analysis.

for i=1:NumDirs
    close all;
    clear StatsStruct
    
    % === RESET PARAMS
    Params=Params_glob;
    
    % === go to dir
    cd([BaseDir '/' MetadataStruct_filtered(i).dirname]);
    
    % === WHAT BATCH OF SONGS?
    % since all dirs are either autolabeled or labeled, just take all labeled
    % songs
    lt_v2_db_transfer_calls(2);
    Params.batch='batch.labeled.all';
    
    % === STIM or NO STIM? (important to determine whether use 'All'; or
    % 'StimCatch' and 'StimNotCatch')
    if strcmp(MetadataStruct_filtered(i).notes, 'STIMoff')
        Params.ExptID='All';
    elseif strcmp(MetadataStruct_filtered(i).notes, 'STIMon')
        Params.ExptID='Stim';
    elseif strcmp(MetadataStruct_filtered(i).condition, 'STIMoff')
        Params.ExptID='All';
    elseif strcmp(MetadataStruct_filtered(i).condition, 'STIMon')
        Params.ExptID='Stim';
    else
        disp('PROBLEM - for one dir not sure if is stim or not (i.e. STIMoff or STIMon not in name)');
        return
    end
  
    
    % ===== DELETE OLD OPTO ANALYSIS FOLDER IF IT EXISTS (AND IF DESIRED)
    if Params_glob.Delete_old_analysis_folder==1
        if exist('lt_Opto_Stim_analy_v2','dir')==7 % then folder exists
            !rm -r lt_Opto_Stim_analy_v2
        end
    end
    
    
    % ==== EXTRACT DATA
    savefigs = 0;
    plotfigs = 0;
    [DatStructCompiled, Params]=lt_Opto_Stim_analy_EXTRACTDATA_v2(Params, ...
        SepBySyl, savefigs, plotfigs);
    
    
    % ====== CREATE A NEW FIELD CALLED NotStim_StimCatch, which combines
    % those two into one. 
    % important if on a given day there is both catch song (one epoch) and
    % catch trials (another epoch).
    DatStructCompiled.NotStim_StimCatch = [DatStructCompiled.NotStim DatStructCompiled.StimCatch];
    % -- sort in temporal order of trials
    ttmp = [DatStructCompiled.NotStim_StimCatch.datenum];
    [~, indtmp] = sort(ttmp);
    DatStructCompiled.NotStim_StimCatch = DatStructCompiled.NotStim_StimCatch(indtmp);
    
    
    % ===== PROCESS DATA
    if strcmp(Params.ExptID, 'All')
        Params.FieldsToCheck{1}='All';
    elseif strcmp(Params.ExptID, 'Stim')
%             Params.FieldsToCheck{1}='NotStim';
%             Params.FieldsToCheck{2}='StimCatch';
%             Params.FieldsToCheck{3}='StimNotCatch';
% modified: now already combines everything before runnign analyses below.
            Params.FieldsToCheck{1}='NotStim_StimCatch';
            Params.FieldsToCheck{2}='StimNotCatch';
    end
    
    % RUN
    plotfigs=0;
    [StatsStruct, Params]=lt_Opto_Stim_analy_PLOT_Compare2_v2(DatStructCompiled,Params, ...
        plotfigs);
    
    
    % ======= EXTRACT TIME WINDOW DATA
    close all
    
    for j=1:length(Params.TimeWindowList)
        TimeWind = Params.TimeWindowList{j};
        if ceil(Params.TimeWindowList{j})~=Params.TimeWindowList{j}; % means that there is fraction
            Params.TimeField{j}=['Wind' num2str(floor(TimeWind(1))) 'fr_' num2str(floor(TimeWind(2))) 'frms']; % put "fr" at end to signify fraction
        else
            Params.TimeField{j}=['Wind' num2str(TimeWind(1)) '_' num2str(TimeWind(2)) 'ms'];
        end
    end
    
    % RUN
    RunStats=1;
    savefigs = 0;
    [StatsStruct, Params]=lt_Opto_Stim_analy_PLOT_TimeWindow_v2(StatsStruct,Params, ...
        RunStats,KeepOutliers, savefigs);
    
    
    % ===== PLOT OVER TIME
    Params.SmthBin=15; % smooth # rends
    savefigs = 0;
    [StatsStruct, Params]=lt_Opto_Stim_analy_PLOT_TimeWindow_OverTime_v2(StatsStruct,Params,1,KeepOutliers, savefigs);
    close all;
    
end


