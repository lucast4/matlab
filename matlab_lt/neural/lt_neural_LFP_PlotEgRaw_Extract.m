function [DatAll, t_onoff, fs, bregionlist, chanlist_toget, i, ii, mm] = ...
    lt_neural_LFP_PlotEgRaw_Extract(PARAMS, SwitchStruct, MOTIFSTATS_pop, SummaryStruct, ...
    SwitchCohStruct, birdplot, exptplot, swplot, motifplot, ...
    extrapad)
%%  lt 11/14/18 - plots example raw traces, + filtered, indiviodual trials

% NOTE: this is diffeerent from other raw plot functions in that this also
% gets raw data (unfiltered ata ll). Useful if want to confirm is not due
% to artifacts, etc.
% NOTE THIs might take some time to run because needs to load raw data.

% ===== FOR A GIVEN EXPERIMENT, PLOT BASE AND WN RAW DATA
% birdplot = 'pu69wh78';
% exptplot = 'RALMANOvernightLearn1';
% swplot = 1;
% motifplot = []; % [string] leave blank for target
% extrapad = 0.05; % seconds, pre and post...


% OUTPUTS:
% DatAll, cell of raw dats. size chan x1. chans are in sorted order (the
% same as the chans in LFPstruct...
% t_onoff onset and offset of data (relative to syl onset)
% fs, for raw data...
% bregionlist
% chanlist_toget
%% sandbox. IGNORE.
if (0) % here just developing code to extract raw dat...
    
    % ==== uses segextract
    i=1;
    ii=1;
    set = 4;
    mm=1;
    motifpredur = PARAMS.motif_predur;
    extrapad = 0.05; % extra on flank, to match extraction for LFP
    
    % ================ GET DATA STRUCTURES
    segextract = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(set).motif(mm).SegExtr_neurfakeID(1).SegmentsExtract;
    neurlist = MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons{set};
    chanlist = [SummaryStruct.birds(i).neurons(neurlist).channel];
    dirmain = cellfun(@fileparts, {SummaryStruct.birds(i).neurons(neurlist).dirname}, 'UniformOutput', 0);
    dirmain = unique(dirmain); dirmain = dirmain{1};
    
    % =============== RUN
    if (0)
        % TOO SLOW!! - since has to load for aech song...
        [DatAll, Tall, Chanlist] = lt_neural_QUICK_Segextr_GetRawNeural(segextract, ...
            dirmain, chanlist, motifpredur, extrapad);
        
    else
        % ================
        segextract = lt_neural_QUICK_GetRawNeural(segextract, SummaryStruct, i, neurlist(1), ...
            extrapad, motifpredur);
    end
end


%% #########################################################

% =============== GET PARAMS
i = find(strcmp({MOTIFSTATS_pop.birds.birdname}, birdplot));
ii = find(strcmp({MOTIFSTATS_pop.birds(i).exptnum.exptname}, exptplot));
if isempty(motifplot)
    % get target
    motifplot = SwitchStruct.bird(i).exptnum(ii).switchlist(swplot).learningDirs{1};
    if length(SwitchStruct.bird(i).exptnum(ii).switchlist(swplot).learningDirs)>2
        disp('NOTE: more than one target, will choose first one...!!');
        pause
    end
end
mm = find(strcmp({SwitchCohStruct.bird(i).exptnum(ii).switchlist(swplot).motifnum.motifname}, motifplot));
setnum = SwitchCohStruct.bird(i).exptnum(ii).switchlist(swplot).motifnum(mm).neursetused;

%% =============== for all channels in this set of data, collect all trial
% raw data
segextract = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(setnum).motif(mm).SegExtr_neurfakeID(1).SegmentsExtract;
neurlist = MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons{setnum};
chanlist_neurons = [SummaryStruct.birds(i).neurons(neurlist).channel];
bregionlist = {SummaryStruct.birds(i).neurons(neurlist).NOTE_Location};

% ---- sort to match order of previous extraction data
[~, indsort] = sort(chanlist_neurons);
chanlist_neurons = chanlist_neurons(indsort);
neurlist = neurlist(indsort);
bregionlist = bregionlist(indsort);

motifpredur = PARAMS.motif_predur;
fs = segextract.fs;

% ---- only relevant channels (e.g. LMAN-RA)
chanlist_toget = chanlist_neurons; % using these maintains order, so that bregions are accurate.
% chanlist_toget = SwitchCohStruct.bird(i).exptnum(ii).switchlist(swplot).motifnum(mm).lfpall_chans;

%% =================== GO THRU EACH NEURON. IF ITS CHAN IS DESEIRED, THEN
% GET
DatAll = cell(length(chanlist_toget),1); % nchans
t_onoff = cell(length(chanlist_toget),1); % nchans
for j=1:length(chanlist_toget)
    
    neurthis = neurlist(chanlist_neurons==chanlist_toget(j));
    assert(length(neurthis)==1);
    
    tmp = lt_neural_QUICK_GetRawNeural(segextract, SummaryStruct, i, neurthis, ...
        extrapad, motifpredur, 1);
    
    % --- extract all trials neural from tmp
    dum = [tmp.neural_rawdat];
    dum2 = [tmp.neural_rawdat_tOnOff_reltok];
    
    DatAll{j} = dum;
    t_onoff{j} = median(dum2,2);
end

