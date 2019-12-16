%% you will need to find the following file in the backup HD:
% '/bluejay0/bluejay2/lucas/analyses/neural/SummaryStruct.mat'
% then copy it to:
%'/bluejay5/lucas/analyses/neural/SummaryStruct.mat' in your local machine

% This is because the code in the next section will look for this. it saves
% all the metadata for all my experiments.
% alternatively save it anywhere but change the code below to look for it
% there. you can find the place in the code by running and seeing where it
% runs into an error.



%% 1) extract a summary structure that contains metadata about all neurons.
clear all; close all; fclose all;
BirdsToKeep = {'gr48bu5', 'wh72pk12', 'wh44wh39', 'pu69wh78'}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
% BirdsToKeep = {'wh72pk12'}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
% BrainArea = {'LMAN', 'RA', 'LMANoutside', 'RAoutside'}; % IF DOING NEGATIVE CONTROLS.
BrainArea = {'LMAN', 'RA'}; % 
% BrainArea = {}; % if want Sam/Mel data, must include "RAmel"
ExptToKeep = {};
RecordingDepth = [];
LearningOnly = 0;
BatchesDesired = {};
ChannelsDesired = [];
extractpreDatenums = 1;
onlySpikes = 1; % [default =1]if 1, then only keeps if have spiking data; if 0, then gets anything.
% if 2, then is for LFP analysis - i.e. for wh72 gets only if not SU or MU
% set to 2 if want to get data for cohere/xcov.
% NOTE: use 1 if want to look at spiking experiments across birds...
[NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database(BirdsToKeep, ...
    BrainArea, ExptToKeep, RecordingDepth, LearningOnly, BatchesDesired, ChannelsDesired, ...
    extractpreDatenums, onlySpikes);

% --- load all neurons
if (0)
    [NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database;
end



%% NOTE: IMPORTANT
% The following functions extract data by finding data by looking at the
% directory saved in 
% SummaryStruct.birds(1).neurons(1).dirname, where 1 and 1 are the bird and
% neuron you want. You need to change dirname to make sense for wherever
% you saved the data. You can do this my making some sort of for loop that
% modifies the string in dirname so that it points to wherever you have my
% raw data saved.

%% note:
% to get metadata for bird 1, neuron 1, do this:
SummaryStruct.birds(1).neurons(1)
% The fields should be somewhat self-explanatory. Ask me if you for help.


%% 2) For a given neuron you run the following to extract all data
% NOTE: You will have to run this for each neuron you want to analyze (e.g. in a
% for loop)

i = 1; % bird of interest
ii= 1; % neuron of this bird.
try
[SongDat, NeurDat, Params] = lt_neural_ExtractDat2(SummaryStruct, i, ii);
catch err
    disp('This may be an error that the directory is incorrect. you need to change "dirname" to that in the backup harddrive');
end

% SongDat is useful. It has all labels across all song files concatenated
% into only long string. All the fields are aligned perfectly and refer to
% that one long concatenated file. Most relevant to you are "AllLabels".

% NeurDat is useful. It has neural data, aligned to the song data in
% SongDat, so it is also the long concatenated file (neural). useful to you
% is NeurDat.spikes_cat.cluster_class, which is a T x 2 matrix. The first
% col are cluster ids, second column are spike times (within this long
% file). cluster classes are almost all 1, since I combined all into one MU
% site usually).

%% extracting motif data.
% You can choose to work directly with SongDat and NeurDat. Alternativeyl I
% have written a lot of helper functions to do things with that. Here is
% one that may be useful. Give it a regexp for a sequence you want to
% extract and it will extract all cases where it finds that sequence, along
% with neural data and bunch of other things.

regexp_str = 'n(a)abh'; % this means get all cases of abc, and align at onset of b. the parantheses
% area always interpreted as the syl to align to. You must have one and
% only one syl with parantheses. 
motif_predur = 0.2; % seconds to extract preceding syl
motif_postdur = 0.5;
alignByOnset = 1; % if 0, then algins by syl offset.
WHOLEBOUTS_edgedur = ''; % see within.
FFparams = '';
keepRawSongDat = 0; % if 1, then keeps song dat. avoid since large file unless need.
suppressout = 1;
collectWNhit = 0; % only useful to look at leanring.
collectWholeBoutPosition = 0; 
LearnKeepOnlyBase = 1; % 1 then excludes any songs after WN begins, if it is a trainign expt.
preAndPostDurRelSameTimept = 1; % 1: motif_predur will be algiend to the syl in parantehses. keep 1 to have all trials same dur.
RemoveIfTooLongGapDur = 0; % see within.
clustnum = [];
extractDirSong = 0; % to get directed song.
keepRawNeuralDat = 0; % NOTE: will be bandpass filtered in spike range [large size, avoid unless needed.]

[SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
    regexp_str, motif_predur, motif_postdur, alignByOnset, WHOLEBOUTS_edgedur, FFparams, ...
    keepRawSongDat, suppressout, collectWNhit, collectWholeBoutPosition, LearnKeepOnlyBase, ...
    preAndPostDurRelSameTimept, RemoveIfTooLongGapDur, clustnum, extractDirSong, ...
    keepRawNeuralDat);


%% example using SegmentsExtract

SegmentsExtract(1) % all data for trial 1.
SegmentsExtract(1).spk_Times % spike times relative to the onset of data for this trial.
SegmentsExtract(1).global_tokenind_DatAlignedToOnsetOfThis; % anything "global" refers to the concatenated file like in NeurDat
SegmentsExtract(1).motifsylOnsets % onset times of the syls in this motif, relative to the onset of data.

%% an example of stuff you can do with SegmentsExtract
figure; hold on;

% will use this in a bit
onsets = SegmentsExtract(1).motifsylOnsets;
offsets = SegmentsExtract(1).motifsylOffsets;

% get smoothed firing rates and plot summary
subplot(2,1,1); hold on;
title('smoothed fr');

SegmentsExtract = lt_neural_SmoothFR(SegmentsExtract);
frmat = [SegmentsExtract.FRsmooth_rate_CommonTrialDur];
x = SegmentsExtract(1).FRsmooth_xbin_CommonTrialDur;
plot(x, frmat, '-k'); 
shadedErrorBar(x', mean(frmat,2), std(frmat,[],2), {'Color', 'r'}, 1);

% overlay syl timings using jsut the first trial (not best, should maybe
% take average over trials).
lt_neural_QUICK_PlotSylPatches(onsets, offsets, [0 max(frmat(:))], [], 'r')

% plot rasters
subplot(2,1,2); hold on;
title('rasters');
for i=1:length(SegmentsExtract)
    seg = SegmentsExtract(i);
    lt_neural_PLOT_rasterline(seg.spk_Times, i, 'k', 0);
end
lt_neural_QUICK_PlotSylPatches(onsets, offsets, [0 i], [], 'r')


