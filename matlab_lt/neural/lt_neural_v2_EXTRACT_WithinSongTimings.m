function lt_neural_v2_EXTRACT_WithinSongTimings(SummaryStruct, birdnum, neurnum)
%% lt 10/1/18 - to go from Segextract to song file data, need to map rends to within song timing

% output:
% saves to batch/channel dir, the mapping. (see readme below)


%%

i= birdnum;
ii = neurnum;
% SegExtract = segextract;


%% ============= extract metaparams
[SongDat, NeurDat, ~] = lt_neural_ExtractDat2(SummaryStruct, i, ...
    ii, 0);

fs = NeurDat.metaDat.fs;

%% ============ save metadata linking each syllable (in labels) to songfile
% and time in file

% ------- 1) get mapping between song num and the cumulative seconds up to the
% start of the song.
CumTimeAtSongOnset = [0 cumsum([NeurDat.metaDat.numSamps])];
CumTimeAtSongOnset(end) = [];
CumTimeAtSongOnset = CumTimeAtSongOnset./fs; % in seconds

% ------- 2) For each ind in labels, figure out how far from start of song.
nlabels = length(SongDat.AllLabels);
Songnames_all = cell(nlabels, 1);
Ons_all = nan(nlabels,1);
Offs_all = nan(nlabels,1);
for j=1:nlabels
    songnum = SongDat.AllSongNum(j);
    ons = SongDat.AllOnsets(j);
    offs = SongDat.AllOffsets(j);
    
    % ---- time at song onset?
    time_songonset = CumTimeAtSongOnset(songnum);
    
    % ---- Convert to time relative to start of song.
    ons = ons - time_songonset;
    offs = offs - time_songonset;
    
    % ============= OUTPUT
    fname = NeurDat.metaDat(songnum).filename;
    
    Songnames_all{j} = fname;
    Ons_all(j) = ons;
    Offs_all(j) = offs;
end


%% =============== SAVE OUTPUT
fnamesave = [SummaryStruct.birds(i).neurons(ii).dirname '/extract_WithinSongTimes.mat'];
timestruct = struct;
timestruct.Songnames_all = Songnames_all;
timestruct.Onsets = single(Ons_all);
timestruct.Offsets = single(Offs_all);
timestruct.readme = '[use this to map from regexp to song files] each rend is one label in SongDat. these onsets and offsets are relative to song onset';

save(fnamesave, 'timestruct')
disp(['DONE! saved at: ' fnamesave]);