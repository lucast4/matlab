function segextract = lt_neural_QUICK_GetWithinSongTime(SummaryStruct, bnum, ...
    neurnum, segextract, motif_predur, motif_postdur,PrePostRelSameTime)
%% lt 10/2/18 - convert global time (i.e. concat) to within song time
% segextract only have inds and times wrt to the concatenated songs. This
% converts that to within song times, to allow to extract stuff from
% filenames

% NOTE: runs a bunch of sanity checks so I am confident that there is match
% between segextract, within timing, and old notmat files.

%% params

% bnum = birdnum;
% neurnum = neuron id for this bird
% segextract. all is really needed is the gloval onsetsa nd offsets,.
% however, I include all of segextact becuase it is used for sanity checks
% (making sure what segextract thinks is the same as what I extract here).
% could modify to not need segextract if don't care about sanity checks.


%% RUn

% --- LOAD struct that tells you how to map trial to within song
% time
dirname_main = fileparts(SummaryStruct.birds(bnum).neurons(neurnum).dirname);
dirthis = SummaryStruct.birds(bnum).neurons(neurnum).dirname;
dat_withinsongtime = load([dirthis '/extract_WithinSongTimes.mat']);


%% ========= FOR EACH TRIAL, EXTRACT SONG FILE AND TIME IN FILE
ntrials = length(segextract);

for tt = 1:ntrials
    
    % ------ info from segextract regarding this trial
    segextract(tt).song_filename;
    segextract(tt).song_ind_in_batch;
    segextract(tt).global_ontime_motifInclFlank;
    segextract(tt).matchlabel;
    
    indtoken_glob = segextract(tt).global_tokenind_DatAlignedToOnsetOfThis;
    indstart_glob = segextract(tt).global_startind_motifonly;
    indend_glob = segextract(tt).global_endind_motifonly;
    
    
    % ==================== FIND TIME WITHIN SONG
    % ------ ind within song
    onsmotif_withinsong = dat_withinsongtime.timestruct.Onsets(indstart_glob);
    ons_withinsong = dat_withinsongtime.timestruct.Onsets(indtoken_glob);
    off_withinsong = dat_withinsongtime.timestruct.Offsets(indend_glob);
    
    
    % ==================== sanity checks
    % 1) filenames correct?
    assert(strcmp(dat_withinsongtime.timestruct.Songnames_all{indstart_glob}, segextract(tt).song_filename));
    assert(strcmp(dat_withinsongtime.timestruct.Songnames_all{indend_glob}, segextract(tt).song_filename));
    
    % 2) segextract motif matches that from the raw not.mat? (using
    % the within song time to extract from raw not.mat)
    % NOTE: this could fail if I have relabeled songs in this notmat...
    datnotmat = load([dirname_main '/' segextract(tt).song_filename '.not.mat']);
    [~, indon_tmp] = min(abs(datnotmat.onsets/1000 - onsmotif_withinsong));
    [~, indoff_tmp] = min(abs(datnotmat.offsets/1000 - off_withinsong));
    assert((indoff_tmp-indon_tmp)==(indend_glob-indstart_glob));
    assert(strcmp(datnotmat.labels(indon_tmp:indoff_tmp), segextract(tt).matchlabel), 'then motifs dont match... [if you relabeled or changed segmentation that culd explain, then is OK.]');
%     assert(strcmp(datnotmat.offsets(indoff_tmp)-datnotmat.onsets(indon_tmp), segextract(tt).matchlabel), 'then motifs dont match... [if you relabeled or changed segmentation that culd explain, then is OK.]');
    
    % 3) duration of motif must match 2 ways: 1) extracted
    % within-song time, and 2) segextract
    motifdur1 = off_withinsong - onsmotif_withinsong;
    
    motifdur2 = segextract(tt).motifsylOffsets(end) - segextract(tt).motifsylOnsets(1);
    assert(length(segextract(tt).motifsylOnsets)==(indend_glob-indstart_glob+1), 'assume motifsylonsets and offsets contains all syls in the motif, camped by instarta nd end)');
    assert(length(segextract(tt).motifsylOnsets)==length(segextract(tt).matchlabel), 'my assumption that motifsylonsets is length from on of first to end of last is wrong...');
    assert(motifdur2-motifdur1<0.001, 'why dur diff greater than one ms?');
    
    %         % 4) based on extracted motif dur
    %         if PrePostRelSameTime==1
    %            motifdur1 =
    %         elseif PrePostRelSameTime==0
    %
    %         end
    
    % 4) Make sure token durations are the same.
    tmp = indtoken_glob-indstart_glob; % to add to start to get token
    tokendur1 = segextract(tt).motifsylOffsets(1+tmp)-segextract(tt).motifsylOnsets(1+tmp);
    tokendur2 = datnotmat.offsets(indon_tmp+tmp) - datnotmat.onsets(indon_tmp+tmp);
    assert((tokendur1-tokendur2/1000)<0.001, 'why diff so large?');
    
    % =========================== OUTPUT
    segextract(tt).WithinSong_OnsetTokenNoflank = ons_withinsong;
    segextract(tt).WithinSong_OffsetMotifNoflank = off_withinsong;
    
    if isfield(segextract(tt), 'WithinSong_TokenOns')
        assert((segextract(tt).WithinSong_TokenOns-ons_withinsong)<0.001)
    end
end


