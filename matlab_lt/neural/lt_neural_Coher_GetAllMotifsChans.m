function [CohAllTrials, Chanpairs, t_relons, ffbins] = lt_neural_Coher_GetAllMotifsChans(SummaryStruct, birdnum, neurnum, ...
    segextract, Chanlist, motif_predur, motif_postdur,PrePostRelSameTime)
%% LT 10/2/18 - for a given motif (i.e. segextract) collects coh across all desired chans

% === output
% CohAllTrials; cell, (chanpair x trial);
% Chanpairs = list of paired channels
% t_relons = timebins at middle of bins, relative to onset of entire motif
% ff = ffbins.


%%
dircoh = [fileparts(SummaryStruct.birds(birdnum).neurons(neurnum).dirname) '/COHERENCE'];

    
%%
% ================= FOR EACH TRIAL IN SEGEXTRACT, GET WITHIN SONG
% TIMES.

segextract = lt_neural_QUICK_GetWithinSongTime(SummaryStruct, birdnum, ...
    neurnum, segextract, motif_predur, motif_postdur,PrePostRelSameTime);

% ================== FOR EACH TRIAL EXTRACT COHERENCE BETWEEN ALL
% CHAN PAIRS
ntrials = length(segextract);
%         CohAllTrials = cell(ntrials,1); % trial, then chan pair
CohAllTrials = cell(nchoosek(length(Chanlist),2), ntrials); % trial, then chan pair

for tt = 1:ntrials
    
    fnthis = segextract(tt).song_filename;
    ons_withinsong = segextract(tt).WithinSong_OnsetTokenNoflank;
    off_withinsong = segextract(tt).WithinSong_OffsetMotifNoflank;
    
    % =============== COLLECT COHERENCE FOR THIS TRIAL ACROSS CHANS
    [Cohcell, Chanpairs, t, ffbins] = lt_neural_Coher_GetDat(...
        [dircoh '/' fnthis], Chanlist);
    
    
    % ================= trim to desired time window
    onthis = ons_withinsong - motif_predur; % to get predur
    
    if PrePostRelSameTime==1
        % then both on and off locked to the motif onset
        offthis = ons_withinsong + motif_postdur;
    elseif PrePostRelSameTime==0
        % then off locked to end of motif + postdur
        offthis = off_withinsong + motif_postdur;
    end
    assert((offthis - onthis)-(segextract(tt).global_offtime_motifInclFlank ...
        -segextract(tt).global_ontime_motifInclFlank)<0.001);
    
    tbins = t>=onthis & t<=offthis;
    
    
    for j=1:length(Cohcell)
        Cohcell{j} = Cohcell{j}(tbins,:);
    end
    
    % ================== OUTPUT
    CohAllTrials(:, tt) = Cohcell;
    
end
try
    % ---- convert tbins to time relative to onset of aligned syl
    t_relons = t(tbins);
    t_relons = t_relons - t_relons(1); % first bin is 0
    binsize = t(2)-t(1); % make first bin 0 plus half of binsize...
    t_relons = t_relons+binsize/2;
    t_relons = t_relons - motif_predur; % convert to relative to alignemnt.
catch err
    t_relons = [];
end
% ==
assert(size(Chanpairs,1) == size(CohAllTrials,1));
