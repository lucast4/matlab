function [CohAllTrials, Chanpairs, t_relons, ffbins] = lt_neural_Coher_GetAllMotifsChans(SummaryStruct, birdnum, neurnum, ...
    segextract, Chanlist, motif_predur, motif_postdur,PrePostRelSameTime)
%% LT 10/2/18 - for a given motif (i.e. segextract) collects coh across all desired chans

% === output
% CohAllTrials; cell, (chanpair x trial);
% ohAllTrials =
%
%   3×60 cell array
%
%   Columns 1 through 6
%
%     [20×19 single]    [20×19 single]    [20×19 single]    [20×19 single]    [20×19 single]    [20×19 single]
%     [20×19 single]    [20×19 single]    [20×19 single]    [20×19 single]    [20×19 single]    [20×19 single]
%     [20×19 single]    [20×19 single]    [20×19 single]    [20×19 single]    [20×19 single]    [20×19 single]
% Chanpairs = list of paired channels (numpairs x 2)
% Chanpairs =
%
%      9    17
%      9    21
%     17    21

% t_relons = timebins at middle of bins, relative to onset of entire motif
% ff = ffbins.

% % ========= TWO METHODS
% preextracted; then precalcualted coherence values
% reextract; then extract lfp and then calculate coherence
cohmethod = 'reextract';


% =============== CHRONUX PARAMS
lt_switch_chronux(1);
movingwin = [0.1 0.01];

params = struct;
params.fpass = [1/movingwin(1) 150];
w = 30; % in hz, for desired frequency resolution of tapers. % note, t is set to movingwin(1)
tw = movingwin(1)*w;
params.tapers = [tw 2*tw-1];
params.Fs = 1500; % hard coded fs for LFP;

extrapad = movingwin(1)/2; % in sec, how much additioanl time to add to account for window size of coherogram.


%%
dircoh = [fileparts(SummaryStruct.birds(birdnum).neurons(neurnum).dirname) '/COHERENCE'];

Chanlist = sort(Chanlist);
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



if strcmp(cohmethod, 'preextracted')
    %% ============= METHOD 1 - previoupsly extracted coh
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
    
    
elseif strcmp(cohmethod, 'reextract')
    %% ============= METHOD 2 - extracts LFP and recalcualtes
    
    % ============== 1) EXTRACT ALL LFP FOR SEGEXTRACT
    fnamebase = fileparts(SummaryStruct.birds(birdnum).neurons(neurnum).dirname);
    [LFPall, Tall] = lt_neural_QUICK_Segextr_GetLFP(segextract, fnamebase, Chanlist, ...
        motif_predur, extrapad);
    
    % ============== CALCULATE COHERENCE FOR ALL DESIRED PAIRS
    % --- go thru all chan pairs and calculate coherence on each trials
    rowcount = 1;
    Chanpairs = [];
    for c = 1:length(Chanlist)
        for cc = c+1:length(Chanlist)
            chan1 = Chanlist(c);
            chan2 = Chanlist(cc);
            Chanpairs = [Chanpairs; [chan1 chan2]];
            % =========== go thru all trials getting coherence
            dat1 = cell2mat(LFPall(c,:));
            dat2 = cell2mat(LFPall(cc,:));
            
            % ===========
            %                 [C,phi,S12,S1,S2,t,f] = cohgramc(dat1, dat2, movingwin, params);
            [C,~,~,~,~,t,ffbins] = cohgramc(dat1, dat2, movingwin, params);
            C = single(C);
            
            % ======== SAVE INTO CELL ARRAY
            tmp = squeeze(mat2cell(C, size(C,1), size(C,2), ones(1, size(C,3))))';
            CohAllTrials(rowcount, :) = tmp;
            rowcount = rowcount+1;
        end
    end
    
    t_relons = t-(motif_predur+extrapad);
    %     t_relons = Tall{1};
    
end
% ==
% if floor(now)~=737362
assert(size(Chanpairs,1) == size(CohAllTrials,1));
% end


lt_switch_chronux(0)

