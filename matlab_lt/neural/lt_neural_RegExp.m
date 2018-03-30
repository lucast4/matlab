%% lt 8/18/16 - input regexpt, outputs segments,each of which is the regexp +/- duration.

function [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
    regexpr_str, predur, postdur, alignByOnset, WHOLEBOUTS_edgedur, FFparams, ...
    keepRawSongDat, suppressout, collectWNhit, collectWholeBoutPosition, LearnKeepOnlyBase, ...
    preAndPostDurRelSameTimept, RemoveIfTooLongGapDur, clustnum, extractDirSong)
%% lt 2/19/18 - labels rendition as in DIR or UNDIR song
% does this by looking at suffix in song file name

% MADE default to ignore DIR

if ~exist('extractDirSong', 'var')
    extractDirSong = 0;
end

%% lt 1/14/18 - throwing out all spikes with clustnum=0 (is noise) [AUTO]
% ALSO - if specify clustnum then will only keep spikes that are that
% number

if ~exist('clustnum', 'var')
    clustnum = [];
end

%% lt 12/20/17 - NOTE, made WN detect independent of FFparams ... (i.e. moved "end")

%% lt 12/20/17 - outputs smoothed rectified amplitude

if ~exist('keepSmAmpl', 'var')
    keepSmAmpl = 0;
end


%% lt 11/29/17 - made RemoveIfTooLongGapDur default to 1 [back to 0]
%%  lt 8/15/17 - added default gap duration (between all syls in a given motif)

if ~exist('RemoveIfTooLongGapDur', 'var')
    RemoveIfTooLongGapDur=0;
end
maxgapdur = 1; % sec

%% note on FF stuff
% WNhit collection is only performed if FF collection is performed - can
% modify easily to make them independent. Needs raw audio to extract WNhit
% data...

%% note 3/21/17
% Might miss some whole bouts. reason is potential motifs are first detected by
% the onset_pre_threshold criterion (which is about 1s). then I ask how many of
% those are preceded and followed by gap larger than WHOLEBOUTS_edgedur, which
% is larger than onset_pre_threshold. So some will fail. If want all that pass
% onset_pre_threshold to be included, then make  onset_pre_threshold = WHOLEBOUTS_edgedur;
% (i.e. make WHOLEBOUTS_edgedur shorter). problem with that is that then might
% get some bouts with some calls ppreceding and following motif.

%% keepRawSongDat

% will keep raw dat if ask for FF.
% will discard if don't ask, (unless set keepRawSongDat = 1 [default = 0]);

if ~exist('keepRawSongDat', 'var')
    keepRawSongDat = 0;
end

if ~exist('suppressout', 'var')
    suppressout = 1;
end

if ~exist('collectWNhit', 'var')
    collectWNhit=1;
end

if ~exist('collectWholeBoutPosition', 'var')
    %     collectWholeBoutPosition=0; % to get position of a given datapoint within its bout
    collectWholeBoutPosition=1; % to get position of a given datapoint within its bout
end

if ~exist('LearnKeepOnlyBase', 'var')
    LearnKeepOnlyBase=0;
end

if ~exist('preAndPostDurRelSameTimept', 'var')
    preAndPostDurRelSameTimept=0; % if 1, then will get pre and postdur aligned to either onset or
    % offset of token (depending on "alignByOnset"). if 0, then predur is
    % rel to token, and postdur is relative to end of matched motif
    % (default)
end

if ~exist('FFparams' ,'var')
    FFparams = [];
end

if isempty(FFparams)
    FFparams.collectFF=1; % note, will try to collect FF for each motif inputed in the cell array. will
    FFparams.FF_PosRelToken=0; % syl to get FF of, relative to token (i.e. -1 is 1 before token;
    % +1 is 1 after token
    FFparams.FF_sylName=''; % Optional: what syl do you expect this to be? if incompatible will raise error
end


%%

if collectWNhit==1
    keepSmAmpl=1; % if collecting WN, then not much extra work to keep smoothe dampltide.
end


%% --
% FFparams.collectFF=1;
% FFparams.cell_of_FFtimebins={'h', [0.042 0.058], 'b', [0.053 0.07], ...
%             'v', [0.052 0.07]}; % in sec, relative to onset (i.e. see vector T)
% FFparams.cell_of_freqwinds={'h', [1100 2600], 'b', [2400 3500], ...
%             'v', [2450 4300]};
% FFparams.FF_PosRelToken=1; % syl to get FF of, relative to token (i.e. -1 is 1 before token;
%     % +1 is 1 after token
% FFparams.FF_sylName='b'; % Optional: what syl do you expect this to be? if incompatible will raise error
%     % not required (can set as []);



%%
% regexpr_str='ghh'; % if 'WHOLEBOUTS' then will automatically extract
% bouts by onset and offsets.
% predur=4; % sec, time before token onset (or offset, depending on
% alignByOnset)
% postdur=4; % sec
% SongDat, NeurDat, Params, shoudl correspond to same dataset, extracted using lt_neural_ExtractDat
% alignByOnset=1, then aligns all by the onset of token ind. if 0, then by offset of token ind.

%     WHOLEBOUTS_edgedur=6; % OPTIONAL (only works if regexpr_str='WHOLEBOUTS', only keeps
%     % those motifs that have long enough pre and post - LEAVE EMPTY TO GET ALL BOUTS


if ~exist('FFparams', 'var')
    FFparams.collectFF=0;
end

if isempty(FFparams)
    FFparams.collectFF=0;
end


if ~exist('WHOLEBOUTS_edgedur', 'var');
    WHOLEBOUTS_edgedur=[];
end

if ~exist('alignByOnset', 'var');
    alignByOnset=1;
end
if isempty('alignByOnset')
    alignByOnset=1;
end

UseLastSylAsToken=0;

%% === extract stuff from inputs

AllLabels=SongDat.AllLabels;
AllOnsets=SongDat.AllOnsets;
AllOffsets=SongDat.AllOffsets;
% spikes_cat=NeurDat.spikes_cat;
% metaDat=NeurDat.metaDat;

fs=NeurDat.metaDat(1).fs;

Params.REGEXP.predur=predur;
Params.REGEXP.postdur=postdur;
Params.REGEXP.alignByOnset = alignByOnset;

% PARAMS ONLY USED FOR WHOLEBOUTS EXTRACTION
%      onset_pre_threshold=2; % OLD (changed on 3/21/17, based on wh6pk36)
onset_pre_threshold=1; % in seconds, minimum time preceding and following bout [DO NOT CHANGE! THIS IS FOR "DEFINING" SONG BOUT]
% If you only want to keep those with a long enough quiet time, then do
% following

%     offset_post_threshold=3; % in seconds
min_syls_in_bout=10;

if isfield(NeurDat, 'TotalSamps')
    TotalSamps = NeurDat.TotalSamps;
else
    TotalSamps = sum([NeurDat.metaDat.numSamps]);
end

TotalDurSec = TotalSamps/fs;

%% params for smoothing ampl

if keepSmAmpl==1
        sm_win = 0.005; % 
        len = round(fs*sm_win);
        h_smooth   = ones(1,len)/len;
end


%% ------ find motifs

SegmentsExtract=struct;

if strcmp(regexpr_str, 'WHOLEBOUTS')
    gapdurs=[AllOnsets(1) AllOnsets(2:end)-AllOffsets(1:end-1) TotalDurSec-AllOffsets(end)]; % first gap is just dur from start of file to first onset
    
    % - a motif is
    potentialOnsets=find(gapdurs>onset_pre_threshold); % ind of onset syl
    
    inds_tmp=find(diff(potentialOnsets)>=min_syls_in_bout); % indexes thos bouts that pass min syls criterion
    
    bout_firstsyls=potentialOnsets(inds_tmp);
    bout_lastsyls=potentialOnsets(inds_tmp+1)-1; % minus 1 because want to end of the current bout, not the start of the next bout.
    
    bout_firstsyls_beforeflankfilter = bout_firstsyls;
    bout_lastsyls_beforeflankfilter = bout_lastsyls;
    
    
    % ------------
    % === ONLY KEEP THOSE WHICH HAVE FLANKING GAP DURS LONGER THAN
    % CRITERION
    if ~isempty(WHOLEBOUTS_edgedur)
        gapdurs_subset=gapdurs(gapdurs>onset_pre_threshold);
        inds_tmp2=gapdurs_subset(inds_tmp)>WHOLEBOUTS_edgedur ...
            & gapdurs_subset(inds_tmp+1)>WHOLEBOUTS_edgedur;
        
        bout_firstsyls=bout_firstsyls(inds_tmp2);
        bout_lastsyls=bout_lastsyls(inds_tmp2);
    else
        disp('PROBLEM - WHOLEBOUTS_edgedur is empty!!!');
    end
    
    % ============================== TROUBLESHOOTING - PLOT ALL ONSETS IN A
    % LINE PLOT
    if (0)
        % 1) song
        hsplots = [];
        lt_figure; hold on;
        
        hsplot= lt_subplot(2,1,1); hold on;
        plot([1:length(SongDat.AllSongs)]./NeurDat.metaDat(1).fs, ...
            SongDat.AllSongs, 'k');
        title('song');
        hsplots = [hsplots hsplot];
        
        % 2) syls
        hsplot = lt_subplot(2,1,2); hold on;
        for i=1:length(AllOnsets)
            line([AllOnsets(i) AllOffsets(i)], [0 0], 'LineWidth', 2);
        end
        
        for i=1:length(bout_firstsyls)
            line([AllOnsets(bout_firstsyls(i)) AllOnsets(bout_firstsyls(i))]-WHOLEBOUTS_edgedur,...
                ylim, 'Color', 'g');
            line([AllOffsets(bout_lastsyls(i)) AllOffsets(bout_lastsyls(i))]+WHOLEBOUTS_edgedur, ...
                ylim , 'Color', 'r');
        end
        
        for i=1:length(bout_firstsyls_beforeflankfilter)
            line([AllOnsets(bout_firstsyls_beforeflankfilter(i)) AllOnsets(bout_firstsyls_beforeflankfilter(i))],...
                ylim, 'Color', 'b');
            line([AllOffsets(bout_lastsyls_beforeflankfilter(i)) AllOffsets(bout_lastsyls_beforeflankfilter(i))], ...
                ylim , 'Color', 'b');
        end
        hsplots = [hsplots hsplot];
        
        title('blue: potential motifs...redgreen: those that pass edge_dur test');
        xlabel('time of syl (s)');
        
        linkaxes(hsplots, 'x');
    end
    
    % ----------------------
    
    
    % ---- convert to variables used below
    tokenExtents=bout_firstsyls;
    startinds=bout_firstsyls;
    endinds=bout_lastsyls;
    matchlabs={};
    for j=1:length(bout_firstsyls)
        matchlabs{j}=AllLabels(bout_firstsyls(j):bout_lastsyls(j));
    end
    
    % --- if want to use last syl as token.
    if UseLastSylAsToken==1;
        tokenExtents=bout_lastsyls;
        startinds=bout_lastsyls;
        endinds=bout_lastsyls;
        matchlabs={};
        for j=1:length(bout_lastsyls)
            matchlabs{j}=AllLabels(bout_lastsyls(j):bout_lastsyls(j));
        end
    end
else
    % then is a regexp string
    % MODIFIED TO ALIGN TO ONSET OF TOKEN
    
    % --------- method 1
    
    if (0)
        [startinds1, endinds1, matchlabs1, tokenExtents1]=regexp(AllLabels, regexpr_str, 'start', 'end', ...
            'match', 'tokenExtents');
        functmp = @(X) X(1);
        tokenExtents1=cellfun(functmp, tokenExtents1); % convert from cell to vector.
    end
    
    % ------- method 2
    [startinds, tokenExtents, endinds, matchlabs] = lt_neural_QUICK_regexp(AllLabels, regexpr_str);

%     indstmp = regexp(regexpr_str, '[\(\)]'); % find left and right parantheses
%     
%     % confirm that there is exactly one token syl.
%     assert(length(indstmp)==2,'sdafasd');
%     assert(indstmp(2)-indstmp(1) ==2, 'sdfsdaf');
%     
%     % -
%     tmptmp = [regexpr_str(1:indstmp(1)-1) ...
%         regexpr_str(indstmp(1)+1) regexpr_str(indstmp(2)+1:end)];
%     
%     %     regexpr_str2 = [regexpr_str(1:indstmp(1)-1) '(?=' ...
%     %         regexpr_str(indstmp(1)+1) regexpr_str(indstmp(2)+1:end) ')']; % for lookahead assertion
%     
%     regexpr_str2 = [tmptmp(1) '(?=' tmptmp(2:end) ')']; % for lookahead assertion
%     
%     [startinds, ~, ~]=regexp(AllLabels, regexpr_str2, 'start', 'end', ...
%         'match');
%     
%     strlength = length(regexpr_str)-2;
%     
%     tokenExtents = startinds + indstmp(1)-1; % i.e. where token was
%     endinds = startinds + strlength -1;
%     
%     % - get match syls
%     if size(startinds,2)==1
%         startinds = startinds';
%     end
%     
%     indmat = [];
%     for j=1:strlength
%         indmat = [indmat startinds'+j-1];
%     end
%     
%     %     if size(AllLabels,1)==1
%     %     matchlabs = mat2cell(AllLabels(indmat)', ones(size(indmat,1),1))';
%     %     else
%     if size(indmat,2)==1
%         % then is one col (i.e. only one syl in motif) so need to make sure
%         % output is col
%         matchlabs = mat2cell(AllLabels(indmat)', ones(size(indmat,1),1))';
%     else
%         matchlabs = mat2cell(AllLabels(indmat), ones(size(indmat,1),1))';
%     end
%     
    % ---- compare methods
    if (0)
        if length(startinds) == length(startinds1)
            assert(isempty(setxor(endinds, endinds1)), 'asfasd');
            assert(isempty(setxor(startinds1, startinds)), 'asfasd');
            assert(isempty(setxor(matchlabs1, matchlabs)), 'asdfsd');
            assert(all(tokenExtents1 == tokenExtents), 'asdfsd');
            
        else
            assert(isempty(setdiff(startinds1, startinds)), 'asfasdf'); % i.e. old version is proper subset of new version
        end
    end
end


%% =================== IF WANT TO GET WHOLEBOUTS INFO TO GET DATAPOINT
% POSITION IN BOUT (BUT NOT ACTUALLY EXTRACT WHOLEBOUTS AS MAIN DATAPOINT)
if collectWholeBoutPosition==1 & strcmp(regexpr_str, 'WHOLEBOUTS') ==0
    gapdurs=[AllOnsets(1) AllOnsets(2:end)-AllOffsets(1:end-1) TotalDurSec-AllOffsets(end)]; % first gap is just dur from start of file to first onset; % last gap is just offset to end of file
    
    % - a motif is
    potentialOnsets=find(gapdurs>onset_pre_threshold); % ind of onset syl
    
    inds_tmp=find(diff(potentialOnsets)>=min_syls_in_bout); % indexes thos bouts that pass min syls criterion
    
    wholebout_firstsyls=potentialOnsets(inds_tmp);
    wholebout_lastsyls=potentialOnsets(inds_tmp+1)-1; % minus 1 because want to end of the current bout, not the start of the next bout.
    
    % ---- get match labels
    wholebout_matchlabs={};
    for j=1:length(wholebout_firstsyls)
        wholebout_matchlabs{j}=AllLabels(wholebout_firstsyls(j):wholebout_lastsyls(j));
    end
    
    % TROUBLESHOOT
    if (0)
        % 1) song
        hsplots = [];
        lt_figure; hold on;
        
        %         hsplot= lt_subplot(2,1,1); hold on;
        %         plot([1:length(SongDat.AllSongs)]./NeurDat.metaDat(1).fs, ...
        %             SongDat.AllSongs, 'k');
        %         title('song');
        %         hsplots = [hsplots hsplot];
        
        % 2) syls
        hsplot = lt_subplot(2,1,2); hold on;
        for i=1:length(AllOnsets)
            line([AllOnsets(i) AllOffsets(i)], [0 0], 'LineWidth', 2);
            %             lt_plot_text(AllOnsets(i), 0.01, AllLabels(i));
        end
        
        %         for i=tokenExtents
        %             lt_plot_text(AllOnsets(i), 0.01, AllLabels(i));
        %         end
        for i=1:length(wholebout_firstsyls)
            line([AllOnsets(wholebout_firstsyls(i)) AllOnsets(wholebout_firstsyls(i))],...
                ylim, 'Color', 'g');
            line([AllOffsets(wholebout_lastsyls(i)) AllOffsets(wholebout_lastsyls(i))], ...
                ylim , 'Color', 'r');
        end
        
        xlabel('time of syl (s)');
        
        plot(AllOnsets(638), 0, 'ok')
    end
    
end

%%
HitSyls_TEMP={};
MisSyls_TEMP={};

BoutNumsAll = [];
PosInBoutAll = [];
RendInBoutAll = [];

segstoremove = []; % those that have too long gap dur

% -- for each match ind, extract audio + spikes
for i=1:length(tokenExtents)
    
    % --- make sure all gap durations are shorted than threshold
    allgapdurs = AllOnsets(startinds(i)+1:endinds(i)) - AllOffsets(startinds(i):endinds(i)-1);
    if any(allgapdurs > maxgapdur)
        segstoremove = [segstoremove i];
    end
    
    
    % ------------------- EXTRACT SYL DUR AND FLANKING GAP DURS
    ind=tokenExtents(i);
    syldur = AllOffsets(ind) - AllOnsets(ind);
    if ind>1
        gapdur_pre = AllOnsets(ind) - AllOffsets(ind-1);
    else
        gapdur_pre = nan;
    end
    
    if ind == length(AllOnsets);
        gapdur_post = nan;
    else
        gapdur_post = AllOnsets(ind+1) - AllOffsets(ind);
    end
    
    SegmentsExtract(i).Dur_syl = syldur;
    SegmentsExtract(i).Dur_gappre = gapdur_pre;
    SegmentsExtract(i).Dur_gappost = gapdur_post;
    
    
    
    % on time
    ind=tokenExtents(i);
    if alignByOnset==1
        aligntime=AllOnsets(ind); % sec
    elseif alignByOnset==0
        % align by offset of token syl
        aligntime=AllOffsets(ind); % sec
    end
    
    
    ontime=aligntime-predur; % - adjust on time
    onsamp=round(ontime*fs);
    
    % off time
    if preAndPostDurRelSameTimept==0
        % default, align to end of motif
        ind=endinds(i);
        offtime=AllOffsets(ind);
        offtime=offtime+postdur;
    elseif preAndPostDurRelSameTimept==1
        offtime = aligntime+postdur;
    end
    offsamp=round(offtime*fs);
    
    
    % --------------------- EXTRACT ONTIMES AND OFFTIMES OF ALL SYLS IN
    % MOTIF
    sylpositions = startinds(i):endinds(i);
    motifsylOnsets = AllOnsets(sylpositions) - ontime;
    motifsylOffsets = AllOffsets(sylpositions) - ontime;
    SegmentsExtract(i).motifsylOnsets = motifsylOnsets;
    SegmentsExtract(i).motifsylOffsets = motifsylOffsets;
    
    
    %     spkinds=(NeurDat.spikes_cat.cluster_class(:,2) > ontime*1000) & ...
    %         (NeurDat.spikes_cat.cluster_class(:,2) < offtime*1000);
    if isfield(NeurDat, 'spikes_cat')
        spk_ClustTimes = NeurDat.spikes_cat.cluster_class((NeurDat.spikes_cat.cluster_class(:,2) > ontime*1000) & ...
            (NeurDat.spikes_cat.cluster_class(:,2) < offtime*1000), :); % in sec, relative to onset of the segment
        
        % -------- throw out anything wtih clust num 0 (is noise)
        spk_ClustTimes(spk_ClustTimes(:,1) == 0, :) = [];
        
        % -------- if desired, then only keep correct cluster
        if ~isempty(clustnum)
            spk_ClustTimes(spk_ClustTimes(:,1) ~= clustnum, :) = [];
        end
        
    else
        % then is RA data from Sober/Mel
        spiketimes = NeurDat.spiketimes(NeurDat.spiketimes>ontime & NeurDat.spiketimes<offtime);
    end
    
    if keepRawSongDat ==1
        assert(isfield(SongDat, 'AllSongs'), 'PROBLEM - need to extract songdat before running this');
        
        % this effectively does nothing if also collecting FF.
        %         AllSongs=SongDat.AllSongs;
        songseg=SongDat.AllSongs(onsamp:offsamp);
        SegmentsExtract(i).songdat=songseg;
    end
    
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % Extract FF of a specific syllable [OPTIONAL]
    if FFparams.collectFF==1
        
        PC=nan;
        FF=nan;
        T=nan;
        
        if isfield(SongDat, 'FFvals');
            %         if (1)
            % -- take old vals
            FF = SongDat.FFvals(tokenExtents(i));
            
        else
            %% === OLD METHOD (CALCULATE FROM RAW AUDIO HERE)
            disp('HAVE TO DO FFVALS MANUALLY :( ');
            % - collect
            %             AllSongs=SongDat.AllSongs;
            songseg=SongDat.AllSongs(onsamp:offsamp);
            
            FF_PosRelToken=FFparams.FF_PosRelToken;
            FF_sylName=FFparams.FF_sylName;
            cell_of_freqwinds=FFparams.cell_of_freqwinds;
            cell_of_FFtimebins=FFparams.cell_of_FFtimebins;
            
            prepad_forFF=0.015; postpad_forFF=0.015; % don't change, as this affects temporal window for FF
            
            indForFF=tokenExtents(i)+FF_PosRelToken;
            sylForFF=AllLabels(indForFF);
            
            ontime_forFF=AllOnsets(indForFF);
            ontime_forFF=ontime_forFF-prepad_forFF;
            onsamp_forFF=round(ontime_forFF*fs);
            
            offtime_forFF=AllOffsets(indForFF);
            offtime_forFF=offtime_forFF+postpad_forFF;
            offsamp_forFF=round(offtime_forFF*fs);
            
            songseg_forFF=SongDat.AllSongs(onsamp_forFF:offsamp_forFF);
            
            % -- debug
            if (0)
                lt_figure; hold on;
                lt_plot_spectrogram(songseg_forFF, fs, 1,0);
                pause
                close all;
            end
            
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EXTRACT FF
            % -- defaults, sets output (e.g., ff) to this if don't collect
            % actual FF
            collectFF=1;
            
            % ---------------- GO THROUGH POTENTIAL REASONS TO NOT COLLECT FF
            if ~isempty(FF_sylName)
                %             assert(strcmp(FF_sylName, sylForFF)==1, 'Syl is not the desired syl!! - stopping');
                if strcmp(FF_sylName, sylForFF)==0
                    collectFF=0;
                end
            end
            
            indtmp=find(strcmp(sylForFF, cell_of_freqwinds));
            if isempty(indtmp)
                % tel;l use
                disp('EMPTY FREQ WINDOW!! - giving nan for FF');
                collectFF=0;
            elseif isempty(cell_of_freqwinds{indtmp+1})
                disp('EMPTY TIME WINDOW!! - giving nan for FF');
                collectFF=0;
            end
            
            % ------ COLLECT
            if collectFF==1
                indtmp=find(strcmp(sylForFF, cell_of_freqwinds));
                F_high=cell_of_freqwinds{indtmp+1}(2);
                F_low=cell_of_freqwinds{indtmp+1}(1);
                
                indtmp=find(strcmp(cell_of_FFtimebins, sylForFF));
                mintime=cell_of_FFtimebins{indtmp+1}(1); % sec
                maxtime=cell_of_FFtimebins{indtmp+1}(2);
                
                [FF, PC, T]= lt_calc_FF(songseg_forFF, fs, [F_low F_high], [mintime maxtime]);
            end
        end
    
    SegmentsExtract(i).FF_val=FF;
    
    end
    %% ===================================================
    % FIGURE OUT IF WN HIT ON THIS TRIAL (based on clipping of sound)
    % - collect
    if collectWNhit==1

        if (0)
            % -- OLD, used FF posotion to determine where to check for WN
        FF_PosRelToken=FFparams.FF_PosRelToken;
        indForFF=tokenExtents(i)+FF_PosRelToken; % can use this to look for WN in flanking syls.
        else
            % -- NEW - look for WN overlayed on token syl.
        indForFF=tokenExtents(i); % can use this to look for WN in flanking syls.
        end
        prepad_forFF=0.015; postpad_forFF=0.015; % don't change, as this affects temporal window for FF
        
        ontime_forFF=AllOnsets(indForFF);
        ontime_forFF=ontime_forFF-prepad_forFF;
        onsamp_forFF=round(ontime_forFF*fs);
        
        offtime_forFF=AllOffsets(indForFF);
        offtime_forFF=offtime_forFF+postpad_forFF;
        offsamp_forFF=round(offtime_forFF*fs);
        
        songseg_forFF=SongDat.AllSongs(onsamp_forFF:offsamp_forFF);
        
        
        wasTrialHit=[];
        WNonset=[];
        WNoffset=[];
        if max(songseg_forFF)>3.2
            % --- debug, plot spectrogram and sound file for all hits
            %             figure; hold on;
            %             lt_subplot(1,2,1); hold on; plot(songseg_forFF);
            %             lt_subplot(1,2,2); hold on; lt_plot_spectrogram(songseg_forFF, fs, 1, 0);
            %             pause
            %             close all;
            % ----
            wasTrialHit=1;
            
%             WNonset=find(songseg>3.2, 1, 'first'); % timepoint of hit (for entire segment)
%             WNoffset=find(songseg>3.2, 1, 'last');
            
            % --- convert to start of syl (adding motifpredur)
            WNonset=find(songseg_forFF>3.2, 1, 'first'); % timepoint of hit (for entire segment)
            WNoffset=find(songseg_forFF>3.2, 1, 'last');
            
            WNonset = WNonset-(prepad_forFF*fs)+predur*fs;
            WNoffset = WNoffset-(prepad_forFF*fs)+predur*fs;
            
            HitSyls_TEMP=[HitSyls_TEMP songseg_forFF];
        else
            wasTrialHit=0;
            MisSyls_TEMP=[MisSyls_TEMP songseg_forFF];
        end
        
        SegmentsExtract(i).hit_WN=wasTrialHit;
        SegmentsExtract(i).WNonset_sec=WNonset/fs;
        SegmentsExtract(i).WNoffset_sec=WNoffset/fs;
        
    end
    
    %% get smoothed amplitude
    if keepSmAmpl==1

        songseg=SongDat.AllSongs(onsamp:offsamp);
        
        % ------- smooth and rectify...
        % 1) FILTER
        if (1)
        filter_type='hanningfirff';
        F_low  = 500;
        F_high = 8000;
        songseg_sm=bandpass(double(songseg),fs,F_low,F_high,filter_type);
        else
        songseg_sm = songseg - mean(songseg);
        end
        
        % 2) SQUARE
        songseg_sm = songseg_sm.^2; % -- square
        
        % 3) SMOOTH
        songseg_sm = conv(h_smooth, songseg_sm);
        offsetTMP = round((length(songseg_sm)-length(songseg))/2);
        songseg_sm=songseg_sm(1+offsetTMP:length(songseg)+offsetTMP);
% 
%         figure; hold on; plot(songseg_sm,'k')
    
% ---- convert to single and downsample
% songseg_smD = downsample(songseg_sm, fs/1000, round(fs/2000));
songseg_sm = downsample(songseg_sm, fs/1000);
t_songseg = 0:0.001:0.001*(length(songseg_sm)-1);
% figure; hold on; plot(xsm, songseg_smD, 'k');
% plot(1/fs:1/fs:(1/fs)*length(songseg_sm), songseg_sm, 'r');

SegmentsExtract(i).songseg_sm = single(songseg_sm);
SegmentsExtract(i).songseg_tOnOff = [t_songseg(1) t_songseg(end)];


    end
    
    
    
    %% =============== FIGURE OUT POSITION OF MOTIF WITHIN ITS BOUT
    if collectWholeBoutPosition==1
        ind; % current syl posotion
        
        boutnum = find(wholebout_firstsyls<=ind & wholebout_lastsyls>=ind);
        if length(boutnum)==0
            SegmentsExtract(i).BOUT_boutnum = nan;
            SegmentsExtract(i).BOUT_posOfTokenInBout = nan;
            SegmentsExtract(i).BOUT_RendInBout = nan;
            
        else
            % assert(length(boutnum)==1, 'sdafasd');
            positionOfTokenInBout = ind - wholebout_firstsyls(boutnum) + 1;
            RendInBout = sum(BoutNumsAll == boutnum)+1;
            
            BoutNumsAll = [BoutNumsAll boutnum];
            PosInBoutAll = [PosInBoutAll positionOfTokenInBout];
            RendInBoutAll = [RendInBoutAll RendInBout];
            
            SegmentsExtract(i).BOUT_boutnum = boutnum;
            SegmentsExtract(i).BOUT_posOfTokenInBout = positionOfTokenInBout;
            SegmentsExtract(i).BOUT_RendInBout = RendInBout;
        end
    end
    
    
    % =================== GET ONSETS/OFFSETS OF ALL SYLS IN MOTIF, RELATIVE
    % TO TIMING OF EACH MOTIF.
    
    % get only those that are within bounds of this segment's data
    inds = AllOnsets>ontime & AllOnsets<offtime;
    ontimesWithinData = AllOnsets(inds);
    ontimesRelStartDat = ontimesWithinData - ontime; % relative to start of segment
    
    inds = AllOffsets>ontime & AllOffsets<offtime;
    offtimesWithinData = AllOffsets(inds);
    offtimesWithinData = offtimesWithinData - ontime; % relative to start of segment
    
    % --- compare onsets to acoustic data
    if (0)
        lt_figure; hold on;
        lt_plot(ontimesRelStartDat, 0, {'Color', 'g'});
        lt_plot(offtimesWithinData, 0, {'Color', 'r'});
        plot([1:length(songseg)]./fs, (songseg-mean(songseg)).^2, 'k');
        line([predur predur], ylim, 'Color', 'b');
    end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % ========= TIME OF SONG FILE FOR THIS SEGMENT
    if isfield(NeurDat.metaDat, 'numSamps')
        % then is my data (i.e. not RA data from mel/sam)
        globalOnsetTime=AllOnsets(tokenExtents(i)); % sec
        globalOnsetSamp=globalOnsetTime*fs;
        cumOnsSamps=cumsum([0 NeurDat.metaDat.numSamps]);
        songind=find((globalOnsetSamp-cumOnsSamps)>0, 1, 'last');
        %     disp(['songind = ' num2str(songind)]);
        songfname=NeurDat.metaDat(songind).filename;
        
        SegmentsExtract(i).song_filename=songfname;
        %     SegmentsExtract(i).song_datenum=lt_neural_fn2datenum(songfname);
        SegmentsExtract(i).song_datenum=NeurDat.metaDat(songind).song_datenum;
        SegmentsExtract(i).song_ind_in_batch=songind;
    end
    
    
    % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if isfield(NeurDat, 'spikes_cat')
        % then is my data
        SegmentsExtract(i).spk_Clust=spk_ClustTimes(:,1)';
        SegmentsExtract(i).spk_Times=(spk_ClustTimes(:,2)/1000)'-ontime;
    else
        % then is RA dat from sam/mel
        SegmentsExtract(i).spk_Times = spiketimes - ontime;
    end
    SegmentsExtract(i).global_ontime_motifInclFlank=ontime;
    SegmentsExtract(i).global_offtime_motifInclFlank=offtime;
    SegmentsExtract(i).matchlabel=matchlabs{i};
    SegmentsExtract(i).fs=fs;
    SegmentsExtract(i).global_startind_motifonly=startinds(i);
    SegmentsExtract(i).global_endind_motifonly=endinds(i);
    SegmentsExtract(i).global_tokenind_DatAlignedToOnsetOfThis=tokenExtents(i);
    
    actualmotifdur=(offtime - postdur) - (ontime + predur);
    SegmentsExtract(i).actualmotifdur=actualmotifdur;
    
    SegmentsExtract(i).sylOnTimes_RelDataOnset=ontimesRelStartDat;
    SegmentsExtract(i).sylOffTimes_RelDataOnset=offtimesWithinData;
    
end

if suppressout==0
    disp(['DONE! EXTRACTED ' num2str(length(SegmentsExtract)) ' segments']);
end


%% ============== remove any with too long gap dur

if RemoveIfTooLongGapDur ==1
    disp(['Removed ' num2str(length(segstoremove)) ' (gap dur too long)']);
    SegmentsExtract(segstoremove) = [];
end


%% ====== REMOVE WN TRIALS IF IS LEARNING AND IF WANT TO REMOVE
if LearnKeepOnlyBase ==1
    
    % is this learning?
    birdname = Params.birdname;
    exptname = Params.exptname;
    [islearning, LearnSummary, switchtime] = lt_neural_v2_QUICK_islearning(birdname, exptname, 1);
    if isfield(SegmentsExtract, 'song_datenum')
        
        
        if islearning ==1
            assert(~isempty(switchtime), 'asdfasdf')
            
            tvals = [SegmentsExtract.song_datenum];
            indstoremove = tvals > switchtime;
            
            SegmentsExtract(indstoremove) = [];
            
            disp([' === removed WN files from (learning) : ' birdname '-' exptname '-' regexpr_str ' - ' num2str(sum(indstoremove)) '/' num2str(length(indstoremove))]);
            
            
            % --- if all data removed, make it an empty structure
            if all(indstoremove)
                SegmentsExtract = struct;
            end
            
        end
    end
end

%% ===== make structure empty if there is no data. LT added 8/18/17 (problem is length is 1 if no data but has fields)
% so force to have length 0;

if ~isfield(SegmentsExtract, 'fs');
    SegmentsExtract = [];
end


%% ##################### DIR, UNDIR? USE FILENAME SUFFIX

for j=1:length(SegmentsExtract)

   ind =  strfind(SegmentsExtract(j).song_filename, '_DIR_');
   
   if isempty(ind)
       SegmentsExtract(j).DirSong=0;
   elseif length(ind)==1
       % -- make sure is at expected location
       SegmentsExtract(j).DirSong=1;
   else
       asdfsdaf;
   end
end

%% =============== remove directed song if desired
if extractDirSong==0
    if ~isempty(SegmentsExtract)
    % -- remove
    SegmentsExtract([SegmentsExtract.DirSong]==1) = [];
    end
end

%% ==== DEBUG - CHECK HIT DETECTION
if suppressout==0
    if ~isempty(HitSyls_TEMP) | ~isempty(MisSyls_TEMP)
        lt_figure; hold on;
        lt_subplot(1,2,1); hold on; title('hits');
        for j=1:length(HitSyls_TEMP)
            plot(HitSyls_TEMP{j}, ':', 'Color', [rand rand rand]);
        end
        lt_subplot(1,2,2); hold on; title('misses');
        for j=1:length(MisSyls_TEMP)
            plot(MisSyls_TEMP{j}, ':', 'Color', [rand rand rand]);
        end
        
        lt_subtitle(regexpr_str);
    end
end

