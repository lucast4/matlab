%% ========== DATABASE OF BOS EXPERIMENTS
% --- ONE expt for each BOS presentation (at a specific depth)
clear all; close all;
ind = 0;


% % ========= NOTE, THIS IS A SMALL BATCH JUST FOR TESTING.
% ind = ind+1;
% SummaryBOS.expt(ind).dirname = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/BOS_110917Night';
% SummaryBOS.expt(ind).batchname = 'Batch2353to2357';
% SummaryBOS.expt(ind).channels = [9 14 14];
% SummaryBOS.expt(ind).clusters = [1 1 2];
% SummaryBOS.expt(ind).bregions= {'test', 'test', 'test'};

% ================================== 11/9 - Night 0
% ind = ind+1; OLD VERSION. NOW EXCLUDING CHAN 8 SINCE IS NOT SONG MOD
% (DURING SINGING), SO NOT IN LMAN.
% SummaryBOS.expt(ind).dirname = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/BOS_110917Night';
% SummaryBOS.expt(ind).batchname = 'Batch2327to2349good';
% SummaryBOS.expt(ind).channels = [8 9 9 14 21]; % note: 21 maybe coule be 2 clusteres?
% SummaryBOS.expt(ind).clusters = [1 1 2 1 1]; 
% SummaryBOS.expt(ind).isSU = [0 0 0 0 0];
% SummaryBOS.expt(ind).bregions= {'LMAN', 'LMAN', 'LMAN', 'LMAN', 'RA'};
ind = ind+1;
SummaryBOS.expt(ind).dirname = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/BOS_110917Night';
SummaryBOS.expt(ind).batchname = 'Batch2327to2349good';
SummaryBOS.expt(ind).channels = [9 9 14 21]; % note: 21 maybe coule be 2 clusteres?
SummaryBOS.expt(ind).clusters = [1 2 1 1]; 
SummaryBOS.expt(ind).isSU = [0 0 0 0];
SummaryBOS.expt(ind).bregions= {'LMAN', 'LMAN', 'LMAN', 'RA'};

ind = ind+1;
SummaryBOS.expt(ind).dirname = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/BOS_110917Night';
SummaryBOS.expt(ind).batchname = 'Batch0004to0032';
SummaryBOS.expt(ind).channels = [9 14 14 17 21 21]; % 
SummaryBOS.expt(ind).clusters = [1 1 2 1 1 2]; 
SummaryBOS.expt(ind).isSU = [0 0 1 0 0 1];
SummaryBOS.expt(ind).bregions= {'LMAN', 'LMAN', 'LMAN', 'RA', 'RA', 'RA'};


ind = ind+1;
SummaryBOS.expt(ind).dirname = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/BOS_110917Night';
SummaryBOS.expt(ind).batchname = 'Batch2353to2356';
SummaryBOS.expt(ind).channels = [9 14 14 21]; % 
SummaryBOS.expt(ind).clusters = [1 1 2 1]; 
SummaryBOS.expt(ind).isSU = [0 0 0 0];
SummaryBOS.expt(ind).bregions= {'LMAN', 'LMAN', 'LMAN', 'RA'};


% ================================== 11/10 - Night 1
ind = ind+1;
SummaryBOS.expt(ind).dirname = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/BOS_111017Night';
SummaryBOS.expt(ind).batchname = 'Batch0223to0251';
SummaryBOS.expt(ind).channels = [14 21]; % 
SummaryBOS.expt(ind).clusters = [1 1]; 
SummaryBOS.expt(ind).isSU = [0 0];
SummaryBOS.expt(ind).bregions= {'LMAN', 'RA'};


ind = ind+1;
SummaryBOS.expt(ind).dirname = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/BOS_111017Night';
SummaryBOS.expt(ind).batchname = 'Batch0158to0213';
SummaryBOS.expt(ind).channels = [14 21]; % 
SummaryBOS.expt(ind).clusters = [1 1]; 
SummaryBOS.expt(ind).isSU = [0 0];
SummaryBOS.expt(ind).bregions= {'LMAN', 'RA'};


% ================================== 11/11 - Night 2
ind = ind+1;
SummaryBOS.expt(ind).dirname = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/BOS_111117Night';
SummaryBOS.expt(ind).batchname = 'Batch0202to0226';
SummaryBOS.expt(ind).channels = [14 21 21]; % 
SummaryBOS.expt(ind).clusters = [1 1 2]; 
SummaryBOS.expt(ind).isSU = [0 0 0];
SummaryBOS.expt(ind).bregions= {'LMAN', 'RA', 'RA'};



% ####################### AUTO EXTRACTION OF BIRD NAME
for i=1:length(SummaryBOS.expt)
   tmp1 = strfind(SummaryBOS.expt(i).dirname, 'birds/');
   tmp2 = strfind(SummaryBOS.expt(i).dirname, '/NEURAL');
   bname = SummaryBOS.expt(i).dirname(tmp1+6:tmp2-1);
   SummaryBOS.expt(i).birdname = bname;
end

% ======== sanity check (entered one thing for each unit..)
for i=1:length(SummaryBOS.expt)
    assert(length(unique([length(SummaryBOS.expt(i).channels); ...
        length(SummaryBOS.expt(i).clusters); ...
        length(SummaryBOS.expt(i).isSU); ...
        length(SummaryBOS.expt(i).bregions); ...
        ]))==1, 'did not enter correcly, all shoudl have same number units');
end

%% ##################### [MAIN EXTRACTION]
close all; 

expttoget = 1;
assert(length(SummaryBOS.expt) ==1, 'currently only sure this works well if only have one expt loaded..');
lt_neural_BOS_Extraction(SummaryBOS, expttoget)


% ============= NOTE: following this, cheeck using wave_clus.

%% ##############################################################
%% ##############################################################
%% ============ 4) PLOT EACH SONG FILE OVERLAID WITH EXTRACTED SPIKES
i = 1;
fs = 30000;
batchname = SummaryBOS.expt(i).batchname;
chanstoget = 14;

% ======================= GO THRU EACH FILENAME
filenames = textread(batchname, '%s');
figcount=1;
subplotrows=3;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for j=1:length(filenames)
    fname = filenames{j};
    
    % ============= load and extract each channel
    % ---- 1) raw neural
    [amplifier_data,board_dig_in_data,frequency_parameters, ...
        board_adc_data, board_adc_channels, amplifier_channels, ...
        board_dig_in_channels, t_amplifier] = pj_readIntanNoGui(fname);
        
    [a, b] = fileparts(fname);
    
    % === go thru each channle
    for jj=1:length(chanstoget)
        chan = chanstoget(jj);
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title([fname '-chan' num2str(chan)]);
        
        neurdat = load([batchname '-Chan' num2str(chan) '/' b '.ch' num2str(chan) '.mat']);
        neurdatold = amplifier_data([amplifier_channels.chip_channel] == chan,:);
        assert(all(neurdat.data==single(neurdatold)));
        
        % --- filter neural dat
        neurdat = lt_neural_filter(double(neurdat.data), fs, 1);
        
        % ------ EXTRACT SPIKES
%         spksold = load([batchname '-Chan' num2str(chan) '/' b '.ch' num2str(chan) '_spikes.mat']);
        spksnew = load([batchname '-Chan' num2str(chan) '/allsongs_spikes.mat']);
        timesdat = load([batchname '-Chan' num2str(chan) '/times_allsongs.mat']);
%         spksnew.filenuminbatch_all==j;
        assert(size(spksnew.index,2) == size(timesdat.cluster_class,1));
        
        clustclassthis = timesdat.cluster_class(spksnew.filenuminbatch_all==j, :);
        
        % ================================ PLOT STUFF
        t = 1:length(neurdat);
        t = (t./fs)*1000;
        plot(t, neurdat, 'k');
        numclust = max(clustclassthis(:,1));
        plotcols = lt_make_plot_colors(numclust, 0, 0);
        for jjj=1:numclust % go thru all clusters (ignore 0, which is noise)
            spkthis = clustclassthis(clustclassthis(:,1)==jjj,2);
            lt_plot(spkthis, jjj, 'o', 'Color', plotcols{jjj});
        end
    end
end




%% ############################ [ANALYSES]
%% ========= PLOT OVERVIEW OF DATA
i = 1;

ls(SummaryBOS.expt(i).dirname, '-l');
disp('LIST OF ALL BOS SONGS');


%% ======================= [EXTRACTION]
% EXTRACT:
% 1) onsets/offsets
% 2) spike times
close all;

throwOutClippedSyls = 1; % do not collect rise/fall times for any syls that are clipped (edge of data)
ignoreIfFileEndsDuringSong=1; % if 1, then will make sure that for song extractio, will ignore any song that ends
% during file (assumes will be duplicated in next file...). Will not
% extract a song that ends a file UNLESS it ends at least 200ms (maxgapdur)
% befoer the end of the file.
ignoreIfFileStartsBeforeTimeZero=1; % if 1, then useful if a file starts and includes
% a song in pre-zero buffer. this occurs if the last file ends before
% finishing the song, and the beginning of a song is carried over to buffer
% on next file (e.g. becuase 1min time max for song file is done).

disp('NOTE: extractions assume that you have some buffer time (reasonably a few seconds)')


% ======= PARAMS
% --- channel for dig pulses and audio
chantype = 'ana'; % dig or ana
pulsechan = 1; % 0, 1, ... % for pulse\
%     threshold = 1; % voltage, or rise and fall detection.
% --- channel for song.
audchan = 0;

plotON = 0; % then plots each file, raw audio, syl onsets, etc.

% ---- for figure out motifs/bouts.
gapdurall = [];
maxgapdur = 0.2;
minedgetime = 2; % seconds of silence at edges - otherwise aborts.

% ===== RUN
for i=1:length(SummaryBOS.expt)
    
    % ==========================
    batchfile = [SummaryBOS.expt(i).dirname '/' SummaryBOS.expt(i).batchname];
    
    
    % ####################### EXTRACT SPIKETIMES FOR EACH CHANNEL
    % - Gets spktimes across all songs. later on will use this stuff to
    % extract song by song spikes.
    chans_spk = SummaryBOS.expt(i).channels;
    clust_spk = SummaryBOS.expt(i).clusters;
    Spkholder = struct; % NOTE: will contain all clusters for a given channel. this is fixed below.
    fs_all = [];
    for j=1:length(chans_spk)
        chanthis = chans_spk(j);
        
        tmp1 = load([SummaryBOS.expt(i).dirname '/' SummaryBOS.expt(i).batchname '-Chan' num2str(chanthis) '/times_allsongs.mat']);
        tmp2 = load([SummaryBOS.expt(i).dirname '/' SummaryBOS.expt(i).batchname '-Chan' num2str(chanthis) '/allsongs_spikes.mat']);
        Spkholder.unitnum(j).filenuminbatch_all = tmp2.filenuminbatch_all;
        Spkholder.unitnum(j).cluster_class = tmp1.cluster_class;
        fs_all = tmp1.par.sr;
    end
    
    
    
    % ############################# GO THRU ALL FILES
    % --- 1) extract timestamps of digital signals indicating bos playback
    fid = fopen(batchfile);
    fline = fgetl(fid);
    
    % --------- for debugging, plotting all pulse data
    if plotON==1
        figcount=1;
        subplotrows=6;
        subplotcols=1;
        fignums_alreadyused=[];
        hfigs=[];
    end
    
    % ==================== COLLECT DATA
    AllTotalSamps = [];
    Allt_amplifier = {};
    AllOnsets_samp = {}; % each index is one song.
    AllOffsets_samp = {};
    AllSpktimesByUnits = {};
    AllSongOnsets_inds = {};
    AllSongOffsets_inds = {};
    AllSongOnsets_ActualTime = {};
    
    songnum = 1;
    % ================= RUN
    while ischar(fline)
        
        % ---------- load file
        [amplifier_data,board_dig_in_data,frequency_parameters, ...
            board_adc_data, board_adc_channels, amplifier_channels, ...
            board_dig_in_channels, t_amplifier] = pj_readIntanNoGui([SummaryBOS.expt(i).dirname '/' fline]);
        fs = frequency_parameters.amplifier_sample_rate;
        fs_all = [fs_all; fs];
        
        flanktime = -t_amplifier(1);
        
        % ================== extract syl pulses
        if strcmp(chantype, 'dig')
            ind = [board_dig_in_channels.chip_channel]==pulsechan;
            pulsedat = board_dig_in_data(ind, :);
        elseif strcmp(chantype, 'ana')
            ind = [board_adc_channels.chip_channel]==pulsechan;
            pulsedat = board_adc_data(ind, :);
        end
        
        % ---------------------
        if plotON==1
            % --- plot pulse dat
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title(fline);
            xlabel('file time (sec');
            ylabel('V');
            plot(t_amplifier, pulsedat);
            % line(xlim, [threshold threshold], 'Color','r');
        end
        
        
        % ================================ IS THERE BOS IN THIS FILE? HOW MANY?
        % GET ONSETS/OFFSETS
        %         C = midcross(pulsedat, 'Tolerance', 15);
        [~, LT, UT, LL, UL] = risetime(pulsedat, 'PercentReferenceLevels', [35 65]);
        rtimes = round(mean([LT; UT],1));
        
        [~, LT, UT] = falltime(pulsedat, 'PercentReferenceLevels', [35 65]);
        ftimes = round(mean([LT; UT],1));
        
        if plotON==1
            line([t_amplifier(1) t_amplifier(end)], [LL LL]);
            line([t_amplifier(1) t_amplifier(end)], [UL UL]);
            plot(t_amplifier(rtimes), 1, 'ob');
            plot(t_amplifier(ftimes), 1, 'sr');
        end
        
        % ============== if syl is on edge of data, then will have only
        % rise or fall.
        if throwOutClippedSyls==1
            if ftimes(1)<rtimes(1)
                ftimes = ftimes(2:end);
            end
            if rtimes(end)>ftimes(end)
                rtimes=rtimes(1:end-1);
            end            
        end
        
        % ====================== SANITY CHECKS ABOUT RISE AND FALL TIMES
        assert(length(rtimes) == length(ftimes), 'syllable clipped?');
        assert(all((ftimes - rtimes)./fs > 0), 'all fall shoudl be after rise');
        assert(max(ftimes-rtimes)/fs<0.5, 'all syls should be <500ms...');
        assert(all((ftimes - rtimes)./fs <0.5), 'all syls shorter than 500ms');
        
        
        % ====================== if short syls, then are actually markers
        % for BOS it
        if min(ftimes-rtimes)/fs<0.002
            keyboard
        end
        
        
        % ===== CURRENTLY ASSUMING THAT NO SONGS WRAP AROUND FILES i.e. for a
        % given song, will be eintrely within one file
        if ignoreIfFileStartsBeforeTimeZero==0
        assert(rtimes(1)./fs>=flanktime, 'this song wraps around multiple files?');
        end
        if ignoreIfFileEndsDuringSong==0
        assert((t_amplifier(end)-t_amplifier(1) - (ftimes(end)./fs))>minedgetime, 'file ends too early?');
        % note: divide flanktime by 2 becuase sometimes I end file early...
        end
        
        % ===================== SANITY CHECK - MAKE SURE AUDIO POWER IS
        % HIGHER DURING SYLS (I.E. AUDIO ACTUALYL PLAYED)
        songdat = board_adc_data([board_adc_channels.chip_channel]==audchan, :);
        songdat_centered = songdat-mean(songdat);
        assert(mean(songdat_centered(rtimes(1):ftimes(end)).^2) > 5*mean(songdat_centered([1:rtimes(1) ftimes(end):length(rtimes)]).^2), 'is there no song played?');

        
        % ====================== SPIKES
        % --- go thru all units, for each unit extracrt spikes correspnding
        % to this song file
        numunits = length(Spkholder.unitnum);
        SpksAllUnits = cell(1,numunits);
        for nn=1:numunits
            clustthis = SummaryBOS.expt(i).clusters(nn);
            indstmp = Spkholder.unitnum(nn).filenuminbatch_all'==songnum ...
                & Spkholder.unitnum(nn).cluster_class(:,1)==clustthis;
            % make sure is correct song and cluster.
            
            spkthis = Spkholder.unitnum(nn).cluster_class(indstmp, 2);
            assert(all(diff(spkthis))>0, 'then not in order, problem with extraction ...');
            SpksAllUnits{nn} = spkthis;
            
            if plotON==1
                % -- plot this unit
                plot(spkthis/1000-flanktime, 1.5+nn*0.5/numunits, 'xk');
                lt_plot_text(spkthis(end)/1000-flanktime, 1.5+nn*0.5/numunits, ...
                    ['unit' num2str(nn)]);
            end
        end
        
        % ====================== extract sound
        if plotON==1
            songdat = board_adc_data([board_adc_channels.chip_channel]==audchan, :);
            %     t= [1:length(songdat)]./fs;
            %     lt_plot_spectrogram(songdat, fs, 0, 0);
            %     plot(songdat, 'm');
            plot(t_amplifier, (songdat-median(songdat)).^2);
            axis tight
        end
        
        
        % ======================== what is actual time of this file?
        [dtnum dtstring] = ...
            lt_neural_fn2datenum([SummaryBOS.expt(i).dirname '/' fline]); % onset of file
        
        
        
        %% ############################3 FIGURE OUT START/STOP OF SONG BOUTS
        onsets = rtimes;
        offsets = ftimes;
        totalsamps = length(pulsedat);
%         fs = fs;
        
        gapdurs = [onsets totalsamps] - [1 offsets]; % append with start and end of song
        gapdurs = gapdurs./fs;
        gapdurall = [gapdurall gapdurs];
        
        
        % =================== FIND LONG DURATIONS --> EDGES OF FILES
        edgeinds = find(gapdurs>maxgapdur); % correspond to onsets.
        
        if ignoreIfFileEndsDuringSong==0
            % then there shoudld be long gap between last syl and the end
            % of the song. otherwise don't care.
            assert(edgeinds(end) == length(onsets)+1, 'last extracted gap should always corerpond to end of song...');
        end
        songonsets = edgeinds(1:end-1); % throw out end as this is end of song.
        songoffsets = edgeinds(2:end)-1; % assume the offset right before each onset is an end.
        
        
        if plotON==1
            lt_plot(t_amplifier(rtimes(songonsets)), 2, {'Color', 'b'});
            lt_plot(t_amplifier(ftimes(songoffsets)), 2, {'Color', 'r'});
        end
        
        % ---
        if (0)
            lt_figure; hold on;
            title('gap durations');
            lt_plot_histogram(gapdurs, 0:0.01:1);
        end
        
        
        % ============= actual time for start of song
        songonsets_actualtime = dtnum+(t_amplifier(rtimes(songonsets)))./(60*60*24); % convert from seconds to days
        
        
        %% ############################### SAVE OUTPUT FOR THIS FILE
        AllTotalSamps = [AllTotalSamps; length(pulsedat)];
        Allt_amplifier = [Allt_amplifier; t_amplifier];
        AllOnsets_samp = [AllOnsets_samp; rtimes]; % each index is one song.
        AllOffsets_samp = [AllOffsets_samp; ftimes];
        AllSpktimesByUnits = [AllSpktimesByUnits; SpksAllUnits];
            AllSongOnsets_inds = [AllSongOnsets_inds; songonsets];
        AllSongOffsets_inds = [AllSongOffsets_inds; songoffsets];
        AllSongOnsets_ActualTime = [AllSongOnsets_ActualTime; songonsets_actualtime];
        
        
        % ========= load next song
        fline = fgetl(fid);
        songnum = songnum+1;
        
    end
    assert(size(AllSpktimesByUnits,1) == songnum-1);
    assert(songnum == max(Spkholder.unitnum(1).filenuminbatch_all)+1, 'make sure while looop counted correct number of songfiles');
    fs = unique(fs_all);
    assert(length(fs)==1, 'diff fs for diff files?');
    
    
    % ################## SAVE OUTPUT FOR THIS EXPERIMENT
    SummaryBOS.expt(i).DAT.AllTotalSamps = AllTotalSamps;
    SummaryBOS.expt(i).DAT.Allt_amplifier = Allt_amplifier;
    SummaryBOS.expt(i).DAT.AllOnsets_samp = AllOnsets_samp;
    SummaryBOS.expt(i).DAT.AllOffsets_samp = AllOffsets_samp;
    SummaryBOS.expt(i).DAT.AllSpktimesByUnits = AllSpktimesByUnits;
    SummaryBOS.expt(i).DAT.AllSongOnsets_inds = AllSongOnsets_inds;
    SummaryBOS.expt(i).DAT.AllSongOffsets_inds = AllSongOffsets_inds;
    SummaryBOS.expt(i).DAT.AllSongOnsets_ActualTime = AllSongOnsets_ActualTime;
    SummaryBOS.expt(i).fs = fs;
end

lt_figure; hold on;
title('gap durations (and line used as threshold for getting bouts)');
lt_plot_histogram(gapdurall, 0:0.01:1);
line([maxgapdur maxgapdur], ylim);




%% ===================== [EXTRACT] time relative data onset

for i=1:length(SummaryBOS.expt)
    
    if isfield(SummaryBOS.expt(i).DAT, 'Allt_fromDatOnset')
        SummaryBOS.expt(i).DAT = rmfield(SummaryBOS.expt(i).DAT, 'Allt_fromDatOnset');
    end
    
    numsongs = length(SummaryBOS.expt(i).DAT.AllOffsets_samp);
    
    for ss=1:numsongs
        
        % ========== collect all song bouts
        t = 1:SummaryBOS.expt(i).DAT.AllTotalSamps(ss);
        fs = SummaryBOS.expt(i).fs;
        
        t = t./fs;
        SummaryBOS.expt(i).DAT.Allt_fromDatOnset{ss} = t;
        
    end
end


%% ====================== [EXTRACTION] - SONG BY SONG BOUT
clear PARAMS
flanktotake = [-2 2]; % sec from onset and offset.
minsongdur = 8; % seconds

% ===============

for i=1:length(SummaryBOS.expt)

numsongs = length(SummaryBOS.expt(i).DAT.AllOffsets_samp);

% ----------- to collect over all songs
DAT_songextract = struct;
SpkTime_RelSongOnset = {};
SongDur = [];
SongOnset_ActualDatenum = [];
SylOnsets = {};
SylOffsets = {};
TEdges = [];

for ss=1:numsongs
    
    % ===== stats for this file
    onsets_samps = SummaryBOS.expt(i).DAT.AllOnsets_samp{ss};
    offsets_samps = SummaryBOS.expt(i).DAT.AllOffsets_samp{ss};
    
    t = SummaryBOS.expt(i).DAT.Allt_fromDatOnset{ss};
    SpkTimes = SummaryBOS.expt(i).DAT.AllSpktimesByUnits(ss,:);
    nunits = length(SpkTimes);
    
    % ========== collect all song bouts
    nbouts = length(SummaryBOS.expt(i).DAT.AllSongOnsets_inds{ss});
    for nn=1:nbouts
        
        on_inds = SummaryBOS.expt(i).DAT.AllSongOnsets_inds{ss}(nn);
        off_inds = SummaryBOS.expt(i).DAT.AllSongOffsets_inds{ss}(nn);
        dtnum = SummaryBOS.expt(i).DAT.AllSongOnsets_ActualTime{ss}(nn);

        % ---------- COLLECT SPIKES ALIGNED TO INDS
        % - time of onset/offset(sec);
        tonset_sec = t(onsets_samps(on_inds));
        toffset_sec = t(offsets_samps(off_inds));
        
        
        % ============== FILTER - IGNORE IF TOO SHORT
        songdur = toffset_sec - tonset_sec;
        if songdur < minsongdur
            disp('SKIPPING SONG RENDITION - TOO SHORT');
        continue
        end
        
        % - time endpoints
        tedges = [flanktotake(1) songdur+flanktotake(2)];
            
        % -- get spikes flanking this - recenter to be relative to onset of
        % song
        % ----- GO THRU ALL UNITS
        spktimes_tmp = cell(1,nunits);
        for j=1:nunits
            
            spkthis = SpkTimes{j}./1000; % convert to sec.
            spkthis = spkthis(spkthis>tonset_sec+flanktotake(1) ...
                & spkthis<toffset_sec+flanktotake(2)); % get desired window
            % -- get time relative to onset of song
            spkthis_relsongonset = spkthis - tonset_sec;
            spktimes_tmp{j} = spkthis_relsongonset;
        end
        
        
        % ============= GET ONSETS and offsets
        allonsets = t(onsets_samps(on_inds:off_inds));
        alloffsets = t(offsets_samps(on_inds:off_inds));
        % -- rel to onset of first syl
        allonsets = allonsets-tonset_sec;
        alloffsets = alloffsets-tonset_sec;
        
        % ============ GET ACTUAL TIMES
        
        
        % #####################3 COLLECT THINGS ABOUT THIS SONG
        SylOnsets = [SylOnsets; allonsets];
        SylOffsets = [SylOffsets; alloffsets];
        SongDur = [SongDur; toffset_sec-tonset_sec];
        SpkTime_RelSongOnset = [SpkTime_RelSongOnset; spktimes_tmp];
        TEdges = [TEdges; tedges];
        SongOnset_ActualDatenum = [SongOnset_ActualDatenum; dtnum];
    end
end

SummaryBOS.expt(i).DAT_bysongrend.TEdges = TEdges;
SummaryBOS.expt(i).DAT_bysongrend.SylOnsets = SylOnsets;
SummaryBOS.expt(i).DAT_bysongrend.SylOffsets = SylOffsets;
SummaryBOS.expt(i).DAT_bysongrend.SongDur = SongDur;
SummaryBOS.expt(i).DAT_bysongrend.SpkTime_RelSongOnset = SpkTime_RelSongOnset;
SummaryBOS.expt(i).DAT_bysongrend.SongOnset_ActualDatenum = SongOnset_ActualDatenum;
PARAMS.flanktotake = flanktotake;



% ==================== CONVERT TO FORMAT OF SEGEXTRACT TO USE OLD CODE
for j=1:nunits
    % -- go thru all trials
    spks = SummaryBOS.expt(i).DAT_bysongrend.SpkTime_RelSongOnset(:, j);
    ntrials = size(SummaryBOS.expt(i).DAT_bysongrend.SpkTime_RelSongOnset,1);
    for k=1:ntrials
        % -- spk times
        SummaryBOS.expt(i).DAT_bysongrend.SegextractFormat.unitnum(j).segextract(k).spk_Times = ...
            (spks{k} - PARAMS.flanktotake(1))'; % convert to start from 0
        % -- trial durations
        datdur = SummaryBOS.expt(i).DAT_bysongrend.SongDur(k) + (PARAMS.flanktotake(2) - PARAMS.flanktotake(1));
        SummaryBOS.expt(i).DAT_bysongrend.SegextractFormat.unitnum(j).segextract(k).datdur = datdur;
    end
end

end


%% ===================== [OPTIONAL] PLOTS OF BOS FILE RELATED STUFF.
% ######################## PLOTS TO HELP FIGURE OUT WHAT BOS ARE USED.
if (1)
    exptnum =3;
    i=exptnum;
    
    % ============== TO VISUALIZE THE ACTUAL TIMES OF ALL BOS PRESENTATIONS
    lt_figure; hold on;
    subplot(2,1,1);
    title('actual times, BOS onsets');
    bostimes = SummaryBOS.expt(i).DAT_bysongrend.SongOnset_ActualDatenum;
    plot(bostimes, '-ok');
    xlabel('BOS #');
    datetick('y', 'HH:MM:SS');
    
    subplot(2,1,2);
    bostimes = (bostimes-bostimes(1))*(24*60*60); % convert to sec
    plot(bostimes, '-ok');
    ylabel('sec from first BOS');
    
    
    % ============== TO VISUALIZE STRUCTURE OF BOS OVER TRIALS... (SYL, GAPS)
    if (1)
        overlayspks = 0;
        % ------------
        numrends = length(SummaryBOS.expt(i).DAT_bysongrend.SylOnsets);
        lt_figure; hold on;
        ylabel('trial (1, 2, ...)');
        
        
        for j=1:numrends
            ons = SummaryBOS.expt(i).DAT_bysongrend.SylOnsets{j};
            offs = SummaryBOS.expt(i).DAT_bysongrend.SylOffsets{j};
            
            lt_neural_QUICK_PlotSylPatches(ons, offs, j,1);
            
            % --- plot subset of spikes
            if overlayspks==1
                spks = SummaryBOS.expt(i).DAT_bysongrend.SpkTime_RelSongOnset{j}(1:100);
                lt_neural_PLOT_rasterline(spks, j, 'r', 0);
            end
            
            try
                bostype = SummaryBOS.expt(i).DAT_bysongrend.BOStype(j);
                lt_plot_text(0, j+0.2, ['bos#' num2str(bostype)], 'm');
            catch err
            end
        end
    end
    
    
    % ================= VISUALIZE ACTUAL BOS FILES
    % ---------- PARAMS
    dirforBOS = '/bluejay5/egret_data/lucas/Test_Songs/BOS/'; % where .wav files are stored
    dirforBOS = '/run/user/1197/gvfs/smb-share:server=egret.cin.ucsf.edu,share=data/lucas/Test_Songs/BOS/';
    BOSfiles = {'pu69wh78_031117_102806.9552.cbin_DigOnsOff.wav', ...
        'pu69wh78_031117_102806.9552.cbin_DigOnsOff_REV.wav'}; % filenames for all BOS files; left = sound, right = ons/offset pulses
    BOScodes = []; % code used for each BOS file, can be empty
    % batchfile = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/BOS_111017Night/batchtest';
    
    % =============== plot all bos files
    datall = {};
    for j=1:length(BOSfiles)
        dat = audioread([dirforBOS BOSfiles{j}]);
        datall = [datall; dat];
        
    end
    
    figcount=1;
    subplotrows=5;
    subplotcols=1;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    % === plot
    for j=1:length(datall)
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(BOSfiles{j});
        plot(datall{j}(:,2));
        ylim([-1 2]);
    end
end

lt_figure; lt_plot_text(0,1, 'not general! fix where looks for BOS files');

%% ===================== [EXTRACTION] Determine type of BOS on each trial.
% ################ CLASSIFY SONG BASED ON CROSSCORRELATION WITH
% TEMPLATES FROM ACTUAL FILES
% exptnum = 1;
% ---------- PARAMS
clear BOSfiles
clear dirforBOS
clear BOSnames
clear BOSbirdname

% ---------------- GENERAL DIRECTORY WHERE ALL BIRDS' BOS ARE STORED.
% dirforBOS = '/bluejay5/egret_data/lucas/Test_Songs/BOS/'; % where .wav files are stored
PARAMS.dirforBOS = '/bluejay5/lucas/analyses/BOS/Songs/BOS/'; % where .wav files are stored
% dirforBOS = '/run/user/1197/gvfs/smb-share:server=egret.cin.ucsf.edu,share=data/lucas/Test_Songs/BOS/';

% ----------------- NAME OF BIRD - must correspond to the directory BOS
% files are saved in.
PARAMS.BOSbirdname{1} = 'pu69wh78';
PARAMS.BOSbirdname{2} = 'wh72pk12';

% ---------------- SONG FILES. ORDER MUST CORRESPOND TO ORDER OF BOSnames
PARAMS.BOSfiles{1} = {'pu69wh78_031117_102806.9552.cbin_DigOnsOff.wav', ...
    'pu69wh78_031117_102806.9552.cbin_DigOnsOff_REV.wav'}; % filenames for all BOS files; left = sound, right = ons/offset pulses
PARAMS.BOSfiles{2} = {};

% ------------------ FOR LABELING - order must correspond to BOSfiles
PARAMS.BOSnames{1} = {'fwd', 'rev'}; % code used for each BOS file, can be empty
PARAMS.BOSnames{2} = {};
% batchfile = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/BOS_111017Night/batchtest';



for i=1:length(SummaryBOS.expt)
    
    % ----------- COLLECT ACTUAL BOS FILES
    % --- GET BOS PARAMS FOR THIS BIRD
    birdthis = SummaryBOS.expt(i).birdname;
    
    indbird = strcmp(PARAMS.BOSbirdname, birdthis);
    dirforBOS = [PARAMS.dirforBOS birdthis '/'];
    BOSfiles = PARAMS.BOSfiles{indbird};
    
    
    % =============== plot all bos files
    datall = {};
    for j=1:length(BOSfiles)
        dat = audioread([dirforBOS BOSfiles{j}]);
        datall = [datall; dat];
    end
    
    % ============== get rise and fall times for each actual songdat
    datall_risetimes = {};
    datall_falltimes = {};
    for j=1:length(datall)
        
        [~, b] = risetime(datall{j}(:,2));
        [~, c] = falltime(datall{j}(:,2));
        
        datall_risetimes{j} = b;
        datall_falltimes{j} = c;
    end
    
    % ============== GO THRU ALL SONG RENDS AND COMPARE TO ALL BOS
    nrends = length(SummaryBOS.expt(i).DAT_bysongrend.SongDur);
    nBosSongs = length(datall);
    BosIndAll = [];
    for j=1:nrends
        ons = SummaryBOS.expt(i).DAT_bysongrend.SylOnsets{j};
        offs = SummaryBOS.expt(i).DAT_bysongrend.SylOffsets{j};
        
        % ===  get corr vs each candidate song
        rhoall = [];
        for jj=1:nBosSongs
            rho = corr(offs'-ons', datall_falltimes{jj}-datall_risetimes{jj});
            rhoall = [rhoall; rho];
        end
        bosindthis = find(rhoall>0.99);
        assert(length(bosindthis)==1, 'should be exactly one bos this corresponds to...');
        
        BosIndAll = [BosIndAll; bosindthis];
    end
    SummaryBOS.expt(i).DAT_bysongrend.BOStype = BosIndAll;
end

% PARAMS.dirforBOS = dirforBOS;
% PARAMS.BOSfiles = BOSfiles;
% PARAMS.BOSnames = BOSnames;

%% ===================== [LOAD LABEL BOS FILES]
if (0)
    dirforBOS = '/bluejay5/egret_data/lucas/Test_Songs/BOS/'; % where .wav files are stored
    % dirforBOS = '/run/user/1197/gvfs/smb-share:server=egret.cin.ucsf.edu,share=data/lucas/Test_Songs/BOS/';
    BOSfiles = {'pu69wh78_031117_102806.9552.cbin_DigOnsOff.wav', ...
        'pu69wh78_031117_102806.9552.cbin_DigOnsOff_REV.wav'}; % filenames for all BOS files; left = sound, right = ons/offset pulses
    BOSnames = {'fwd', 'rev'}; % code used for each BOS file, can be empty
end
try
    PARAMS = rmfield(PARAMS, 'BOSlabels');
catch err
end

for i=1:length(PARAMS.BOSfiles)
    nbos = length(PARAMS.BOSfiles{i});
    BOSlabels = {};
    for j=1:nbos
        
        fthis = [PARAMS.dirforBOS PARAMS.BOSbirdname{i} '/' PARAMS.BOSfiles{i}{j} '.labels'];
        tmp = load(fthis, '-mat');
        
        BOSlabels{j} = tmp.label;
    end
    PARAMS.BOSlabels{i} = BOSlabels;
end
%% ===================== [PLOT] SUMMARY FOR RESPONSE [BOSTYPES]
% EACH UNIT PLOT MEAN RESPONSE TO EACH BOS TYPE
close all;
i=1;
birdthis = SummaryBOS.expt(i).birdname;

% =============
nunits = size(SummaryBOS.expt(i).DAT_bysongrend.SpkTime_RelSongOnset,2);
% nsongs = size(SummaryBOS.expt(i).DAT_bysongrend.SpkTime_RelSongOnset, 1);

figcount=1;
subplotrows=4;
subplotcols = 1;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for nn=1:nunits
   
    chan = SummaryBOS.expt(i).channels(nn);
    bregion = SummaryBOS.expt(i).bregions{nn};
    spks = SummaryBOS.expt(i).DAT_bysongrend.SpkTime_RelSongOnset(:,nn);
    onsets = SummaryBOS.expt(i).DAT_bysongrend.SylOnsets;
    offsets = SummaryBOS.expt(i).DAT_bysongrend.SylOffsets;
    segextract = SummaryBOS.expt(i).DAT_bysongrend.SegextractFormat.unitnum(nn).segextract;
    bostype = SummaryBOS.expt(i).DAT_bysongrend.BOStype;
    BOSnames = PARAMS.BOSnames{strcmp(PARAMS.BOSbirdname, birdthis)};
    BOSlabels = PARAMS.BOSlabels{strcmp(PARAMS.BOSbirdname, birdthis)};
    
    % ======= collect smoothed fr for each trial
    segextract = lt_neural_SmoothFR(segextract, [], [], [], [], [], PARAMS.flanktotake(1));    
    frmat = [segextract.FRsmooth_rate_CommonTrialDur];
    x = segextract(1).FRsmooth_xbin_CommonTrialDur;
    
%     % ====== 2) rasters
%     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%     title('rasters [syls shaded]');
%     ylabel(['unit' num2str(nn)]);
%     for ss=1:nsongs
%        
%        % -- overlay onsets and offsets
%         ons = onsets{ss};
%         offs = offsets{ss};
%         lt_neural_QUICK_PlotSylPatches(ons, offs, ss);
%        lt_neural_PLOT_rasterline(spks{ss}, ss, 'r', 0);
% %           
% %                        X=[segextract(j).WNonset_sec  segextract(j).WNoffset_sec ...
% %                     segextract(j).WNoffset_sec  segextract(j).WNonset_sec];
% %                 Y=[-j-0.4 -j-0.4 -j+0.4 -j+0.4];
% 
%     end
    
    % ================= ONE PLOT FOR EACH BOS TYPE
    nbostypes = max(bostype);
    pcols_bos = lt_make_plot_colors(nbostypes, 0,0);
    
    for bb=1:nbostypes
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['unit' num2str(nn) ' [ch' num2str(chan) '-' bregion '] -' BOSnames{bb}]);
        hsplots = [hsplots hsplot];
        
        indsthis = bostype==bb;
        
        frthis = frmat(:, indsthis);
        
        frmean = mean(frthis,2);
        frsem = lt_sem(frthis');
        if size(frthis,2)==1
           plot(x, frthis, 'Color', pcols_bos{bb}); 
        else
        shadedErrorBar(x, frmean, frsem, {'Color', pcols_bos{bb}},1);
        end
        lt_plot_zeroline;
        
        % ======== overlay song
        ons_mat = cell2mat(onsets(indsthis));
        if ~all(max(ons_mat)-min(ons_mat)<0.001)
            tmp = max(max(ons_mat)-min(ons_mat));
            lt_plot_annotation(1, ['some onsets up to: ' num2str(tmp) ' sec diff (across tirals)...'], 'm');
        end
        ons = median(ons_mat,1);
        
        offs_mat = cell2mat(offsets(indsthis));
        if ~all(max(offs_mat)-min(offs_mat)<0.001)
            tmp = max(max(offs_mat)-min(offs_mat));
            lt_plot_annotation(1, ['some oiffests up to: ' num2str(tmp) ' sec diff (across tirals)...'], 'm');
        end
        offs = median(offs_mat,1);
        
        lt_neural_QUICK_PlotSylPatches(ons, offs, median(frmean), 1);
        
        % ----------- PLOT LABELS
        labels = BOSlabels{bb};
        assert(length(labels)==length(ons), 'wierd');
        for k=1:length(labels)
           lt_plot_text(ons(k)+0.01, median(frmean)+5, labels(k), 'b'); 
        end
    end
    
    
    
    % sanity check
    if (0)
       figure; hold on;
       % ----
       fr = segextract(1).FRsmooth_rate_CommonTrialDur;
       t = segextract(1).FRsmooth_xbin_CommonTrialDur;
%        t = t+SummaryBOS.expt(1).DAT_bysongrend.TEdges(1);
       plot(t, fr, '-k');
       spks = SummaryBOS.expt(1).DAT_bysongrend.SpkTime_RelSongOnset{1,nn};
       lt_neural_PLOT_rasterline(spks, 100, 'r', 0);
       
       % ----
       fr = segextract(1).FRsmooth_rate_CommonTrialDur;
       lt_neural_PLOT_rasterline(spks, 100, 'r', 0);
    end
    
end

linkaxes(hsplots, 'xy');



%% ===================== [PLOT] EXTRACT MOTIFS (FROM LABEL);




%% ====================== [PLOT] RASTERS OVER TRIALS, ONE FOR EACH UNIT
close all;
i = 2;

% =============
nunits = size(SummaryBOS.expt(i).DAT_bysongrend.SpkTime_RelSongOnset,2);
nsongs = size(SummaryBOS.expt(i).DAT_bysongrend.SpkTime_RelSongOnset, 1);

figcount=1;
subplotrows=4;
subplotcols = 1;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for nn=1:nunits
   
    
    spks = SummaryBOS.expt(i).DAT_bysongrend.SpkTime_RelSongOnset(:,nn);
    onsets = SummaryBOS.expt(i).DAT_bysongrend.SylOnsets;
    offsets = SummaryBOS.expt(i).DAT_bysongrend.SylOffsets;
    segextract = SummaryBOS.expt(i).DAT_bysongrend.SegextractFormat.unitnum(nn).segextract;
    
    % ====== 2) rasters
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('rasters [syls shaded]');
    ylabel(['unit' num2str(nn)]);
    for ss=1:nsongs
       
       % -- overlay onsets and offsets
        ons = onsets{ss};
        offs = offsets{ss};
        lt_neural_QUICK_PlotSylPatches(ons, offs, ss);
       lt_neural_PLOT_rasterline(spks{ss}, ss, 'r', 0);
%           
%                        X=[segextract(j).WNonset_sec  segextract(j).WNoffset_sec ...
%                     segextract(j).WNoffset_sec  segextract(j).WNonset_sec];
%                 Y=[-j-0.4 -j-0.4 -j+0.4 -j+0.4];

    end
    
    % ====== 3) mean fr
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('rasters [syls shaded]');
%     segextract = segextract';
    segextract = lt_neural_SmoothFR(segextract, [], [], [], [], [], PARAMS.flanktotake(1));    
    frmat = [segextract.FRsmooth_rate_CommonTrialDur];
    x = segextract(1).FRsmooth_xbin_CommonTrialDur;
    plot(x, frmat, 'Color', [0.7 0.7 0.7]);
    frmean = mean(frmat,2);
    frsem = lt_sem(frmat');
    shadedErrorBar(x, frmean, frsem, {'Color', 'k'},1);
    
    % sanity check
    if (0)
       figure; hold on;
       % ----
       fr = segextract(1).FRsmooth_rate_CommonTrialDur;
       t = segextract(1).FRsmooth_xbin_CommonTrialDur;
%        t = t+SummaryBOS.expt(1).DAT_bysongrend.TEdges(1);
       plot(t, fr, '-k');
       spks = SummaryBOS.expt(1).DAT_bysongrend.SpkTime_RelSongOnset{1,nn};
       lt_neural_PLOT_rasterline(spks, 100, 'r', 0);
       
       % ----
       fr = segextract(1).FRsmooth_rate_CommonTrialDur;
       lt_neural_PLOT_rasterline(spks, 100, 'r', 0);
    end
    
end

%% ====================== [OLD] EXTRACT UNCLUSTERED NEURAL DATA
% -- neural
chanstoplot = [9 14 17];


% --- 1) extract timestamps of digital signals indicating bos playback
fid = fopen(batchfile);

fline = fgetl(fid);

% --------- for debugging, plotting all pulse data
figcount=1;
subplotrows=6;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

while ischar(fline)
    
    % ---------- load file 
    [amplifier_data,board_dig_in_data,frequency_parameters, ...
        board_adc_data, board_adc_channels, amplifier_channels, ...
        board_dig_in_channels, t_amplifier] = pj_readIntanNoGui(fline);
    
    fs = frequency_parameters.amplifier_sample_rate;
    
    % ================== extract syl pulses
    if strcmp(chantype, 'dig')
        ind = [board_dig_in_channels.chip_channel]==pulsechan;
        pulsedat = board_dig_in_data(ind, :);
    elseif strcmp(chantype, 'ana')
        ind = [board_adc_channels.chip_channel]==pulsechan;
        pulsedat = board_adc_data(ind, :);
    end

    % ---------------------
    if (1)
       % --- plot pulse dat       
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    
plot(pulsedat);
line(xlim, [threshold threshold], 'Color','r');
    end
    
    
    % ================================ IS THERE BOS IN THIS FILE? HOW MANY?
    % GET ONSETS/OFFSETS
    risetime
    
    
    
    % ====================== extract sound    
    songdat = board_adc_data([board_adc_channels.chip_channel]==audchan, :);
%     t= [1:length(songdat)]./fs;
%     lt_plot_spectrogram(songdat, fs, 0, 0);
%     plot(songdat, 'm');
    plot((songdat-median(songdat)).^2);
    
    % ================ extract neural
    if isempty(chanstoplot)
        chans = [amplifier_channels.chip_channel];
    else
        chans = chanstoplot;
    end
    
    
    for i=chans
       neurdat = amplifier_data([amplifier_channels.chip_channel]==i, :);
       
       
        
    end
        
    
    
    fline = fgetl(fid);
end

%% ============================= EXTRACT PRE-CLUSTERED NEURAL DATA
% ---

% -- defaults
extractsound = 0;

NeurDat = struct;
numsampsAll = []; % to confirm that all data are aligned (samme num samps);
fsAll = [];
filenamesAll = {};

for i = 1:length(ChansToGet)
    cd(dirname);
    chan = ChansToGet(i);
    clust = ClustToGet(i);
    
    [SongDatTMP, NeurDatTMP, ParamsTMP] = lt_neural_ExtractDat(Batchname, chan, ...
        extractsound, clust);
    
    % - convert from ms to sec
    NeurDatTMP.spikes_cat.cluster_class(:,2) = ...
    NeurDatTMP.spikes_cat.cluster_class(:,2)./1000;
        
    
    % -
    NeurDat.metaDat = NeurDatTMP.metaDat;
    NeurDat.unit(i).spikes_cat = NeurDatTMP.spikes_cat;
    NeurDat.unit(i).chan = chan;
    NeurDat.unit(i).clust = clust;
    NeurDat.dirname = dirname;
    NeurDat.unit(i).BrainRegion = BrainRegion{i};
    
    % ====== sanity checks - all data are from same files
    numsampsAll = [numsampsAll; [NeurDatTMP.metaDat.numSamps]];
    fsAll = [fsAll; [NeurDatTMP.metaDat.fs]];
    filenamesAll = [filenamesAll; [NeurDatTMP.metaDat.filename]];
end

% === sanity checks
assert(size(unique(numsampsAll, 'rows'),1) ==1, 'problem');
assert(size(unique(fsAll, 'rows'),1)==1, 'problem');
assert(numel(unique(filenamesAll)) ==1, 'problem');


%% +++++++++++++++++++++++++++++++++++++++++ OLD STUFF 
%% ========== TO MAKE BOS, ETC.
%% given batch of BOS, plot responses for all channels
clear all; close all;

% ---------- PARAMS
dirforBOS = '/bluejay5/egret_data/lucas/Test_Songs/BOS/'; % where .wav files are stored
BOSfiles = {'pu69wh78_031117_102806.9552.cbin_DigOnsOff.wav', ...
    'pu69wh78_031117_102806.9552.cbin_DigOnsOff_REV.wav'}; % filenames for all BOS files; left = sound, right = ons/offset pulses
BOScodes = []; % code used for each BOS file, can be empty
batchfile = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/BOS_111017Night/batchtest';

% --- channel for dig pulses and audio
chantype = 'ana'; % dig or ana
pulsechan = 1; % 0, 1, ... % for pulse\
threshold = 1; % voltage, or rise and fall detection.

% --- channel for song.
audchan = 0;

% -- neural
chanstoplot = [9 14 17];


% --- 1) extract timestamps of digital signals indicating bos playback
fid = fopen(batchfile);

fline = fgetl(fid);

% --------- for debugging, plotting all pulse data
figcount=1;
subplotrows=6;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

while ischar(fline)
    
    % ---------- load file 
    [amplifier_data,board_dig_in_data,frequency_parameters, ...
        board_adc_data, board_adc_channels, amplifier_channels, ...
        board_dig_in_channels, t_amplifier] = pj_readIntanNoGui(fline);
    
    fs = frequency_parameters.amplifier_sample_rate;
    
    % ================== extract syl pulses
    if strcmp(chantype, 'dig')
        ind = [board_dig_in_channels.chip_channel]==pulsechan;
        pulsedat = board_dig_in_data(ind, :);
    elseif strcmp(chantype, 'ana')
        ind = [board_adc_channels.chip_channel]==pulsechan;
        pulsedat = board_adc_data(ind, :);
    end

    % ---------------------
    if (1)
       % --- plot pulse dat       
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    
plot(pulsedat);
line(xlim, [threshold threshold], 'Color','r');
    end
    
    
    % ================================ IS THERE BOS IN THIS FILE? HOW MANY?
    % GET ONSETS/OFFSETS
    risetime
    
    
    
    % ====================== extract sound    
    songdat = board_adc_data([board_adc_channels.chip_channel]==audchan, :);
%     t= [1:length(songdat)]./fs;
%     lt_plot_spectrogram(songdat, fs, 0, 0);

    
    % ================ extract neural
    if isempty(chanstoplot)
        chans = [amplifier_channels.chip_channel];
    else
        chans = chanstoplot;
    end
    
    
    for i=chans
       neurdat = amplifier_data([amplifier_channels.chip_channel]==i, :);
       
       
        
    end
        
    
    
    fline = fgetl(fid);
end



%% ######################################################################
%% ######################################################################
%%  TO MAKE BOS FILES (with digitual pulse for syl onset/offset

% ================== CAN CHOOSE TO PUT DIGITAL MARKER for each song (will
% correspond to number put in file name)
useRandID = 1; % number between 1 and 100, will put thi


% ================ pu69
% songfname = 'pu69wh78_031117_102806.9552.cbin'; % shoudl have notmat also
% OutDir = '/bluejay5/egret_data/lucas/Test_Songs/BOS';

% ================ wh72
songfname = '/bluejay0/bluejay2/lucas/birds/wh72pk12/112518_TetheredCMOS/wh72pk12_251118_111424.3033.cbin'; % use this as is better quality (louder)
% songfname = '/bluejay0/bluejay2/lucas/birds/wh72pk12/112818_RALMANLearn2/wh72pk12_281118_132419.445.cbin'; % this is longer though./..
OutDir = '/bluejay5/lucas/analyses/BOS/Songs/BOS/wh72pk12';


% ######################## convert to wavfile
% -- chan 1 (song)
[dat, fs] = ReadCbinFile(songfname);

dat = dat/max(dat);

% --- bandpass filter
filter_type='hanningfirff';
F_low  = 500;
F_high = 8000;
dat=bandpass(dat,fs,F_low,F_high,filter_type);
datbeforesmth = dat;


% -- chan 2 (onsets and offsets)
notdat = load([songfname '.not.mat']);

datdig = zeros(length(dat),1);

for i=1:length(notdat.onsets)
    ons = notdat.onsets(i);
    off = notdat.offsets(i);
    
    % convert from ms to samps
    
    ons = round(fs*(ons/1000));
    off = round(fs*(off/1000));
    
    datdig(ons:off) = 1;
end



% ######################## exponential falloff from sound onset and offset
rolltime = 0.7;
expsize = round(rolltime*fs);

% -- roll on
smth = exp(-((expsize:-1:1)-1)/(expsize/5));
dat(1:expsize) = dat(1:expsize).*smth';

% --- roll off
smth = exp(-((1:expsize)-1)/(expsize/5));
dat(end-expsize+1:end) = dat(end-expsize+1:end).*smth';

figure; hold on;
plot(datbeforesmth, 'k');
plot(dat, 'c');
plot(datdig, 'r');


% ############# combine song and dig signal
[tmpa, tmpb, tmpc] = fileparts(songfname);
if useRandID==1
    [datdig_ID, idthis] = lt_neural_BOS_StimID(datdig, fs);
    fnameout = [tmpb tmpc '_DigOnsOff_randID' num2str(idthis) '.wav'];
    fnameout = [OutDir '/' fnameout];
    audiowrite(fnameout, [dat datdig_ID], fs, 'BitsPerSample', 16)
else
    fnameout = [tmpb tmpc '_DigOnsOff.wav'];
    fnameout = [OutDir '/' fnameout];
    audiowrite(fnameout, [dat datdig], fs, 'BitsPerSample', 16)
end
% wavwrite([dat datdig], fs, 16, fnameout);

% ############# MAKE REVERSE VERSION
dat_rev = flipud(dat);
datdig_rev = flipud(datdig);

if useRandID==1
    [datdig_ID, idthis] = lt_neural_BOS_StimID(datdig_rev, fs);
    fnameout = [tmpb tmpc '_DigOnsOff_REV_randID' num2str(idthis) '.wav'];
    fnameout = [OutDir '/' fnameout];
%     fnameout = [tmpb tmpc '_DigOnsOff_randID' num2str(idthis) '.wav'];
    audiowrite(fnameout, [dat_rev datdig_ID], fs, 'BitsPerSample', 16)
else
%     fnameout = [tmpb tmpc '_DigOnsOff.wav'];
    fnameout = [tmpb tmpc '_DigOnsOff_REV.wav'];
    % fnameout = [songfname '_DigOnsOff_REV.wav'];
    fnameout = [OutDir '/' fnameout];
    audiowrite(fnameout, [dat_rev datdig_rev], fs, 'BitsPerSample', 16)
    % wavwrite([dat_rev datdig_rev], fs, 16, fnameout);
end

% ########### plot
figure; hold on;
subplot(211); hold on;
title('forward');
plot(datbeforesmth, 'k');
plot(dat, 'c');
plot(datdig, 'r');

subplot(212); hold on;
title('rev');
plot(dat_rev, 'c');
plot(datdig_rev, 'r');


%%   script for getting song names and syl positions, given only part of song name and the motif number
% songlist = {'9.124', '46.203', '06.290'}; % partially enter strings for song names
% motifnums = [4 2 4]; % rendition of the motif within the song
% motifname = 'jbh'; % motif
% sylnuminmotif = 3; % which syl do you want in the motif?
songlist = {'46.203'}; % partially enter strings for song names
motifnums = [4]; % rendition of the motif within the song
motifname = 'abh'; % motif
sylnuminmotif = 3; % which syl do you want in the motif?

disp('--')

songsout = {};
sylposout = [];
for i=1:length(songlist)
    tmp = dir(['*' songlist{i} '*.cbin']);
    disp(tmp(1).name);
    assert(length(tmp)==1, 'asdfasd');
    
    % -- find location of desired syl
    load([tmp(1).name '.not.mat']);
    
    inds = strfind(labels, motifname);
    
    sylpos = inds(motifnums(i))+sylnuminmotif-1;
    
    disp(['POSITION: ' num2str(sylpos)])
    disp(labels(sylpos))
    disp(labels(sylpos-2:sylpos+2))
    
    
    songsout = [songsout tmp(1).name];
    sylposout = [sylposout sylpos];
end

songsout
sylposout

%% ================== splicing syllables into a BOS file
clear all; close all;
songfname = '/bluejay5/lucas/birds/pu69wh78/110317_songs/pu69wh78_031117_102806.9552.cbin'; % --- load cbin file, splice syllables, then save as wav file
inds_abh = [18 34 49]; % inds within the song (not.mat)
inds_jbh = [28 43 58];
sylstoremove = [inds_abh inds_jbh];

% --- abh
replfiles_abh_hi = {'pu69wh78_101117_223206.290.cbin', ...
    'pu69wh78_101117_222646.203.cbin', 'pu69wh78_101117_203411.123.cbin'};
notes_abh_hi = [38 48 17];

replfiles_abh_lo = { 'pu69wh78_101117_114745.3.cbin', ...
    'pu69wh78_101117_114457.188.cbin', 'pu69wh78_101117_114745.3.cbin'};
notes_abh_lo = [36 83 36];

% --- jbh
replfiles_jbh_hi = {'pu69wh78_101117_212709.124.cbin', ...
    'pu69wh78_101117_222646.203.cbin', 'pu69wh78_101117_223206.290.cbin'};
notes_jbh_hi = [50 27 64];

replfiles_jbh_lo = {'pu69wh78_101117_104946.142.cbin', ...
    'pu69wh78_101117_114457.188.cbin', 'pu69wh78_101117_114745.3.cbin'};
notes_jbh_lo = [27 51 71];


% ReplacementFiles = {...
%     '/bluejay5/lucas/birds/pu69wh78/110317_songs/pu69wh78_031117_102806.9552.cbin', ...
%     '/bluejay5/lucas/birds/pu69wh78/110317_songs/pu69wh78_031117_102806.9552.cbin', ...
%     '/bluejay5/lucas/birds/pu69wh78/110317_songs/pu69wh78_031117_102806.9552.cbin'}; % list of files - extract replacement syls (provide in order, of sylstoremove)
% NotesToTake = [16 32 47]; % array, corresponding to Replacement files



% ========================== 1) ab(h) high; jb(h) low
ReplacementFiles = [replfiles_abh_hi replfiles_jbh_lo];
NotesToTake = [notes_abh_hi notes_jbh_lo];
suffix = 'aHI_jLO'; % for naming

UseOriginalDigPulse = 0; % if 0, then resegements with new onsets, offsets.
OutDir = '/bluejay5/egret_data/lucas/Test_Songs/BOS';

lt_songtools_splicesyls(songfname, sylstoremove, ReplacementFiles, NotesToTake, ...
    UseOriginalDigPulse, suffix, OutDir)



% ========================== 1) ab(h) high; jb(h) high
ReplacementFiles = [fliplr(replfiles_abh_hi) replfiles_jbh_hi];
NotesToTake = [fliplr(notes_abh_hi) notes_jbh_hi];
suffix = 'aHI_jHI'; % for naming

UseOriginalDigPulse = 0; % if 0, then resegements with new onsets, offsets.
OutDir = '/bluejay5/egret_data/lucas/Test_Songs/BOS';

lt_songtools_splicesyls(songfname, sylstoremove, ReplacementFiles, NotesToTake, ...
    UseOriginalDigPulse, suffix, OutDir)


% ========================== 1) ab(h) lo; jb(h) high
ReplacementFiles = [replfiles_abh_lo fliplr(replfiles_jbh_hi)];
NotesToTake = [notes_abh_lo fliplr(notes_jbh_hi)];
suffix = 'aLO_jHI'; % for naming

UseOriginalDigPulse = 0; % if 0, then resegements with new onsets, offsets.
OutDir = '/bluejay5/egret_data/lucas/Test_Songs/BOS';

lt_songtools_splicesyls(songfname, sylstoremove, ReplacementFiles, NotesToTake, ...
    UseOriginalDigPulse, suffix, OutDir)



% ========================== 1) ab(h) lo; jb(h) lo
ReplacementFiles = [fliplr(replfiles_abh_lo) fliplr(replfiles_jbh_lo)];
NotesToTake = [fliplr(notes_abh_lo) fliplr(notes_jbh_lo)];
suffix = 'aLO_jLO'; % for naming

UseOriginalDigPulse = 0; % if 0, then resegements with new onsets, offsets.
OutDir = '/bluejay5/egret_data/lucas/Test_Songs/BOS';

lt_songtools_splicesyls(songfname, sylstoremove, ReplacementFiles, NotesToTake, ...
    UseOriginalDigPulse, suffix, OutDir)









