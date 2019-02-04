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
disp('NOTE: will extract song bouts that are entirely within a file, with at least maxgapdur time from edges of file');

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
    
AllLFP_datByUnits = {};
AllLFP_tedges = [];
AllLFP_chanlist = {};

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
            SpksAllUnits{nn} = single(spkthis);
            
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
        
        
        %% =============== LFP
        
        [LFP_dat, LFP_tedges, LFP_chanlist] = ...
            lt_neural_BOS_Extraction_sub(SummaryBOS, fline, i, amplifier_data, ...
            amplifier_channels, fs, t_amplifier);
        
        
        
        
        %% ############################### SAVE OUTPUT FOR THIS FILE
        AllTotalSamps = [AllTotalSamps; length(pulsedat)];
        Allt_amplifier = [Allt_amplifier; single(t_amplifier)];
        AllOnsets_samp = [AllOnsets_samp; single(rtimes)]; % each index is one song.
        AllOffsets_samp = [AllOffsets_samp; single(ftimes)];
        AllSpktimesByUnits = [AllSpktimesByUnits; SpksAllUnits];
        AllSongOnsets_inds = [AllSongOnsets_inds; single(songonsets)];
        AllSongOffsets_inds = [AllSongOffsets_inds; single(songoffsets)];
        AllSongOnsets_ActualTime = [AllSongOnsets_ActualTime; songonsets_actualtime];
        
        AllLFP_datByUnits = [AllLFP_datByUnits; LFP_dat];
        AllLFP_tedges = [AllLFP_tedges; LFP_tedges];
        AllLFP_chanlist = [AllLFP_chanlist; LFP_chanlist];
            
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
    SummaryBOS.expt(i).DAT.AllLFP_datByUnits = AllLFP_datByUnits;
    SummaryBOS.expt(i).DAT.AllLFP_tedges = AllLFP_tedges;
    SummaryBOS.expt(i).DAT.AllLFP_chanlist = AllLFP_chanlist;
    SummaryBOS.expt(i).fs = fs;
end

lt_figure; hold on;
title('gap durations (and line used as threshold for getting bouts)');
lt_plot_histogram(gapdurall, 0:0.01:1);
line([maxgapdur maxgapdur], ylim);



%% ==================== [EXTRACT - SANITY PLOT] all files, with extracted songs
% [USEFUL!!!] quickly plots data and extracted song files.
if (0)
    i=1;
    
    figcount=1;
    subplotrows=8;
    subplotcols=1;
    fignums_alreadyused=[];
    hfigs=[];
    
    nfiles = length(SummaryBOS.expt(i).DAT.AllOffsets_samp);
    for j=1:nfiles
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        t = SummaryBOS.expt(i).DAT.Allt_amplifier{j};
        ons = SummaryBOS.expt(i).DAT.AllOnsets_samp{j};
        offs = SummaryBOS.expt(i).DAT.AllOffsets_samp{j};
        
        
        lt_neural_QUICK_PlotSylPatches(t(ons), t(offs), 1, 1);
        xlim([t(1) t(end)]);
        
        % ===== song extractions
        sons = SummaryBOS.expt(i).DAT.AllSongOnsets_inds{j};
        soffs = SummaryBOS.expt(i).DAT.AllSongOffsets_inds{j};
        
        lt_neural_QUICK_PlotSylPatches(t(ons(sons)), t(offs(soffs)), 2, 1);
        ylim([0 3]);
        ylabel('syls - songs extracted');
        xlabel('range of entire file (sec)');
    end
end
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
flanktotake = [-3 3]; % sec from onset and offset.
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
    
    LFPall = {};
    LFPall_t = {};
    LFPall_chanlist = {};
    for ss=1:numsongs
        
        % ===== stats for this file
        onsets_samps = SummaryBOS.expt(i).DAT.AllOnsets_samp{ss};
        offsets_samps = SummaryBOS.expt(i).DAT.AllOffsets_samp{ss};
        
        t = SummaryBOS.expt(i).DAT.Allt_fromDatOnset{ss};
        SpkTimes = SummaryBOS.expt(i).DAT.AllSpktimesByUnits(ss,:);
        nunits = length(SpkTimes);
        
        LFP_dat = SummaryBOS.expt(i).DAT.AllLFP_datByUnits(ss,:);
        LFP_t = SummaryBOS.expt(i).DAT.AllLFP_tedges(ss,:);
        LFP_chanlist = SummaryBOS.expt(i).DAT.AllLFP_chanlist{ss};
        % -- time bins for LFP
        t_LFP = linspace(LFP_t(1), LFP_t(2), length(LFP_dat{1}));
        t_LFP = (t_LFP-t_LFP(1)) + (t_LFP(2)-t_LFP(1))/2; % get relative to dat onset (to match other stuff)
        
        
        
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
            
            
            
            % ============= GET LFP 
            indsgood_LFP = t_LFP>tonset_sec+flanktotake(1) ...
                    & t_LFP<toffset_sec+flanktotake(2);
            lfpall = cellfun(@(x)x(indsgood_LFP), LFP_dat, 'UniformOutput', 0);
            t_LFP_short = t_LFP(indsgood_LFP)-tonset_sec;
            
            
            % #####################3 COLLECT THINGS ABOUT THIS SONG
            SylOnsets = [SylOnsets; allonsets];
            SylOffsets = [SylOffsets; alloffsets];
            SongDur = [SongDur; toffset_sec-tonset_sec];
            SpkTime_RelSongOnset = [SpkTime_RelSongOnset; spktimes_tmp];
            TEdges = [TEdges; tedges];
            SongOnset_ActualDatenum = [SongOnset_ActualDatenum; dtnum];
            
            LFPall = [LFPall; lfpall];
            LFPall_chanlist = [LFPall_chanlist; LFP_chanlist];
            LFPall_t = [LFPall_t; t_LFP_short];
        end
    end
    
        SummaryBOS.expt(i).DAT_bysongrend.TEdges = TEdges;
        SummaryBOS.expt(i).DAT_bysongrend.SylOnsets = SylOnsets;
        SummaryBOS.expt(i).DAT_bysongrend.SylOffsets = SylOffsets;
        SummaryBOS.expt(i).DAT_bysongrend.SongDur = SongDur;
        SummaryBOS.expt(i).DAT_bysongrend.SpkTime_RelSongOnset = SpkTime_RelSongOnset;
        SummaryBOS.expt(i).DAT_bysongrend.SongOnset_ActualDatenum = SongOnset_ActualDatenum;
        SummaryBOS.expt(i).DAT_bysongrend.LFPall = LFPall;
        SummaryBOS.expt(i).DAT_bysongrend.LFPall_t = LFPall_t;
        SummaryBOS.expt(i).DAT_bysongrend.LFPall_chanlist = LFPall_chanlist;
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
if (0)
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