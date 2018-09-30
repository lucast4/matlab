function lt_neural_Coher_Extract(SummaryStruct)
%% goes thru all songs and channels in SUmmary struct and extracts coherence and spec
% extracts all pairs and chans, but excludes those not in Summarystruct
% does not overwrite old extractions

%% lt 9/20/18 - 


%%

lt_switch_chronux(1);
movingwin = [0.1 0.01];
params = struct;
params.fpass = [1/movingwin(1) 150];
w = 30; % in hz, for desired frequency resolution of tapers. % note, t is set to movingwin(1)
tw = movingwin(1)*w;
params.tapers = [tw 2*tw-1];


% ==================== 1) Get list of all song files and correspond chans.
[Allbird_Fnames, Allbird_chanlist, Allbird_birdnum] = lt_neural_tools_allsongs(SummaryStruct);


% ===================== 2) GET COHEROGRAM FOR ALL SONG FILES AND CHAN PAIRS
for j=1:length(Allbird_Fnames)
    
    fnamethis = Allbird_Fnames{j};
    channellist = Allbird_chanlist{j};
    [filepart1, filepart2] = fileparts(fnamethis);
    
    if ~exist(fnamethis, 'file')
        disp(['skipping [doesnt exist]: ' fnamethis]);
        continue
    end
    
    % ---- make COHERENCE dir if doesn't already exist
    if ~exist([filepart1 '/COHERENCE'])
        mkdir([filepart1 '/COHERENCE'])
    end
    
    % ---- save params if doesn't already exist
    fname_params = [filepart1 '/COHERENCE/params.mat'];
    if ~exist(fname_params, 'file')
        save(fname_params, 'params');
        save([filepart1 '/COHERENCE/movingwin.mat'], 'movingwin');
    end
    
    if length(channellist)==1
        disp('skip, only one chan for this song ...');
        continue
    end
    
    % ========= CHECK IF ALL FILES ALREADY EXIST, IF SO THEN CONTINUE
    doanalysis = 0;
    for c=1:length(channellist)
        chan1 = channellist(c);
        fnamespec1 = [filepart1 '/COHERENCE/' filepart2 '.rhd.ch' num2str(chan1) '.spec'];
        
        % ###################### CHECK WHETHER THESE CHANS ALREADY DONE
        if exist(fnamespec1, 'file')
            continue
        end
        doanalysis = 1;
    end
    % --------------- COHEROGRAMS - GO THRU ALL PAIRS OF CHANNELS
    for c=1:length(channellist)
        for cc=c+1:length(channellist)
            chan1 = channellist(c);
            chan2 = channellist(cc);
            
            % ###################### CHECK WHETHER THESE CHANS ALREADY DONE
            fnamecoh = [filepart1 '/COHERENCE/' filepart2 '.rhd.ch' num2str(chan1) '-' num2str(chan2) '.coh'];
            if exist(fnamecoh, 'file')
                continue
            end
        doanalysis =1;
        end
    end
    
    if doanalysis==0
        disp('skip, alrady done');
        continue
    end
    
    % ========= EXTRACT SONG DATA FOR THIS FILE
    [amplifier_data, ~, frequency_parameters, board_adc_data, ~, amplifier_channels, ...
        ~, t_amplifier] = pj_readIntanNoGui(fnamethis);
    params.Fs = frequency_parameters.amplifier_sample_rate;
    
    
    % --------------- SPECTROGRAMS - GO THRU ALL CHANNELS
    tall = [];
    fall = [];
    for c=1:length(channellist)
        chan1 = channellist(c);
        
        fnamespec1 = [filepart1 '/COHERENCE/' filepart2 '.rhd.ch' num2str(chan1) '.spec'];
        
        % ###################### CHECK WHETHER THESE CHANS ALREADY DONE
        if exist(fnamespec1, 'file')
            continue
        end
        
        
        ind1 = [amplifier_channels.chip_channel] == chan1;
        dat1 = amplifier_data(ind1,:);
        
        % ================ SPETROGRAMS
        [S1,t,f]= mtspecgramc(dat1', movingwin, params);
        S1 = single(S1);
        tall = [tall; t];
        fall = [fall; f];
        
        % ================ SAVE TO DISK
        disp(fnamespec1)
        save(fnamespec1, 'S1');
        % ----- save t and f information [should be identical 
        fnamet = [filepart1 '/COHERENCE/' filepart2 '.rhd.tbins'];
        if ~exist(fnamet, 'file')
            save(fnamet, 't')
        end
        fnameff = [filepart1 '/COHERENCE/' filepart2 '.rhd.ffbins'];
        if ~exist(fnameff, 'file')
            save(fnameff, 'f')
        end
        
    end
    
    
    % --------------- COHEROGRAMS - GO THRU ALL PAIRS OF CHANNELS
    for c=1:length(channellist)
        for cc=c+1:length(channellist)
            chan1 = channellist(c);
            chan2 = channellist(cc);
            
            % ###################### CHECK WHETHER THESE CHANS ALREADY DONE
            fnamecoh = [filepart1 '/COHERENCE/' filepart2 '.rhd.ch' num2str(chan1) '-' num2str(chan2) '.coh'];
            if exist(fnamecoh, 'file')
                continue
            end
            
            % ============================
            ind1 = [amplifier_channels.chip_channel] == chan1;
            ind2 = [amplifier_channels.chip_channel] == chan2;
            
            dat1 = amplifier_data(ind1,:);
            dat2 = amplifier_data(ind2,:);
            
            % ================ COHEREHENCE
            [C,phi,S12,S1,S2,t,f] = cohgramc(dat1', dat2', movingwin, params);
            C = single(C);
            tall = [tall; t];
            fall = [fall; f];
            
            % ================ SAVE TO DISK
            disp(fnamecoh);
            save(fnamecoh, 'C');
            
            %% ========= sanity check
            % plot song spec, raw lfp, spectra, and coherogram\
            if (0)
                lt_figure; hold on;
                hsplots = [];
                
                % -------------------- SONG SPECTROGRAM
                hsplot = lt_subplot(6,1,1); hold on;
                hsplots = [hsplots; hsplot];
                lt_plot_spectrogram(board_adc_data, params.Fs, 1, 0);
                
                
                % -------------------- RAW LFP
                t = [1:length(dat1)]./params.Fs;
                hsplot = lt_subplot(6,1,2);
                hsplots = [hsplots; hsplot];
                plot(t, dat1, '-k');
                axis tight;
                
                hsplot = lt_subplot(6,1,3);
                hsplots = [hsplots; hsplot];
                plot(t, dat2, '-b');
                axis tight;
                
                % -------------------- INDIVIDUAL SPECTROGRAM
                [S1,t,f]= mtspecgramc(dat1', movingwin, params);
                [S2,t,f]= mtspecgramc(dat2', movingwin, params);
                
                hsplot = lt_subplot(6,1,4); hold on;
                hsplots = [hsplots; hsplot];
                imagesc(t, f, 10*log10(S1)');
                axis tight;
                ylim([1/(movingwin(1)) params.fpass(2)]);
                
                hsplot = lt_subplot(6,1,5); hold on;
                hsplots = [hsplots; hsplot];
                imagesc(t, f, 10*log10(S2)');
                axis tight;
                ylim([1/(movingwin(1)) params.fpass(2)]);
                
                
                % -------------------- COHEROGRAM
                hsplot= lt_subplot(6,1,6); hold on;
                hsplots = [hsplots; hsplot];
                % title(['ch' num2str(DatAll(j).chan_amp) '[' DatAll(j).bregion ']']);
                Cdb = 10*log10(C); % convert to db
                Cdb = C; % convert to db
                imagesc(t, f, Cdb');
                axis tight;
                ylim([1/(movingwin(1)) params.fpass(2)]);
                colormap('jet');
                
                % ----------------------------------
                linkaxes(hsplots, 'x');
                
                
            end
            
            
        end
    end
    
    % === sdanity check, make sure akll have same t and f.
    if ~isempty(tall)
    assert(all(all((diff(tall)<0.001))));
    assert(all(all((diff(fall)<0.001))));
    end
end

lt_switch_chronux(0);