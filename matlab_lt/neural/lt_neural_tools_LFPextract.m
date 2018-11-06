function lt_neural_tools_LFPextract(SummaryStruct)
%% lt 10/19/18 - extracts and saves LFP for all chanesl
% DEFAULTS:
% 400hz lowpass, 4th order butterworth, downsampled to 1500. Based on
% plotting pre vs post, and comparing diff frequencies and Ns, determined
% this to be appropriate. downsampling is fine.

[Allbird_Fnames, Allbird_chanlist, Allbird_birdnum] = lt_neural_tools_allsongs(SummaryStruct);

for j=1:length(Allbird_Fnames)
    
    fnamethis = Allbird_Fnames{j};
    channellist = Allbird_chanlist{j};
    [filepart1, filepart2] = fileparts(fnamethis);
    
    if ~exist(fnamethis, 'file')
        disp(['skipping [doesnt exist]: ' fnamethis]);
        continue
    end
    
    % =================== SAVE LFP FOR ALL RELEVANT CHANNELS
    [amplifier_data,~,frequency_parameters, ~, ...
        ~, amplifier_channels, ~, ~] =...
        pj_readIntanNoGui(fnamethis);
    fs = frequency_parameters.amplifier_sample_rate;
    
    indschan = ismember([amplifier_channels.chip_channel], channellist);
    dat = amplifier_data(indschan, :);
    %     dat = amplifier_data([amplifier_channels.chip_channel] == chanamp, :);
    
    
    
    % ##################### GET LFP
    if (1)
        % ================= LOWPASS
        N=4;
        fpass = 400;
        [filt_b,filt_a]=butter(N, [fpass*2/fs]);
        datfilt = filtfilt(filt_b, filt_a, dat'); % do multiple channels at once
        datfilt = datfilt';
        if (0)
            lt_figure; hold on;
            plot(dat, 'k');
            plot(datfilt, 'm');
        end
        if (0)
            % NOTE: from iterating over some pass freque4mncies and N
            % orders, I determined that good pass frequency is 400 and good
            % N is 4.
            passflist = [100 200 300 400 500];
            Nlist = [1 2 3 5 8];
            lt_figure; hold on;
            plot(dat, 'k');
            plotcols = lt_make_plot_colors(length(passflist), 1, [1 0 0]);
            %             for kkk=1:length(passflist)
            kkk=4;
            N=4;
            for nnn=1:length(Nlist)
                N = Nlist(nnn);
                
                fpass = passflist(kkk);
                [filt_b,filt_a]=butter(N, [fpass*2/fs]);
                datfilt = filtfilt(filt_b, filt_a, dat);
                %                 plot(datfilt, 'Color', plotcols{kkk});
                plot(datfilt, 'Color', plotcols{nnn});
            end
        end
    else
        % ---- FILTER (0.5 to 400hz) (2pole bessel?)
        [datfilt,neuralFiltLow,neuralFiltHi] =lt_neural_filter(dat, ...
            frequency_parameters, 0, 5, 400);
        if (0)
            figure; hold on;
            t = [1:length(datthis)]/30000;
            plot(t, datthis, '-k'); plot(t, datfilt, 'r');
        end
    end
    
    
    % ---- DOWNSAMPLE
    fs_new = 1500;
    t = [1:size(datfilt,2)]/fs;
    factor = floor(fs/fs_new);
    datfilt = downsample(datfilt', factor);
    t = downsample(t, factor);
    
    % ########################### FINAL STUFF.
    datfilt = single(datfilt);
    t = single(t);
    
    if (0)
        lt_figure; hold on;
        plot(t, datfilt, 'b');
        plot(t(1:factor:end), datfilt_dn, 'r');
    end
    
    % ########################### SAVE LFP
    lfpstruct = struct;
    lfpstruct.dat = datfilt;
    lfpstruct.t = t;
    lfpstruct.chanlist = channellist;
    
    % 1) lfp
    fnamesave = [fnamethis '.lfp'];
    save(fnamesave, 'lfpstruct');
    
    % ===== final troubleshoot (compare original dat to processed)
    if (0)
        indtoplot = 2;
        lt_figure; hold on;
        plot([1:size(dat(indtoplot,:),2)]/fs, dat(indtoplot,:), '-k');
        plot(t, datfilt(:,indtoplot), 'b');
    end
    
    
    disp(['DONE: ' fnamethis]);
end
