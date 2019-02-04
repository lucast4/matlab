function lt_neural_tools_LFPextract(SummaryStruct, freqvals, skipifdone)
%% LT 12/13/18 - uses filters from Jai Yu (frank lab), freq bands are (f-1, f+2).

%% ==== get list of all files across dataset
[Allbird_Fnames, Allbird_chanlist, Allbird_birdnum] = lt_neural_tools_allsongs(SummaryStruct);

%% go thru all files and extract LFP
for j=1:length(Allbird_Fnames)
    
    fnamethis = Allbird_Fnames{j};
%     channellist = Allbird_chanlist{j};
%     [filepart1, filepart2] = fileparts(fnamethis);
    
    if ~exist(fnamethis, 'file')
        disp(['skipping [doesnt exist]: ' fnamethis]);
        continue
    end
    
    
    
    % ########################### SAVE LFP
    fnamesave = [fnamethis '.lfp'];
    load(fnamesave, '-mat');
    
        
    % --=-- skuip if already done
    if skipifdone ==1
        if exist([fnamethis '.filt'], 'file')
            % === cjheck that all chans are dopne
            
            tmp = load([fnamethis '.filt'], '-mat');
            if all(ismember(lfpstruct.chanlist, tmp.filtstruct.chanlist))
                            disp(['skip ' fnamethis '[alraedy done, including desired chans]']);
            continue
            end
        end
    end
    
    
    % ################# ACROSS ALL CHANNELS, DO FILTER
    nchans = length(lfpstruct.chanlist);
    datfilt_chans = cell(1,nchans);
    for cc=1:nchans
        
        datthis = lfpstruct.dat(:,cc);
        
        % ===== run data thru filter bank
        datthis_filt = jai_waveletfilt_FiltDat(datthis, freqvals);
        
        % ----- convert to single
        datthis_filt = cellfun(@single, datthis_filt, 'UniformOutput', 0);
        
        if (0)
           figure; hold on;
           plot(lfpstruct.t, datthis, '-k');
           plot(lfpstruct.t, real(datthis_filt{freqvals==10}), '-r');
           plot(lfpstruct.t, real(datthis_filt{freqvals==10}), '-g');
        end
        
        datfilt_chans{cc} = datthis_filt;
    end
    
    % ======= SAVE
    filtstruct = struct;
    filtstruct.datfilt_chans = datfilt_chans;
    filtstruct.chanlist = lfpstruct.chanlist;
    filtstruct.freqvals = freqvals;
    filtstruct.t_edges = [lfpstruct.t(1) lfpstruct.t(end)];
    
    fnamesave = [fnamethis '.filt'];
    save(fnamesave, 'filtstruct');
    disp(['DONE: ' fnamesave]);
end
