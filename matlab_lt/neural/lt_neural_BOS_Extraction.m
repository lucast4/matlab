function lt_neural_BOS_Extraction(SummaryBOS, expttoget)
%% lt 12/1/18 - Main extraction, from raw songs, of BOS data


i=expttoget; % CAN ONLY DO ONE AT A TIME

%%
cd(SummaryBOS.expt(i).dirname);
%% =============== PLOT SOME RAW SONG TRIALS FOR EACH BATCH/CHANNEL


%% =============== EXTRACT RAW DATA IN MAT FORMAT
batchname = SummaryBOS.expt(i).batchname;
chanstoget = SummaryBOS.expt(i).channels;

% ==== make directories, one for each batch/channel combination
for chan=chanstoget
    dirname = [batchname '-Chan' num2str(chan)];
    if ~exist(dirname)
        mkdir([batchname '-Chan' num2str(chan)]);
    end
end

% ---- for each file, extrat neural data for desired channels and save in
% folder.
filenames = textread(batchname, '%s');

for fname=filenames'
    fname = fname{1};
    % ============= load and extract each channel
    [amplifier_data,board_dig_in_data,frequency_parameters, ...
        board_adc_data, board_adc_channels, amplifier_channels, ...
        board_dig_in_channels, t_amplifier] = pj_readIntanNoGui(fname);
    
    [a, b] = fileparts(fname);
    
    % === go thru each channle
    for chan = chanstoget
        
        indthis = [amplifier_channels.chip_channel]==chan;
        assert(sum(indthis)==1, 'bad chans?');
        data = single(amplifier_data(indthis,:));
        
        % === save
        
        fnamethis = [batchname '-Chan' num2str(chan) '/' b '.ch' num2str(chan) '.mat'];
        save(fnamethis, 'data');
        
        disp(['saved ' fnamethis]);
    end
end

%% === save a new batch file, but with updated names
for chan=chanstoget
    fname = [batchname '-Chan' num2str(chan) '/' batchname '_ch' num2str(chan) '.txt'];
    fid = fopen(fname, 'w');
    for j=1:length(filenames)
        fname = filenames{j};
        [~, b] = fileparts(fname);
        fnamewrite = [b '.ch' num2str(chan) '.mat'];
        fprintf(fid, '%s\n', fnamewrite);
        %        fwrite(fid, fnamewrite);
    end
    fclose(fid);
    
    %     % ====== also save one with "spikes"
    %     fname = [batchname '_ch' num2str(chan) '_spikes.txt'];
    %     fid = fopen(fname, 'w');
    %     for j=1:length(filenames)
    %        fname = filenames{j};
    %        [~, b] = fileparts(fname);
    %        fnamewrite = [b '.ch' num2str(chan) '_spikes.mat'];
    %        fprintf(fid, '%s\n', fnamewrite);
    % %        fwrite(fid, fnamewrite);
    %     end
    %     fclose(fid);
end




%% =============== LFP




%% =============== SPIKES
% ============ 1) EXTRACT THRESHOLD CROSSES FOR EACH DATA FILE
for chan=chanstoget
    cd([batchname '-Chan' num2str(chan)])
    Get_spikes([batchname '_ch' num2str(chan) '.txt']);
    disp(['===========' batchname '-Chan' num2str(chan)]);
    cd ..
end

%
% Get_spikes('Batch2353to2357_ch9.txt');


%% ============ 2) combine all threshold crosses into one file
% --- in that file, indicate where each spike originated from
for chan=chanstoget
    cd([batchname '-Chan' num2str(chan)])
    disp([batchname '-Chan' num2str(chan)]);
    % ================== go thru all spike files and concatenate all spike
    % times and shapes into large matrix
    filenames = textread([batchname '_ch' num2str(chan) '.txt'], '%s');
    
    spikes_all = [];
    index_all = [];
    threshold_all = [];
    %     filename_all = {};
    filenuminbatch_all = [];
    spikenumberInSong_all = [];
    spiketimeInSong_all = [];
    for j=1:length(filenames)
        % --- load spikes
        [a, b] = fileparts(filenames{j});
        filethis = [b '_spikes.mat'];
        
        tmp = load(filethis);
        spikes_all = [spikes_all; tmp.spikes];
        index_all = [index_all tmp.index];
        threshold_all = [threshold_all tmp.threshold];
        
        filenuminbatch_all = [filenuminbatch_all j*ones(size(tmp.index))];
        spikenumberInSong_all = [spikenumberInSong_all 1:length(tmp.index)];
    end
    
    % ============ save a summary spikes mat
    spikes = spikes_all;
    index = index_all;
    par = tmp.par;
    %     save(['allsongs_spikes.mat'], 'spikes', 'index', 'par','psegment','sr_psegment')
    save(['allsongs_spikes.mat'], 'spikes', 'index', 'par', 'filenuminbatch_all', 'spikenumberInSong_all', 'threshold_all')
    %     spikes_struct = struct;
    %     spikes_struct.spikes = spikes_all;
    %     spikes_struct.index = index_all;
    %     spikes_struct.threshold_all = threshold_all;
    %     spikes_struct.filenuminbatch_all = filenuminbatch_all;
    %     spikes_struct.spikenumberInSong_all = spikenumberInSong_all;
    %
    %     save('spikes_struct.mat');
    cd ..
end


%% ============ 3) Cluster the singel file
for chan=chanstoget
    cd([batchname '-Chan' num2str(chan)])
    Do_clustering('allsongs_spikes.mat')
    cd ..
end

% Do_clustering('Batch2353to2357_ch9_spikes.txt', 'make_plots', true);