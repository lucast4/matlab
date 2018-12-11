function lt_batchsong_extractWNhit(ListOfDirs, ListOfBatch, SongFnames)

%% lt 8/4/18 - saves file indication whether there was WN hit. This is useful mainly
% for neural analysis (.rhd files), since evtaf already saaves this info.
% for rhd files will determine based on audio power (WN presumably has
% high power)

%%
% ListOfDirs = {...
%     '/bluejay5/lucas/birds/pu69wh78/NEURAL/110117_RALMANlearn1'};
% ListOfBatch = {...
%     'BatchAll'};

% SongFnames =
%   336×1 cell array
%
%     '/bluejay5/lucas/birds/wh44wh39/NEURAL/032118_RALMANlearn3/wh44wh39_180321_112614.rhd'
%     '/bluejay5/lucas/birds/wh44wh39/NEURAL/032118_RALMANlearn3/wh44wh39_180321_112740.rhd'
%     '/bluejay5/lucas/birds/wh44wh39/NEURAL/032118_RALMANlearn3/wh44wh39_180321_112830.rhd'
%     '/bluejay5/lucas/birds/wh44wh39/NEURAL/032118_RALMANlearn3/wh44wh39_180321_112850.rhd'
%     '/bluejay5/lucas/birds/wh44wh39/NEURAL/032118_RALMANlearn3/wh44wh39_180321_112956.rhd'
%     '/bluejay5/lucas/birds/wh44wh39/NEURAL/032118_RALMANlearn3/wh44wh39_180321_113034.rhd'
% Will use SongFnames if the first 2 arguments are empty ([]).


%% run
% NOTE: currently has ugly duplication (if, else...), combine into one
% loop.

WNdurations = [];
NumMaxes = [];

% ========= go thru all dirs and batches. for each file read audio and
% determine timepoitns of WN

if ~isempty(ListOfDirs)
    for i=1:length(ListOfDirs)
        
        cd(ListOfDirs{i});
        
        % =======
        fid = fopen(ListOfBatch{i});
        fline = fgetl(fid);
        
        while ischar(fline)
            disp(fline);
            
            % ============== check extension. if is cbin then skip. i.e.
            % continue if is rhd
            [~, ~, fe] = fileparts(fline);
            
            
            if strcmp(fe, '.rhd')
                
                
                % ============== figure out if already done. if yes, then skip
                % (optional)
                fname_save = [fline '.wntime.mat'];
                if exist(fname_save, 'file')
                    disp(['skipping (wn extract already done): ' fline]);
                    fline = fgetl(fid);
                    continue
                else
                    disp(['EXTRACTING (wn): ' fline]);
                end
                
                
                % ========= load file
                [fp, songdat] = ...
                    pj_readIntanNoGui_AudioOnly(fline);
                fs=fp.amplifier_sample_rate;
                
                % ========= call WN any time that reaches threshold
                WNtimes = find(songdat>3.2);
                
                if ~isempty(WNtimes)
                    % - pad to make extraction of onsets and offsets easier
                    WNtimes = [0 WNtimes length(songdat)];
                end
                WNedges = find(diff(WNtimes)>0.08*fs);
                
                % --- samples that are to the left of gaps
                WNoffsets = WNtimes(WNedges(2:end))+0.0015*fs; % add 2ms to flank to account for erise and fall.
                WNonsets = WNtimes(WNedges(1:end-1)+1)-0.0015*fs;
                
                
                % ========================== COLLECT THINGS FOR CLULSTERING
                % LATER.
                WNdurations = [WNdurations WNoffsets-WNonsets];
                
                for k = 1:length(WNonsets)
                    
                    nm = sum(WNtimes>=WNonsets(k) & WNtimes<=WNoffsets(k)); %
                    NumMaxes = [NumMaxes nm];
                    % number of times there are maxes within the duration of
                    % this putative WN.
                end
                
                
                % ======= sanity cehck - plot audio
                if (0)
                    lt_figure; hold on;
                    
                    % --- 1) audito
                    lt_subplot(3,2,1); hold on;
                    plot(songdat, '-k');
                    
                    for k=1:length(WNoffsets)
                        line([WNonsets(k) WNonsets(k)], ylim, 'Color', 'b');
                        line([WNoffsets(k) WNoffsets(k)], ylim, 'Color', 'r');
                    end
                    
                    % --- histo
                    lt_subplot(3,2,2); hold on;
                    title('ampl^2');
                    lt_plot_histogram(songdat.^2);
                    
                end
                
                
                
                % ========= SAVE
                wnstruct = struct;
                wnstruct.WNonsets = single(WNonsets./fs); % convert to sec
                wnstruct.WNoffsets = single(WNoffsets./fs);
                wnstruct.Nmaxes_withindur = single(NumMaxes(end-length(WNonsets)+1:end));
                
                save(fname_save, 'wnstruct');
            else
                disp(['skipping wn extract, not rhd file: ' fline]);
            end
            
            fline = fgetl(fid);
        end
        
    end
else
    % Use SongFnames
    for i=1:length(SongFnames)
        
        fline = SongFnames{i};
        
        disp(fline);
        
        % ============== check extension. if is cbin then skip. i.e.
        % continue if is rhd
        [~, ~, fe] = fileparts(fline);
        
        
        if strcmp(fe, '.rhd')
            
            
            % ============== figure out if already done. if yes, then skip
            % (optional)
            fname_save = [fline '.wntime.mat'];
            if exist(fname_save, 'file')
                disp(['skipping (wn extract already done): ' fline]);
%                 fline = fgetl(fid);
                continue
            else
                disp(['EXTRACTING (wn): ' fline]);
            end
            
            % ==== ognore if file doesn't exist
            if ~exist(fline)
            disp(['skipping (audio file doesnt exist): ' fline]);
            continue
            end
            
            % ========= load file
            [fp, songdat] = ...
                pj_readIntanNoGui_AudioOnly(fline);
            fs=fp.amplifier_sample_rate;
            
            % ========= call WN any time that reaches threshold
            WNtimes = find(songdat>3.2);
            
            if ~isempty(WNtimes)
                % - pad to make extraction of onsets and offsets easier
                WNtimes = [0 WNtimes length(songdat)];
            end
            WNedges = find(diff(WNtimes)>0.08*fs); % gaps between offsets and onsets
            
            % --- samples that are to the left of gaps
            WNoffsets = WNtimes(WNedges(2:end))+0.0015*fs; % add 2ms to flank to account for erise and fall.
            WNonsets = WNtimes(WNedges(1:end-1)+1)-0.0015*fs;
            
            
            % ========================== COLLECT THINGS FOR CLULSTERING
            % LATER.
            WNdurations = [WNdurations WNoffsets-WNonsets];
            
            for k = 1:length(WNonsets)
                
                nm = sum(WNtimes>=WNonsets(k) & WNtimes<=WNoffsets(k)); %
                NumMaxes = [NumMaxes nm];
                % number of times there are maxes within the duration of
                % this putative WN.
            end
            
            
            % ======= sanity cehck - plot audio
            if (0)
                lt_figure; hold on;
                
                % --- 1) audito
                lt_subplot(3,2,1); hold on;
                plot(songdat, '-k');
                
                for k=1:length(WNoffsets)
                    line([WNonsets(k) WNonsets(k)], ylim, 'Color', 'b');
                    line([WNoffsets(k) WNoffsets(k)], ylim, 'Color', 'r');
                end
                
                % --- histo
                lt_subplot(3,2,2); hold on;
                title('ampl^2');
                lt_plot_histogram(songdat.^2);
            end
            
            
            
            % ========= SAVE
            wnstruct = struct;
            wnstruct.WNonsets = single(WNonsets./fs); % convert to sec
            wnstruct.WNoffsets = single(WNoffsets./fs);
            wnstruct.Nmaxes_withindur = single(NumMaxes(end-length(WNonsets)+1:end));
            
            save(fname_save, 'wnstruct');
        else
            disp(['skipping wn extract, not rhd file: ' fline]);
        end
        
    end
end
%% ===== DIAGNOSIS, plot in 2d space all WN detects
if ~isempty(WNdurations)
    lt_figure; hold on;
    xlabel('duration (s)');
    ylabel('n maxes within dur');
    
    plot(WNdurations./fs, NumMaxes, 'kx');
end
