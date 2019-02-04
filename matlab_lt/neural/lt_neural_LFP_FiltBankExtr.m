function FILTLFPSTRUCT = lt_neural_LFP_FiltBankExtr(SummaryStruct, SwitchCohStruct, MOTIFSTATS_pop, ...
    PARAMS, BirdsToPlot, extrapad, saveON, skipifdone)

%% lt 12/13/18 - filter bacnk data, extracts singasl.
%
% NOTE: only extracts for neuron sets that are used in at leat one of the
% the actual switcvhes in SwitchCohStruct.[in this setsne differs from
% LFPSTRUCT, since latter extracts all swtiche.]
freqvals = [20:5:40];
skipifOnlyOneChan = 1;
firstExtrAllSongs = 1; % keep at 1, much faster. if 0 then extracts each song once for each trial and motif.
% if 1, then extracts once for each song file...

% ========= AUTO PARAMS
motif_predur = PARAMS.motif_predur;
FILTLFPSTRUCT = struct;

%% ======= save dfir
savedir = ['/bluejay5/lucas/analyses/neural/FILTER/MOTIFEXTRACT/' PARAMS.savemarker];
if ~exist(savedir, 'dir')
    mkdir(savedir)
end

%% ============== RUN

numbirds = length(SummaryStruct.birds);
for i=1:numbirds
    birdname = SummaryStruct.birds(i).birdname;
    disp(birdname);
    
    if ~any(strcmp(birdname, BirdsToPlot))
        continue
    end
    
    numexpts = length(SwitchCohStruct.bird(i).exptnum);
    for ee = 1:numexpts
        disp(['expt ' num2str(ee)]);
        
        numsets = length(MOTIFSTATS_pop.birds(i).exptnum(ee).Sets_neurons);
        
        % === get list of sets to get
        setstoget = [];
        for j=1:length(SwitchCohStruct.bird(i).exptnum(ee).switchlist)
            if isempty(SwitchCohStruct.bird(i).exptnum(ee).switchlist(j).motifnum)
                continue
            end
            setstoget = [setstoget ...
                unique([SwitchCohStruct.bird(i).exptnum(ee).switchlist(j).motifnum.neursetused])];
        end
        if size(setstoget,1)>1
            setstoget = setstoget';
        end
        if isempty(setstoget)
            continue
        end
        
        % =============== FOR THIS BIRD, COLLECT ALL COHERENCE
        for ss = setstoget
            disp(['neuron set ' num2str(ss)]);
            
            dat = MOTIFSTATS_pop.birds(i).exptnum(ee).DAT.setnum(ss);
            nummotifs = length(dat.motif);
            
            % ================== FOR THIS SET, COLLECT ALL LFP DATA
            Neurlist = MOTIFSTATS_pop.birds(i).exptnum(ee).Sets_neurons{ss};
            %             Fnamelist = MOTIFSTATS_pop.birds(i).exptnum(ee).Sets_songfiles{ss};
            %             %     dirname = SummaryStruct.birds(i).neurons(Neurlist(1)).dirname;
            Chanlist = [SummaryStruct.birds(i).neurons(Neurlist).channel];
            %             Chanlist_sort = sort(Chanlist);
            Bregionlist = {SummaryStruct.birds(i).neurons(Neurlist).NOTE_Location};
            [~, indsort] = sort(Chanlist);
            Chanlist = Chanlist(indsort);
            Bregionlist = Bregionlist(indsort);
            
            % ----- only one chan, then skip
            if skipifOnlyOneChan==1
                if length(unique(Chanlist))==1
                    continue
                end
            end
            
            % --- get dirname
            dirname_main = {};
            for j=1:length(Neurlist)
                dname = fileparts(SummaryStruct.birds(i).neurons(Neurlist(j)).dirname);
                dirname_main = [dirname_main; dname];
            end
            dirname_main = unique(dirname_main); assert(length(dirname_main)==1);
            dirname_main = dirname_main{1};
            
            % =============== FIRST LOAD ALL SONG FILES THAT WILL NEED -
            % WILL THEN EXTRACT SPECIFIC MOTIFS AND TRIALS
            if firstExtrAllSongs==1
                songfiles = {};
                for mm=1:nummotifs
                    if isempty(dat.motif(mm).SegExtr_neurfakeID)
                        continue
                    end
                    
                    % --- segextract is shared across channels
                    segextract = dat.motif(mm).SegExtr_neurfakeID(1).SegmentsExtract;
                    songfiles = [songfiles {segextract.song_filename}];
                end
                songfiles = unique(songfiles);
                
                % ========== channels needed for analysis
                FILTMATall = cell(1, length(songfiles));
                tall = {};
                for ff =1:length(songfiles)
%                     disp(ff);
                    fnamethis = [dirname_main '/' songfiles{ff} '.filt'];
                    filtstruct = load(fnamethis, '-mat');
                    filtstruct = filtstruct.filtstruct;
                    
                    indschans = ismember(filtstruct.chanlist, Chanlist);
                    if size(Chanlist',1)~=size([filtstruct.chanlist(indschans)], 1)
                        Chanlist = Chanlist';
                    end
%                         assert(all([filtstruct.chanlist(indschans)]==Chanlist'), 'not matcjed..');    
                    assert(all([filtstruct.chanlist(indschans)]==Chanlist'), 'not matcjed..');
                    indsfreq = ismember(filtstruct.freqvals, freqvals);
                    assert(sum(indsfreq) == length(freqvals), 'have not previusly extracted all desired freqeunces...');
                    assert(all([filtstruct.freqvals(indsfreq)]==freqvals), 'not matcjed..');
                    
                    tmp = cellfun(@cell2mat, cellfun(@transpose, filtstruct.datfilt_chans(indschans), 'UniformOutput', 0), ...
                        'UniformOutput', 0);
                    
                    filtmat = nan(size(tmp{1},1), sum(indsfreq), length(tmp), 'single');
                    for j=1:length(tmp)
                        filtmat(:,:, j) = single(tmp{j}(:, indsfreq));
                    end
                    FILTMATall{ff} = filtmat;
                    tall = [tall ...
                        linspace(filtstruct.t_edges(1), filtstruct.t_edges(2), size(filtstruct.datfilt_chans{1}{1},1))];
                end
            end
            
            % =============== FOR EACH MOTIF, COLLECT ALL TRIALS FOR ALL
            % CHANNELS
            for mm=1:nummotifs
                disp(['motif ' num2str(mm)]);
                if isempty(dat.motif(mm).SegExtr_neurfakeID)
                    continue
                end
                
                sname = [savedir '/filtdat_bird' num2str(i) '_expt' num2str(ee) '_set' num2str(ss) '_mot' num2str(mm) '.mat'];
                if skipifdone ==1
                    if exist(sname, 'file')
                        disp(['SKIUPPING [aklrady done] - ' sname]);
                        continue
                    end
                end
                
                % --- segextract is shared across channels
                segextract = dat.motif(mm).SegExtr_neurfakeID(1).SegmentsExtract;
                
                if firstExtrAllSongs==0
                    % ================ COLLECT FILTERED DATA FOR ALL CHANNELS
                    [filtdat, t, freqsout, chansout] = lt_neural_QUICK_Segextr_GetFiltBank(...
                        segextract, dirname_main, Chanlist, motif_predur, extrapad, freqvals);
                elseif firstExtrAllSongs==1
                    % ================ USE PREVIUSLY EXTRACTED DATA (ALL SONGS)
                    % TO EXTRACT SPECIFIC TIMEWINDOWS
                    [filtdat, t] = lt_neural_QUICK_Segextr_GetFiltBankv2(...
                        segextract, motif_predur, extrapad, FILTMATall, songfiles, tall);
                end
                
                % ================ CONVERT TO A UNIVERSAL FORMAT
                
                
                % ============== COLLECT
                filtdatstruct = struct;
                filtdatstruct.filtdat_t_f_chan = filtdat;
                filtdatstruct.Chanlist = Chanlist;
                filtdatstruct.freqvals = freqvals;
                filtdatstruct.Bregionlist = Bregionlist;
                filtdatstruct.t_relons = t{1};
                
                if saveON==1
                    save(sname, 'filtdatstruct');
                end
                
            end
        end
        
        % ============= SAVE
    end
end