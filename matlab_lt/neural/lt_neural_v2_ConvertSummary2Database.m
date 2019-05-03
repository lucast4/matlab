function [NeuronDatabase, SummaryStruct_filtered] = ...
    lt_neural_v2_ConvertSummary2Database(BirdsToKeep, BrainArea, ExptToKeep, ...
    RecordingDepth, LearningOnly, BatchesDesired, ChannelsDesired, ...
    extractpreDatenums, onlySpikes)
%%
if (0)
if onlySpikes==0
    disp('not sure what to do with wh72 - some chans multiple extractions, LFP vs spikes');
    pause;
end
disp('NOTE: if OnlySpikes=1, then this will not work for any birds that added after gr48bu5, since assumes that any bird other than gr48bu5 and wh72pk12 are all spike data (which is true)');
pause
end

% NOTE:
% For all birds did not enter what cell type it is, except following
% situation:
% wh72
%     - empty ==> is LFP
%     - "LF" ==> is LFP
%     - "MU"/"SU" ==> is spikes.
%     Note that there are some cases where same data extracted both spikes and LFP (e.g. spikes fewer trials).
%     To get only LFP in these cases, use onlySpikes=2;
% gr48 and after, explicitly noted down either LF, MU, or SU.
%     except RALMANLearn2 [is spikes, but forgot to note down MU or SU]...

    
%% lt 12/1/17 - for RA, added this, so extracts both SAm/Mel and my RA dat

%% TO DO:
% 1) Extract "notes" from summary struct


%% lt 3/6/17 - converts from summary structure to neuron database to use in analyses
% in conversion can filter to only extract certain experiments/birds, etc

% BirdsToKeep = {'bk7'}; % optional - if empty, then gets all birds.
% BirdsToKeep = {'bk7', []}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
% BrainArea = {'LMAN', 'X'}; % empty for all.
% ExptToKeep = {'LMANlearn2'}; % emopty for all;'
% RecordingDepth = [1800 1950] % in microns
% LearningOnly = 1; % then only if expt is in LearningMetastruct; 2: then
% excludes elarning expts [previous - only if summarysturct has WN on
% date][

%%
SummaryStruct =  lt_neural_v2_LoadSummary;
LearnStruct = lt_neural_v2_LoadLearnMetadat;


if ~exist('BirdsToKeep', 'var')
    BirdsToKeep = {};
end

if isempty(BirdsToKeep)
    BirdsToKeep = unique({SummaryStruct.birds.birdname}); % keep all birds
end

if isempty(BrainArea);
    %     BrainArea = {};
    %  DEFAULT IS TO ONLY KEEP IF IS SONG SYSTEMS
    BrainArea = {'HVC', 'RA', 'LMAN', 'X'};
    
end

if ~exist('BatchesDesired', 'var')
    BatchesDesired = {};
end


if ~exist('ExptToKeep', 'var')
    ExptToKeep = {};
end

if ~exist('RecordingDepth', 'var')
    RecordingDepth = [];
end

if ~exist('LearningOnly', 'var')
    LearningOnly = 0;
end

if ~exist('ChannelsDesired', 'var')
    ChannelsDesired = [];
end



%% =====

numbirds = length(SummaryStruct.birds);
NeuronDatabase = struct;
count = 1;

SummaryStruct_filtered = struct; % pulls out only the desired neurons
birdcounter = 1;
for i=1:numbirds
    
    if ~any(strcmp(SummaryStruct.birds(i).birdname, BirdsToKeep))
        disp(['skipped ' SummaryStruct.birds(i).birdname]);
        continue
    end
    
    
    numneurons = length(SummaryStruct.birds(i).neurons);
    birdname = SummaryStruct.birds(i).birdname;
    neuroncounter = 1;
    
    for ii=1:numneurons
        
        % ======== note that is not Sam/mel data
        %         SummaryStruct.birds(i).neurons(ii).isRAsobermel = 0;
        
        if ~isempty(BrainArea)
            if ~any(strcmp(SummaryStruct.birds(i).neurons(ii).NOTE_Location, BrainArea))
                disp(['skipped ' SummaryStruct.birds(i).birdname '; nueron ' num2str(ii) ' wrong brain area (' ...
                    SummaryStruct.birds(i).neurons(ii).NOTE_Location ')']);
                continue
            end
        end
        
        
        
        if ~isempty(ExptToKeep)
            if ~any(strcmp(SummaryStruct.birds(i).neurons(ii).exptID, ExptToKeep))
                disp(['skipped ' SummaryStruct.birds(i).birdname '; nueron ' num2str(ii) ' wrong experiment (' ...
                    SummaryStruct.birds(i).neurons(ii).exptID ')']);
                continue
            end
        end
        
        if ~isempty(RecordingDepth)
            if ~any(SummaryStruct.birds(i).neurons(ii).electrode_depth == RecordingDepth)
                disp(['skipped ' SummaryStruct.birds(i).birdname '; nueron ' num2str(ii) ' not correct depth (' ...
                    num2str(SummaryStruct.birds(i).neurons(ii).electrode_depth) ')']);
                continue
            end
        end
        
        if ~isempty(ChannelsDesired)
            if ~any(SummaryStruct.birds(i).neurons(ii).channel == ChannelsDesired)
                disp(['skipped ' SummaryStruct.birds(i).birdname '; nueron ' num2str(ii) ' not correct chan (' ...
                    num2str(SummaryStruct.birds(i).neurons(ii).channel) ')']);
                continue
            end
        end
        
        if ~isempty(BatchesDesired)
            if ~any(strcmp(SummaryStruct.birds(i).neurons(ii).batchfilename, BatchesDesired))
                disp(['skipped ' SummaryStruct.birds(i).birdname '; nueron ' num2str(ii) ' wrong batchname (' ...
                    SummaryStruct.birds(i).neurons(ii).batchfilename ')']);
                continue
            end
        end
        
        
        if (0) % OLD METHOD
            if LearningOnly==1
                if isempty(SummaryStruct.birds(i).neurons(ii).LEARN_WNonDatestr)
                    disp(['skipped ' SummaryStruct.birds(i).birdname '; nueron ' num2str(ii) ' (not learning)']);
                    continue
                end
            end
            
            if LearningOnly==2
                if ~isempty(SummaryStruct.birds(i).neurons(ii).LEARN_WNonDatestr)
                    disp(['skipped ' SummaryStruct.birds(i).birdname '; nueron ' num2str(ii) ' (is learning)']);
                    continue
                end
            end
            
        else
            % ==== is this learning expt/
            exptname = SummaryStruct.birds(i).neurons(ii).exptID;
            [islearning, LearnSummary, switchtime] = lt_neural_v2_QUICK_islearning(birdname, exptname, 1);
            
            %                 if max(SummaryStruct.birds(i).neurons(ii).Filedatenum_unsorted) < switchtime
            %                     % then latest song is before onset of WN
            %
            %                 end
            
            %                 birdindtmp = strcmp({LearnStruct.bird.birdname}, birdname);
            %                 if any(strcmp([LearnStruct.bird(birdindtmp).info(1,:)], exptname))
            %                     % then is learning
            %                     islearning =1;
            %                 else
            %                     islearning = 0;
            %                 end
            
            SummaryStruct.birds(i).neurons(ii).INFO_islearning = islearning;
            
            if LearningOnly==1
                % then only keeps if learniong
                if islearning==0
                    disp(['skipped ' SummaryStruct.birds(i).birdname '; nueron ' num2str(ii) ' (not learning)']);
                    continue
                end
            elseif LearningOnly==2
                % onl keeps if not learning
                if islearning==1
                    disp(['skipped ' SummaryStruct.birds(i).birdname '; nueron ' num2str(ii) ' (is learning)']);
                    continue
                end
            end
        end
        
        % ============== ONLY SPIKES?
        if onlySpikes==0 % THEN all data
            
        elseif onlySpikes==1 % THEN only spiking data
            
            if strcmp(SummaryStruct.birds(i).birdname, 'wh72pk12')
                % ---- if is "LFP, tjhen skip
                if strcmp(SummaryStruct.birds(i).neurons(ii).NOTE_PutativeCellType, 'LF')
                    continue
                    
                    % ---- if doesn't have entyr, then skip (since is LFP, I
                    % did not give entry earlier on, and those were all LFP)
                elseif isempty(SummaryStruct.birds(i).neurons(ii).NOTE_PutativeCellType)
                    continue
                    
                    
                    % --- if is MU or SU, then continue
                elseif ismember(SummaryStruct.birds(i).neurons(ii).NOTE_PutativeCellType, {'MU', 'SU'})
                    disp('is spikes, keeping');
                    % --- good
                else
                    disp('PROBLEM - what cell type is this?')
                    pause;
                end
                
            elseif strcmp(SummaryStruct.birds(i).birdname, 'gr48bu5')
                if strcmp(SummaryStruct.birds(i).neurons(ii).NOTE_PutativeCellType, 'LF')
                    continue
                elseif ismember(SummaryStruct.birds(i).neurons(ii).NOTE_PutativeCellType, {'MU', 'SU'})
                    disp('is spikes, keeping');
                elseif isempty(SummaryStruct.birds(i).neurons(ii).NOTE_PutativeCellType)
                    assert(strcmp(SummaryStruct.birds(i).neurons(ii).exptID, 'RALMANLearn2'));
                else
                    disp('PROBLEM - what cell type is this?')
                    asdfasdfasdfasdf;
                end
            else
                % ---- if is any bird other than wh72, then they are
                % allspikes..
            end
        elseif onlySpikes==2
            % This is special code that applies only for birds that have
            % multiple extractions for same channel, some with and some without
            % LFP. This will make sure that only keeps if is explicitly LFP
            % USE THIS IF WANT TO DO COHERE/XCOV ANALYSES.
            if strcmp(birdname, 'wh72pk12') & ...
                    ismember(SummaryStruct.birds(i).neurons(ii).NOTE_PutativeCellType, {'MU', 'SU'})
                continue
            end
        end
        
            
        tmpstruct = SummaryStruct.birds(i).neurons(ii);
        disp(['---- EXTRACTING TO NEURON DATABASE ... ' SummaryStruct.birds(i).birdname ...
            ' neuron ' num2str(ii)]);
        
        NeuronDatabase.neurons(count).birdname = birdname; %
        NeuronDatabase.neurons(count).exptID=tmpstruct.exptID; %
        NeuronDatabase.neurons(count).date=tmpstruct.date; % date
        NeuronDatabase.neurons(count).batchfile=tmpstruct.batchfilename; % batchfile (songs)
        NeuronDatabase.neurons(count).chan=tmpstruct.channel; % channel
        NeuronDatabase.neurons(count).clustnum=tmpstruct.clustnum; % cluster
        NeuronDatabase.neurons(count).electrode_depth=tmpstruct.electrode_depth; % cluster
        NeuronDatabase.neurons(count).NOTE_is_single_unit=tmpstruct.NOTE_is_single_unit; % is_single_unit
        NeuronDatabase.neurons(count).LEARN_WNonDatestr= tmpstruct.LEARN_WNonDatestr;
        NeuronDatabase.neurons(count).LEARN_WNotherImportantDates=tmpstruct.LEARN_WNotherImportantDates; % leave empty if nothing.
        NeuronDatabase.neurons(count).NOTE_random=tmpstruct.Notes{end}; % random note
        try
            NeuronDatabase.neurons(count).NOTE_Location=tmpstruct.NOTE_Location; % random note
        catch err
            NeuronDatabase.neurons(count).NOTE_Location=''; % random note
        end
        
        NeuronDatabase.neurons(count).NOTE_clust_qual_confirmed=''; % cluster_quality_confirmed
        NeuronDatabase.neurons(count).NOTE_all_songs_gotten=''; % all_songs_gotten
        NeuronDatabase.neurons(count).NOTE_all_labeled=''; %
        
        slashes = strfind(tmpstruct.dirname, '/');
        basedir = tmpstruct.dirname(1:slashes(end-1)-1);
        NeuronDatabase.neurons(count).basedir = basedir;
        
        count = count +1;
        
        % ====== keep in filtered summarystruct
        SummaryStruct_filtered.birds(birdcounter).neurons(neuroncounter) ...
            = SummaryStruct.birds(i).neurons(ii);
        
        SummaryStruct_filtered.birds(birdcounter).birdname = birdname;
        
        neuroncounter = neuroncounter +1;
        
    end
    birdcounter = birdcounter+1;
end

% ========== IF ANY BIRDS EMPTY, REMOVE
if isfield(SummaryStruct_filtered, 'birds')
    birdstoremove = [];
    for i=1:length(SummaryStruct_filtered.birds)
        if isempty(SummaryStruct_filtered.birds(i).neurons)
            birdstoremove = [birdstoremove i];
        end
    end
    SummaryStruct_filtered.birds(birdstoremove) = [];
    
    
    %% =========== POST INFO
    
    SummaryStruct_filtered = lt_neural_v2_PostInfo(SummaryStruct_filtered);
    
    %% ==== confirm all metadat have permanent datenum. if not, extract
    if extractpreDatenums==1;
        lt_neural_v2_PRE_datenums(SummaryStruct);
    end
end

%% ===== save current params to summary struct

SummaryStruct_filtered.loadparams.BirdsToKeep = BirdsToKeep;
SummaryStruct_filtered.loadparams.BrainArea = BrainArea;
SummaryStruct_filtered.loadparams.ExptToKeep = ExptToKeep;
SummaryStruct_filtered.loadparams.RecordingDepth = RecordingDepth;
SummaryStruct_filtered.loadparams.LearningOnly = LearningOnly;
SummaryStruct_filtered.loadparams.BatchesDesired = BatchesDesired;
SummaryStruct_filtered.loadparams.ChannelsDesired = ChannelsDesired;

tstamp = lt_get_timestamp(0);
SummaryStruct_filtered.loadparams.extractiontime = tstamp;


%% ========== if RA, then append sam/mel
% ======= first, extract sam/mel
if any(strcmp(BrainArea, 'RAmel'))
    SummaryStruct_SamMel = lt_neural_RASamMel_SummaryStruct;
    
    % ==== append
    birdcount = length(SummaryStruct_filtered.birds)+1;
    for i=1:length(SummaryStruct_SamMel.birds)
        
        %         if ~any(strcmp(SummaryStruct_SamMel.birds(i).birdname, BirdsToKeep))
        %             disp(['skipped ' SummaryStruct_SamMel.birds(i).birdname]);
        %             continue
        %         end
        
        SummaryStruct_filtered.birds(birdcount) = SummaryStruct_SamMel.birds(i);
        birdcount= birdcount+1;
    end
end

%% ========== find simulataneous neurons

SummaryStruct_filtered = lt_neural_PRE_GetSimultNeur(SummaryStruct_filtered);




