function NGRAMSTRUCT = lt_neural_NGRAMS_Extract(SummaryStruct, Params, saveON)
%% lt 4/23/18 - extract for each neuron all ngrams and summary stats
% NOTE:
% ALSO EXTRACT STATS (E.G. FR TRACE) FOR EACH NGRAM.
% WILL SKIP IF SAMPLE SIZE TOO LOW
% throws out segmentsextract to reduce size of struct.


%%
% LearnKeepOnlyBase = 1;
% strtype = 'xaa';
% Nmin = 7; % num trials minimum
% alignsyl = 2;
%
% Params.regexpr.motifpredur = 0.05;
% Params.regexpr.motifpostdur = 0.1;
% Params.regexpr.alignOnset = 1;
% Params.regexpr.preAndPostDurRelSameTimept=1;
% Params.regexpr.RemoveIfTooLongGapDur = 1;

%% SAVE?
if saveON==1
    tstamp = lt_get_timestamp(0);
    savedir = '/bluejay5/lucas/analyses/neural/NGRAMS';
    savedir = [savedir '/' Params.strtype '_' tstamp];
    mkdir(savedir);
    %
    %     save(savefname, 'NGRAMSTRUCT', '-v7.3');
end

if saveON==1
    % --- save params
    save([savedir '/SummaryStruct.mat'], 'SummaryStruct');
    save([savedir '/Params.mat'], 'Params');
end


%% params

LearnKeepOnlyBase = Params.LearnKeepOnlyBase;
strtype = Params.strtype;
Nmin = Params.Nmin;
alignsyl = Params.alignsyl;


%% what regexp string will be used to find ngrams?

switch length(strtype)
    case 2
        searchstring = '[a-z](?=[a-z])';
    case 3
        searchstring = '[a-z](?=[a-z][a-z])';
    case 4
        searchstring = '[a-z](?=[a-z][a-z][a-z])';
end

%% ============= GET LIST OF ALL N-GRAMS THAT EXIST FOR THIS NEURON
% ALSO EXTRACT STATS (E.G. FR TRACE) FOR EACH NGRAM.
% WILL SKIP IF SAMPLE SIZE TOO LOW


NGRAMSTRUCT = struct;
NGRAMSTRUCT.params.strtype = strtype;
NGRAMSTRUCT.params.searchstring = searchstring;
NGRAMSTRUCT.params.Params = Params;

numbirds = length(SummaryStruct.birds);
for i=1:numbirds
    numneurons = length(SummaryStruct.birds(i).neurons);
    birdname = SummaryStruct.birds(i).birdname;
    
    for ii=1:numneurons
        
        
        % ============== extract metadat
        [SongDat, NeurDat, Prmstmp] = lt_neural_ExtractDat2(SummaryStruct, i, ii);
        clustnum = SummaryStruct.birds(i).neurons(ii).clustnum;
        
        
        % ============== only keep songs before start of learning
        if LearnKeepOnlyBase==1 & ~isfield(SummaryStruct.birds(i).neurons(ii), 'isRAsobermel')
            exptname = SummaryStruct.birds(i).neurons(ii).exptID;
            [islearning, ~, switchtime] = ...
                lt_neural_v2_QUICK_islearning(birdname, exptname, 1);
            
            if islearning==1
                % --- pare down labels to just those before learning start.
                songstokeep = find([NeurDat.metaDat.song_datenum]'<switchtime);
                
                if isempty(songstokeep)
                    % ---- then no data for this neuron, skip
                    continue
                end
                
                rendstokeep = find(SongDat.AllSongNum <= max(songstokeep));
                SongDat = lt_structure_subsample_all_fields(SongDat, rendstokeep);
            end
        end
        
        
        % ============= get list of all ngrams
        [start] = regexp(SongDat.AllLabels, searchstring, 'start');
        %             lt_neural_QUICK_regexp(SongDat.AllLabels, regexpcontrol)
        
        % - get match syls
        strlength = length(strtype);
        indmat = [];
        for j=1:strlength
            indmat = [indmat start'+j-1];
        end
        allmotifs = tabulate(SongDat.AllLabels(indmat));
        
        % ------ only keep if sample size large enough
        allmotifs = allmotifs(cell2mat(allmotifs(:,2))>=Nmin, :);
        
        
        % ============ SAVE ALLMOTIFS
        NGRAMSTRUCT.bird(i).neuron(ii).ngramlist = allmotifs(:,1);
        NGRAMSTRUCT.bird(i).neuron(ii).ngramN = cell2mat(allmotifs(:,2));
        
        
        
        % #############################################################
        % ########## FOR EACH NGRAM EXTRACT SUMMARY DATA
        
        for m=1:size(allmotifs,1)
            
            % ----------- 1) what position to put the token (i.e. what syl
            % to align to)
            regexprstr = allmotifs{m,1};
            regexprstr = [regexprstr(1:alignsyl-1) '(' ...
                regexprstr(alignsyl) ')' regexprstr(alignsyl+1:end)];
            
            % ---------- struct
            NGRAMSTRUCT.bird(i).neuron(ii).ngramnum(m).DAT = [];
            NGRAMSTRUCT.bird(i).neuron(ii).ngramnum(m).regexprstr = regexprstr;
            
            % ---------- 2) EXtRACT DAT
            FFparams.collectFF=1;
            extractDirSong=0;
            %             keepRawNeuralDat=0;
            [SegmentsExtract, Prmstmp]=lt_neural_RegExp(SongDat, NeurDat, ...
                Prmstmp, regexprstr, Params.regexpr.motifpredur, ...
                Params.regexpr.motifpostdur, Params.regexpr.alignOnset, ...
                '', FFparams, 0, 1, 0, 0, LearnKeepOnlyBase, ...
                Params.regexpr.preAndPostDurRelSameTimept, ...
                Params.regexpr.RemoveIfTooLongGapDur, clustnum, extractDirSong);
            

            % =================== throw out this motif if there is not
            % enough trials. % usually mismatch from regexp above occurs
            % becuase of the added restriction of gap durations in
            % lt_neural_RegExp [or throw out if empty]
            if length(SegmentsExtract) < Params.Nmin
                continue
            end
            
            % ============================ EXTRACT THINGS
            % ++++++++++++++++++++++++++++++ ENTIRE SEGETRACT
            if (0) % don't keep as filesize too large. take out things I want instead
                NGRAMSTRUCT.bird(i).neuron(ii).ngramnum(m).SegmentsExtract = SegmentsExtract;
            end
            
            
            % ++++++++++++++++++++++++++++++ take frmat
            SegmentsExtract = lt_neural_SmoothFR(SegmentsExtract, clustnum);
            frmat = [SegmentsExtract.FRsmooth_rate_CommonTrialDur];
            % make duration consistent
            if Params.regexpr.preAndPostDurRelSameTimept==1
                % then expect all syls to have aligned data timebase
                maxsamps = round(1000*(Params.regexpr.motifpredur + Params.regexpr.motifpostdur)-1);
                frmat = frmat(1:maxsamps,:);
            end
            NGRAMSTRUCT.bird(i).neuron(ii).ngramnum(m).DAT.frmat = frmat;
            
            
            % ++++++++++++++++++++++++++++++ syl timing
            motifsylOn = cell2mat({SegmentsExtract.motifsylOnsets}');
            motifsylOff = cell2mat({SegmentsExtract.motifsylOffsets}');
            
            NGRAMSTRUCT.bird(i).neuron(ii).ngramnum(m).DAT.motifsylOn = motifsylOn;
            NGRAMSTRUCT.bird(i).neuron(ii).ngramnum(m).DAT.motifsylOff = motifsylOff;
           
            
            % ++++++++++++++++++++++++++++++ FF
            ffvals = [SegmentsExtract.FF_val];
            if isfield(SegmentsExtract, 'song_datenum')
                tvals = [SegmentsExtract.song_datenum];
            else
                tvals = nan(size(ffvals));
            end
            
            NGRAMSTRUCT.bird(i).neuron(ii).ngramnum(m).DAT.ffvals = ffvals';
            NGRAMSTRUCT.bird(i).neuron(ii).ngramnum(m).DAT.tvals = tvals';
            
            % ++++++++++++++++++++++++++++++ 1) Smoothed FR
            if (0) % don't do this, since now I am extracting frmat, so can get smoothed fr easiyl
                SegmentsExtract = lt_neural_SmoothFR(SegmentsExtract, clustnum);
                frmean = mean([SegmentsExtract.FRsmooth_rate_CommonTrialDur],2);
                frx = SegmentsExtract(1).FRsmooth_xbin_CommonTrialDur;
                
                NGRAMSTRUCT.bird(i).neuron(ii).ngramnum(m).DAT.frmean = frmean;
                NGRAMSTRUCT.bird(i).neuron(ii).ngramnum(m).DAT.frx = frx;
            end
        end
    end
    
    % ====================== SAVE EACH BIRD'S DATA INDIVIDUALLY
    if saveON==1
        savefname = [savedir '/bird' num2str(i) '.mat'];
        birdstruct = NGRAMSTRUCT.bird(i);
        save(savefname, 'birdstruct');
    end
end

%%
if (0) % old version, now don't save struct since too much data.
    if saveON==1
        tstamp = lt_get_timestamp(0);
        savedir = '/bluejay5/lucas/analyses/neural/NGRAMS';
        savefname = [savedir '/NGRAMSTRUCT_' tstamp '.mat'];
        save(savefname, 'NGRAMSTRUCT', '-v7.3');
    end
end

if saveON==1
    % --- save params
    save([savedir '/SummaryStruct.mat'], 'SummaryStruct');
    save([savedir '/Params.mat'], 'Params');
end



