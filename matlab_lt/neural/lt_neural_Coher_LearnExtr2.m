function SwitchCohStruct = lt_neural_Coher_LearnExtr(COHSTRUCT, MOTIFSTATS_pop, SwitchStruct, pairtoget, ...
    LFPSTRUCT, PARAMS)
%% lt 10/12/18 - extract lerning related coherence dataset

% pairtoget = 'LMAN-RA';
removeifnan=0; % if 1, then runs into rpoblem of not saving some motfis. if 0 then incomoplete dat?


if exist('LFPSTRUCT', 'var')
    doLFP = 1;
else
    doLFP = 0;
end

%%

savedir = '/bluejay5/lucas/analyses/neural/LFP/PROCESSED';

%%
SwitchCohStruct = struct;
numbirds = length(SwitchStruct.bird);
for i=1:numbirds
    numexpts = length(SwitchStruct.bird(i).exptnum);
%     birdname = SwitchStruct.bird(i)
    for ii=1:numexpts
        numswitch = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        for ss=1:numswitch
            
            % =========== for this switch, get edge times and switch time
            tstart = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum_previous;
            tswitch = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum;
            tend = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum_next;
            
            % =========== Find channel pairs that have data both pre and post.
            numsets = length(MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons);
            setskept = [];
            
            
            for k=1:numsets
                %    COHSTRUCT.bird(i).experiment(ii).setnum(k).motif(1);
                nummotifs = length(MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(k).motif);
                
                
                for mm=1:nummotifs
                    if isempty(MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(k).motif(mm).SegExtr_neurfakeID)
                        continue
                    end
                    segextract = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(k).motif(mm).SegExtr_neurfakeID(1).SegmentsExtract;
                    
                    tvals = [segextract.song_datenum];
                    
                    % ---- ARE THERE BOTH PRE AND POST TVALS?
                    if ~any(tvals>tstart & tvals<tswitch) | ~any(tvals>tswitch & tvals<tend)
                        % then don't keep...
                        continue
                    end
                    
                    disp([num2str(i) '-' num2str(ii) '-' num2str(mm)]);
                    
                    %% =============== SAVE DATA FOR THIS MOTIF
                    cohdat = COHSTRUCT.bird(i).experiment(ii).setnum(k).motif(mm);
                    motifname = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(k).motif(mm).regexpstr;
                    if isempty(motifname)
                        keyboard
                    end
                    indtmp = strcmp(cohdat.bregionpairs_sorted, pairtoget);
                    indtmp_num = find(indtmp);
                    if ~any(indtmp)
                        continue
                    end
                    
%                     % ---- collect all coherence across all pairs
%                     numchanpairs = length(indtmp_num);
%                     if isfield(cohdat, 'Coh_ChpairByTrial')
%                         Coh_ChpairByTrial = cohdat.Coh_ChpairByTrial;
%                     else
%                         % -- have to load
%                         Coh_ChpairByTrial = load([savedir '/' PARAMS.savemarker '/Coh_bird' num2str(i) '_expt' num2str(ii) '_set' num2str(k) '_mot' num2str(mm) '.mat']);
%                         Coh_ChpairByTrial = Coh_ChpairByTrial.CohAllTrials;
%                     end
%                     tmp = size(Coh_ChpairByTrial{1});
%                     cohmatall = nan(tmp(1), tmp(2), length(Coh_ChpairByTrial), numchanpairs); % t, ff, trials, chanpairs
%                     for j=1:numchanpairs
%                         intmpthis = indtmp_num(j);
%                         cohmat = lt_neural_Coher_Cell2Mat(Coh_ChpairByTrial(intmpthis,:));
%                         cohmatall(:,:,:,j) = cohmat;
%                     end
%                     
                    
                    %                     cohmat = lt_neural_Coher_Cell2Mat(cohdat.Coh_ChpairByTrial(indtmp,:));
                    bregionpair = cohdat.bregionpairs_sorted(indtmp);
                    bregionpair_originalorder = cohdat.bregionpairs_unsorted(indtmp);
                    chanpair = cohdat.Chanpairs(indtmp,:);
                    tvals = [segextract.song_datenum];
                    ffvals = [segextract.FF_val];
                    
                    % ======================= WHAT WIL CALL BASELINE AND WN
                    % INDS?
                    indsbase = find(tvals>tstart & tvals<tswitch);
                    indsWN = find(tvals>tswitch & tvals<tend);
                    
                    % -------------- take second half of WN inds
                    indsbase = indsbase(round(length(indsbase)/2):end);
                    indsWN = indsWN(round(length(indsWN)/2):end);
                    
                    
                    %% ==================== EXTRACT LFP
                    if doLFP==1
                        lfpdat = LFPSTRUCT.bird(i).experiment(ii).setnum(k).motif(mm);
                        
                        % ----- which channels to get?
                        chanstoget = unique(chanpair(:)');
                        
                        % ----- get those chans
                        indstmp = ismember(lfpdat.Chanlist, chanstoget);
                        lfpall = lfpdat.LFP_chanbytrial(indstmp,:)';
                        lfpall_chans = lfpdat.Chanlist(indstmp)';
                        
                        t_lfp = lfpdat.t_relons;
                        
                        % ==================== SAVE OUTPUT
                        SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).lfpall = lfpall;
                        SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).t_lfp = t_lfp;
                        SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).lfpall_chans = lfpall_chans;
                    end
                    
                    %% ====================== SAVE OUTPUT
                    SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).neursetused = k;
                    SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).motifname = motifname;
                    SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).fileprefix = ...
                        [savedir '/' PARAMS.savemarker];
                    SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).chanpairstokeep = indtmp_num;
                    SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).filesuffix = ...
                        ['_bird' num2str(i) '_expt' num2str(ii) '_set' num2str(k) '_mot' num2str(mm) '.mat'];
%                     SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).cohmat = cohmatall;
                    SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).tvals = tvals;
                    SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).ffvals = ffvals;
                    SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).bregionpair = bregionpair;
                    SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).bregionpair_originalorder = ...
                        bregionpair_originalorder;
                    SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).chanpair = chanpair;
                    SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).indsbase_epoch = indsbase;
                    SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).indsWN_epoch = indsWN;
                    setskept = [setskept k];
                    
                end
            end
            %             try
            %             if any(isempty({SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm).motifname}))
            %                 keyboard
            %             end
            %             catch err
            %             end
            % --- sanity check, confirm that for any switch at most only
            % one set of neurons is included.
            if ~isempty(setskept)
                assert(length(unique(setskept))==1, 'if this fails, then must combine across sets...');
            end
        end
    end
end