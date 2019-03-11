function OUTSTRUCT_lfp = lt_neural_POPLEARN_Cons_Extr_LFP(OUTSTRUCT, OUTSTRUCT_XCOV, SwitchStruct, ...
    SwitchCohStruct, MOTIFSTATS_Compiled, MOTIFSTATS_pop, SummaryStruct, ...
    PARAMS, onlykeepgoodcase, onlyadjacentpairs, corrwind, baseallinds, ...
    wnallinds, usehalfwn, lfpfilt)

% % ===================== SKIP? CAN CHOOSE TO ONLY KEEP
% % IF THIS IS MATCHED WITH WHAT IS IN CURRENT OUTSTRUCT.
% onlykeepgoodcase==1

% baseallinds = 0; % 0 is epocjh, 1 is all
% wnallinds = 0;

% === for computing correlations.
% corrwind = [-0.1 0.01]; % rel syl onset, to get pairwise trial by trial correlations.

%% lt 2/5/19 - across trial consistency [both LFP and SPIKES]

%% ################# [SPIKES]
%% [EXTRACT] go thru all switches. only extract if this is a desired neuron/motif.

OUTSTRUCT_lfp.bnum = [];
OUTSTRUCT_lfp.enum = [];
OUTSTRUCT_lfp.motifnum = [];
OUTSTRUCT_lfp.issame= [];
OUTSTRUCT_lfp.istarg= [];
OUTSTRUCT_lfp.learndirTarg= [];
OUTSTRUCT_lfp.switch= [];
OUTSTRUCT_lfp.motifname= {};
OUTSTRUCT_lfp.xtrialFrRho_BaseWn = {};
OUTSTRUCT_lfp.chanthis = [];
OUTSTRUCT_lfp.bregion = {};
                    OUTSTRUCT_lfp.indsWN_epoch = {};
                    OUTSTRUCT_lfp.indsBase_epoch = {};

numbirds = length(SwitchStruct.bird);
for i=1:numbirds
    numexpts = length(SwitchCohStruct.bird(i).exptnum);
    birdname = SwitchStruct.bird(i).birdname;
    for ii=1:numexpts
        numswitches = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
        for iii=1:numswitches
            if length(SwitchCohStruct.bird(i).exptnum(ii).switchlist)<iii
                continue
            end
            
            if isempty(SwitchCohStruct.bird(i).exptnum(ii).switchlist(iii).motifnum)
                continue
            end
            disp(num2str([i ii iii]));
            
            % ======= figure out which dataset to use
            DAT = SwitchCohStruct.bird(i).exptnum(ii).switchlist(iii);
            nummotifs = length(DAT.motifnum);
            
            
            % ==== which chans to get
            indstmp = OUTSTRUCT.bnum==i & OUTSTRUCT.enum==ii & ...
                OUTSTRUCT.switch==iii;
            chanlist = unique((OUTSTRUCT.chanpair(indstmp,:)));
            
            
            
            %             % ########################################## GO THRU ALL MOTIFS
            %             % ====================== 1) get liust of motifst ot plot
            %             % ===== then find out what are valid targ motifs
            %             indstmp = OUTSTRUCT_XCOV.bnum==i & OUTSTRUCT_XCOV.enum==ii ...
            %                 & OUTSTRUCT_XCOV.switch==iii & OUTSTRUCT_XCOV.istarg==1;
            %             motiflist = unique(OUTSTRUCT_XCOV.motifnum(indstmp));
            %
            % ################################ GO THRU EACH NEURON AND EACH MOTIF.
            for mm=1:nummotifs
                
                for cc=1:length(chanlist)
                    chanthis = chanlist(cc);
                    
                    %                motifname = DAT.motif(mm).regexpstr;
                    %                     motifpredur = MOTIFSTATS_Compiled.birds(1).MOTIFSTATS.params.motif_predur;
                    %                     bregionthis = bregionlist{nn};
                    
                    % -- what brain region is this?
                    indstmp = [SummaryStruct.birds(i).neurons.channel]==chanthis;
                    bregionthis = unique({SummaryStruct.birds(i).neurons(indstmp).NOTE_Location});
                    assert(length(bregionthis)==1, 'same number for actualyl different chans across diferernt days?');
                    bregionthis = bregionthis{1};
                    
                    % ===================== SKIP? CAN CHOOSE TO ONLY KEEP
                    % IF THIS IS MATCHED WITH WHAT IS IN CURRENT OUTSTRUCT.
                    if onlykeepgoodcase==1
                        %                         % == check for matching case
                        %                         indstmp = find(OUTSTRUCT_XCOV.bnum==i & ...
                        %                             OUTSTRUCT_XCOV.enum==ii & OUTSTRUCT_XCOV.switch==iii ...
                        %                             & OUTSTRUCT_XCOV.motifnum==mm & any(OUTSTRUCT_XCOV.neurpair'==nID')');
                        indstmp = OUTSTRUCT.bnum==i & OUTSTRUCT.enum==ii & ...
                            OUTSTRUCT.switch==iii & OUTSTRUCT.motifnum==mm & any(OUTSTRUCT.chanpair'==chanthis)';
                        if ~any(indstmp)
                            continue
                        end
                    end
                    
                    %% ===================== WHAT ARE TRIAL INDS TO USE?
                    indstmp = find(OUTSTRUCT.bnum==i & OUTSTRUCT.enum==ii & OUTSTRUCT.switch==iii ...
                        & OUTSTRUCT.motifnum==mm);
                    
                    if isempty(indstmp)
                        continue
                    end
                    
                    
                    assert(length(unique(cellfun(@(x)x(1), OUTSTRUCT.indsWN_epoch(indstmp))))==1, 'they are not identical..');
                    indstmp=indstmp(1); % assumes that multip[el cases are multip[el chans..
                    
                    if baseallinds==1
                        inds_base = find(OUTSTRUCT.indsbase{indstmp});
                    else
                        inds_base = OUTSTRUCT.indsbase_epoch{indstmp};
                    end
                    
                    if wnallinds==1
                        inds_wn = find(OUTSTRUCT.indsWN{indstmp});
                    else % WN epoch that was used for analsyis.
                        inds_wn = OUTSTRUCT.indsWN_epoch{indstmp};
                    end
                    
                    if usehalfwn==1
                        inds_wn = OUTSTRUCT.indsWN_epoch{indstmp};
                        inds_wn = inds_wn(round(length(inds_wn)/2):end);
                    end
                    
                    
                    % ========================== COLLECT OTHER THINGS FROM
                    % OUTSTRTUCT
                    istarg = OUTSTRUCT.istarg(indstmp);
                    issame = OUTSTRUCT.issame(indstmp);
                    learndirTarg = OUTSTRUCT.learndirTarg(indstmp);
                    motifname = OUTSTRUCT.motifname{indstmp};
                    
                    assert(length(DAT.motifnum(mm).tvals) == length(OUTSTRUCT.indsWN{indstmp}), 'not matcjhed!!');
                    
                    %% ====================== COMPUTE DAT
                    % ################### PAIRWISE CORRELATION OF FR
                    Ybasewn = cell(1,2); % base, wn, each trial p[air corrleations.
                    
                    % ================= 1) BASE
                    indsthis = inds_base;
                    yidx = 1;
                    
                    % ----------------- COLLECT
                    lfpthis = DAT.motifnum(mm).lfpall(:, DAT.motifnum(mm).lfpall_chans==chanthis);
                    % -- get correct trials
                    lfpthis = lfpthis(indsthis);
                    lfpthis = squeeze(lt_neural_Coher_Cell2Mat(lfpthis));
                    % -- filter LFP
                    lfpthis = lt_neural_filter(lfpthis, 1500, 0, lfpfilt(1), lfpfilt(2));
                    
                    % -- get correct within-trial time
                    x = DAT.motifnum(mm).t_lfp;
                    indx = x>=corrwind(1) & x<=corrwind(2);
                    lfpthis = lfpthis(indx,:);
                    
                    % --- GET ALL PAIRWISE CORRELATIONS (n vs. n+1);
                    rhoall = [];
                    if onlyadjacentpairs==1
                        for kk=1:size(lfpthis,2)-1
                            fr1 = lfpthis(:,kk);
                            fr2 = lfpthis(:,kk+1);
                            %                     rho =
                            rhoall = [rhoall; corr(fr1, fr2)];
                        end
                    else
                        rhomat = corr(lfpthis);
                        rhomat = triu(rhomat, 1);
                        rhoall = rhomat(rhomat(:)~=0);
                        assert(length(rhoall)==nchoosek(size(lfpthis,2), 2), 'mistake eitehr because some turned out to be 0, or duplicates');
                    end
                    Ybasewn{yidx} = rhoall;
                    
                    
                    % ================= 1) WN
                    indsthis = inds_wn;
                    yidx = 2;
                    
                    % ----------------- COLLECT
                    lfpthis = DAT.motifnum(mm).lfpall(:, DAT.motifnum(mm).lfpall_chans==chanthis);
                    % -- get correct trials
                    lfpthis = lfpthis(indsthis);
                    lfpthis = squeeze(lt_neural_Coher_Cell2Mat(lfpthis));
                    % -- filter LFP
                    lfpthis = lt_neural_filter(lfpthis, 1500, 0, lfpfilt(1), lfpfilt(2));
                    
                    % -- get correct within-trial time
                    x = DAT.motifnum(mm).t_lfp;
                    indx = x>=corrwind(1) & x<=corrwind(2);
                    lfpthis = lfpthis(indx,:);
                    
                    % --- GET ALL PAIRWISE CORRELATIONS (n vs. n+1);
                    rhoall = [];
                    if onlyadjacentpairs==1
                        for kk=1:size(lfpthis,2)-1
                            fr1 = lfpthis(:,kk);
                            fr2 = lfpthis(:,kk+1);
                            %                     rho =
                            rhoall = [rhoall; corr(fr1, fr2)];
                        end
                    else
                        rhomat = corr(lfpthis);
                        rhomat = triu(rhomat, 1);
                        rhoall = rhomat(rhomat(:)~=0);
                        assert(length(rhoall)==nchoosek(size(lfpthis,2), 2), 'mistake eitehr because some turned out to be 0, or duplicates');
                    end
                    Ybasewn{yidx} = rhoall;
                    
                    
                    
                    %% ========================= OUTPUT
                    OUTSTRUCT_lfp.chanthis = [OUTSTRUCT_lfp.chanthis; chanthis];
                    OUTSTRUCT_lfp.bregion = [OUTSTRUCT_lfp.bregion; bregionthis];
                    OUTSTRUCT_lfp.bnum = [OUTSTRUCT_lfp.bnum; i];
                    OUTSTRUCT_lfp.enum = [OUTSTRUCT_lfp.enum; ii];
                    OUTSTRUCT_lfp.motifnum = [OUTSTRUCT_lfp.motifnum; mm];
                    OUTSTRUCT_lfp.issame= [OUTSTRUCT_lfp.issame; issame];
                    OUTSTRUCT_lfp.istarg= [OUTSTRUCT_lfp.istarg; istarg];
                    OUTSTRUCT_lfp.learndirTarg= [OUTSTRUCT_lfp.learndirTarg; learndirTarg];
                    OUTSTRUCT_lfp.switch= [OUTSTRUCT_lfp.switch; iii];
                    OUTSTRUCT_lfp.motifname= [OUTSTRUCT_lfp.motifname; motifname];
                    OUTSTRUCT_lfp.xtrialFrRho_BaseWn = [OUTSTRUCT_lfp.xtrialFrRho_BaseWn; Ybasewn];
                    
                    OUTSTRUCT_lfp.indsWN_epoch = [OUTSTRUCT_lfp.indsWN_epoch; inds_wn];
                    OUTSTRUCT_lfp.indsBase_epoch = [OUTSTRUCT_lfp.indsBase_epoch; inds_base];
                    
                    
                end
                
            end
        end
    end
end