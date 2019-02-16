function OUTSTRUCT_units = lt_neural_POPLEARN_Cons_Extr(OUTSTRUCT, OUTSTRUCT_XCOV, SwitchStruct, ...
    SwitchCohStruct, MOTIFSTATS_Compiled, MOTIFSTATS_pop, SummaryStruct, ...
    PARAMS, onlykeepgoodcase, onlyadjacentpairs, corrwind, baseallinds, ...
    wnallinds, usehalfwn, usesameindsasXcov)

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

OUTSTRUCT_units.bnum = [];
OUTSTRUCT_units.enum = []; 
OUTSTRUCT_units.motifnum = []; 
OUTSTRUCT_units.issame= []; 
OUTSTRUCT_units.istarg= []; 
OUTSTRUCT_units.learndirTarg= []; 
OUTSTRUCT_units.switch= []; 
OUTSTRUCT_units.motifname= {};
OUTSTRUCT_units.xtrialFrRho_BaseWn = {};
OUTSTRUCT_units.neurID = [];
OUTSTRUCT_units.bregion = {};

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
            neurset = SwitchCohStruct.bird(i).exptnum(ii).switchlist(iii).motifnum(1).neursetused;
            
            % ========= get segextract for the neuron pair (using the correct neuron set)
            DAT = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(neurset);
            neurlist = MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons{neurset};
            bregionlist = {SummaryStruct.birds(i).neurons(neurlist).NOTE_Location};
            
            numneur = length(neurlist);
            nummotifs = length(DAT.motif);
            
            
            
            %             % ########################################## GO THRU ALL MOTIFS
            %             % ====================== 1) get liust of motifst ot plot
            %             % ===== then find out what are valid targ motifs
            %             indstmp = OUTSTRUCT_XCOV.bnum==i & OUTSTRUCT_XCOV.enum==ii ...
            %                 & OUTSTRUCT_XCOV.switch==iii & OUTSTRUCT_XCOV.istarg==1;
            %             motiflist = unique(OUTSTRUCT_XCOV.motifnum(indstmp));
            %
            % ################################ GO THRU EACH NEURON AND EACH MOTIF.
            for mm=1:nummotifs
                if isempty(DAT.motif(mm).SegExtr_neurfakeID)
                    continue
                end
                for nn=1:numneur
                    
                    nID = DAT.motif(mm).SegExtr_neurfakeID(nn).neurID_orig;
                    assert(nID==neurlist(nn), 'asdfd');
                    segthis = DAT.motif(mm).SegExtr_neurfakeID(nn).SegmentsExtract;
                    seg_global = DAT.motif(mm).SegExtr_neurfakeID(1).SegmentsExtract;
                    motifname = DAT.motif(mm).regexpstr;
                    motifpredur = MOTIFSTATS_Compiled.birds(1).MOTIFSTATS.params.motif_predur;
                    bregionthis = bregionlist{nn};
                    
                    % ===================== SKIP? CAN CHOOSE TO ONLY KEEP
                    % IF THIS IS MATCHED WITH WHAT IS IN CURRENT OUTSTRUCT.
                    if onlykeepgoodcase==1
                        % == check for matching case
                        indstmp = find(OUTSTRUCT_XCOV.bnum==i & ...
                            OUTSTRUCT_XCOV.enum==ii & OUTSTRUCT_XCOV.switch==iii ...
                            & OUTSTRUCT_XCOV.motifnum==mm & any(OUTSTRUCT_XCOV.neurpair'==nID')');
                        if isempty(indstmp)
                            continue
                        end
                    end
                    
                    %% ===================== WHAT ARE TRIAL INDS TO USE?
                    indstmp = find(OUTSTRUCT.bnum==i & OUTSTRUCT.enum==ii & OUTSTRUCT.switch==iii ...
                        & OUTSTRUCT.motifnum==mm);
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
                    
                    if usesameindsasXcov==1
                        indstmp = find(OUTSTRUCT_XCOV.bnum==i & ...
                            OUTSTRUCT_XCOV.enum==ii & OUTSTRUCT_XCOV.switch==iii ...
                            & OUTSTRUCT_XCOV.motifnum==mm);
                        
                        inds_base = OUTSTRUCT_XCOV.inds_base_epoch{indstmp(1)};
                        inds_wn = OUTSTRUCT_XCOV.inds_WN_epoch{indstmp(1)};
                    end
                    
                    % ========================== COLLECT OTHER THINGS FROM
                    % OUTSTRTUCT
                    indstmp = find(OUTSTRUCT.bnum==i & OUTSTRUCT.enum==ii & OUTSTRUCT.switch==iii ...
                        & OUTSTRUCT.motifnum==mm);
                    assert(length(unique(cellfun(@(x)x(1), OUTSTRUCT.indsWN_epoch(indstmp))))==1, 'they are not identical..');
                    indstmp=indstmp(1); % assumes that multip[el cases are multip[el chans..

                    istarg = OUTSTRUCT.istarg(indstmp);
                    issame = OUTSTRUCT.issame(indstmp);
                    learndirTarg = OUTSTRUCT.learndirTarg(indstmp);
                    motifname = OUTSTRUCT.motifname{indstmp};
                    
                    assert(length(segthis) == length(OUTSTRUCT.indsWN{indstmp}), 'not matcjhed!!');
                    
                    %% ====================== COMPUTE DAT
                    % ==== 1) get smoothed FR
                    segthis = lt_neural_SmoothFR(segthis, [], [], [], 0, seg_global);
                    
                    
                    % ################### PAIRWISE CORRELATION OF FR
                    Ybasewn = cell(1,2); % base, wn, each trial p[air corrleations.
                    
                    
                    % ================= 1) BASE
                    indsthis = inds_base;
                    yidx = 1;
                    
                    % ----------------- COLLECT
                    frmat = [segthis(indsthis).FRsmooth_rate_CommonTrialDur];
                    frx = segthis(indsthis(1)).FRsmooth_xbin_CommonTrialDur - motifpredur;
                    % --- only keep within time bin
                    indx = frx>=corrwind(1) & frx<=corrwind(2);
                    frmat = frmat(indx,:);
                    % --- GET ALL PAIRWISE CORRELATIONS (n vs. n+1);
                    rhoall = [];
                    if onlyadjacentpairs==1
                        for kk=1:size(frmat,2)-1
                            fr1 = frmat(:,kk);
                            fr2 = frmat(:,kk+1);
                            %                     rho =
                            rhoall = [rhoall; corr(fr1, fr2)];
                        end
                    else
                        rhomat = corr(frmat);
                        rhomat = triu(rhomat, 1);
                        rhoall = rhomat(rhomat(:)~=0);
                        assert(length(rhoall)==nchoosek(size(frmat,2), 2), 'mistake eitehr because some turned out to be 0, or duplicates');
                    end
                    Ybasewn{yidx} = rhoall;
                    
                    
                    % ================= 1) WN
                    indsthis = inds_wn;
                    yidx = 2;
                    
                    % ----------------- COLLECT
                    frmat = [segthis(indsthis).FRsmooth_rate_CommonTrialDur];
                    frx = segthis(indsthis(1)).FRsmooth_xbin_CommonTrialDur - motifpredur;
                    % --- only keep within time bin
                    indx = frx>=corrwind(1) & frx<=corrwind(2);
                    frmat = frmat(indx,:);
                    % --- GET ALL PAIRWISE CORRELATIONS (n vs. n+1);
                    rhoall = [];
                    if onlyadjacentpairs==1
                        for kk=1:size(frmat,2)-1
                            fr1 = frmat(:,kk);
                            fr2 = frmat(:,kk+1);
                            %                     rho =
                            rhoall = [rhoall; corr(fr1, fr2)];
                        end
                    else
                        rhomat = corr(frmat);
                        rhomat = triu(rhomat, 1);
                        rhoall = rhomat(rhomat(:)~=0);
                        assert(length(rhoall)==nchoosek(size(frmat,2), 2), 'mistake eitehr because some turned out to be 0, or duplicates');
                    end
                    Ybasewn{yidx} = rhoall;
                    
                    
                    
                    %% ========================= OUTPUT
                    OUTSTRUCT_units.neurID = [OUTSTRUCT_units.neurID; nID];
                    OUTSTRUCT_units.bregion = [OUTSTRUCT_units.bregion; bregionthis];
                    OUTSTRUCT_units.bnum = [OUTSTRUCT_units.bnum; i];
                    OUTSTRUCT_units.enum = [OUTSTRUCT_units.enum; ii]; 
                    OUTSTRUCT_units.motifnum = [OUTSTRUCT_units.motifnum; mm]; 
                    OUTSTRUCT_units.issame= [OUTSTRUCT_units.issame; issame]; 
                    OUTSTRUCT_units.istarg= [OUTSTRUCT_units.istarg; istarg]; 
                    OUTSTRUCT_units.learndirTarg= [OUTSTRUCT_units.learndirTarg; learndirTarg]; 
                    OUTSTRUCT_units.switch= [OUTSTRUCT_units.switch; iii]; 
                    OUTSTRUCT_units.motifname= [OUTSTRUCT_units.motifname; motifname];
                    OUTSTRUCT_units.xtrialFrRho_BaseWn = [OUTSTRUCT_units.xtrialFrRho_BaseWn; Ybasewn];

                    
                end
                
            end
        end
    end
end
