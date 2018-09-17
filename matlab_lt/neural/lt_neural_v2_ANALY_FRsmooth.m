function OUTDAT = lt_neural_v2_ANALY_FRsmooth(MOTIFSTATS_Compiled, SwitchStruct, onlyFirstSwitch, ...
    onlyIfSameTarg, BirdsToPlot, BrainLocation, throwoutlonggap)

%% lt - 9/6/18, smoothed FR, average of low and high FF trials, compare to change during learning

singledayonly = 1; % if 1 then gets only data on the day of the switch...
pretime = 0.15; % how much to collect, relative to alignemnt (sec); leave empty to get all.
posttime = 0.1;
minN = 20; % min num trials, both base and train


%% ================= params to extract
motif_predur = MOTIFSTATS_Compiled.birds(1).exptnum(1).MOTIFSTATS.params.motif_predur;
motif_postdur = MOTIFSTATS_Compiled.birds(1).exptnum(1).MOTIFSTATS.params.motif_postdur;

if isempty(pretime)
    pretime = motif_predur;
end
if isempty(posttime)
    posttime = motif_postdur;
end


%% ============== FOR EACH MOTIF/NEURON, GET SMOOTHED FR IN SEPARATE FF BINS


% ============================== INITIATE OUTPUT HOLDERS
OUTDAT.All_birdnum = [];
OUTDAT.All_exptnum = [];
OUTDAT.All_swnum = [];
OUTDAT.All_istarg = [];
OUTDAT.All_issame = [];
OUTDAT.All_neurnum = [];
OUTDAT.All_motifnum = [];

OUTDAT.All_FRsmooth = {};
OUTDAT.All_FRsmooth_t = {};
OUTDAT.All_FF = {};
OUTDAT.All_FF_t = {};

% ============================= RUN THRU ALL DATA
numbirds = length(SwitchStruct.bird);
for i=1:numbirds
    
    bname = SwitchStruct.bird(i).birdname;
    
    if ~isempty(BirdsToPlot)
        if ~any(strcmp(BirdsToPlot, bname))
            continue
        end
    end
    
    numexpt = length(SwitchStruct.bird(i).exptnum);
    for ii=1:numexpt
        
        
        ename = SwitchStruct.bird(i).exptnum(ii).exptname;
        
        numswitch = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        
        for ss =1:numswitch
            
            if onlyFirstSwitch==1
                if ss>1
                    continue
                end
            end
            
            targsyls = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningDirs(1:2:end);
            
            if onlyIfSameTarg==1
                if SwitchStruct.bird(i).exptnum(ii).switchlist(ss).targsAreSameSyl~=1
                    continue
                end
            end
            
            % =========================== FOR THIS SWITCH, COLLECT SMOOTHED FR
            neuronlist  = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).goodneurons;
            neuronlist_origID = find(MOTIFSTATS_Compiled.birds(i).exptnum(ii).neurIDOriginal_inorder);
            motiflist = {SwitchStruct.bird(i).exptnum(ii).switchlist(ss).STATS_motifsyl.sylname};
            
            sametypelist = lt_neural_QUICK_ExtractSameType(motiflist, targsyls);
            
            for jjj = 1:length(neuronlist)
                nn = neuronlist(jjj);
                nn_orig = neuronlist_origID(nn);
                
                bregionthis = MOTIFSTATS_Compiled.SummaryStruct.birds(i).neurons(nn_orig).NOTE_Location;
                if ~isempty(BrainLocation)
                    if ~any(strcmp(BrainLocation, bregionthis))
                        continue
                    end
                end
                
                for mm=1:length(motiflist)
                    
                    
                    % ----------------- STATS FOR THIS SYL
                    motifthis = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif_regexpr_str{mm};
                    assert(strcmp(motifthis, motiflist{mm}), 'not expected motif, based on switch struct common motifs');
                    istarg = any(strcmp(targsyls, motifthis));
                    issame = any(strcmp(sametypelist, motifthis)) | istarg;
                    
                    % ========================== REMOVE SYLS THAT ARE BAD
                    % (I.E. POST TARGET SYL)
                    sylbad = lt_neural_QUICK_LearnRemoveBadSyl(bname, ename, ss, motifthis);
                    if sylbad==1
                        disp(['Removing syl, bad syl -----' bname '-' ename '-s' num2str(ss) '-' motifthis]);
                    end
                    
                    
                    if length(SwitchStruct.bird(i).exptnum(ii).switchlist(ss).neuron(nn).DATA.motif) < mm
                        continue
                    end
                    
                    % ------- collect the trial inds for this syl/motif combo
                    if singledayonly==1
                        inds_base = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).neuron(nn).DATA.motif(mm).baseInds_WithinDayOfSw;
                        inds_train = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).neuron(nn).DATA.motif(mm).trainInds_WithinDayOfSw;
                    else
                        inds_base = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).neuron(nn).DATA.motif(mm).baseInds;
                        inds_train = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).neuron(nn).DATA.motif(mm).trainInds;
                    end
                    
                    % ---- skip if not enough data
                    if sum(inds_base)<minN | sum(inds_train)<minN
                        disp('skip, not enouguh data');
                        continue
                    end
                    
                    % ------ collect the data
                    segextract = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(mm).SegmentsExtract;
                    assert(length(segextract) == size([inds_base; inds_train],2));
                    if isempty(segextract)
                        continue
                    end
                    assert(all(diff([segextract.song_datenum])>=0), 'need to be in temporal order...');
                    
                    % ====================== COLLECT SMOOTHED FR
                    FRateall = cell(1, 2); % base, train
                    FRtbinall = cell(1,2);
                    FFall = cell(1,2); %
                    Tall = cell(1,2);
                    
                    % --------- BASE
                    indstmp = inds_base;
                    
                    fratemat = [segextract(indstmp).FRsmooth_rate_CommonTrialDur];
                    t = segextract(1).FRsmooth_xbin_CommonTrialDur;
                    [fratemat, t] = lt_neural_QUICK_GetFRmat(fratemat, t, motif_predur, pretime, posttime);
                    
                    ffthis =  [segextract(indstmp).FF_val];
                    tthis = [segextract(indstmp).song_datenum];
                    
                    FRateall{1} = fratemat;
                    FRtbinall{1} = t;
                    FFall{1} = ffthis;
                    Tall{1} = tthis;
                    tlast = tthis(end);
                    
                    % --------- TRAIN
                    indstmp = inds_train;
                    
                    fratemat = [segextract(indstmp).FRsmooth_rate_CommonTrialDur];
                    t = segextract(1).FRsmooth_xbin_CommonTrialDur;
                    [fratemat, t] = lt_neural_QUICK_GetFRmat(fratemat, t, motif_predur, pretime, posttime);
                    
                    ffthis =  [segextract(indstmp).FF_val];
                    tthis = [segextract(indstmp).song_datenum];
                    
                    FRateall{2} = fratemat;
                    FRtbinall{2} = t;
                    FFall{2} = ffthis;
                    Tall{2} = tthis;
                    tfirst = tthis(1);
                    
                    
                    % ========= make sure gap between base and train is not
                    % too big
                    if throwoutlonggap==1
                    if (tfirst - tlast)>1/24 % one hour
                        disp('skip, gap between train and base too large (>1hr)');
                        continue
                    end
                    end
                    
                    % ======== SANITY CEHCK, PLOT FOR THIS NEURON/MOTIF
                    if (0)
                        
                        lt_figure; hold on;
                        
                        % -- base
                        fratethis = mean(FRateall{1}, 2);
                        plot(fratethis, 'r')
                        
                        
                    end
                    
                    
                    % ======================== APPEND TO OUTPUT
                    OUTDAT.All_birdnum = [OUTDAT.All_birdnum; i];
                    OUTDAT.All_exptnum = [OUTDAT.All_exptnum; ii];
                    OUTDAT.All_swnum = [OUTDAT.All_swnum; ss];
                    OUTDAT.All_istarg = [OUTDAT.All_istarg; istarg];
                    OUTDAT.All_issame = [OUTDAT.All_issame; issame];
                    OUTDAT.All_neurnum = [OUTDAT.All_neurnum; nn];
                    OUTDAT.All_motifnum = [OUTDAT.All_motifnum; mm];
                    OUTDAT.All_FRsmooth = [OUTDAT.All_FRsmooth; FRateall];
                    OUTDAT.All_FRsmooth_t = [OUTDAT.All_FRsmooth_t; FRtbinall];
                    OUTDAT.All_FF = [OUTDAT.All_FF; FFall];
                    OUTDAT.All_FF_t = [OUTDAT.All_FF_t; Tall];
                    
                end
                
            end
        end
    end
end


%% make durations similar
% ======================
ntbins = min(cellfun(@numel, OUTDAT.All_FRsmooth_t(:))); % min numbner 0f time bins

for j=1:length(OUTDAT.All_FRsmooth)
    OUTDAT.All_FRsmooth{j,1} = OUTDAT.All_FRsmooth{j,1}(1:ntbins, :);
    OUTDAT.All_FRsmooth{j,2} = OUTDAT.All_FRsmooth{j,2}(1:ntbins, :);
    
    OUTDAT.All_FRsmooth_t{j,1} = OUTDAT.All_FRsmooth_t{j,1}(1:ntbins);
    OUTDAT.All_FRsmooth_t{j,2} = OUTDAT.All_FRsmooth_t{j,2}(1:ntbins);
end





