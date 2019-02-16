function [OUTSTRUCT, PARAMS] = lt_neural_Coher_RecalcMat(OUTSTRUCT, OUTSTRUCT_CohMatOnly, ...
    SwitchCohStruct, PARAMS, twind, fwind, wntouse, useWNtiming, ...
    prewind_relWN, COHSTRUCT)
%% lt 2/15/19 - uses new set of inds to recalcalate coherence matrices
% NOTE: saves old version as backup
Nmin = 15; % n trials... otherwise will ignore... [TO DO]

% == for getting scal rel wn
% - linearly interpolate coherence?
interpol =0; % across time  [NOT YET DONE].
onlyusegoodtargsyls = 1; % default = 1, i.e. gets timing only from that that will actually analyze.


%%

Cohmean_base = {};
Cohmean_wn = {};

cohscal_all = {};
cohscal_diff_all = [];

% ========= params for scalars
indT = PARAMS.tbins>=twind(1) & PARAMS.tbins<=twind(2);
indF = PARAMS.ffbins>=fwind(1) & PARAMS.ffbins<=fwind(2);

if length(PARAMS.tbins)~=size(OUTSTRUCT_CohMatOnly{1},1)
    disp('REEXTRACTING tbine: been overwritten with wrong one (wn realign...)');
    % reextract correct tbins
    % ================
    settmp = 1;
    if isempty(COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif)
        settmp=2;
        assert(~isempty(COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif));
    end
    PARAMS.tbins = COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif(1).t_relons;
end

assert(length(PARAMS.ffbins)==size(OUTSTRUCT_CohMatOnly{1},2), 'fbins not corrent unsure why..');

% ===============
for i=1:length(OUTSTRUCT.bnum)
    disp(i);
    
    indsbase = find(OUTSTRUCT.indsbase{i});
    indswn = find(OUTSTRUCT.indsWN{i});
    
    % ======================
    indsbase_epoch = indsbase(end-round(length(indsbase)/2):end);
    
    if strcmp(wntouse, 'half')
        indswn_epoch = indswn(end-round(length(indswn)/2):end);
    elseif strcmp(wntouse, 'quarter')
        indswn_epoch = indswn(end-round(length(indswn)/4):end);
    elseif strcmp(wntouse, 'third')
        indswn_epoch = indswn(end-round(length(indswn)/3):end);
    elseif strcmp(wntouse, 'firsthalf')
        indswn_epoch = indswn(1:round(length(indswn)/2));
    end
    
    % ====================== GET COHERENCE DATA
    bnum = OUTSTRUCT.bnum(i);
    enum = OUTSTRUCT.enum(i);
    sw = OUTSTRUCT.switch(i);
    mm = OUTSTRUCT.motifnum(i);
    
    %% =============== recalc coherence mat
    % ================== GET RAW COHMAT
    if (1)
        cohmat = OUTSTRUCT_CohMatOnly{i};
    else
        % === olde version, go back and load - tajkes a while since does
        % once for each case.
        chpair = OUTSTRUCT.chanpair(i,:);
        DAT = SwitchCohStruct.bird(bnum).exptnum(enum).switchlist(sw).motifnum(mm);
        fpre = DAT.fileprefix;
        fpost = DAT.filesuffix;
        
        tmp = ismember(DAT.chanpair, chpair, 'rows');
        chind = DAT.chanpairstokeep(tmp);
        cohmat = lt_neural_LFP_loadProcessDat([fpre '/Coh' fpost], chind);
    end
    
    % =================== GET PRE AND POST
    cohbase = mean(cohmat(:,:, indsbase_epoch), 3);
    cohwn = mean(cohmat(:,:, indswn_epoch), 3);
    
    
    % =========== SAVE IN OUTSTRUCT
    Cohmean_base = [Cohmean_base; cohbase];
    Cohmean_wn = [Cohmean_wn; cohwn];
    
    
    
    %% ========== recalc coherence scalar
    if useWNtiming==1
        % -- recalculate indT
        wnontime = OUTSTRUCT.WNonset(i);
        
        twindtmp = wnontime + prewind_relWN;
        indT = PARAMS.tbins>=twindtmp(1) & PARAMS.tbins<=twindtmp(2);
    end
    
    cohscal = squeeze(mean(mean(cohmat(indT, indF, :), 1),2));
    cohscal_base = mean(cohscal(indsbase_epoch));
    cohscal_wn = mean(cohscal(indswn_epoch));
    cohscal_diff = cohscal_wn - cohscal_base;
    
    % ====================== OUTPUT
    cohscal_all = [cohscal_all; cohscal'];
    cohscal_diff_all = [cohscal_diff_all; cohscal_diff];
    
    
    
end

% === save old versions
OUTSTRUCT.CohMean_Base_orig = OUTSTRUCT.CohMean_Base;
OUTSTRUCT.CohMean_WN_orig = OUTSTRUCT.CohMean_WN;

OUTSTRUCT.cohscal_orig = OUTSTRUCT.cohscal;
OUTSTRUCT.cohscal_diff_orig= OUTSTRUCT.cohscal_diff;

% ===== save new versions
OUTSTRUCT.CohMean_Base = Cohmean_base;
OUTSTRUCT.CohMean_WN = Cohmean_wn;

OUTSTRUCT.cohscal = cohscal_all;
OUTSTRUCT.cohscal_diff = cohscal_diff_all;


