function [OUTSTRUCT, PARAMS] = lt_neural_Coher_RecalcScal(OUTSTRUCT, PARAMS, twind, ...
    fwind)

% NOTE: saves old version as backup
Nmin = 15; % n trials... otherwise will ignore... [TO DO]

%%

% Cohmean_base = {};
% Cohmean_wn = {};

for i=1:length(OUTSTRUCT.bnum)
    disp(i);
    
    indsbase = find(OUTSTRUCT.indsbase{i});
    indswn = find(OUTSTRUCT.indsWN{i});
    
    % ======================
    if (1)
        indsbase_epoch = indsbase(end-round(length(indsbase)/2):end);
        indswn_epoch = indswn(end-round(length(indswn)/2):end);
    end
    
    % ====================== GET COHERENCE DATA
    bnum = OUTSTRUCT.bnum(i);
    enum = OUTSTRUCT.enum(i);
    sw = OUTSTRUCT.switch(i);
    mm = OUTSTRUCT.motifnum(i);
    chpair = OUTSTRUCT.chanpair(i,:);
%     
%     DAT = SwitchCohStruct.bird(bnum).exptnum(enum).switchlist(sw).motifnum(mm);
%     fpre = DAT.fileprefix;
%     fpost = DAT.filesuffix;
%     
%     tmp = ismember(DAT.chanpair, chpair, 'rows');
%     chind = DAT.chanpairstokeep(tmp);
%     cohmat = lt_neural_LFP_loadProcessDat([fpre '/Coh' fpost], chind);
%     
%     % =================== GET PRE AND POST
%     cohbase = mean(cohmat(:,:, indsbase_epoch), 3);
%     cohwn = mean(cohmat(:,:, indswn_epoch), 3);
%     
%     
%     % =========== SAVE IN OUTSTRUCT
%     Cohmean_base = [Cohmean_base; cohbase];
%     Cohmean_wn = [Cohmean_wn; cohwn];
end

% % === save old versions
% OUTSTRUCT.CohMean_Base_orig = OUTSTRUCT.CohMean_Base;
% OUTSTRUCT.CohMean_WN_orig = OUTSTRUCT.CohMean_WN;
% 
% % ===== save new versions
% OUTSTRUCT.CohMean_Base = Cohmean_base;
% OUTSTRUCT.CohMean_WN = Cohmean_wn;
% 
% PARAMS.tbins = PARAMS.tbins_BeforeAlignWN;