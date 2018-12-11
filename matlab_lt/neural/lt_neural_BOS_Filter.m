function [AllUnits, birdname] = lt_neural_BOS_Filter(DATSTRUCT, SummaryBOS,PARAMS, exptlist_toget, bostype_toget, ...
    bregionthis_toget, motif_toget)
%% lt 12/2/18 - extracts trial data, grouped by unit, for specific desired dataset (filter inputs)

% 
% 
% exptlist_toget = 1:3; % takes units from these experiments (1, ,2 ...,
% across all birds)

% bostype_toget = 'fwd';
% bregionthis_toget = 'LMAN';
% motif_toget = '(a)abhh'; 
% motifs.

% outputws:
% AllUnits.FRmat = {};
% AllUnits.FRmat_t = {};
% AllUnits.SylOnsets = {};
% AllUnits.SylOffsets = {};
% AllUnits.bostype = {};
% AllUnits.motifnum = [];
% AllUnits.exptnum = [];
% AllUnits.unitnum = [];
% AllUnits.bregions = {};


%% --- iterate over experiments and units.
[indsgrp, indsgrp_unique] = lt_tools_grp2idx({DATSTRUCT.exptnum DATSTRUCT.unitnum});

motifslist = unique(DATSTRUCT.motifnum);
bostypelist = unique(DATSTRUCT.bostype);

AllUnits.FRmat = {};
AllUnits.FRmat_t = {};
AllUnits.SylOnsets = {};
AllUnits.SylOffsets = {};
AllUnits.bostype = {};
AllUnits.motifnum = [];
AllUnits.exptnum = [];
AllUnits.unitnum = [];
AllUnits.bregions = {};

for i=1:length(indsgrp_unique)
    
    for mm=motifslist'
        for bb=1:length(bostypelist)
            bosthis = bostypelist(bb);
            
            indsthis = indsgrp==indsgrp_unique(i) & DATSTRUCT.motifnum==mm & DATSTRUCT.bostype==bosthis;
            if isempty(indsthis)
                continue
            end
            
            rendnum = DATSTRUCT.songrendnum(indsthis);
            spks = DATSTRUCT.Spktimes(indsthis);
            sylonsets = cell2mat(DATSTRUCT.Onstimes(indsthis));
            syloffsets = cell2mat(DATSTRUCT.Offtimes(indsthis));
            
            minmaxtime = DATSTRUCT.minmax_time(indsthis,:);
            
            % === for pplotting.
            birdname = unique(DATSTRUCT.birdname(indsthis));
            assert(length(birdname)==1, 'has to be one bird, since specified motif...'); 
            birdname = birdname{1};
            exptnum = unique(DATSTRUCT.exptnum(indsthis));
            unitnum = unique(DATSTRUCT.unitnum(indsthis));
            bregion = SummaryBOS.expt(exptnum).bregions{unitnum};
            motif = unique(DATSTRUCT.motifregexp(indsthis));
            bostype = PARAMS.BOSnames{strcmp(PARAMS.BOSbirdname, birdname)}{bosthis};
            
            %% ======================== FILTER...
            if ~any(ismember(exptlist_toget, exptnum))
                continue
            end
            if ~strcmp(bostype_toget, bostype)
                continue
            end
            if ~strcmp(bregion, bregionthis_toget)
                continue
            end
            if ~strcmp(motif, motif_toget)
                continue
            end
            
            %% ========== PLOT FOR THIS UNIT, RESPONSES OVER ALL RENDITIONS OF
            % THIS MOTIF
            assert(length(unique(minmaxtime(:,1)))==1, 'assumes all predur from same alignment');
            mintime = minmaxtime(1,1);
            maxtime = min(minmaxtime(:,2));
            
            
            [FRmat, t] = lt_neural_QUICK_Spk2FRmat(spks, mintime, maxtime);
            
            frmean = mean(FRmat,2);
            frsem = lt_sem(FRmat');
            
            
            % --- SYLS ONSET OFFSET
            on = median(sylonsets,1);
            off = median(syloffsets,1);
            
            
            % =========================== COLLECT
            AllUnits.FRmat = [AllUnits.FRmat; FRmat];
            AllUnits.FRmat_t = [AllUnits.FRmat_t; t];
            AllUnits.SylOnsets = [AllUnits.SylOnsets; sylonsets];
            AllUnits.SylOffsets = [AllUnits.SylOffsets; syloffsets];
            
            AllUnits.motifnum = [AllUnits.motifnum; mm];
            AllUnits.bostype = [AllUnits.bostype; bostype];
            AllUnits.exptnum = [AllUnits.exptnum; exptnum];
            AllUnits.unitnum = [AllUnits.unitnum; unitnum];
            
            AllUnits.bregions = [AllUnits.bregions; bregion];

        end
    end
end