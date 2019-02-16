function lt_neural_Coher_SumPlotMotifs2(OUTSTRUCT, SwitchStruct, ...
    averageOverChanPairs, statmethod, indtoget_b_e_s, meanOverExpt, ...
    birdstoplot, expttoplot, swtoplot)
%% lt 12/18/18 - [MOTIFPLOT] FOR EACH MOTIF, PLOT COH CHANGE DEPEDNING ON STATUS
% ########## FOR A GIVEN SYL, PLOT ITS CHANGE IN COHERNECE DEPENDING ON
% WHETHER IT IS TARG, SAME OR DIFF ACROSS EXPERIMENTS


%% ==== filter types of trial.
OUTSTRUCT = lt_structure_RmvEmptyField(OUTSTRUCT); % so followniog code works...

% === filter each dimension independenyl
if ~isempty(birdstoplot)
    indstokeep = ismember(OUTSTRUCT.bnum, birdstoplot);
    OUTSTRUCT = lt_structure_subsample_all_fields(OUTSTRUCT, indstokeep, 1);
end

if ~isempty(expttoplot)
    indstokeep = ismember(OUTSTRUCT.enum, expttoplot);
    OUTSTRUCT = lt_structure_subsample_all_fields(OUTSTRUCT, indstokeep, 1);
end

if ~isempty(swtoplot)
    indstokeep = ismember(OUTSTRUCT.switch, swtoplot);
    OUTSTRUCT = lt_structure_subsample_all_fields(OUTSTRUCT, indstokeep, 1);
end


if ~isempty(indtoget_b_e_s)
    indstokeep = ismember([OUTSTRUCT.bnum OUTSTRUCT.enum OUTSTRUCT.switch], indtoget_b_e_s, 'rows');
    OUTSTRUCT = lt_structure_subsample_all_fields(OUTSTRUCT, indstokeep, 1);
end

%%

% -- switches
[indsgrp_switch, indsgrp_switch_unique] = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch});
% -- chan pairs
% [indsgrp_chanpair, indsgrp_chanpair_unique] = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch, ...
%     OUTSTRUCT.chanpair});


%%
All_bnum = [];
All_istarg = [];
All_issame = [];
All_cohscal = [];
All_motifID = [];
All_motifname = {};

% --- COLLECT DATA FOR EACH SWITCH.
for i=1:length(indsgrp_switch_unique)
    %     swgrpthis = indsgrp_switch_unique(i);
    indsthis = indsgrp_switch==indsgrp_switch_unique(i);
    
    bnum = unique(OUTSTRUCT.bnum(indsthis));
    enum = unique(OUTSTRUCT.enum(indsthis));
    swnum = unique(OUTSTRUCT.switch(indsthis));
    bname = SwitchStruct.bird(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    
    % ====== overall
    
    istarg = OUTSTRUCT.istarg(indsthis);
    issame = OUTSTRUCT.issame(indsthis);
    %     motifs = OUTSTRUCT.motifname(indsthis);
    [motifID, motiflist] = lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)); % ---- get positions within global motif
    %     cohscal = OUTSTRUCT.CohMean_WNminusBase_scalar(indsthis);
    cohscal = OUTSTRUCT.cohscal_diff(indsthis);
    %     chanpair = OUTSTRUCT.chanpair(indsthis,:);
    
    % ---- COLLAPSE ACROSS CHANNEL PAIRS
    if averageOverChanPairs ==1
        cohscal = grpstats(cohscal, motifID, {'mean'});
        istarg = grpstats(istarg, motifID, {'mean'});
        issame = grpstats(issame, motifID, {'mean'});
        motifID = unique(motifID);
    else
        disp('NOT CODED if not averaging over chan pairs...')
        pause;
        return
    end
    
    % === get mean across motifs for this switch
    if strcmp(statmethod, 'zscore')
        ymean = mean(cohscal);
        ystd = std(cohscal);
        cohscal = (cohscal-ymean)/ystd;
    elseif strcmp(statmethod, 'diff')
        ymean = mean(cohscal);
        cohscal = cohscal-ymean;
    elseif strcmp(statmethod, 'minusbase')
        % do nothing
    elseif strcmp(statmethod, 'minDiffType')
        ymean = mean(cohscal(istarg==0 & issame==0));
        cohscal = cohscal-ymean;
    end
    
    % ==================== SAVE ALL OUTPUTS
    All_bnum = [All_bnum; bnum*ones(size(istarg))];
    All_motifname = [All_motifname; motiflist(motifID)'];
    All_istarg = [All_istarg ; istarg];
    All_issame = [All_issame; issame];
    All_cohscal = [All_cohscal; cohscal];
    All_motifID = [All_motifID; motifID];
    
end

%% ======= convert to x values based on type of syl
% i.e. x = 1 2 3 corresponds to targ, same, diff

All_xval = nan(size(All_issame));

All_xval(All_istarg==1) =1;
All_xval(All_istarg==0 & All_issame==1) =2;
All_xval(All_istarg==0 & All_issame==0) =3;


%% =========== [PLOT]
figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


% ================= 1) One plot for each bird.


% =================== ALL IN ONE PLOT, ONE LINE FOR EACH UNIQUE MOTIF
% ============== 1) ONLY CASES WITH TARG, SAME, DIFF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
[indsgrp, indsunique] = lt_tools_grp2idx({All_bnum, All_motifID});
typetoget = [1 2 3]; % i.e. [1 2 3] means must have targ, same and diff.
xlabel('targ - same - diff');
lt_neural_Coher_SumPlotMotifs2_sub;


% ============== 1) ONLY CASES WITH TARG, SAME
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
[indsgrp, indsunique] = lt_tools_grp2idx({All_bnum, All_motifID});
typetoget = [1 2]; % i.e. [1 2 3] means must have targ, same and diff.
xlabel('targ - same');
lt_neural_Coher_SumPlotMotifs2_sub;

% ============== 1) ONLY CASES WITH TARG, DIFF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
[indsgrp, indsunique] = lt_tools_grp2idx({All_bnum, All_motifID});
typetoget = [1 3]; % i.e. [1 2 3] means must have targ, same and diff.
xlabel('targ - diff');
lt_neural_Coher_SumPlotMotifs2_sub;

% ============== 1) ONLY CASES WITH TARG, DIFF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
[indsgrp, indsunique] = lt_tools_grp2idx({All_bnum, All_motifID});
typetoget = [2 3]; % i.e. [1 2 3] means must have targ, same and diff.
xlabel('same - diff');
lt_neural_Coher_SumPlotMotifs2_sub;

% ============== 1) ONLY CASES WITH TARG,
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
[indsgrp, indsunique] = lt_tools_grp2idx({All_bnum, All_motifID});
typetoget = [1]; % i.e. [1 2 3] means must have targ, same and diff.
xlabel('targ');
lt_neural_Coher_SumPlotMotifs2_sub;

% ============== 1) ONLY CASES WITH TARG,
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
[indsgrp, indsunique] = lt_tools_grp2idx({All_bnum, All_motifID});
typetoget = [2]; % i.e. [1 2 3] means must have targ, same and diff.
xlabel('same');
lt_neural_Coher_SumPlotMotifs2_sub;

% ============== 1) ONLY CASES WITH TARG,
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
[indsgrp, indsunique] = lt_tools_grp2idx({All_bnum, All_motifID});
typetoget = [3]; % i.e. [1 2 3] means must have targ, same and diff.
xlabel('diff');
lt_neural_Coher_SumPlotMotifs2_sub;


%% ======= convert to x values based on type of syl
% i.e. x = 1 2 3 corresponds to targ, same, diff

All_xval = nan(size(All_issame));

All_xval(All_istarg==1) =1;
All_xval(All_istarg==0) =2;
% All_xval(All_istarg==0 & All_issame==0) =3;

% All_xval(All_istarg==1) =1;
% All_xval(All_istarg==0 & All_issame==1) =2;
% All_xval(All_istarg==0 & All_issame==0) =3;

%% =========== [PLOT]
figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


% ================= 1) One plot for each bird.
[indsgrp, indsunique] = lt_tools_grp2idx({All_bnum, All_motifID});


% =================== ALL IN ONE PLOT, ONE LINE FOR EACH UNIQUE MOTIF
% ============== 1) ONLY CASES WITH TARG, SAME, DIFF
% ============== 1) ONLY CASES WITH TARG, SAME, DIFF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
typetoget = [1 2]; % i.e. [1 2 3] means must have targ, same and diff.
xlabel('targ - [nontarg]');
lt_neural_Coher_SumPlotMotifs2_sub;




