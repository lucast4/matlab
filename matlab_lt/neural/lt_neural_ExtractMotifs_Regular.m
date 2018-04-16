%% gets motifs, not learning, gets DIR.
% used for 
% - DIR vs UNDIR
% - Pop analysis (not learning)



%%
close all; clear MOTIFSTATS_Compiled;
collectWNhit=0; % NOTE!!! temporary - need to change so don't need to extract audio each time (i.e. do and save)
onlyCollectTargSyl=0;
LearnKeepOnlyBase = 1;
saveOn = 1;
OrganizeByExpt =0;
collectFF=1;

% --- to make sure extracts motifs
% MotifsToCollect = {'pu69wh78', {'(j)jjbhhg', '(a)abhhg'}};
%     Params_regexp.motif_predur = 0.05;
%     Params_regexp.motif_postdur = 0.05;
%     Params_regexp.preAndPostDurRelSameTimept = 0;
%     Params_regexp.RemoveIfTooLongGapDur = 1;

% --- to make sure extracts DIR
MotifsToCollect = [];
Params_regexp.motif_predur = [];
    Params_regexp.motif_postdur = [];
    Params_regexp.preAndPostDurRelSameTimept = 1; % CHANGED to 1 (3/31/18) - change back to 0 if causes problems.
    Params_regexp.RemoveIfTooLongGapDur = [];
    Params_regexp.extractDirSong = 1;


MOTIFSTATS_Compiled = lt_neural_v2_ANALY_MultExtractMotif(SummaryStruct, ...
    collectWNhit, LearnKeepOnlyBase, saveOn, onlyCollectTargSyl, OrganizeByExpt,...
    collectFF, MotifsToCollect, Params_regexp);


%% ============ IMPORTANT, FILTER OUT DIR IF NOT WANTED
