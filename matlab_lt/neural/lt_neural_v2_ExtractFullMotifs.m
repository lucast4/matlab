%% lt - quick script for extracting specific bird and motifs into Motifstruct
% ====== by default does time warping.

close all; clear MOTIFSTATS_Compiled;
collectWNhit=0; % NOTE!!! temporary - need to change so don't need to extract audio each time (i.e. do and save)
onlyCollectTargSyl=0;
LearnKeepOnlyBase = 1;
saveOn = 1;
OrganizeByExpt =0;
collectFF=1;

% --- to make sure extracts motifs
% MotifsToCollect = {'pu69wh78', {'(j)jjbhhg', '(a)abhhg'}, ...
%     'wh44wh39', {'(n)hh', '(d)kccbb'}};
MotifsToCollect = {'pu69wh78', {'(j)bhh', 'j(b)hh', '(a)bhh', 'a(b)hh'}, ...
    'wh44wh39', {'(n)hh', '(d)k', '(c)cbb'}}; % these versions limit the distortion caused
% by time warping --> I think distortion leads to wide cross correlations.
Params_regexp.motif_predur = 0.05;
Params_regexp.motif_postdur = 0.15;
Params_regexp.preAndPostDurRelSameTimept = 0;
Params_regexp.RemoveIfTooLongGapDur = 1;
Params_regexp.extractDirSong = 1;

MOTIFSTATS_Compiled = lt_neural_v2_ANALY_MultExtractMotif(SummaryStruct, ...
    collectWNhit, LearnKeepOnlyBase, saveOn, onlyCollectTargSyl, OrganizeByExpt,...
    collectFF, MotifsToCollect, Params_regexp);



%% TIME WARP ALL DATA (FOR A GIVEN MOTIF, SAME WARP ACROSS ALL NEURONS)
% close all;

MOTIFSTATS_Compiled = lt_neural_MOTIF_TimeWarpAll(MOTIFSTATS_Compiled);