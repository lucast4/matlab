%% 7/3/15 - Notes on all autolabel methods up to today

%% Evtaf amp sim version
% Using modified version of Joanne's, which uses evtaf_amp simulation
% This autolabels, and does post-analysis, looking at all syls and
% replacing any mistakes.

% ================== 1) Input template parameters. This runs evtafsim and autolabels.
% Can choose whether to overwrite old .cbin.not.mat or to add onto old
% ones.

Params.batch='batch.rec_FB.rand';

Params.ampThresh=21000;
Params.min_dur=13;
Params.min_int=1;

Params.syl.pre='';
Params.syl.post='';
Params.syl.targ='b';

Params.overwrite_notmat=1;


% TEMPLATE SETTINGS
Params.TEMPLATE.templatefile = 'autolabel_templ_b1_v2.dat'; % should be one dir up.
% Template params: one array entry for each column. 
Params.TEMPLATE.cntrng(1).MIN=1;
Params.TEMPLATE.cntrng(1).MAX=3;
Params.TEMPLATE.cntrng(1).NOT=0;
Params.TEMPLATE.cntrng(1).MODE=1;
Params.TEMPLATE.cntrng(1).TH=1;
Params.TEMPLATE.cntrng(1).AND=0;
Params.TEMPLATE.cntrng(1).BTMIN=0;

Params.TEMPLATE.cntrng(2).MIN=1;
Params.TEMPLATE.cntrng(2).MAX=3;
Params.TEMPLATE.cntrng(2).NOT=0;
Params.TEMPLATE.cntrng(2).MODE=1;
Params.TEMPLATE.cntrng(2).TH=2.2;
Params.TEMPLATE.cntrng(2).AND=0;
Params.TEMPLATE.cntrng(2).BTMIN=0;

Params.TEMPLATE.refract=0.2;

% RUN
[fnames, sylnum, vlsorfn, vlsorind]=lt_autolabel_function(Params);

% ============ 2) Use evsonganaly manually on the .wav file created above
% (contains only the syls you chose)
% INSTRUCTIONS: 
% 1) open .wav file using evsonganaly
% 2) change threshold to segment all syls indivudally
% 3) any syl labeled "-" (default) will remain unchanged (i.e. will stay autolabeled). 
%     give a new label to any mislabeled syl - that will be its new actual label
evsonganaly


% ============ 3) Replace hand-checekd mislabeld syls 
lt_autolabel_FixHandCheckedSyls(fnames, sylnum, vlsorfn, vlsorind)



%% === Evtafv4 sim version
% After this, can run stuff same as for evtaf amp version to replace false
% positives.

batch = 'batch.rand.keep';
config= '/bluejay4/lucas/birds/pk32/config.evconfig2';

syl.targ='g';
syl.pre='jabcddef';
syl.post=''; 
NoteNum=0; 

ampThresh=15000;
min_dur=30;
min_int=4;

overwrite_notmat=1; % will always make backup folder

[fnames, sylnum, vlsorfn, vlsorind]=lt_autolabel_EvTAFv4(batch, config, syl, NoteNum, ampThresh, min_dur, min_int, overwrite_notmat);

% ============ 2) Use evsonganaly manually on the .wav file created above
% (contains only the syls you chose)
% INSTRUCTIONS: 
% 1) open .wav file using evsonganaly
% 2) change threshold to segment all syls indivudally
% 3) any syl labeled "-" (default) will remain unchanged (i.e. will stay autolabeled). 
%     give a new label to any mislabeled syl - that will be its new actual label
evsonganaly


% ============ 3) Replace hand-checekd mislabeld syls 
lt_autolabel_FixHandCheckedSyls(fnames, sylnum, vlsorfn, vlsorind)



%% using feature vectors

lt_autolabel

%% ======= EVTAFV4 sim version, but good workflow if doing multiple configs
% AND ITERATION OVER DAY
clear all; close all;

% ========================= 0) ECTRACT DIRECTORIRES
basedir = '/bluejay5/lucas/birds/pu53wh88';
date_range_base={'10Feb2015','12Feb2015'};
date_range_WN={'14Feb2015','16Feb2015'};
experiment = 'SeqDepPitchShift3';

% -------- COLLECT METADAT
cd(basedir);
MetadataStruct=lt_metadata_collect;

condition='';
notes='';
only_labeled_dirs=0;

% ----- BASELINE
ListOfDirs1=lt_metadata_find_dirs(MetadataStruct, experiment, condition, ...
    notes, date_range_base, only_labeled_dirs, 2);

% ------ WN
ListOfDirs2=lt_metadata_find_dirs(MetadataStruct, experiment, condition, ...
    notes, date_range_WN, only_labeled_dirs, 2);

% ------- COMBINE
ListOfDirs = [ListOfDirs1 ListOfDirs2];


% ============================== RUN, ITERATE OVER DAYS
for j=1:length(ListOfDirs)

    % ==================== 0) go to day folder
    cd(ListOfDirs(j).dirname);
    
    
    % ==================== 1) extract all s
    lt_make_batch(3);
    lt_cleandirAuto('batch', 1000, 5, 5);
    % randsamp('batch.keep', 0.1);
    % batch = 'batch.keep.rand';
    randsamp('batch.keep', 0.05);
    batch = 'batch.keep.rand';
    
    % ==================== 2) move old .notmat
    if ~exist('OLDNOTMAT_SeqDepPitch', 'dir')
        mkdir OLDNOTMAT_SeqDepPitch
        !cp *.not.mat* OLDNOTMAT_SeqDepPitch
    else
        disp('not making OLDNOTMAT_SeqDepPitch, since already made!!');
    end
    
    
    % ==================== 3) run autolabel
    % ---- GENERAL PARAMS
    ampThresh=30000;
    min_dur=30;
    min_int=3;
    
    % ---- MOTIF 1 [abB] [not labeling abb(b) % since WN target is ab(b)]
    config= '/bluejay5/lucas/birds/pu53wh88/config_AL_aB.evconfig2';
    syl.targ='b';
    syl.pre='ab';
    syl.post='';
    NoteNum=0;
    overwrite_notmat=1; % will always make backup folder
    
    [fnames, sylnum, vlsorfn, vlsorind]=lt_autolabel_EvTAFv4(batch, config, syl, NoteNum, ampThresh, min_dur, min_int, overwrite_notmat);
    
    
    % ------ MOTIF 2 [accBb]
    config= '/bluejay5/lucas/birds/pu53wh88/config_AL_cB.evconfig2';
    syl.targ='b';
    syl.pre='acc';
    syl.post='b';
    NoteNum=0;
    overwrite_notmat=0; % will always make backup folder
    
    [fnames, sylnum, vlsorfn, vlsorind]=lt_autolabel_EvTAFv4(batch, config, syl, NoteNum, ampThresh, min_dur, min_int, overwrite_notmat);
    
    
    % ================ 4) MAKE WAV FILE TO LOOK FOR FALSE POSITIVES
    sylstocheck = {syl.targ};  % cell array of all syls you want to check.
    
    for j=1:length(sylstocheck)
        [fnames, sylnum]=lt_jc_chcklbl(batch, sylstocheck{j}, 0.025,0.025,'','','');
        [vlsorfn vlsorind]=jc_vlsorfn(batch, sylstocheck{j},'','');
    end
    
end


% ============ 2) Use evsonganaly manually on the .wav file created above
% (contains only the syls you chose)
% INSTRUCTIONS: 
% 1) open .wav file using evsonganaly
% 2) change threshold to segment all syls indivudally
% 3) any syl labeled "-" (default) will remain unchanged (i.e. will stay autolabeled). 
%     give a new label to any mislabeled syl - that will be its new actual label
evsonganaly


% ============ 3) Replace hand-checekd mislabeld syls 
lt_autolabel_FixHandCheckedSyls(fnames, sylnum, vlsorfn, vlsorind)


