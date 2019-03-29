function LearningMetaDat = lt_neural_v2_LoadLearnMetadat(ver)
%%

% ver, leave empty. otherwise pick out specific versions.

if ~exist('ver', 'var')
    ver = '';
end

%%
if isempty(ver)
    cd('/bluejay5/lucas/analyses/neural/');
    load('LearningMetaDat.mat');
    eval('!cp LearningMetaDat.mat BACKUP/');  % -- save a backup
    
else
    % for RALMAN learning, correlations. put fake swtiches so that code extracts data at those switches
    % i.e. recordings did not have real switches sometimes
    % STRATEGY: added a new column with a new syllable and entered switch
    % times
    cd('/bluejay5/lucas/analyses/neural/');
    load(['LearningMetaDat_' ver '.mat']);
    % DO NOT SAVE BACKUP - i.e. do not overwrite old backup
end

