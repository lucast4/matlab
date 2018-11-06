function dat = lt_neural_LFP_loadProcessDat(filename, chanpairtoget)
%% lt 11/1/18 - extracts coh, phi, s1, s2, or s12.
% filename =
% 
%     '/bluejay5/lucas/analyses/neural/LFP/PROCESSED/14Oct2018_2147/Coh_bird1_expt1_set4_mot1.mat'
% 
% 
% pairstoget =
% 
%      1

% NOTE: names use Coh but is general (to phi, s1, etc).

%% load
cohtmp = load(filename);
tmp1 = fieldnames(cohtmp); assert(length(tmp1)==1, 'why more than one field?');
Coh_ChpairByTrial = cohtmp.(tmp1{1});

%% collect all coherence across all pairs
numchanpairs = length(chanpairtoget);
tmp = size(Coh_ChpairByTrial{1});
dat = nan(tmp(1), tmp(2), length(Coh_ChpairByTrial), numchanpairs); % t, ff, trials, chanpairs
for j=1:numchanpairs
    intmpthis = chanpairtoget(j);
    cohmatthis = lt_neural_Coher_Cell2Mat(Coh_ChpairByTrial(intmpthis,:));
    dat(:,:,:,j) = cohmatthis;
end

