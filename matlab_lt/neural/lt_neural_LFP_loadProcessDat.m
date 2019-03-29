function dat = lt_neural_LFP_loadProcessDat(filename, chanpairtoget, celloutput)
%% lt 11/1/18 - extracts coh, phi, s1, s2, or s12.
% filename =
%
%     '/bluejay5/lucas/analyses/neural/LFP/PROCESSED/14Oct2018_2147/Coh_bird1_expt1_set4_mot1.mat'
%
%
% pairstoget =
%
%      1
% celloutput = 1, then outputs cohcell. if 0, then cohmat

% NOTE: names use Coh but is general (to phi, s1, etc).

% out:
% cohmat

if ~exist('celloutput', 'var')
    celloutput=0;
end
%% load
cohtmp = load(filename);
tmp1 = fieldnames(cohtmp); assert(length(tmp1)==1, 'why more than one field?');
Coh_ChpairByTrial = cohtmp.(tmp1{1});

%% collect all coherence across all pairs

if celloutput==1
    % then just output the cell data.
    dat = Coh_ChpairByTrial;
else
    % convert to mat.
    if isempty(chanpairtoget)
        % then get all chan pairs
        chanpairtoget = 1:size(Coh_ChpairByTrial,1);
    end
    numchanpairs = length(chanpairtoget);
    
    
    tmp = size(Coh_ChpairByTrial{1});
    dat = nan(tmp(1), tmp(2), size(Coh_ChpairByTrial,2), numchanpairs); % t, ff, trials, chanpairs
    for j=1:numchanpairs
        intmpthis = chanpairtoget(j);
        cohmatthis = lt_neural_Coher_Cell2Mat(Coh_ChpairByTrial(intmpthis,:));
        dat(:,:,:,j) = cohmatthis;
    end
end
