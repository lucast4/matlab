function [spkmeanNonadAdapt_LMAN, spkmeanNonadAdapt_RA] = lt_neural_POPLEARN_Xcov_Epochs_DistrSub(...
    OUTSTRUCT_XCOV, epochtoplot)
%% lt 3/25/19 - extract nspks


%%

assert(all(strcmp(OUTSTRUCT_XCOV.bregionpair, 'LMAN-RA')), 'assumes n1 is lMAN, n2 is RA');

epochthis = epochtoplot +1;
nspksAll= cellfun(@(x)x(epochthis).nspksByNeuron, OUTSTRUCT_XCOV.epochSplitStatsAll, 'UniformOutput', 0);
indshi = cellfun(@(x)x(epochthis).inds_hi, OUTSTRUCT_XCOV.epochSplitStatsAll, 'UniformOutput', 0);
indslo = cellfun(@(x)x(epochthis).inds_lo, OUTSTRUCT_XCOV.epochSplitStatsAll, 'UniformOutput', 0);
learndirTarg = OUTSTRUCT_XCOV.learndirTarg;

% indslo = cellfun(@(x)x(epochthis).inds_lo, OUTSTRUCT_XCOV.epochSplitStatsAll, 'UniformOutput', 0);
% ffall = cellfun(@(x)x(epochthis).ffthis, OUTSTRUCT_XCOV.epochSplitStatsAll, 'UniformOutput', 0);

% === get mean and sem spikes for each case

% -------------------- LMAN
spks = cellfun(@(x)x{1}, nspksAll, 'UniformOutput', 0);
% - get mean for lo and hi pitch trials
spkmean_hi = cellfun(@(x,y)nanmean(x(y)), spks, indshi);
spkmean_lo = cellfun(@(x,y)nanmean(x(y)), spks, indslo);
% - convert to spk on adaptive vs. nonadaptivwe
spkmeanNonadAdapt = [spkmean_lo spkmean_hi];
spkmeanNonadAdapt(learndirTarg==-1, :) = fliplr(spkmeanNonadAdapt(learndirTarg==-1, :));

spkmeanNonadAdapt_LMAN = spkmeanNonadAdapt;


% -------------------- RA
spks = cellfun(@(x)x{2}, nspksAll, 'UniformOutput', 0);
% - get mean for lo and hi pitch trials
spkmean_hi = cellfun(@(x,y)nanmean(x(y)), spks, indshi);
spkmean_lo = cellfun(@(x,y)nanmean(x(y)), spks, indslo);
% - convert to spk on adaptive vs. nonadaptivwe
spkmeanNonadAdapt = [spkmean_lo spkmean_hi];
spkmeanNonadAdapt(learndirTarg==-1, :) = fliplr(spkmeanNonadAdapt(learndirTarg==-1, :));

spkmeanNonadAdapt_RA = spkmeanNonadAdapt;


