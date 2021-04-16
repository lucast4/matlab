function [spkmeanNonadAdapt_LMAN, spkmeanNonadAdapt_RA] = lt_neural_POPLEARN_Xcov_Epochs_DistrSub(...
    OUTSTRUCT_XCOV, epochtoplot)
%% lt 3/25/19 - extract nspks
dohack =1; % if 1, then does hack where I tried to make it work if epoch was more than one number. 
% but ran into issues, so now I just run each epoch separately and average
% the outputs. [def = 1, since otherwise doesn't work if anything has nan]

%%

assert(all(strcmp(OUTSTRUCT_XCOV.bregionpair, 'LMAN-RA')), 'assumes n1 is lMAN, n2 is RA');

epochthis = epochtoplot +1; % since 0 is baseline.

if dohack==0
    % == old version, doesn't work if epochthis is array length>1
    nspksAll= cellfun(@(x)x(epochthis).nspksByNeuron, OUTSTRUCT_XCOV.epochSplitStatsAll, 'UniformOutput', 0);
else
    nspksAll= cell(length(OUTSTRUCT_XCOV.bregionpair),1);
    for j=1:length(OUTSTRUCT_XCOV.epochSplitStatsAll)
        
        nspk1 =[];
        nspk2 =[];
        
        for jj=1:length(epochthis)
            if length(OUTSTRUCT_XCOV.epochSplitStatsAll{j}(epochthis(jj)).nspksByNeuron)==1
                % i.e. is missing data ...
                continue
            end
            nspk1 = [nspk1; OUTSTRUCT_XCOV.epochSplitStatsAll{j}(epochthis(jj)).nspksByNeuron{1}];
            nspk2 = [nspk2; OUTSTRUCT_XCOV.epochSplitStatsAll{j}(epochthis(jj)).nspksByNeuron{2}];
            %         nspk = [nspk; OUTSTRUCT_XCOV.epochSplitStatsAll{j}(epochthis(jj)).nspksByNeuron{1}]
        end
        nspksAll{j} = {nspk1, nspk2};
    end
end

indshi = cellfun(@(x)[x(epochthis).inds_hi], OUTSTRUCT_XCOV.epochSplitStatsAll, 'UniformOutput', 0);
indslo = cellfun(@(x)[x(epochthis).inds_lo], OUTSTRUCT_XCOV.epochSplitStatsAll, 'UniformOutput', 0);
learndirTarg = OUTSTRUCT_XCOV.learndirTarg;


% ============== any nan in inds? REMOVE
if dohack==1
for i=1:length(indshi)
    if sum(isnan(indshi{i}))>2
        disp('why so many nan?');
        keyboard
    elseif sum(isnan(indshi{i}))==1
        % then is because was nan for an epoch was compressed into one nan
        assert(find(isnan(indslo{i}))==find(isnan(indshi{i})))
        
        indshi{i}(isnan(indshi{i})) = [];
        indslo{i}(isnan(indslo{i})) = [];
    end
end
end
% indslo = cellfun(@(x)x(epochthis).inds_lo, OUTSTRUCT_XCOV.epochSplitStatsAll, 'UniformOutput', 0);
% ffall = cellfun(@(x)x(epochthis).ffthis, OUTSTRUCT_XCOV.epochSplitStatsAll, 'UniformOutput', 0);

% === get mean and sem spikes for each case

% -------------------- LMAN
spks = cellfun(@(x)x{1}, nspksAll, 'UniformOutput', 0);
% - get mean for lo and hi pitch trials
spkmean_hi = cellfun(@(x,y)nanmean(x(y)), spks, indshi, 'UniformOutput', 0);
spkmean_lo = cellfun(@(x,y)nanmean(x(y)), spks, indslo, 'UniformOutput', 0);
% - convert to spk on adaptive vs. nonadaptivwe
spkmeanNonadAdapt = [spkmean_lo spkmean_hi];
spkmeanNonadAdapt(learndirTarg==-1, :) = fliplr(spkmeanNonadAdapt(learndirTarg==-1, :));

spkmeanNonadAdapt_LMAN = spkmeanNonadAdapt;


% -------------------- RA
spks = cellfun(@(x)x{2}, nspksAll, 'UniformOutput', 0);
% - get mean for lo and hi pitch trials
spkmean_hi = cellfun(@(x,y)nanmean(x(y)), spks, indshi, 'UniformOutput', 0);
spkmean_lo = cellfun(@(x,y)nanmean(x(y)), spks, indslo, 'UniformOutput', 0);
% - convert to spk on adaptive vs. nonadaptivwe
spkmeanNonadAdapt = [spkmean_lo spkmean_hi];
spkmeanNonadAdapt(learndirTarg==-1, :) = fliplr(spkmeanNonadAdapt(learndirTarg==-1, :));

spkmeanNonadAdapt_RA = spkmeanNonadAdapt;


% ==== convert to mat
tmp = [];
for i=1:size(spkmeanNonadAdapt_LMAN,1)
    if isnan(spkmeanNonadAdapt_LMAN{i,1})
        tmp = [tmp; [nan nan]];
    else
    tmp = [tmp; cell2mat(spkmeanNonadAdapt_LMAN(i,:))];
    end
end
spkmeanNonadAdapt_LMAN = tmp;

tmp = [];
for i=1:size(spkmeanNonadAdapt_RA,1)
    if isnan(spkmeanNonadAdapt_RA{i,1})
        tmp = [tmp; [nan nan]];
    else
    tmp = [tmp; cell2mat(spkmeanNonadAdapt_RA(i,:))];
    end
end
spkmeanNonadAdapt_RA = tmp;
