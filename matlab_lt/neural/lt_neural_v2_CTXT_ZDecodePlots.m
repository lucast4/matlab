function lt_neural_v2_CTXT_ZDecodePlots(ALLBRANCH)
%% lt 4/16/18 - does decode correlate with FR and is this corrected using z-score

plotstat = 'F1';

%% make sure have multiple negative control trials.
apos = 1;

% if ALLBRANCH.alignpos(apos).ParamsFirstIter.ClassSlide.NumNegControls<2
%     disp('CANNOT RUN - not enough negative controls to get shuffle distribution');
%     return
% else
%     disp(['Num shuffles: ' num2str(ALLBRANCH.alignpos(apos).ParamsFirstIter.ClassSlide.NumNegControls)]);
% end

motifpredur = ALLBRANCH.alignpos(apos).ParamsFirstIter.motifpredur;
%%

numbirds = length(ALLBRANCH.alignpos(apos).bird);
for i = 1:numbirds
    numbranch = length(ALLBRANCH.alignpos(apos).bird(i).branch);
    for b=1:numbranch
        numneuron = length(ALLBRANCH.alignpos(apos).bird(i).branch(b).neuron);
        for n=1:numneuron
            
            dat = ALLBRANCH.alignpos(apos).bird(i).branch(b).neuron(n);
            
            [dat.FR.classnum(:).FRsmooth_rate_CommonTrialDur]
            
        end
    end
end