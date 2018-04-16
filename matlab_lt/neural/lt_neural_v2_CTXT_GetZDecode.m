function ALLBRANCH = lt_neural_v2_CTXT_GetZDecode(ALLBRANCH)
%% lt 4/12/18 - for all data, converts decode to z-score rel shuffle
% NOTE: must have multiple neg control shuffle datapoints for this to work
plotstat = 'F1';

%% make sure have multiple negative control trials.
apos = 1;

if ALLBRANCH.alignpos(apos).ParamsFirstIter.ClassSlide.NumNegControls<2
    disp('CANNOT RUN - not enough negative controls to get shuffle distribution');
    return
else
    disp(['Num shuffles: ' num2str(ALLBRANCH.alignpos(apos).ParamsFirstIter.ClassSlide.NumNegControls)]);
end

motifpredur = ALLBRANCH.alignpos(apos).ParamsFirstIter.motifpredur;
%%

numbirds = length(ALLBRANCH.alignpos(apos).bird);
for i = 1:numbirds
    numbranch = length(ALLBRANCH.alignpos(apos).bird(i).branch);
    for b=1:numbranch
        numneuron = length(ALLBRANCH.alignpos(apos).bird(i).branch(b).neuron);
        for n=1:numneuron
            
            dat = ALLBRANCH.alignpos(apos).bird(i).branch(b).neuron(n);
            
            if isempty(dat.yvals)
                continue
            end
            
            disp([num2str(i) '-' num2str(b) '-' num2str(n)]);
            
            % ============
            assert(length(dat.yvals) == size(dat.ConfMatAll_NEG_Mult,4), 'not enought ime bins?');
            
            % ============ GET ZSCORE FOR EACH TIME BIN
            ntimebins = length(dat.yvals);
            figcount=1;
            subplotrows=5;
            subplotcols=3;
            fignums_alreadyused=[];
            hfigs=[];
            hsplots = [];
            
            decode_z = nan(1,ntimebins);
            decodeneg_mean = nan(1,ntimebins);
            decodeneg_std = nan(1,ntimebins);
            
            for t = 1:ntimebins
                confmatall = dat.ConfMatAll_NEG_Mult(:,:,:, t);
                confcell = squeeze(mat2cell(confmatall, size(confmatall,1), size(confmatall,2), ones(1,size(confmatall,3))));
                
                sts = lt_neural_ConfMatStats(confcell);
                decode_neg = [sts.F1];
                
                % === collect mean and std
                mean_neg = mean(decode_neg);
                std_neg = std(decode_neg);
                
                decode_z(t) = (dat.yvals(t) - mean_neg)/std_neg;
                decodeneg_mean(t) = mean_neg;
                decodeneg_std(t) = std_neg;
                
                
                % ======= compare actual decode to decode neg
                if (0)
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    title(num2str(t));
                    lt_plot_histogram(decode_neg);
                end
            end
            
            % =================== OUTPUT
            ALLBRANCH.alignpos(apos).bird(i).branch(b).neuron(n).yvals_z = decode_z;
            ALLBRANCH.alignpos(apos).bird(i).branch(b).neuron(n).yvals_neg_shuffmean = decodeneg_mean;
            ALLBRANCH.alignpos(apos).bird(i).branch(b).neuron(n).yvals_neg_shuffstd = decodeneg_std;
            
                        
            if (0) % sanity check
                lt_figure; hold on;
                
                lt_subplot(2,2,1); hold on;
                title('decode + neg control distrib');
                plot(dat.xtimes, dat.yvals, '-ok');
                shadedErrorBar(dat.xtimes, decodeneg_mean, decodeneg_std, {'Color', 'r'}, 1);
                
                lt_subplot(2,2,2); hold on;
                ylabel('fr (mean across class)');
                frdat = [dat.FR.classnum.FRsmooth_rate_CommonTrialDur];
                frx = dat.FR.classnum(1).FRsmooth_xbin_CommonTrialDur;
                shadedErrorBar(frx-motifpredur, mean(frdat,2), std(frdat, 0,2))
                
                lt_subplot(2,2, 4); hold on;
                ylabel('fr (mean across class)');
                numclasses = length(dat.FR.classnum);
                for cc =1:numclasses
                frdat = [dat.FR.classnum(cc).FRsmooth_rate_CommonTrialDur];
                frx = dat.FR.classnum(1).FRsmooth_xbin_CommonTrialDur;
                shadedErrorBar(frx-motifpredur, mean(frdat,2), std(frdat, 0,2), {'Color', [rand rand rand]}, 1);
                end
                
                lt_subplot(2,2,3); hold on;
                ylabel('zscore decode');
                plot(dat.xtimes, decode_z, '-ob');
            end
            
        end
    end
end