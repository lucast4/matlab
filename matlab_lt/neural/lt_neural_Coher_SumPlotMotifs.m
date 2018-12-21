function lt_neural_Coher_SumPlotMotifs(OUTSTRUCT, SwitchStruct,MOTIFSTATS_Compiled,PARAMS)

%% lt 12/17/18 - plots cioherence scalar change during laerning, relative to posotion of motif
% i.e. for all motifs plot change in coh scalar.
clim = [-0.2 0.2];

% %% for each case get wn minus base, mean coh scalar
% cohscal_diff = [];
% for i=1:length(OUTSTRUCT.bnum)
%     
%     indsbase = OUTSTRUCT.indsbase_epoch{i};
%     indswn = OUTSTRUCT.indsWN_epoch{i};
%     
%     cohscal = OUTSTRUCT.cohscal{i};
%     
%     cohdiff = mean(cohscal(indswn)) - mean(cohscal(indsbase));
%     
%     cohscal_diff = [cohscal_diff; cohdiff];
% end
% 
% OUTSTRUCT.cohscal_diff = cohscal_diff;

%% ############

[indsgrp_switch, indsgrp_switch_unique] = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch});
[indsgrp_chanpair, indsgrp_chanpair_unique] = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch, ...
    OUTSTRUCT.chanpair});

figcount=1;
subplotrows=6;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% indsgrp_switch_unique = unique(indsgrp_switch);
% indsgrp_chanpair_unique = unique(indsgrp_chanpair);
for i=1:length(indsgrp_switch_unique)
    swgrpthis = indsgrp_switch_unique(i);
    bnum = unique(OUTSTRUCT.bnum(indsgrp_switch==swgrpthis));
    enum = unique(OUTSTRUCT.enum(indsgrp_switch==swgrpthis));
    swnum = unique(OUTSTRUCT.switch(indsgrp_switch==swgrpthis));
    bname = MOTIFSTATS_Compiled.birds(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(swnum)]);
    ylabel(['WN - base, coh']);
    
    % ====== plot each channel pair its own line
    for chanpair = indsgrp_chanpair_unique'
        indsthis = indsgrp_switch==swgrpthis & indsgrp_chanpair==chanpair;
        if ~any(indsthis)
            continue
        end
        motifID = lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)); % ---- get positions within global motif
        
        cohscal = OUTSTRUCT.cohscal_diff(indsthis);
        %         cohscal = OUTSTRUCT.CohMean_WNminusBase_scalar(indsthis);
        
        assert(length(unique(motifID)) == length(motifID), 'each motif shoudl have at most 1 datapoint');
        % ---- sort in order of motifs
        [~, indsort] = sort(motifID);
        motifID = motifID(indsort);
        cohscal = cohscal(indsort);
        plot(motifID, cohscal, 'o-k');
    end
    lt_plot_zeroline;
    
    % ====== overall
    indsthis = indsgrp_switch==swgrpthis;
    
    istarg = OUTSTRUCT.istarg(indsthis);
    motifs = OUTSTRUCT.motifname(indsthis);
    motifID = lt_neural_QUICK_MotifID(bname, motifs); % ---- get positions within global motif
    %     cohscal = OUTSTRUCT.CohMean_WNminusBase_scalar(indsthis);
    cohscal = OUTSTRUCT.cohscal_diff(indsthis);
    
    [ymean, ysem] = grpstats(cohscal, motifID, {'mean', 'sem'});
    x = unique(motifID);
    lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
    
    % -------- NOTE DOWN POSITION OF TARGET SYSL
    indsthis = indsgrp_switch==swgrpthis & OUTSTRUCT.istarg==1;
    xtarg = unique(lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)));
    plot(xtarg, clim(1)+0.02, '^r');
    
    % ------- NOTE POSITION OF SAME_TYPES
    indsthis = indsgrp_switch==swgrpthis & OUTSTRUCT.istarg==0 & OUTSTRUCT.issame==1;
    xtarg = unique(lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)));
    if ~isempty(xtarg)
        plot(xtarg, clim(1)+0.02, '^b');
    end
    
    % ----- labels
    [~, motiflabels] = lt_neural_QUICK_MotifID(bname);
    set(gca, 'XTick', 1:length(motiflabels));
    set(gca, 'XTickLabel', motiflabels);
    rotateXLabels(gca, 90);
    ylim(clim);
    
end



