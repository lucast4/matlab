function lt_neural_Coher_CohScalMotifPlot(OUTSTRUCT, SwitchStruct, MOTIFSTATS_Compiled, PARAMS, ...
    onlyplotgood, useshuffpval)
%% lt 1/6/19 -  OVERVIEW PLOT - === each switch, each motif, each chan pair, pre and post WN

%% filter to only training onset expts?
if onlyplotgood==1
    OUTSTRUCT = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT, SwitchStruct);
end

%% get inds for switchs and channel pairs.

[indsgrp_switch, indsgrp_switch_unique] = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch});
[indsgrp_chanpair, indsgrp_chanpair_unique] = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch, ...
    OUTSTRUCT.chanpair});


%% RUN

figcount=1;
subplotrows=6;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


Yall = []; % targ, same, diff, nontarg (means across chans, motifs)
All_bname ={};
All_bnum = [];
All_ename = {};
All_swnum = [];
for i=1:length(indsgrp_switch_unique)
    swgrpthis = indsgrp_switch_unique(i);
    bnum = unique(OUTSTRUCT.bnum(indsgrp_switch==swgrpthis));
    enum = unique(OUTSTRUCT.enum(indsgrp_switch==swgrpthis));
    swnum = unique(OUTSTRUCT.switch(indsgrp_switch==swgrpthis));
    bname = MOTIFSTATS_Compiled.birds(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    
    
    %% ====== plot each channel pair its own line
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(swnum)]);
    ylabel(['rho (ff vs coh) [bk=base; rd=WN, bu=all]']);
    for chanpair = indsgrp_chanpair_unique'
        indsthis = indsgrp_switch==swgrpthis & indsgrp_chanpair==chanpair;
        if ~any(indsthis)
            continue
        end
        
        %% 1) baseline
        XSHIFT = -0.2;
        pcol = 'k';
        fieldtoget = 'cohscal_ff_rho_rhoSE_base';
        
        % =============
        motifID = lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)); % ---- get positions within global motif
        chnums = OUTSTRUCT.chanpair(indsthis,:);
        chnums = unique(chnums)'; assert(length(chnums)==2);
        rho = OUTSTRUCT.(fieldtoget)(indsthis, 1);
        rho_sem = OUTSTRUCT.(fieldtoget)(indsthis, 2);
        learndir = unique(OUTSTRUCT.learndirTarg(indsthis));
        
%         cohscal = OUTSTRUCT.cohscal_diff(indsthis);
        
        assert(length(unique(motifID)) == length(motifID), 'each motif shoudl have at most 1 datapoint');
        
        % ---- sort in order of motifs
        [~, indsort] = sort(motifID);
        motifID = motifID(indsort);
%         cohscal = cohscal(indsort);
        rho = rho(indsort);
        rho_sem = rho_sem(indsort);
        
        indstmp = ~isnan(rho);
        x = motifID(indstmp) + XSHIFT;
        y = rho(indstmp);
        ysem = rho_sem(indstmp);
        if all(isnan(ysem))
        plot(x, y, '-o', 'Color', pcol);    
        else
        errorbar(x, y, ysem, '-o', 'Color', pcol);
        end
        lt_plot_text(motifID(end)+0.3, rho(end), num2str(chnums), 'k', 9);
        lt_plot_annotation(1, ['targdir ' num2str(learndir)], 'm');
        %% 2) WN
        XSHIFT = 0;
        pcol = 'r';
        fieldtoget = 'cohscal_ff_rho_rhoSE_WN';
        
        % =============
        motifID = lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)); % ---- get positions within global motif
        chnums = OUTSTRUCT.chanpair(indsthis,:);
        chnums = unique(chnums)'; assert(length(chnums)==2);
        rho = OUTSTRUCT.(fieldtoget)(indsthis, 1);
        rho_sem = OUTSTRUCT.(fieldtoget)(indsthis, 2);
%         cohscal = OUTSTRUCT.cohscal_diff(indsthis);
        
        assert(length(unique(motifID)) == length(motifID), 'each motif shoudl have at most 1 datapoint');
        
        % ---- sort in order of motifs
        [~, indsort] = sort(motifID);
        motifID = motifID(indsort);
%         cohscal = cohscal(indsort);
        rho = rho(indsort);
        rho_sem = rho_sem(indsort);
        
        indstmp = ~isnan(rho);
        x = motifID(indstmp) + XSHIFT;
        y = rho(indstmp);
        ysem = rho_sem(indstmp);
        if all(isnan(ysem))
        plot(x, y, '-o', 'Color', pcol);    
        else
        errorbar(x, y, ysem, '-o', 'Color', pcol);
        end
        lt_plot_text(motifID(end)+0.3, rho(end), num2str(chnums), 'k', 9);
        
        %% 3) ALL TRIOALS
                
        XSHIFT = 0.2;
        pcol = 'b';
        fieldtoget = 'cohscal_ff_rho_rhoSE_All';
        
        % =============
        motifID = lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)); % ---- get positions within global motif
        chnums = OUTSTRUCT.chanpair(indsthis,:);
        chnums = unique(chnums)'; assert(length(chnums)==2);
        rho = OUTSTRUCT.(fieldtoget)(indsthis, 1);
        rho_sem = OUTSTRUCT.(fieldtoget)(indsthis, 2);
%         cohscal = OUTSTRUCT.cohscal_diff(indsthis);
        
        assert(length(unique(motifID)) == length(motifID), 'each motif shoudl have at most 1 datapoint');
        
        % ---- sort in order of motifs
        [~, indsort] = sort(motifID);
        motifID = motifID(indsort);
%         cohscal = cohscal(indsort);
        rho = rho(indsort);
        rho_sem = rho_sem(indsort);
        
        indstmp = ~isnan(rho);
        x = motifID(indstmp) + XSHIFT;
        y = rho(indstmp);
        ysem = rho_sem(indstmp);
        if all(isnan(ysem))
        plot(x, y, '-o', 'Color', pcol);    
        else
        errorbar(x, y, ysem, '-o', 'Color', pcol);
        end
        lt_plot_text(motifID(end)+0.3, rho(end), num2str(chnums), 'k', 9);
        

        
        
    end
    
        % ----- labels
        lt_plot_zeroline;
[~, motiflabels] = lt_neural_QUICK_MotifID(bname);
    set(gca, 'XTick', 1:length(motiflabels));
    set(gca, 'XTickLabel', motiflabels);
    rotateXLabels(gca, 90);

    
    %% ======= PLOT MEANS OVER CHANNELS
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(swnum)]);
    ylabel(['rho (ff vs coh) [bk=base; rd=WN]']);

    %% ====== BASELINE
        XSHIFT = -0.2;
        pcol = 'k';
        fieldtoget = 'cohscal_ff_rho_rhoSE_base';
        
        indsthis = indsgrp_switch==swgrpthis;
    
    istarg = OUTSTRUCT.istarg(indsthis);
    issame = OUTSTRUCT.issame(indsthis);
    motifs = OUTSTRUCT.motifname(indsthis);
    motifID = lt_neural_QUICK_MotifID(bname, motifs); % ---- get positions within global motif
    %     cohscal = OUTSTRUCT.CohMean_WNminusBase_scalar(indsthis);
%     cohscal = OUTSTRUCT.cohscal_diff(indsthis);
    rho = OUTSTRUCT.(fieldtoget)(indsthis, 1);
    rho_sem = OUTSTRUCT.(fieldtoget)(indsthis, 2);
    
    % --- only keep those that have defined ff
    indstmp = ~isnan(rho);
    istarg = istarg(indstmp);
    issame = issame(indstmp);
    motifs = motifs(indstmp);
    motifID = motifID(indstmp);
    rho = rho(indstmp);
    rho_sem = rho_sem(indstmp); 

    % ==== collapse across channels
    [ymean, ysem] = grpstats(rho, motifID, {'mean', 'sem'});
    
%     [ymean, ysem] = grpstats(cohscal, motifID, {'mean', 'sem'});
    [istarg] = grpstats(istarg, motifID, {'mean'});
    [issame] = grpstats(issame, motifID, {'mean'});
    
    % ===== PLOT MEAN ACROSS CHANNEL PAIRS
    x = unique(motifID) + XSHIFT;
    lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', pcol, 'LineStyle', '-'});
    
    
    
    %% ====== WN
        XSHIFT = 0;
        pcol = 'r';
        fieldtoget = 'cohscal_ff_rho_rhoSE_WN';
        
        indsthis = indsgrp_switch==swgrpthis;
    
    istarg = OUTSTRUCT.istarg(indsthis);
    issame = OUTSTRUCT.issame(indsthis);
    motifs = OUTSTRUCT.motifname(indsthis);
    motifID = lt_neural_QUICK_MotifID(bname, motifs); % ---- get positions within global motif
    %     cohscal = OUTSTRUCT.CohMean_WNminusBase_scalar(indsthis);
%     cohscal = OUTSTRUCT.cohscal_diff(indsthis);
    rho = OUTSTRUCT.(fieldtoget)(indsthis, 1);
    rho_sem = OUTSTRUCT.(fieldtoget)(indsthis, 2);
    
    % --- only keep those that have defined ff
    indstmp = ~isnan(rho);
    istarg = istarg(indstmp);
    issame = issame(indstmp);
    motifs = motifs(indstmp);
    motifID = motifID(indstmp);
    rho = rho(indstmp);
    rho_sem = rho_sem(indstmp); 

    % ==== collapse across channels
    [ymean, ysem] = grpstats(rho, motifID, {'mean', 'sem'});
    
%     [ymean, ysem] = grpstats(cohscal, motifID, {'mean', 'sem'});
    [istarg] = grpstats(istarg, motifID, {'mean'});
    [issame] = grpstats(issame, motifID, {'mean'});
    
    % ===== PLOT MEAN ACROSS CHANNEL PAIRS
    x = unique(motifID) + XSHIFT;
    lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', pcol, 'LineStyle', '-'});
    

        %% ====== All
        XSHIFT = 0.2;
        pcol = 'b';
        fieldtoget = 'cohscal_ff_rho_rhoSE_All';
        
        indsthis = indsgrp_switch==swgrpthis;
    
    istarg = OUTSTRUCT.istarg(indsthis);
    issame = OUTSTRUCT.issame(indsthis);
    motifs = OUTSTRUCT.motifname(indsthis);
    motifID = lt_neural_QUICK_MotifID(bname, motifs); % ---- get positions within global motif
    %     cohscal = OUTSTRUCT.CohMean_WNminusBase_scalar(indsthis);
%     cohscal = OUTSTRUCT.cohscal_diff(indsthis);
    rho = OUTSTRUCT.(fieldtoget)(indsthis, 1);
    rho_sem = OUTSTRUCT.(fieldtoget)(indsthis, 2);
    
    % --- only keep those that have defined ff
    indstmp = ~isnan(rho);
    istarg = istarg(indstmp);
    issame = issame(indstmp);
    motifs = motifs(indstmp);
    motifID = motifID(indstmp);
    rho = rho(indstmp);
    rho_sem = rho_sem(indstmp); 

    % ==== collapse across channels
    [ymean, ysem] = grpstats(rho, motifID, {'mean', 'sem'});
    
%     [ymean, ysem] = grpstats(cohscal, motifID, {'mean', 'sem'});
    [istarg] = grpstats(istarg, motifID, {'mean'});
    [issame] = grpstats(issame, motifID, {'mean'});
    
    % ===== PLOT MEAN ACROSS CHANNEL PAIRS
    x = unique(motifID) + XSHIFT;
    lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', pcol, 'LineStyle', '-'});

    %% ##################### other things
    lt_plot_zeroline;
    YLIM = ylim;
    
    % -------- NOTE DOWN POSITION OF TARGET SYSL
    indsthis = indsgrp_switch==swgrpthis & OUTSTRUCT.istarg==1;
    xtarg = unique(lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)));
    plot(xtarg, YLIM(1)+0.02, '^r');
    
    % ------- NOTE POSITION OF SAME_TYPES
    indsthis = indsgrp_switch==swgrpthis & OUTSTRUCT.istarg==0 & OUTSTRUCT.issame==1;
    xtarg = unique(lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)));
    if ~isempty(xtarg)
        plot(xtarg, YLIM(1)+0.02, '^b');
    end
    
    
    
    % ----- labels
    [~, motiflabels] = lt_neural_QUICK_MotifID(bname);
    set(gca, 'XTick', 1:length(motiflabels));
    set(gca, 'XTickLabel', motiflabels);
    rotateXLabels(gca, 90);
%     ylim();
    
    
%     % ========= NOTE DOWN IF THIS IS AN EXPERIMENT TO GET
%     if ismember([bnum enum swnum], indtoget_b_e_s, 'rows')
%         lt_plot_annotation(3, 'expttoget', 'm');
%         
%         % ============= ALSO COLLECT DATA FOR SUMMARY PLOT
%         Y = nan(1,4);
%         
%         % - targ
%         Y(1) = nanmean(ymean(istarg==1));
%         
%         % - same
%         Y(2) = nanmean(ymean(istarg==0 & issame==1));
%         
%         % - diff
%         Y(3) = nanmean(ymean(istarg==0 & issame==0));
%         
%         % -- nontarg
%         Y(4) = nanmean(ymean(istarg==0));
%         
%         Yall = [Yall; Y];
%         
%         All_bname = [All_bname; bname];
%         All_ename = [All_ename; ename];
%         All_swnum = [All_swnum; swnum];
%         All_bnum = [All_bnum; bnum];
%     end
        % ============= ALSO COLLECT DATA FOR SUMMARY PLOT
        Y = nan(1,4);
        
        % - targ
        Y(1) = nanmean(ymean(istarg==1));
        
        % - same
        Y(2) = nanmean(ymean(istarg==0 & issame==1));
        
        % - diff
        Y(3) = nanmean(ymean(istarg==0 & issame==0));
        
        % -- nontarg
        Y(4) = nanmean(ymean(istarg==0));
        
        Yall = [Yall; Y];
        
        All_bname = [All_bname; bname];
        All_ename = [All_ename; ename];
        All_swnum = [All_swnum; swnum];
        All_bnum = [All_bnum; bnum];
end
