function lt_neural_Coher_COHSCALAR_Overview(OUTSTRUCT, FFcorrCoh, FFcorrCoh_pval, ...
    FFcorrCoh_pctileVsShuff, FFcorrCoh_zscoreVsShuff, FFcorrCoh_shuffCI, ...
    SwitchStruct, onlygood, PARAMS, plotEachSylChan, useshuffpval)
%% lt 1/8/19 - overview plots of coherence scalar (systematic, for each t, ff bin)

clim = [-0.3 0.3];
combineChans=1; % 1, combines across chan pairs (simult).

% === for distance matrix
tlims = [-0.1 0]; % limit matric to this.
flims = [20 125];


%% ========= recalculatuion correlatio p val using shuflfed?
if useshuffpval==1
    if (0) % NOTE: This preextract p val seems to be incorrect..
%    FFcorrCoh_pval = FFcorrCoh_pctileVsShuff<0.025 | FFcorrCoh_pctileVsShuff>0.975;
    
    else % kluge, here p wil be either 0 or 1.
   nrows = size(FFcorrCoh,1);
   ncols = size(FFcorrCoh, 2);
   for rr=1:nrows
       for cc=1:ncols
        
           tmp = squeeze(FFcorrCoh_shuffCI(rr, cc, :, :));
           corrthis = squeeze(FFcorrCoh(rr, cc, :))';
           
           FFcorrCoh_pval(rr, cc, :) = ~(corrthis<tmp(1,:) | corrthis>tmp(2,:)); % p = 0 or 1
       end
   end
    end
end
%% filter data?

if onlygood==1
    [OUTSTRUCT, indsgood] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT, SwitchStruct);
end


FFcorrCoh = FFcorrCoh(:,:, indsgood);
FFcorrCoh_pval = FFcorrCoh_pval(:,:, indsgood);
FFcorrCoh_pctileVsShuff = FFcorrCoh_pctileVsShuff(:,:, indsgood);
FFcorrCoh_zscoreVsShuff = FFcorrCoh_zscoreVsShuff(:,:, indsgood);
FFcorrCoh_shuffCI = FFcorrCoh_shuffCI(:,:,:, indsgood);


%% only keep those that have FF defined
indsgood = [];
for i=1:size(FFcorrCoh, 3)
    tmp = FFcorrCoh(:,:, i);
    if ~all(isnan(tmp(:)))
        indsgood = [indsgood; i];
    end
end

OUTSTRUCT = lt_structure_subsample_all_fields(OUTSTRUCT, indsgood, 1);
FFcorrCoh = FFcorrCoh(:,:, indsgood);
FFcorrCoh_pval = FFcorrCoh_pval(:,:, indsgood);
FFcorrCoh_pctileVsShuff = FFcorrCoh_pctileVsShuff(:,:, indsgood);
FFcorrCoh_zscoreVsShuff = FFcorrCoh_zscoreVsShuff(:,:, indsgood);
FFcorrCoh_shuffCI = FFcorrCoh_shuffCI(:,:,:, indsgood);



%% ======= ONE PLOT FOR EACH EXPERIMENT/SYLLABLES/CHANPAIR.

if plotEachSylChan ==1
    [indsgrp, indsgrpunique] = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch});
    % [indsgrp, indsgrpunique] = lt_tools_grp2idx(});
    
    figcount=1;
    subplotrows=6;
    subplotcols=3;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    for i=1:length(indsgrpunique)
        
        indsthis = indsgrp == indsgrpunique(i);
        cormat = FFcorrCoh(:,:, indsthis);
        pvalmat = FFcorrCoh_pval(:,:, indsthis);
        
        bnum = unique(OUTSTRUCT.bnum(indsthis));
        bname = SwitchStruct.bird(bnum).birdname;
        enum = unique(OUTSTRUCT.enum(indsthis));
        ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
        swnum = unique(OUTSTRUCT.switch(indsthis));
        motiflist = OUTSTRUCT.motifname(indsthis);
        chanpairlist = OUTSTRUCT.chanpair(indsthis,:);
        
        % ==== one plot for each syl/chanpair
        for j=1:size(cormat,3)
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title({[bname '-' ename '-sw' num2str(swnum)], [motiflist{j} ' - ' num2str(chanpairlist(j,:))]});
            
            corthis = cormat(:,:, j);
            pthis = pvalmat(:,:, j);
            
            lt_neural_Coher_Plot(corthis, PARAMS.tbins, PARAMS.ffbins, 1, '', clim, 0, 0, PARAMS.ffbinsedges);
            
            % --- mark places with pval <0.05
            for rr=1:length(PARAMS.tbins)
                for cc=1:length(PARAMS.ffbins)
                    if pthis(rr, cc)<0.05
                        plot(PARAMS.tbins(rr), PARAMS.ffbins(cc), 'xy');
                    end
                    
                end
            end
            
            
            ylim([15 130]);
        end
        
        tmp = input('close? (y or n)', 's');
        %     pause;
        if strcmp(tmp, 'y')
            close all;
        else
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        end
        
    end
    
end



%% ====== MEAN OVER CHANNEL PAIRS - THEY LOOK VERY SIMILAR.

%% ====== FLATTEN ALL - THEN GET STATS FOR EACH BIN

% ============= COMBINE BIRDS
birdthis = [];
lt_neural_Coher_COHSCALAR_Overview_flat(birdthis, OUTSTRUCT, FFcorrCoh, ...
    FFcorrCoh_pval, FFcorrCoh_pctileVsShuff, FFcorrCoh_zscoreVsShuff, FFcorrCoh_shuffCI, ...
    PARAMS, tlims, flims, combineChans);
lt_subtitle('all bird combined');


% ============= DO SEPARATELY FOR EACH BIRD
for i=1:max(OUTSTRUCT.bnum)
    birdthis = i;
    lt_neural_Coher_COHSCALAR_Overview_flat(birdthis, OUTSTRUCT, FFcorrCoh, ...
        FFcorrCoh_pval, FFcorrCoh_pctileVsShuff, FFcorrCoh_zscoreVsShuff, FFcorrCoh_shuffCI, ...
        PARAMS, tlims, flims, combineChans);
    birdname = SwitchStruct.bird(i).birdname;
    lt_subtitle(birdname);
end


%% ======= GET DISTANCE MATRIX BETWEEN ALL CASES

% ============= DO SEPARATELY FOR EACH BIRD
birdthis = 1;

lt_neural_Coher_COHSCALAR_Overview_sub(birdthis, OUTSTRUCT, FFcorrCoh, ...
    PARAMS, tlims, flims);


%% =============== PICK ONE SYL, PLOT ITS DISTANCE TO EVERY OTHER SYL
seedmotif = 'aa(b)hh';
locthis = 'RA';
distmetric = 'correlation';

indneur = strcmp(AllNeurLocation, locthis);
indmotif = strcmp(AllMotifRegexp, seedmotif);

frmat = FRmatMotifByNeur(:,indneur);

if strcmp(distmetric, 'correlation')
    
    % === get correlation between all motifs
    rhomat = corr(frmat');
    corrvalues = rhomat(indmotif, :);
    
    lt_figure; hold on;
    lt_plot_bar(1:length(corrvalues), corrvalues);
    
    set(gca, 'XTick', 1:length(AllMotifRegexp), 'XTickLabel', AllMotifRegexp);
    title([locthis ', pop corr to ' seedmotif]);
end





