pcols = lt_make_plot_colors(length(indsunique), 0, 0);

ytmp = nan(length(indsunique), length(typetoget));

for i=1:length(indsunique)
    indsthis = indsgrp==indsunique(i); % a unique motif.
    
    coh = All_cohscal(indsthis);
    istarg = All_istarg(indsthis);
    issame = All_issame(indsthis);
    xval = All_xval(indsthis);
    
    motifthis = unique(All_motifname(indsthis)); assert(length(motifthis)==1);
    birdname = SwitchStruct.bird(unique(All_bnum(indsthis))).birdname;
    
    pcol = pcols{i};
    
    if meanOverExpt==1
        [coh, cohsem] = grpstats(coh, xval, {'mean', 'sem'});
        xval = unique(xval);
    else
        cohsem = [];
    end
    
    %    disp(xval)
    %    if ~all(ismember(xval, typetoget'))
    %        continue
    %    end
    if ~isempty(setxor(xval, typetoget'))
        continue
    end
    if isempty(cohsem)
        plot(xval+0.4*rand-0.2, coh, '-o', 'Color', pcol);
    else
        lt_plot(xval+0.4*rand-0.2, coh, {'Color', pcol, 'Errors', cohsem, 'LineStyle', '-'});
    end
    
    
    % ---- plot text id
    [~, tmp] = max(xval);
    lt_plot_text(max(xval)+0.3, coh(tmp), [birdname '-' motifthis{1}], pcol);
    
    % --- collect
    ytmp(i,xval) = coh;
end
xlim([0 5]);
lt_plot_zeroline;
if strcmp(statmethod, 'zscore')
    ylim([-2 2]);
elseif strcmp(statmethod, 'diff')
    ylim([-0.15 0.15]);
end

if size(ytmp,2)==2
   [~, p] = ttest(ytmp(:,1), ytmp(:,2))
   lt_plot_pvalue(p, 'ttest',1);
end
