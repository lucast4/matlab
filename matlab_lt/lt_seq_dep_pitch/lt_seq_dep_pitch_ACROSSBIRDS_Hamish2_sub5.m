function [fignums_alreadyused, hfigs, figcount, hsplot] = ...
    lt_seq_dep_pitch_ACROSSBIRDS_Hamish2_sub5(default_type, DATSTRUCT, ...
    indthis, Nrends, dd, ...
    subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount)
hsplots = [];

if strcmp(default_type, 'MUSC')
    PC_default = DATSTRUCT.All_PitchCont_BASE_MUSC(indthis).All_PCmat{dd};
elseif strcmp(default_type, 'PBS')
    PC_default = DATSTRUCT.All_PitchCont_BASE_PBS(indthis).All_PCmat{dd};
end


%%

PC_dat = DATSTRUCT.All_PitchCont_BASE_PBS(indthis).All_PCmat{dd};
ffmat = DATSTRUCT.All_PitchCont_BASE_PBS(indthis).All_ffvals{dd};
tthis = DATSTRUCT.All_PitchCont_BASE_PBS(indthis).All_tbins{dd};
twind = DATSTRUCT.All_PitchCont_BASE_PBS(indthis).All_twind(dd,:);
twindsamp = twind;

if ~any(isnan(twind))

    twind = tthis(twind);

% ==================== EXTRACT HI AND LOW PITCH VARIANTS
[~, indsort] = sort(ffmat);
% ---- hi pitch
pcthis_PBS_hi = PC_dat(indsort(end-Nrends+1:end), :);
% ---- lo pitch
pcthis_PBS_lo = PC_dat(indsort(1:Nrends), :);


% ------------ 1) overlay means and std
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('hi/lo/default');


% ---- hi pitch
shadedErrorBar(tthis, mean(pcthis_PBS_hi,1), lt_sem(pcthis_PBS_hi), ...
    {'Color', [0.7 0.7 0.7]}, 1);
% ---- lo pitch
shadedErrorBar(tthis, mean(pcthis_PBS_lo,1), lt_sem(pcthis_PBS_lo), ...
    {'Color', [0.7 0.7 0.7]}, 1);
% ---- default
shadedErrorBar(tthis, mean(PC_default,1), lt_sem(PC_default), ...
    {'Color', 'r'}, 1);

axis tight;
line([twind(1) twind(1)], ylim);
line([twind(2) twind(2)], ylim);
xlim([twind(1)-0.01 twind(2)+0.01]);


% ----------- 2) subtract default, plot trials

% -- LO
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('lo (minus default)');

pctmp = pcthis_PBS_lo - mean(PC_default,1);
plot(tthis, pctmp, '-k');
% axis tight;
line([twind(1) twind(1)], ylim);
line([twind(2) twind(2)], ylim);
xlim([twind(1)-0.01 twind(2)+0.01]);

% --- HI
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('hi (minus default)');

pctmp = pcthis_PBS_hi - mean(PC_default,1);
plot(tthis, pctmp, '-k');
% axis tight;
line([twind(1) twind(1)], ylim);
line([twind(2) twind(2)], ylim);
xlim([twind(1)-0.01 twind(2)+0.01]);

% ------------ 3) THEN, subtract within trial means. (to get wiggle)

% --- LO
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('LO (minus default, centered)');

pctmp = pcthis_PBS_lo - mean(PC_default,1);
pctmp = pctmp - mean(pctmp(:, twindsamp(1):twindsamp(2)),2);

plot(tthis, pctmp, '-k');
% axis tight;
line([twind(1) twind(1)], ylim);
line([twind(2) twind(2)], ylim);
xlim([twind(1)-0.01 twind(2)+0.01]);

% --- HI
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('hi (minus default, centered)');

pctmp = pcthis_PBS_hi - mean(PC_default,1);
pctmp = pctmp - mean(pctmp(:, twindsamp(1):twindsamp(2)),2);

plot(tthis, pctmp, '-k');
% axis tight;
line([twind(1) twind(1)], ylim);
line([twind(2) twind(2)], ylim);
xlim([twind(1)-0.01 twind(2)+0.01]);




% #################### OVERLAY EXTRACTED STATS FOR WIGGLE, ETC
if strcmp(default_type, 'PBS')
    % ------------------------- FORMAT FIGURES
linkaxes(hsplots, 'xy');
ylim([-200 200]);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    xlabel('AFPBIAS - LO(wiggle) - HI(wiggle)');
    ylabel('hz');
    assert(length(DATSTRUCT.BaseBiasAll) == length(DATSTRUCT.All_Birdnum), 'need to run wiggle analyses on all data');
    
    Y = [DATSTRUCT.BaseBiasAll(indthis) cell2mat(DATSTRUCT.WiggleAll(indthis, :))];
    lt_plot_bar([1 2 3], Y);
    lt_plot_zeroline;
    xlim([0 4]);
    % --- put marker over bar expected to be larger
    if Y(1)<0
        % then expect lo renditions to have larger wiggle
        plot(2, 0, 'xr');
        lt_plot(2, 0, {'Color', 'r'});
    elseif Y(1)>0
        plot(3, 0, 'xr');
        lt_plot(3, 0, {'Color', 'r'});
    end
    
elseif strcmp(default_type, 'MUSC')
    
    % ============ overlay wiggle of MUSC
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('MUSC');
    
    % --- take random subset
    if size(PC_default,1)>=Nrends
        
        indrand = randperm(size(PC_default,1), Nrends);
        pctmp = PC_default(indrand, :);
        
        pctmp = pctmp - mean(PC_default,1);
        pctmp = pctmp - mean(pctmp(:, twindsamp(1):twindsamp(2)),2);
        
        plot(tthis, pctmp, '-r');
        % axis tight;
        line([twind(1) twind(1)], ylim);
        line([twind(2) twind(2)], ylim);
        xlim([twind(1)-0.01 twind(2)+0.01]);
    end
    % ------------------------- FORMAT FIGURES
linkaxes(hsplots, 'xy');
ylim([-200 200]);
end
end


