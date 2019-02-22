function [fignums_alreadyused, hfigs, figcount, hsplot] = ...
    lt_neural_POPLEARN_SylLocked_Over_sub3(DATSTRUCT, SummaryStruct, ...
    subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount, ...
    i, xtoplot, indsthis,hsplots, indsgrp, indsgrpU, plotWN)
%%

% ============= PLOT FR SMOOTH
if plotWN==1
    frmat = DATSTRUCT.FRsmooth_dat_WN(indsthis);
else
    frmat = DATSTRUCT.FRsmooth_dat(indsthis);
end
frmat = cell2mat(cellfun(@(x)x', frmat, 'UniformOutput', 0));
neurlist = DATSTRUCT.FRsmooth_neurID(indsthis);
neurlist = cell2mat(cellfun(@(x)x(:), neurlist, 'UniformOutput', 0));

motifID = unique(DATSTRUCT.motifID(indsthis));
bnum = unique(DATSTRUCT.bnum(indsthis));
bname = SummaryStruct.birds(bnum).birdname;

[~, motiflist_out] = ...
    lt_neural_QUICK_MotifID(bname);
motifname_common = motiflist_out{motifID};


bnum = unique(DATSTRUCT.bnum(indsthis));
fr_bregion = {SummaryStruct.birds(bnum).neurons(neurlist).NOTE_Location};

frx = DATSTRUCT.FRsmooth_t(indsthis);
frx = frx{1};


% ================= EXTRACT LFP
if plotWN==1
    lfpmat = DATSTRUCT.LFP_dat_WN(indsthis);
else
    lfpmat = DATSTRUCT.LFP_dat(indsthis);
end
xlen = min(cellfun(@(x)size(x,1), lfpmat)); % since some are not same length...
lfpmat = cellfun(@(x)x(1:xlen, :), lfpmat, 'UniformOutput', 0);
lfpmat = cell2mat(cellfun(@(x)x', lfpmat, 'UniformOutput', 0));

lfpx = DATSTRUCT.LFP_t(indsthis);
lfpx = lfpx{1}(1:xlen);

lfp_bregion = DATSTRUCT.LFP_bregions(indsthis);
lfp_bregion = [lfp_bregion{:}]';


% ############################## PLOT [LMAN]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
regionthis = 'LMAN';
pcol = [0.3 0.7 0.3];
title([bname  '-' motifname_common '-' regionthis]);
if plotWN==1
    ylabel('WN');
end

% ========================== FR SMOOTH
indsthis = strcmp(fr_bregion, regionthis);
frthis = frmat(indsthis,:);

plot(frx, frthis', 'Color', pcol);
frmean = mean(frthis);
frsem = lt_sem(frthis);
%     if sum(indsthis)>1
%         %         shadedErrorBar(frx, frmean, frsem, {'Color', pcol}, 1);
%         %         plot(frx, frmean, 'Color', pcol, 'LineWidth', 2);
%     end
% ---- PLOT LOWER, FOR OVERLAY
if sum(indsthis)>1
    shadedErrorBar(frx, frmean-4, frsem, {'Color', pcol}, 1);
else
    plot(frx, frmean-4, 'Color', pcol, 'LineWidth', 2);
end


% ============================= LFP
indsthis = strcmp(lfp_bregion, regionthis);
lfpthis = lfpmat(indsthis,:);

% --- flip sign, and add to mov eup plot
lfpthis = -lfpthis + 4;

plot(lfpx, lfpthis', 'Color', 'k');
ymean = mean(lfpthis);
ysem= lt_sem(lfpthis);
%     if sum(indsthis)>1
%         %         shadedErrorBar(lfpx, ymean, ysem, {'Color', 'k'}, 1);
%         %         plot(lfpx, ymean, 'Color', 'k', 'LineWidth', 2);
%     end
% --- PLOT LOWER, FOR OVERLAY
if sum(indsthis)>1
    shadedErrorBar(lfpx, ymean-8, ysem, {'Color', 'k'}, 1);
else
    plot(lfpx, ymean-8, 'Color', 'k', 'LineWidth', 2);
end
% plot(lfpx, ymean-8, 'Color', 'k', 'LineWidth', 2);


% ================ OVERLYA SYL PATCHES
indsthis = indsgrp==indsgrpU(i);
ons = DATSTRUCT.sylonsets(indsthis);
offs = DATSTRUCT.syloffsets(indsthis);

% -- only keep last ons at t=0 and one preceding
tmp = cell2mat(cellfun(@(x)find(abs(x)<0.001), ons, 'UniformOutput', 0));

onstmp = [];
for j=1:length(ons)
    try
        onstmp = [onstmp; ons{j}(tmp(j)-1:tmp(j))];
    catch err
        onstmp = [onstmp; ons{j}(tmp(j))];
    end
end

offstmp = [];
for j=1:length(ons)
    try
        offstmp = [offstmp; offs{j}(tmp(j)-1:tmp(j))];
    catch err
        offstmp = [offstmp; offs{j}(tmp(j))];
    end
end

ons = mean(onstmp,1);
offs = mean(offstmp, 1);

YLIM = ylim;
lt_neural_QUICK_PlotSylPatches(ons, offs, [YLIM(2)-0.5 YLIM(2)]);


% =========== format
axis tight;
xlim(xtoplot);
lt_plot_zeroline;
lt_plot_zeroline_vert




% ############################## PLOT [LMAN]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
regionthis = 'RA';
pcol = [0.8 0.2 0.2];

title([bname  '-' motifname_common '-' regionthis]);
if plotWN==1
    ylabel('WN');
end

% ========================== FR SMOOTH
indsthis = strcmp(fr_bregion, regionthis);
frthis = frmat(indsthis,:);

plot(frx, frthis', 'Color', pcol);
frmean = mean(frthis);
frsem = lt_sem(frthis);
%     if sum(indsthis)>1
%         %         shadedErrorBar(frx, frmean, frsem, {'Color', pcol}, 1);
%         %         plot(frx, frmean, 'Color', pcol, 'LineWidth', 2);
%     end
% ---- PLOT LOWER, FOR OVERLAY
if sum(indsthis)>1
    shadedErrorBar(frx, frmean-4, frsem, {'Color', pcol}, 1);
else
    plot(frx, frmean-4, 'Color', pcol, 'LineWidth', 2);
end


% ============================= LFP
indsthis = strcmp(lfp_bregion, regionthis);
lfpthis = lfpmat(indsthis,:);

% --- flip sign, and add to mov eup plot
lfpthis = -lfpthis + 4;

plot(lfpx, lfpthis', 'Color', 'k');
ymean = mean(lfpthis);
ysem= lt_sem(lfpthis);
%     if sum(indsthis)>1
%         %         shadedErrorBar(lfpx, ymean, ysem, {'Color', 'k'}, 1);
%         %         plot(lfpx, ymean, 'Color', 'k', 'LineWidth', 2);
%     end
% --- PLOT LOWER, FOR OVERLAY
if sum(indsthis)>1
    shadedErrorBar(lfpx, ymean-8, ysem, {'Color', 'k'}, 1);
else
    plot(lfpx, ymean-8, 'Color', 'k', 'LineWidth', 2);
end
% plot(lfpx, ymean-8, 'Color', 'k', 'LineWidth', 2);


% ================ OVERLYA SYL PATCHES
indsthis = indsgrp==indsgrpU(i);
ons = DATSTRUCT.sylonsets(indsthis);
offs = DATSTRUCT.syloffsets(indsthis);

% -- only keep last ons at t=0 and one preceding
tmp = cell2mat(cellfun(@(x)find(abs(x)<0.001), ons, 'UniformOutput', 0));

onstmp = [];
for j=1:length(ons)
    try
        onstmp = [onstmp; ons{j}(tmp(j)-1:tmp(j))];
    catch err
        onstmp = [onstmp; ons{j}(tmp(j))];
    end
end

offstmp = [];
for j=1:length(ons)
    try
        offstmp = [offstmp; offs{j}(tmp(j)-1:tmp(j))];
    catch err
        offstmp = [offstmp; offs{j}(tmp(j))];
    end
end

ons = mean(onstmp,1);
offs = mean(offstmp, 1);

YLIM = ylim;
lt_neural_QUICK_PlotSylPatches(ons, offs, [YLIM(2)-0.5 YLIM(2)]);

% =========== format
axis tight;
%     xlim([xtoplot(1)-0.02 xtoplot(2)+0.02]);
xlim(xtoplot);
lt_plot_zeroline;
lt_plot_zeroline_vert