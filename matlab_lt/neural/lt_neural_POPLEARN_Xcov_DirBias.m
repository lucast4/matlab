function lt_neural_POPLEARN_Xcov_DirBias(OUTSTRUCT_XCOV, PARAMS, ...
    onlygoodexpt, SwitchStruct, dirstruct, useff_zscore, targsylonly, ...
    cvdiffmethod, giveFakeAFPbias)
%% lt 3/5/19 - divides up training into epochs -- here PLOTS

if giveFakeAFPbias==1
    cvdiffmethod = 'diff'; % because std will not be defined
    useff_zscore = 0;
end
%% ====== [PREPROCESS] only plot good experiments

if onlygoodexpt==1
    
    % ===== filter outstruct
    [OUTSTRUCT_XCOV] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT_XCOV, ...
        SwitchStruct, 'xcov_spikes');
    %
    %     % ===== filter outstruct
    %     [OUTSTRUCT indstokeep] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT, ...
    %         SwitchStruct, 'xcov_spikes');
    
end



%% ========= SUMMARY PLOT OF ALL DIRECTED SONG
nbirds = length(dirstruct.bird);
for i=1:nbirds
    
    nmotifs = length(dirstruct.bird(i).DAT.motifID);
end


%% ========= FOR EACH MOTIF (UNIQUE ID) COLLECT AFP BIAS EFFECT
[indsgrp, indsgrpU] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, OUTSTRUCT_XCOV.motifID_unique});

figcount=1;
subplotrows=5;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


All_neuralHiMinusLo = {};
All_neuralXcovBase = {};

All_ffUndirMinusDir = {};
All_ffUndir_cv = {};
All_ffDir_cv = {};
All_ffUndirOverDir_cv = {};
All_bnum = [];
All_motifID = [];

for i=1:length(indsgrpU)
    
    if targsylonly==0
        indsthis = indsgrp==indsgrpU(i);
    elseif targsylonly==1
        indsthis = indsgrp==indsgrpU(i) & OUTSTRUCT_XCOV.istarg==1;
    end
    
    if ~any(indsthis)
        continue
    end
    
    % ===== NEURAL
    assert(all(size(OUTSTRUCT_XCOV.xcovscalBase_window_FFsplit{1})==[2 1]), 'have not confinred that dimensions are correct if is more than 2x1')
    yneur = cellfun(@(x)x(2)-x(1), OUTSTRUCT_XCOV.xcovscalBase_window_FFsplit(indsthis)); % each case, get all neuraldifferences
    
    % -------------- magnitude of xcov
    if (1)
        yxcov_base = cellfun(@(x)x(1), OUTSTRUCT_XCOV.Xcovscal_window_BaseWN(indsthis));
    else
        yxcov_base = cellfun(@(x)mean(x), OUTSTRUCT_XCOV.xcovscalBase_window_FFsplit(indsthis));
    end
    
    % ======= DIRECTED SONG
    bnum = unique(OUTSTRUCT_XCOV.bnum(indsthis));
    motifID = unique(OUTSTRUCT_XCOV.motifID_unique(indsthis));
    motifname = unique(OUTSTRUCT_XCOV.motifname(indsthis));
    
    bname = SwitchStruct.bird(bnum).birdname;
    ind1 = strcmp({dirstruct.bird.birdname}, bname);
    
    ff = [dirstruct.bird(ind1).DAT.motifID(motifID).rendnum.ff];
    t = [dirstruct.bird(ind1).DAT.motifID(motifID).rendnum.datenum_song_SecRes];
    isDir = [dirstruct.bird(ind1).DAT.motifID(motifID).rendnum.isDIR];
    
    
    if all(isnan(ff))
    if giveFakeAFPbias==0
        continue
    elseif giveFakeAFPbias==1
        ff = zeros(size(ff));
        
        isDir(1:2)=1; % need at least some dir data to get afp bias.
    end
    end
    
    % === sort by t (for fun)
    [~, indsort] = sort(t);
    t = t(indsort);
    ff = ff(indsort);
    isDir = isDir(indsort);
    motifname = motifname{1};
    assert(length(motifID)==1);
%     motifname = motifname{1};
    
    if rand<1.1
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%         title([bname '-mot' num2str(motifID)]);
        title([bname '-' motifname]);
        plot(t(isDir==0), ff(isDir==0), 'ok');
        plot(t(isDir==1), ff(isDir==1), 'ob');
        datetick('x', 'mm/dd');
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%         title([bname '-mot' num2str(motifID)]);
        title([bname '-' motifname]);
        
        x = 1:length(t);
        plot(x(isDir==0), ff(isDir==0), 'ok');
        plot(x(isDir==1), ff(isDir==1), 'ob');
    end
    
    
    % ========== for each day, get mean DIR and UNDIR (difference)
    t_day = floor(t);
    
    [ff_day_UNDIR, ff_day_UNDIR_std] = grpstats(ff(isDir==0), t_day(isDir==0), {'mean', 'std'});
    [ff_day_DIR, ff_day_DIR_std] = grpstats(ff(isDir==1), t_day(isDir==1), {'mean', 'std'});
    %     ff_day_UNDIR = grpstats(ff(isDir==0), t_day(isDir==0), {'median'});
    %     ff_day_DIR = grpstats(ff(isDir==1), t_day(isDir==1), {'median'});
    assert(length(ff_day_UNDIR)==length(ff_day_DIR));
    
    ff_UndirMinusDir = ff_day_UNDIR - ff_day_DIR;
    
    %     ff_UndirMinusDir_z = [];
    ff_UndirMinusDir_z = -(ff_day_DIR - ff_day_UNDIR)./ff_day_UNDIR_std;
    
    
    % ==== CV UNDIR AND DIR
    cv_UNDIR = ff_day_UNDIR_std./ff_day_UNDIR;
    cv_DIR = ff_day_DIR_std./ff_day_DIR;
    
    All_ffUndir_cv = [All_ffUndir_cv; cv_UNDIR];
    All_ffDir_cv = [All_ffDir_cv; cv_DIR];

    if strcmp(cvdiffmethod, 'diff')
        cvdiff_UndirOverDir = cv_UNDIR-cv_DIR;
    elseif strcmp(cvdiffmethod, 'percent')
        cvdiff_UndirOverDir = (cv_UNDIR-cv_DIR)./cv_UNDIR;
    end
    All_ffUndirOverDir_cv = [All_ffUndirOverDir_cv; cvdiff_UndirOverDir];
    
    % =============== SAVE
    if useff_zscore==2
        All_ffUndirMinusDir = [All_ffUndirMinusDir; (ff_day_UNDIR - ff_day_DIR)./ff_day_UNDIR];
    elseif useff_zscore==1
        All_ffUndirMinusDir = [All_ffUndirMinusDir; ff_UndirMinusDir_z];
    elseif useff_zscore==0
        All_ffUndirMinusDir = [All_ffUndirMinusDir; ff_UndirMinusDir];
    end
    All_neuralHiMinusLo = [All_neuralHiMinusLo; yneur];
    All_bnum = [All_bnum; bnum];
    All_motifID = [All_motifID; motifID];
    
    All_neuralXcovBase = [All_neuralXcovBase; yxcov_base];
end


%% ============ CORRELATION BETWEEN NERUAL XCOV AND CHANGE IN PITCH CV?

lt_figure; hold on;

lt_subplot(3,2,1); hold on;
title('mean over channels/days');
xlabel('ff CV (Undir - Dir)');
ylabel('mean xcov (baseline)');

x = cellfun(@nanmean, All_ffUndirOverDir_cv);
y = cellfun(@nanmean, All_neuralXcovBase);
% x = cellfun(@nanmedian, All_ffUndirMinusDir);
% y = cellfun(@nanmedian, All_neuralHiMinusLo);
% plot(x, y, 'ok');

% -- diff color per bird
pcols = lt_make_plot_colors(max(All_bnum), 0,0);
for i=1:max(All_bnum)
    inds = All_bnum==i;
    lt_plot(x(inds), y(inds), {'Color', pcols{i}});
end
lt_regress(y, x, 1, 0, 0, 1);
lt_plot_zeroline;
lt_plot_zeroline_vert
lt_plot_annotation(1, ['method: ' cvdiffmethod]);


lt_subplot(3,2,2); hold on;
title('mean over channels/days');
xlabel('ff CV (Undir)');
ylabel('mean xcov (baseline)');

x = cellfun(@nanmean, All_ffUndir_cv);
y = cellfun(@nanmean, All_neuralXcovBase);
% x = cellfun(@nanmedian, All_ffUndirMinusDir);
% y = cellfun(@nanmedian, All_neuralHiMinusLo);
% plot(x, y, 'ok');

% -- diff color per bird
pcols = lt_make_plot_colors(max(All_bnum), 0,0);
for i=1:max(All_bnum)
    inds = All_bnum==i;
    lt_plot(x(inds), y(inds), {'Color', pcols{i}});
end
lt_regress(y, x, 1, 0, 0, 1);
lt_plot_zeroline;
lt_plot_zeroline_vert


%% =========== CORRELATION BETWEN FF AND MAGNTIUDE OF XCOV?
lt_figure; hold on;

lt_subplot(3,2,1); hold on;
title('mean over channels/days');
xlabel('ff (Undir - Dir)');
ylabel('mean xcov (baseline)');

x = cellfun(@nanmean, All_ffUndirMinusDir);
y = cellfun(@nanmean, All_neuralXcovBase);
% x = cellfun(@nanmedian, All_ffUndirMinusDir);
% y = cellfun(@nanmedian, All_neuralHiMinusLo);
% plot(x, y, 'ok');

% -- diff color per bird
pcols = lt_make_plot_colors(max(All_bnum), 0,0);
for i=1:max(All_bnum)
    inds = All_bnum==i;
    lt_plot(x(inds), y(inds), {'Color', pcols{i}});
end
lt_regress(y, x, 1, 0, 0, 1);
lt_plot_zeroline;
lt_plot_zeroline_vert

lt_subplot(3,2,2); hold on;
title('mean over channels/days');
xlabel('ff (Undir - Dir) [abs val]');
ylabel('mean xcov (baseline)');

x = abs(cellfun(@nanmean, All_ffUndirMinusDir));
y = cellfun(@nanmean, All_neuralXcovBase);
% x = cellfun(@nanmedian, All_ffUndirMinusDir);
% y = cellfun(@nanmedian, All_neuralHiMinusLo);
% plot(x, y, 'ok');

% -- diff color per bird
pcols = lt_make_plot_colors(max(All_bnum), 0,0);
for i=1:max(All_bnum)
    inds = All_bnum==i;
    lt_plot(x(inds), y(inds), {'Color', pcols{i}});
end
lt_regress(y, x, 1, 0, 0, 1);
lt_plot_zeroline;
lt_plot_zeroline_vert


%% ========= CORRELATION BETWEEN FF AND BIAS?

lt_figure; hold on;
title('mean over channels/days');
xlabel('ff (Undir - Dir)');
ylabel('neural bias (xcov(hiFF) - xcov(loFF) [baseline]');

x = cellfun(@nanmean, All_ffUndirMinusDir);
y = cellfun(@nanmean, All_neuralHiMinusLo);
% x = cellfun(@nanmedian, All_ffUndirMinusDir);
% y = cellfun(@nanmedian, All_neuralHiMinusLo);
% plot(x, y, 'ok');

% -- diff color per bird
pcols = lt_make_plot_colors(max(All_bnum), 0,0);
for i=1:max(All_bnum)
    inds = All_bnum==i;
    lt_plot(x(inds), y(inds), {'Color', pcols{i}});
end
lt_regress(y, x, 1, 0, 0, 1);
lt_plot_zeroline;
lt_plot_zeroline_vert

%% ========= SAME SYLLABLE ACROSS EXPERIMENTS, SIMILAR?

[indsgrp, indsgrpU] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, OUTSTRUCT_XCOV.motifID_unique});
indsgrp_Expt = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, OUTSTRUCT_XCOV.enum, ...
    OUTSTRUCT_XCOV.switch});

figcount=1;
subplotrows=5;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


for i=1:length(indsgrpU)
    indsthis = indsgrp==indsgrpU(i);
    
    yneur = cellfun(@(x)x(2)-x(1), OUTSTRUCT_XCOV.xcovscalBase_window_FFsplit(indsthis));
    bnum = unique(OUTSTRUCT_XCOV.bnum(indsthis));
    mID = unique(OUTSTRUCT_XCOV.motifID_unique(indsthis));
    motifname = OUTSTRUCT_XCOV.motifname(indsthis);
    motifname = motifname{1};
    %    istarg = OUTSTRUCT_XCOV.istarg(indsthis);
    exptID = indsgrp_Expt(indsthis);
    
    % ==== collect ff
    ff = nan(size(mID));
    for j=1:length(mID)
        
        indstmp = find(All_bnum==bnum & All_motifID==mID(j));
        if isempty(indstmp)
            continue
        end
        assert(length(indstmp)==1);
        ff(j) = mean(All_ffUndirMinusDir{indstmp});
    end
    if all(isnan(ff))
        ff = nan;
    else
        ff = unique(ff);
    end
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(['bird' num2str(bnum) '-' motifname]);
    xlabel('recording session');
    lt_plot_annotation(1, ['UndirMinusDir=' num2str(unique(ff))], 'm');
    ylabel('neural bias (xcov(hiFF) - xcov(loFF) [baseline]');
    plot(exptID, yneur, 'ok');
    yneur2 = grpstats(yneur, exptID);
    exptID2 = unique(exptID);
    lt_plot(exptID2, yneur2, {'Color', 'r'});
    xlim([0 max(exptID)+1]);
    lt_plot_zeroline;
    
end
linkaxes(hsplots, 'y');

%% ======== WITHIN EACH EXPERIMENT, SEE RELATINOSHIP?
figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


[indsgrp, indsgrpU] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, OUTSTRUCT_XCOV.enum, OUTSTRUCT_XCOV.switch});
for i=1:length(indsgrpU)
   indsthis = indsgrp==indsgrpU(i);
   
   yneur = cellfun(@(x)x(2)-x(1), OUTSTRUCT_XCOV.xcovscalBase_window_FFsplit(indsthis));
   bnum = unique(OUTSTRUCT_XCOV.bnum(indsthis));
   enum = unique(OUTSTRUCT_XCOV.enum(indsthis));
   mID = OUTSTRUCT_XCOV.motifID_unique(indsthis);
   motifname = OUTSTRUCT_XCOV.motifname(indsthis);
   istarg = OUTSTRUCT_XCOV.istarg(indsthis);
   
   
   ff = nan(size(mID));
   for j=1:length(mID)
      
       indstmp = find(All_bnum==bnum & All_motifID==mID(j));
       if isempty(indstmp)
           continue
       end
       assert(length(indstmp)==1);
       ff(j) = mean(All_ffUndirMinusDir{indstmp});
   end
   
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title(['bird' num2str(bnum) '-expt' num2str(enum)]);
xlabel('ff (Undir - Dir)');
ylabel('neural bias (xcov(hiFF) - xcov(loFF) [baseline]');
plot(ff, yneur, 'or');
lt_plot_zeroline;
lt_plot_zeroline_vert;

% -- label motif
for j=1:length(mID)
    if istarg(j)==1
    lt_plot_text(ff(j)+0.01, yneur(j), motifname{j}, 'r', 8);
    else
        lt_plot_text(ff(j)+0.01, yneur(j), motifname{j}, 'm', 8);
    end
end

end
linkaxes(hsplots, 'xy');


