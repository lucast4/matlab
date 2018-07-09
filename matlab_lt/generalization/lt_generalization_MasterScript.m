%% lt 4/4/18 - looking at timecourse of generalization
%% ===========================

%% ####################### PREPROCESS
close all;
GenStruct = lt_generalization_PreProc;


%% ############################# PLOT RAW DAT
%% PLOT EACH EXPERIMENT (OVERLAY TARGET VS. NONTARGETS)
% OLD VERSION: USES RENDITION BY RENDITION, AND NOT COMBINED SYLS
j = 3; % expt to plot

figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% ===== extract learning window times
tOn_learn = GenStruct.expt(j).Params.Datenum_LearnWind_On;
tOff_learn = GenStruct.expt(j).Params.Datenum_LearnWind_Off;

% ============================================== targ
titlename = 'targ';
plotcol = 'k';

% -------
indmotall = find(GenStruct.expt(j).Params.IsTarg==1);
for k=1:length(indmotall)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots, hsplot];
    indmot = indmotall(k);
    
    tvals = [GenStruct.expt(j).DAT.motif(indmot).rendnum.datenum_song_SecRes];
    ffvals = [GenStruct.expt(j).DAT.motif(indmot).rendnum.ff];
    motifname = GenStruct.expt(j).DAT.motif(indmot).motif;
    
    title([titlename '-' motifname]);
    plot(tvals, ffvals, 'o', 'Color', plotcol);
    % ----- overlay song by song data
    ffsong = GenStruct.expt(j).DAT.SongByMotifMat(:, indmot) + ...
        GenStruct.expt(j).DAT.SongByMotif_baseMeans(indmot);
    tsong = datenum(GenStruct.expt(j).DAT.SongByMotif_songnames, 'yymmdd_HHMMSS');
    lt_plot(tsong, ffsong, {'Color', plotcol});
    % ----- overlay onset and offset of learning window
    line([tOn_learn tOn_learn], ylim, 'Color', 'b', 'LineStyle', '--');
    line([tOff_learn tOff_learn], ylim, 'Color', 'b', 'LineStyle', '--');
end


% ========================================= same type
titlename = 'same';
plotcol = 'b';

% -------
indmotall = find(GenStruct.expt(j).Params.IsTarg==0 & ...
    GenStruct.expt(j).Params.IsSame==1);
for k=1:length(indmotall)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots, hsplot];
    indmot = indmotall(k);
    
    tvals = [GenStruct.expt(j).DAT.motif(indmot).rendnum.datenum_song_SecRes];
    ffvals = [GenStruct.expt(j).DAT.motif(indmot).rendnum.ff];
    motifname = GenStruct.expt(j).DAT.motif(indmot).motif;
    
    title([titlename '-' motifname]);
    plot(tvals, ffvals, 'o', 'Color', plotcol);
    % ----- overlay song by song data
    ffsong = GenStruct.expt(j).DAT.SongByMotifMat(:, indmot) + ...
        GenStruct.expt(j).DAT.SongByMotif_baseMeans(indmot);
    tsong = datenum(GenStruct.expt(j).DAT.SongByMotif_songnames, 'yymmdd_HHMMSS');
    lt_plot(tsong, ffsong, {'Color', plotcol});
    % ----- overlay onset and offset of learning window
    line([tOn_learn tOn_learn], ylim, 'Color', 'b', 'LineStyle', '--');
    line([tOff_learn tOff_learn], ylim, 'Color', 'b', 'LineStyle', '--');
end

% ========================================= diff type
titlename = 'diff';
plotcol = 'r';

% -------
indmotall = find(GenStruct.expt(j).Params.IsTarg==0 & ...
    GenStruct.expt(j).Params.IsSame==0);
for k=1:length(indmotall)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots, hsplot];
    indmot = indmotall(k);
    
    tvals = [GenStruct.expt(j).DAT.motif(indmot).rendnum.datenum_song_SecRes];
    ffvals = [GenStruct.expt(j).DAT.motif(indmot).rendnum.ff];
    motifname = GenStruct.expt(j).DAT.motif(indmot).motif;
    
    title([titlename '-' motifname]);
    plot(tvals, ffvals, 'o', 'Color', plotcol);
    % ----- overlay song by song data
    ffsong = GenStruct.expt(j).DAT.SongByMotifMat(:, indmot) + ...
        GenStruct.expt(j).DAT.SongByMotif_baseMeans(indmot);
    tsong = datenum(GenStruct.expt(j).DAT.SongByMotif_songnames, 'yymmdd_HHMMSS');
    lt_plot(tsong, ffsong, {'Color', plotcol});
    % ----- overlay onset and offset of learning window
    line([tOn_learn tOn_learn], ylim, 'Color', 'b', 'LineStyle', '--');
    line([tOff_learn tOff_learn], ylim, 'Color', 'b', 'LineStyle', '--');
end


% =====
linkaxes(hsplots, 'x');


%% PLOT RAW DAT - NEW VERSION, USING COMBINED SYLS AND SONG BY SONG

j = 3; % expt to plot

figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% ===== extract learning window times
tOn_learn = GenStruct.expt(j).Params.Datenum_LearnWind_On;
tOff_learn = GenStruct.expt(j).Params.Datenum_LearnWind_Off;
motiflist = GenStruct.expt(j).DAT_MotifRenamed.SongByMotif_motifnames;

for k=1:length(motiflist)
    motifthis = motiflist{k};
   
    % =============== is this targ? same?
    istarg = GenStruct.expt(j).DAT_MotifRenamed.Motifs_IsTarg(k);
    issame = GenStruct.expt(j).DAT_MotifRenamed.Motifs_IsSame(k);
    
    if istarg==1
        % then is targ
        titlename = 'targ';
        plotcol = 'k';
    elseif istarg==0 & issame==1
        titlename = 'same';
        plotcol = 'b';
    elseif istarg==0 & issame==0
        titlename = 'diff';
        plotcol = 'r';
    else
        dfafasdfasfasdf;
    end
    
    % ================= DATA
    ffvals = GenStruct.expt(j).DAT_MotifRenamed.SongByMotifMat(:, k);
    tsong = datenum(GenStruct.expt(j).DAT.SongByMotif_songnames, 'yymmdd_HHMMSS');
        
    % =============== PLOT
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots, hsplot];
    title([titlename '-' motifthis]);
    plot(tsong, ffvals, 'o', 'Color', plotcol);
    
    % ----- overlay onset and offset of learning window
    line([tOn_learn tOn_learn], ylim, 'Color', 'b', 'LineStyle', '--');
    line([tOff_learn tOff_learn], ylim, 'Color', 'b', 'LineStyle', '--');
end

linkaxes(hsplots, 'xy');


%% =============== QUICK XCOV
All_Xcov = [];
maxwind = 20;
for j=1:numexpts
    
    % ==== figure out which syls to compare
%     sylinds_targ = find(GenStruct.expt(j).Params.IsTarg==1);
%     sylinds_same = find(GenStruct.expt(j).Params.IsSame==1 & ...
%         GenStruct.expt(j).Params.IsTarg==0);
    sylinds_targ = find(GenStruct.expt(j).DAT_MotifRenamed.Motifs_IsTarg==1);
    sylinds_same = find(GenStruct.expt(j).DAT_MotifRenamed.Motifs_IsSame==1 & ...
        GenStruct.expt(j).DAT_MotifRenamed.Motifs_IsTarg==0);
    
    % ======== figure out which song inds are within bounds of leraning window
    tsongs = GenStruct.expt(j).DAT.SongByMotif_songnames;
    tsongs = datenum(tsongs, 'yymmdd_HHMMSS');
    
    tOn = GenStruct.expt(j).Params.Datenum_LearnWind_On;
    tOff = GenStruct.expt(j).Params.Datenum_LearnWind_Off;
    
    indsongs = tsongs>tOn & tsongs<tOff;
    
    
    % ================= GO THRU ALL PAIRS OF SYLS
    figcount=1;
    subplotrows=5;
    subplotcols=3;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    for i=1:length(sylinds_targ)
        ind1 = sylinds_targ(i);
%         motif1 = GenStruct.expt(j).Params.MotifList{ind1};
        motif1 = GenStruct.expt(j).DAT_MotifRenamed.SongByMotif_motifnames{ind1};
        for ii=1:length(sylinds_same)
            ind2 = sylinds_same(ii);
            motif2 = GenStruct.expt(j).DAT_MotifRenamed.SongByMotif_motifnames{ind2};
            
            % ====================== EXTRACT ALL DATA
%             y1 = GenStruct.expt(j).DAT.SongByMotifMat(:, ind1);
%             y2 = GenStruct.expt(j).DAT.SongByMotifMat(:, ind2);
            y1 = GenStruct.expt(j).DAT_MotifRenamed.SongByMotifMat(:, ind1);
            y2 = GenStruct.expt(j).DAT_MotifRenamed.SongByMotifMat(:, ind2);
            
            % --- only keep songs within learning range
            y1 = y1(indsongs);
            y2 = y2(indsongs);
            
            % --- remove nans
            indsnan = unique([find(isnan(y1)); find(isnan(y2))]);
            y1(indsnan) = [];
            y2(indsnan) = [];
            
            [cc, lags] = xcov(y1, y2, maxwind, 'unbiased');
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title(['[targ]' motif1 '-[same]' motif2]);
            plot(lags, cc, '-k');
            
            lt_plot_zeroline;
            lt_plot_zeroline_vert;
        end
    end
end


