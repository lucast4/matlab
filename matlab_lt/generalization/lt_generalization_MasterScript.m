%% lt 4/4/18 - looking at timecourse of generalization
%% ################## COLLECT ALL LEARNING DATA
GenStruct = struct;
ind = 0;

% ================================= SEQ DEP PITCH EXPERIMENTS



% ================================= NEURAL RECORDING EXPERIMENTS


% ======= pu69wh78 - RALMANlearn1
birdname = 'pu69wh78';
exptname = 'RALMANlearn1';
[DATSTRUCT, Params] = lt_generalization_database(birdname, exptname);
% --- add to overall data structure
ind = ind+1;
GenStruct.expt(ind).birdname = birdname;
GenStruct.expt(ind).exptname = exptname;
GenStruct.expt(ind).DAT = DATSTRUCT;
GenStruct.expt(ind).Params = Params;


% ======= pu69wh78 - RALMANlearn2
%  NOTE: many switches, but maybe use up to first switch if enough data

% ======= pu69wh78 - RALMANOvernightLearn1
%  NOTE: added target in middle of first day ...


% ======= wh44wh39 - RALMANlearn1
% NOTE: added target middle of first day, and no clear learning before
% then.

% ======= wh44wh39 - RALMANlearn2
birdname = 'wh44wh39';
exptname = 'RALMANlearn2';
[DATSTRUCT, Params] = lt_generalization_database(birdname, exptname);
% --- add to overall data structure
ind = ind+1;
GenStruct.expt(ind).birdname = birdname;
GenStruct.expt(ind).exptname = exptname;
GenStruct.expt(ind).DAT = DATSTRUCT;
GenStruct.expt(ind).Params = Params;


% =======
birdname = 'wh44wh39';
exptname = 'RALMANlearn3';
[DATSTRUCT, Params] = lt_generalization_database(birdname, exptname);
% --- add to overall data structure
ind = ind+1;
GenStruct.expt(ind).birdname = birdname;
GenStruct.expt(ind).exptname = exptname;
GenStruct.expt(ind).DAT = DATSTRUCT;
GenStruct.expt(ind).Params = Params;


% ====== wh44wh39 - RALMANlearn4 
% NOTE: no learning ...


%% ############### FIGURE OUT WHAT IS TRAINING WINDOW
% ----- DEFAULT, from WN onset to min([end of first day, switch time]);

for j=1:numexpts
    
    tOn = datenum(GenStruct.expt(j).Params.Date_WNon, 'ddmmmyyyy-HHMM');
    if ~isempty(GenStruct.expt(j).Params.Date_SwitchTimes)
        tSwitch = min(datenum(GenStruct.expt(j).Params.Date_SwitchTimes, 'ddmmmyyyy-HHMM'));
    else
        tSwitch = [];
    end
    
    tOff = min([ceil(tOn) tSwitch]); % the minimum of end of day or switch time.
    
    % ================ save
    GenStruct.expt(j).Params.Datenum_LearnWind_On = tOn;
    GenStruct.expt(j).Params.Datenum_LearnWind_Off = tOff;
    
end


%% ############### GET SONG BY SONG DATA by averaging

figcount=1;
subplotrows=3;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


numexpts = length(GenStruct.expt);
for j=1:numexpts
    
    nummotifs = length(GenStruct.expt(j).DAT.motif);
    
    % =========== 1) get list of all songs that exist across all motifs
    SongDnumAll = {};
    for m=1:nummotifs
        songdnum = {GenStruct.expt(j).DAT.motif(m).rendnum.datestr};
        SongDnumAll = [SongDnumAll songdnum];
    end
    SongDnumUnique = unique(SongDnumAll);
    numsongs = length(SongDnumUnique);
    
    % ---------- what are baseline song inds?
            tsongs_dnum = datenum(SongDnumUnique, 'yymmdd_HHMMSS');
        tOn = GenStruct.expt(j).Params.Datenum_LearnWind_On;
        
        indsbase = tsongs_dnum<tOn;

    % =========== 2) make matrix of songbout x motif (mean of motif)
    SongByMotifMat = nan(numsongs, nummotifs);
    BaseMeans = nan(1, nummotifs);
    for m = 1:nummotifs
        
        tsongs = {GenStruct.expt(j).DAT.motif(m).rendnum.datestr};
        ffvals = [GenStruct.expt(j).DAT.motif(m).rendnum.ff];
        
        
        songcount = 0;
        for s = 1:numsongs
            songthis = SongDnumUnique{s};
            indtmp = strcmp(tsongs, songthis);
            
            y = mean(ffvals(indtmp));
            
            SongByMotifMat(s, m) = y;
            
            songcount = songcount+sum(indtmp);
        end
        assert(songcount == length(ffvals), 'did not account for all data ...');
        
        % ============ normalize by baseline songs
        %  --- 2) get mean FF across baseline
        base_mean = nanmean(SongByMotifMat(indsbase, m));
%         base_std = std(SongByMotifMat(indsbase, m));
      
        % ------- DO NORMLAIZAZTION.
        SongByMotifMat(:, m) = SongByMotifMat(:, m) - base_mean;
        BaseMeans(m) = base_mean;
    end
    
        % ===== plot to get sense of structure of nans
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([GenStruct.expt(j).birdname '-' GenStruct.expt(j).exptname]);
    spy(isnan(SongByMotifMat));
    xlabel('motif');
    ylabel('songbout (dot = isnan)');
    
    % ====== put back into structure
    GenStruct.expt(j).DAT.SongByMotifMat = SongByMotifMat;
    GenStruct.expt(j).DAT.SongByMotif_songnames = SongDnumUnique;
    GenStruct.expt(j).DAT.SongByMotif_motifnames = GenStruct.expt(j).Params.MotifList;
    GenStruct.expt(j).DAT.SongByMotif_baseInds = indsbase;
    GenStruct.expt(j).DAT.SongByMotif_baseMeans = BaseMeans;
    
end



%% ####### RENAME ANY MOTIFS (E.G. ACROSS CONTEXT, BUT BOTH TARGETED)
for j=1:numexpts
   % =========
   motiflist_old = GenStruct.expt(j).Params.MotifList;
   numcombine = length(GenStruct.expt(j).Params.MotifsToRename)/2;
   
   % ========== make a copy of song by motif matrix. will keep the old one
   % as is and update this one (i.e. remove combined, add new ones)
   GenStruct.expt(j).DAT_MotifRenamed.SongByMotifMat = ...
       GenStruct.expt(j).DAT.SongByMotifMat;
   GenStruct.expt(j).DAT_MotifRenamed.SongByMotif_motifnames = ...
       GenStruct.expt(j).DAT.SongByMotif_motifnames;
   GenStruct.expt(j).DAT_MotifRenamed.SongByMotif_baseMeans = ...
       GenStruct.expt(j).DAT.SongByMotif_baseMeans;
   
   if isempty(GenStruct.expt(j).Params.MotifsToRename)
       continue
   end
    
   % ============ GO thru all the sets I want to combine
   Inds_MotToRemove = [];
   for cc =1:numcombine
      
       oldsyls = GenStruct.expt(j).Params.MotifsToRename{2*cc-1};
       newsyl = GenStruct.expt(j).Params.MotifsToRename{2*cc};
       
       % ========= Take average of all old syls
       [~, indtmp] = intersect(motiflist_old, oldsyls);
       ff_old = GenStruct.expt(j).DAT.SongByMotifMat(:, indtmp');
       % ----- subtract baseline mean
       indsbase = GenStruct.expt(j).DAT.SongByMotif_baseInds;
       ff_basemeans = nanmean(ff_old(indsbase, :),1);
       ff_old = ff_old - ff_basemeans;
       % ----- get average across syls;
       ff_new = mean(ff_old,2);
       
       
       % =========== keep track of which motifs I will remove
       Inds_MotToRemove = [Inds_MotToRemove indtmp'];
       
       % =========== add this onto the matrix, and update the list of syls
       GenStruct.expt(j).DAT_MotifRenamed.SongByMotifMat = ...
           [GenStruct.expt(j).DAT_MotifRenamed.SongByMotifMat ff_new];
       GenStruct.expt(j).DAT_MotifRenamed.SongByMotif_motifnames = ...
           [GenStruct.expt(j).DAT_MotifRenamed.SongByMotif_motifnames newsyl];
   end
   
   % =========== 1) Remove the renamed syls
   GenStruct.expt(j).DAT_MotifRenamed.SongByMotifMat(:,Inds_MotToRemove) = [];
   GenStruct.expt(j).DAT_MotifRenamed.SongByMotif_motifnames(Inds_MotToRemove) = [];
   GenStruct.expt(j).DAT_MotifRenamed.SongByMotif_baseMeans(Inds_MotToRemove) = [];
end


%% ================= DETERMINE TARG/NONTARG BASED ON ALL NEW MOTIFS

for j=1:numexpts
    
    motiflist = GenStruct.expt(j).DAT_MotifRenamed.SongByMotif_motifnames;
    targsyl = GenStruct.expt(j).Params.TargSyl;
    
    IsTarg = [];
    IsSame = [];
    for k=1:length(motiflist)
        motifthis = motiflist{k};
        
        % -------- target?
        if any(strcmp(motifthis, targsyl))
            IsTarg = [IsTarg 1];
        else
            IsTarg = [IsTarg 0];
        end
        
        % --- is this motif same type to any of the targets?
        issamevec = [];
        for kk=1:length(targsyl)
            tsylthis = targsyl{kk};
            [issame] = lt_neural_QUICK_SylPairType(motifthis, tsylthis);
            issamevec = [issamevec issame];
        end
        issame = any(issamevec); % if is same type as any target, can call it same
        IsSame = [IsSame issame];
        
    end
    
    % ========= output
    GenStruct.expt(j).DAT_MotifRenamed.Motifs_IsTarg = IsTarg;
    GenStruct.expt(j).DAT_MotifRenamed.Motifs_IsSame = IsSame;
end


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


