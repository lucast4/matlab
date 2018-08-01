function GenStruct = lt_generalization_PreProc
%% lt - extracts genearlization structure. combines syl names if required.

%% ################## COLLECT ALL LEARNING DATA
GenStruct = struct;
ind = 0;

% ############################# SEQ DEP PITCH EXPERIMENTS



% ####################################### NEURAL RECORDING EXPERIMENTS
% ====== bk7 - LearnLMAN1 
% laerning about 30-50hz over a day to next morning. Not great learning but 
% could include.

% ====== bk7 - LearnLMAN2
% NO LEARNING - about 10hz learning, just one day...



% ====== wh6pk36 - LMANlearn2
% NOTE: ALL LABELED
birdname = 'wh6pk36';
exptname = 'LMANlearn2';
[DATSTRUCT, Params] = lt_generalization_database(birdname, exptname);
% --- add to overall data structure
ind = ind+1;
GenStruct.expt(ind).birdname = birdname;
GenStruct.expt(ind).exptname = exptname;
GenStruct.expt(ind).DAT = DATSTRUCT;
GenStruct.expt(ind).Params = Params;



% ====== bu77wh13 - LMANlearn1
% NOTE: ALL LABELED (up to first switch)
birdname = 'bu77wh13';
exptname = 'LMANlearn1';
[DATSTRUCT, Params] = lt_generalization_database(birdname, exptname);
% --- add to overall data structure
ind = ind+1;
GenStruct.expt(ind).birdname = birdname;
GenStruct.expt(ind).exptname = exptname;
GenStruct.expt(ind).DAT = DATSTRUCT;
GenStruct.expt(ind).Params = Params;


% ======= br92br54 - LMANlearn2
% WEAK/NO learning (don't use) [like 0-20 hz over entire day]

% ======= br92br54 - LMANlearn3
% NO LEARNING

% ======= br92br54 - LMANlearn4
% NO LEARNING (at least on first day) - some weak learning later. Also
% problem is starts from UP, so no clear baseline.
% POTENTIALLY can look at prelearning day. but I am not sure, doesn't look
% like substantial learning there either.

% ======= br92br54 - LMANlearn5
% NO LEARNING - also 2 different targets ...

% ======= br92br54 - LMANlearn6
% NO LEARNING - also dirving 2 contexts in opposite directions

% ======= br92br54 - LMANlearn7
% HAVE NOT CHECKED, but likely NO LEARNING 
% also, washout before learning was only today (i.e. ended previous
% experiment same day that learning began).


% ======= or74bk35 - LMANnearal2
% NOTE: ALL LABELED.
birdname = 'or74bk35';
exptname = 'LMANneural2';
[DATSTRUCT, Params] = lt_generalization_database(birdname, exptname);
% --- add to overall data structure
ind = ind+1;
GenStruct.expt(ind).birdname = birdname;
GenStruct.expt(ind).exptname = exptname;
GenStruct.expt(ind).DAT = DATSTRUCT;
GenStruct.expt(ind).Params = Params;


% ======= or74bk35 - LMANlearn3
% NO, dirving 2 contexts in opposite directions, no clear learning in one
% day...


% ======= pu59wh78 - RAlearn1
% NOTE: could add (there is learning for target of (b), but not useful
% becuase there is no nontarget (i.e. even if use (h), all h occur in same
% motif as b, so confounded by acute WN effects
if (0)
birdname = 'pu69wh78';
exptname = 'RAlearn1';
[DATSTRUCT, Params] = lt_generalization_database(birdname, exptname);
% --- add to overall data structure
ind = ind+1;
GenStruct.expt(ind).birdname = birdname;
GenStruct.expt(ind).exptname = exptname;
GenStruct.expt(ind).DAT = DATSTRUCT;
GenStruct.expt(ind).Params = Params;
end


% ======= pu69wh78 - RALMANlearn1
% ALL LABELED!!
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
%  NOTE: weird rapid jump in pitch, doesn't seem to be like normal learning
%  ..


% ======= pu69wh78 - RALMANOvernightLearn1
%  NOTE: added target in middle of first day ... only 3.5 hours for
%  learning for first target, so do not use.


% ======= wh44wh39 - RALMANlearn1
% NOTE: added target middle of first day, and no clear learning before
% then.


% ======= wh44wh39 - RALMANlearn2
% LABELED ALL!
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
% ALL LABELED!
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
% ----- DEFAULT, first siwtch time. if ~exist, then ceiling of last day.
% ----- OLD, from WN onset to min([end of first day, switch time]);
numexpts = length(GenStruct.expt);

for j=1:numexpts
    
    tOn = datenum(GenStruct.expt(j).Params.Date_WNon, 'ddmmmyyyy-HHMM');
    if ~isempty(GenStruct.expt(j).Params.Date_SwitchTimes)
        tSwitch = min(datenum(GenStruct.expt(j).Params.Date_SwitchTimes, 'ddmmmyyyy-HHMM'));
    else
        tSwitch = [];
    end
    
    % ================= what is time of last rendition across syls?
    nsyls = length(GenStruct.expt(j).DAT.motif);
    tmaxall =[];
    for jj=1:nsyls
       tmax = max([GenStruct.expt(j).DAT.motif(jj).rendnum.datenum_song_SecRes]);
       tmaxall = [tmax tmaxall];
    end
    tmax = max(tmaxall);
    
    % ================= if switch exists, then use that. otherwise use end
    % of last day with data.
    if (1)
        if ~isempty(tSwitch)
            tOff = tSwitch;
        else
            tOff = ceil(tmax);
        end
    else
        tOff = min([ceil(tOn) tSwitch]); % the minimum of end of day or switch time.
    end
    
    % === sanity check, where to put end of learning
    % plots WNonset, time of switch (if exist), both dashed red
    % plot tOffset, which is decided above.
    if (0)
    disp([num2str(ceil(tOn)) '--' num2str(tSwitch)]);
        lt_figure; hold on;
        targind = strcmp({GenStruct.expt(j).DAT.motif.motif}, ...
            GenStruct.expt(j).Params.TargSyl);
        if ~any(targind)
            targind = 1;
        end
        t = [GenStruct.expt(j).DAT.motif(targind).rendnum.datenum_song_SecRes];
        ff = [GenStruct.expt(j).DAT.motif(targind).rendnum.ff];
        plot(t, ff, 'ok');
        
        line([tOff tOff], ylim, 'Color','k');
        line([tOn tOn], ylim, 'Color', 'r', 'LineStyle', '--');
        try
            line([tSwitch tSwitch], ylim, 'Color', 'r','LineStyle', '--');
        catch err
        end
    end
    
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
    
    % ======== adding individual renditions
    GenStruct.expt(j).DAT_MotifRenamed.motif= ...
        GenStruct.expt(j).DAT.motif;
    
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
        ff_new = nanmean(ff_old,2);
        
        
        % =========== keep track of which motifs I will remove
        Inds_MotToRemove = [Inds_MotToRemove indtmp'];
        
        % =========== add this onto the matrix, and update the list of syls
        GenStruct.expt(j).DAT_MotifRenamed.SongByMotifMat = ...
            [GenStruct.expt(j).DAT_MotifRenamed.SongByMotifMat ff_new];
        GenStruct.expt(j).DAT_MotifRenamed.SongByMotif_motifnames = ...
            [GenStruct.expt(j).DAT_MotifRenamed.SongByMotif_motifnames newsyl];
        
        
        % ##################### do same for individual rends (just combine
        % all trials)
        % ---------- 1) collect all rends across these 2 syls
        allrends = [];
        for k=1:length(indtmp)
            oldsylthis = indtmp(k);
            assert(strcmp(GenStruct.expt(j).DAT_MotifRenamed.motif(oldsylthis).motif, oldsyls{k}), 'asdf')
            
            % ---------- subtract baseline mean for this syllable
            tvals = [GenStruct.expt(j).DAT_MotifRenamed.motif(oldsylthis).rendnum.datenum_song_SecRes];
            indsbasetmp = tvals<GenStruct.expt(j).Params.Datenum_LearnWind_On;
            
            ffvals = [GenStruct.expt(j).DAT_MotifRenamed.motif(oldsylthis).rendnum.ff];
            ffbase = nanmean(ffvals(indsbasetmp));
            ffvals = ffvals-ffbase;
            
            % ----------- put back into structure
            nrends = length(GenStruct.expt(j).DAT_MotifRenamed.motif(oldsylthis).rendnum);
            for kk=1:nrends
                % --- save actual ff (not minus base...)
                GenStruct.expt(j).DAT_MotifRenamed.motif(oldsylthis).rendnum(kk).ff_notminusbase = ...
                    GenStruct.expt(j).DAT_MotifRenamed.motif(oldsylthis).rendnum(kk).ff;
                
                GenStruct.expt(j).DAT_MotifRenamed.motif(oldsylthis).rendnum(kk).ff = ...
                    ffvals(kk);
            end
            
            allrends = [allrends GenStruct.expt(j).DAT_MotifRenamed.motif(oldsylthis).rendnum];
        end
        
        % ----------- 2) sort all rends by time of rendition
        t_songs = [allrends.datenum_song_SecRes];
        t_withinsongs = [allrends.time_withinsong];
        t_withinsongs = t_withinsongs./(60*60*24); % convert from seconds to units of days
        t_actual = t_songs + t_withinsongs;
        
        % ----------- 3) sort trials
        [~, indssort] = sort(t_actual);
        allrends = allrends(indssort);
        
        % ----------- 4) return to struct
        motiftmp.motif = newsyl;
        motiftmp.rendnum = allrends;
        GenStruct.expt(j).DAT_MotifRenamed.motif = [GenStruct.expt(j).DAT_MotifRenamed.motif ...
            motiftmp];
        
    end
    
    % =========== 1) Remove the renamed syls
    GenStruct.expt(j).DAT_MotifRenamed.SongByMotifMat(:,Inds_MotToRemove) = [];
    GenStruct.expt(j).DAT_MotifRenamed.SongByMotif_motifnames(Inds_MotToRemove) = [];
    GenStruct.expt(j).DAT_MotifRenamed.SongByMotif_baseMeans(Inds_MotToRemove) = [];
    
    GenStruct.expt(j).DAT_MotifRenamed.motif(Inds_MotToRemove) = [];
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

