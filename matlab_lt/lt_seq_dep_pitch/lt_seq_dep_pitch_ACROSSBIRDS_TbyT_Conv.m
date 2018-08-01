function TrialStruct = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Conv(TrialStruct, ParamsTrial, ...
    GenStruct, addWithinSongTime, sortbytime)

%% lt 5/30/18 - Combines generalizationstruct (for neural learning) with seq dep pitch
% Note: see lt_generalization_MasterScript for GenStruct

%%

% addWithinSongTime = 1; % if 1, then resolution is to within song. otherwise 
% time is time of song. NOTE: this currently only works for experiments
% from Genstruct (i.e. neural expts.)

% sortbytime = 1; % then sorts all trials before including

throwOutRendsPostWNOff = 1; % 

%% ======= first confgirm that format for base dates already converted

assert(isfield(TrialStruct.birds(1).exptnum(1), 'WNontime'), 'need to convert first ...');

%% GO thru all Genstruct expts and place into trialstruct

numexpt_gen = length(GenStruct.expt);

for i=1:numexpt_gen
    
    bname = GenStruct.expt(i).birdname;
    ename = GenStruct.expt(i).exptname;
    motiflist = {GenStruct.expt(i).DAT_MotifRenamed.motif.motif};
    %     motiflist = {GenStruct.expt(i).DAT.motif.motif};
    Targsyls = GenStruct.expt(i).Params.TargSyl;
    Targlearndirs = GenStruct.expt(i).Params.TargLearnDirs;
    nummotifs = length(motiflist);
    WNon_datenum = GenStruct.expt(i).Params.Datenum_LearnWind_On;
    WNoff_datenum = GenStruct.expt(i).Params.Datenum_LearnWind_Off;
    
    % ====================================== only continue if there is one
    % target
    assert(length(Targsyls)==1, 'multipel targets learn same dir? then combine in gen preproess. if diff dir then dont keeop...');
    assert(length(Targlearndirs)==1, 'see above ...');
    
    
    
    % ######################## INSERT INTO TRIALSTRUCT
    indbird = find(strcmp({TrialStruct.birds.birdname}, bname));
    if isempty(indbird)
        % new bird and expt.
        indbird = length(TrialStruct.birds)+1;
        indexpt = 1;
    else
        % old bird, find where to slide expt.
        indbird = indbird;
        indexpt = length(TrialStruct.birds(indbird).exptnum)+1;
    end
    
    
    TrialStruct.birds(indbird).birdname = bname;
    TrialStruct.birds(indbird).exptnum(indexpt).exptname = ename;
    TrialStruct.birds(indbird).exptnum(indexpt).SylsUnique = motiflist;
    %     TrialStruct.birds(indbird).exptnum(indexpt).numEmptyDays = '';
    %     TrialStruct.birds(indbird).exptnum(indexpt).LMANinactivated = '';
    TrialStruct.birds(indbird).exptnum(indexpt).targlearndir = Targlearndirs;
    TrialStruct.birds(indbird).exptnum(indexpt).targsyl = Targsyls{1};
    
    
    
    % +++++++++++++++++++++++++++++++++++++++++ COLLECT EACH SYL DATA
    tvalsall = [];
    for mm=1:nummotifs
        
        motifname = GenStruct.expt(i).DAT_MotifRenamed.motif(mm).motif;
        ff = [GenStruct.expt(i).DAT_MotifRenamed.motif(mm).rendnum.ff];
        tvals_dnum = [GenStruct.expt(i).DAT_MotifRenamed.motif(mm).rendnum.datenum_song_SecRes];
        
        % ------- ADD ON WITHIN SONG TIMING
        if addWithinSongTime==1
        tvals_dnum_withinsong = [GenStruct.expt(i).DAT_MotifRenamed.motif(mm).rendnum.time_withinsong];
        % convert within song time to day
        tvals_dnum_withinsong = tvals_dnum_withinsong./(60*60*24);
        tvals_dnum = tvals_dnum + tvals_dnum_withinsong;
        end
        
        
        % ============== sort by time?
        if sortbytime==1
            [~, indsort] = sort(tvals_dnum);
            
            tvals_dnum = tvals_dnum(indsort);
            ff = ff(indsort);
            
        end
        
        
        % =============== only keep trials before end of learni
        if throwOutRendsPostWNOff==1
            indstoremove = tvals_dnum>=WNoff_datenum;
            
            tvals_dnum(indstoremove) = [];
            ff(indstoremove) = [];
        end
        
        tvalsall = [tvalsall, tvals_dnum];
        
        % ---- is this target? same type?
        istarg = GenStruct.expt(i).DAT_MotifRenamed.Motifs_IsTarg(mm);
        issame = GenStruct.expt(i).DAT_MotifRenamed.Motifs_IsSame(mm);
        
        
        % =========== PUT INTO OUTPUT STRUCT
        TrialStruct.birds(indbird).exptnum(indexpt).sylnum(mm).syl = motifname;
        TrialStruct.birds(indbird).exptnum(indexpt).sylnum(mm).Tvals_datenum = tvals_dnum';
        TrialStruct.birds(indbird).exptnum(indexpt).sylnum(mm).FFvals = ff';
        TrialStruct.birds(indbird).exptnum(indexpt).sylnum(mm).INFO_istarget = istarg;
        TrialStruct.birds(indbird).exptnum(indexpt).sylnum(mm).INFO_similar = issame;
        %         TrialStruct.birds(indbird).exptnum(indexpt).sylnum(mm).INFO_missingsomedat = '';
        %         TrialStruct.birds(indbird).exptnum(indexpt).sylnum(mm).INFO_SylDimensions = '';
    end
    
    
    % ############ THINGS THAT REQUIRE CONVERTING TIME TO DAYS REL EXPT START
    firstday = datestr(floor(min(tvalsall)), 'ddmmmyyyy');
    
    % ----------- go thru each motif and get updated timestamp
    for mm=1:nummotifs
        tvals_dnum = [TrialStruct.birds(indbird).exptnum(indexpt).sylnum(mm).Tvals_datenum];
        
        tmp = lt_convert_EventTimes_to_RelTimes(firstday, tvals_dnum); % convert to units of days
        TrialStruct.birds(indbird).exptnum(indexpt).sylnum(mm).Tvals = tmp.FinalValue;
    end
    
    % ----------------- WN ONSET
    tmp = lt_convert_EventTimes_to_RelTimes(firstday, WNon_datenum);
    TrialStruct.birds(indbird).exptnum(indexpt).WNontime = tmp.FinalValue;
    
    % ============== 
    disp(['added ' bname '-' ename '!']);
end

%% NOTE: do not include trials past WNoff period...
