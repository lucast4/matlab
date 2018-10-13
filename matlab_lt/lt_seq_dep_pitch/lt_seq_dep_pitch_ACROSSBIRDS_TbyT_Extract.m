function [DATBYREND, TrialStruct] = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Extract(TrialStruct, ParamsTrial, ...
    ignoreLMANexpt, songasrend, singleRendOnly)

%% 7/31/18 - lt diverged from lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse, here
% focusing on key extractions and analyses.

%% 7/31/18 - lt diverged from lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse, here
% focusing on key extractions and analyses.


%% ##################################################################
%% ############################# PREPROCESSING

Numbirds = length(TrialStruct.birds);

%% ======= CONVERT TO USING SONG AS REND?

if songasrend==1
    TrialStruct = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_SongRend(TrialStruct);
end

%% ====== CONVERT ALL FF TO FF DEVIATION FROM BASELINE (IN LEARN DIR)
% AND DECIDE IF USE SONG AS REND

if ~isfield(TrialStruct, 'FFalreadyFlippedLearnDir')
    TrialStruct.FFalreadyFlippedLearnDir = 0;
end

if TrialStruct.FFalreadyFlippedLearnDir==0
    disp('NOTE TO SELF: ffvals will be automatically flipped to be dir of learning (so no need to flip again later)');
    pause;
    
    
    for i=1:Numbirds
        Numexpt = length(TrialStruct.birds(i).exptnum);
        
        for ii=1:Numexpt
            
            % ---------- SKIP IF NO DATA
            if isempty(TrialStruct.birds(i).exptnum(ii).sylnum)
                disp(['[no DATA] skipping ' TrialStruct.birds(i).exptnum(ii).exptname]);
                continue
            end
            
            % --------- ignore if lMAN?
            if ignoreLMANexpt==1
                if isfield(TrialStruct.birds(i).exptnum(ii), 'LMANinactivated')
                    % if no field then is not LMAN inact experiemts...
                    isLMAN = TrialStruct.birds(i).exptnum(ii).LMANinactivated;
                    if isLMAN==1
                        disp(['[is LMAN] skipping ' TrialStruct.birds(i).exptnum(ii).exptname]);
                        continue
                    end
                end
            end
            
            Numsyls = length(TrialStruct.birds(i).exptnum(ii).sylnum);
            %         birdname =  TrialStruct.birds(i).birdname;
            %         exptname = TrialStruct.birds(i).exptnum(ii).exptname;
            targlearndir = TrialStruct.birds(i).exptnum(ii).targlearndir;
            
            % =========== collect syls for this experiment
            %         ffedges_allsyls =[];
            %         tedges_allsyls = [];
            %         istarg_allsyls = [];
            %         issame_allsyls =[];
            %         sylnames_allsyls = {};
            
            % ===================== FOR ALL SYLS, each rendition collect time and
            % deviation from running avg.
            for ss =1:Numsyls
                
                % ============== subplot for this syl
                t = TrialStruct.birds(i).exptnum(ii).sylnum(ss).Tvals;
                ff = TrialStruct.birds(i).exptnum(ii).sylnum(ss).FFvals;
                
                % --------------- lines for base mean and 1std
                WNontime = TrialStruct.birds(i).exptnum(ii).WNontime;
                
                %             basedays = TrialStruct.birds(i).exptnum(ii).BaseDays;
                %             wndays = TrialStruct.birds(i).exptnum(ii).WNDays;
                indsbase = t<WNontime;
                %             indsbase = t<basedays(end)+1;
                ffmean_base = nanmean(ff(indsbase));
                
                assert(abs(targlearndir)==1);
                
                % --------------- subtract mean and flip if negative larning
                ff = (ff-ffmean_base).*targlearndir;
                TrialStruct.birds(i).exptnum(ii).sylnum(ss).FFvals = ff;
                
                % ------------- for all renditions calculate deviation from
                % recent trials
                % === method1 - fit regression line to one hour of data
                % (directly preceding this rendition...) record deviation from
                % that hour's prediction
                
            end
        end
    end
    TrialStruct.FFalreadyFlippedLearnDir = 1;
    
end

%% ========== FOR ALL RENDS, COMPUTE FLANKING RENDITIONS (FF AND T DEV)

% =============== BINS
% xedges = [-60 -31:5:-1 1:5:31 60]; % minutes
% xedges = [-30:10:-10 -0.1 0.1 10:10:30];
xedges = [-30 -0.5 0.5 30]; % GOOD, for comparison to baseline
% xedges = [-30:2.5:-2.5 -0.1 0.1 2.5:2.5:30];
% xedges = [-60 -30:5:30 60];
% xedges = [-60 -30:5:-5 -0.1 0.1 5:5:30 60];
% xedges = [-40 -20:5:20 40];
% xedges = [-60 -15 -3 0 3 15 60];
% xedges = [-1440 0 1440];
% % NOTE: the edge of xedges will be used as flank. so will only keep rends
% with at least that much time flanking.

xcenters = xedges(1:end-1)+diff(xedges)/2;

flipTargNontarg=0; % if 1, then for each targ/nontarg pair, flips the dataset
% so that is target dev relative to target pitch, as function of nontarget
% density [default =0];

useMeanFFRef = 0; % 1: uses mean in window; 0; uses value for each rend
flanktime_targ = 2.5; % minutes (will do +/- this number). This only used if
% note: this is used to decide how many rends of targ/nontarg are present
% during reference period.

tflankplot = min(abs([xedges(1) xedges(end)]))./60;
tflankplot = 0.5;
% tflankplot = 0;

twind_plot = [xedges(1) xedges(end)]./60;

TrialStruct = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S5(TrialStruct, ...
    ignoreLMANexpt, tflankplot, twind_plot, flanktime_targ, flipTargNontarg, ...
    useMeanFFRef);




%% ====== TO COLLECT DATA TO PLOT (NEW VERSION, COLLECTS FIRST, THEN PLOTS)
% Instead of collecitng and plotting at same time

% ======== what data to take?
% dotrain = 1; % during training?
% singleRendOnly=1; % 1: only takes one datapt (first nontarg after the referen)
% if 0, then takes all rends in window


% ====== what data to use to compute density?
densitymethod = 'refperiod';
% refperiod: will use reference period
% entirebin: will use time from rendtion to maxtime_fromref (see maxtime_fromref below)
% beforefirstnontarg = will get number that is >= ref time(i.e.0) and < time of first nontarget post-reference.


% ============ method for decideing hi and lo density trials
cutoffmethod = 'medianslice';
% medianslice: for each nontarg bin, finds median for targ
% medianoverall: overall median of targ
% medianslice_rand: gets median by first addaing random jitter. useful for
% i) when even num bins. if only one y bin, then forces that to be in low
% density bin. if odd num bins, then will split favoring putting
% more data into low density category.
% targffdev: then splits based on mean FF dev of target syllables


% ======== summary plot: what time period to plot (locked to ref period)
% NOTE: This is also used to define period for estimating density if the
% choice above for densitymethod is entirebin.
% NOTE ONLY INPORTANT FOR DENSITY MEASURE.
mintime_fromref = 5;
maxtime_fromref = 30; % minutes


% =============== DEFAULTS
% singleRendOnly=0; % only take data up to first nontarg after the referen
% densitymethod = 'refperiod';
% cutoffmethod = 'medianslice';
% mintime_fromref = 5;
% maxtime_fromref = 30; % minutes
% minrends_inbin = 20;

% ================== include target syl?
includeTarg=0;

% ================ RUN
DATBYREND = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S2(TrialStruct, ...
    singleRendOnly, densitymethod, cutoffmethod, mintime_fromref, ...
    maxtime_fromref, ignoreLMANexpt, includeTarg);


assert(all(isnan(DATBYREND.IsDurTrain) == isnan(DATBYREND.Density_isHigh)), 'then nan likely matches up');
assert(all(isnan(DATBYREND.Density_targ) == isnan(DATBYREND.Density_nontarg)), 'then nan likely matches up');



%% ################################################
%% ######################## PLOTS

%% TO DO:
if (0)
    %% =========== [PLOT DIAGNOSTIC] FOR INDIVIDUAL EXPERIMENTS
    
    if (1)
        birdtoplot = 13;
        expttotplot = 1;
        
        lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S7(DATBYREND, TrialStruct, ...
            birdtoplot, expttotplot, mintime_fromref, maxtime_fromref, xedges);
    end
    
    % ========== PLOT MULTIPLE
    birdtoplot_list = [13:17];
    % birdtoplot_list = 1:3:17;
    expttoplot_list = [1:10];
    for bb=birdtoplot_list
        for ee=expttoplot_list
            lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S7(DATBYREND, TrialStruct, ...
                bb, ee, mintime_fromref, maxtime_fromref, xedges);
        end
    end
end


%% ###########################################################
%% ###########################################################
%% [IMPORTANT]
%% COLLECT ALL DATA. ANALYSIS WITHOUT SEPARATING INTO HI AND LO DENS
% [ALTHOUGH BECAUSE GOT ALL DATA, CAN EASILY PARSE BY THAT]
%% ==================== [COLLECT - 1] collect deviations

% =============== BINS
% xedges = [-65 -45 -25 -5 5 25 45 65]; % minutes
xedges = [-80 -55 -30 -5 5 30 55 80]; % minutes
xedges = [-60 -31:5:-1 1:5:31 60]; % minutes
% xedges = [-50:15:-5 -0.5 0.5 5:15:50]; % minutes
xedges = [-60 -30 -7 -0.5 0.5 7 30 60]; % minutes
xcenters = xedges(1:end-1)+diff(xedges)/2;
% xcenters = xedges(1:end-1)+binsize/2;

flipTargNontarg=0; % if 1, then for each targ/nontarg pair, flips the dataset
% so that is target dev relative to target pitch, as function of nontarget
% density [default =0];

useMeanFFRef = 0; % 1: uses mean in window; 0; uses value for each rend

collectTarg =1;

TrialStruct = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S5(TrialStruct, ...
    ignoreLMANexpt, tflankplot, twind_plot, flanktime_targ, flipTargNontarg, ...
    useMeanFFRef, collectTarg, songasrend);


%% ====== TO COLLECT DATA TO PLOT (NEW VERSION, COLLECTS FIRST, THEN PLOTS)
% Instead of collecitng and plotting at same time

% ======== what data to take?
% dotrain = 1; % during training?
% singleRendOnly=1;
densitymethod = 'refperiod';
cutoffmethod = 'medianslice';
mintime_fromref = 5;
maxtime_fromref = 30; % minutes
includeTarg=1; % IMPORTANT: will need to parse by syl type in all foture analyses

% ================ RUN

DATBYREND = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S2(TrialStruct, ...
    singleRendOnly, densitymethod, cutoffmethod, mintime_fromref, ...
    maxtime_fromref, ignoreLMANexpt, includeTarg);


assert(all(isnan(DATBYREND.IsDurTrain) == isnan(DATBYREND.Density_isHigh)), 'then nan likely matches up');
assert(all(isnan(DATBYREND.Density_targ) == isnan(DATBYREND.Density_nontarg)), 'then nan likely matches up');


