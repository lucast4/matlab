function DATBYREND = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S2(TrialStruct, ...
    singleRendOnly, densitymethod, cutoffmethod, mintime_fromref, ...
    maxtime_fromref, ignoreLMANexpt, includeTarg)

if ~exist('includeTarg', 'var')
    includeTarg = [];
end

if isempty(includeTarg)
    includeTarg = 0;
end

%% NOTE: checked on 7/31 everything good, except did not look closely at hihgh/lo dens stuff (previously did look at)..

% % ======== what data to take?
% % dotrain = 1; % during training?
% singleRendOnly=1; % 1: only takes one datapt (first nontarg after the referen)
% % if 0, then takes all rends in window
%
%
% % ====== what data to use to compute density?
% densitymethod = 'refperiod';
% % refperiod: will use reference period
% % entirebin: will use time from rendtion to maxtime_fromref (see maxtime_fromref below)
% % beforefirstnontarg = will get number that is >= ref time(i.e.0) and < time of first nontarget post-reference.
%
%
% % ============ method for decideing hi and lo density trials
% cutoffmethod = 'medianslice';
% % medianslice: for each nontarg bin, finds median for targ
% % medianoverall: overall median of targ
% % medianslice_rand: gets median by first addaing random jitter. useful for
% % i) when even num bins. if only one y bin, then forces that to be in low
% % density bin. if odd num bins, then will split favoring putting
% % more data into low density category.
% % targffdev: then splits based on mean FF dev of target syllables
%
%
% % ======== summary plot: what time period to plot (locked to ref period)
% % NOTE: This is also used to define period for estimating density if the
% % choice above for densitymethod is entirebin.
% % NOTE ONLY INPORTANT FOR DENSITY MEASURE.
% mintime_fromref = 5;
% maxtime_fromref = 30; % minutes
%
%
% % =============== DEFAULTS
% % singleRendOnly=0; % only take data up to first nontarg after the referen
% % densitymethod = 'refperiod';
% % cutoffmethod = 'medianslice';
% % mintime_fromref = 5;
% % maxtime_fromref = 30; % minutes
% % minrends_inbin = 20;
%
% % ================== include target syl?
% includeTarg=0;


%% lt 6/13/18 - pulls out all data by rend (for all nontarget syls)

% if singleRendOnly==1
%     % needs to be this to make sense ...
%     densitymethod = 'beforefirstnontarg';
% end
% changed because this does not need to bet he case.


%%

Numbirds = length(TrialStruct.birds);


% ========== NEW, COLLECT EACH RENDITION
DATBYREND.Density_targ = [];
DATBYREND.Density_nontarg = [];
DATBYREND.Density_isHigh = [];
DATBYREND.FF_dev = {};
DATBYREND.Time_dev = {};
DATBYREND.Time_dev_targ = {};
DATBYREND.IsDurTrain = [];
DATBYREND.SigLearn  = [];

DATBYREND.Birdnum= [];
DATBYREND.Exptnum= [];
DATBYREND.Sylnum= [];
DATBYREND.Sylcounter = [];
DATBYREND.IsSame = [];
DATBYREND.IsTarg = [];
DATBYREND.LearnMag_regr = [];
DATBYREND.Isfrom_SDP = [];

DATBYREND.IsWN = [];
DATBYREND.IsCatch = [];

DATBYREND.LearnLocal = [];
DATBYREND.LearnLocal_targ = [];

          DATBYREND.WN_hits = [];
          DATBYREND.WN_miss = [];

          
sylcount = 1;

for i=1:Numbirds
    Numexpt = length(TrialStruct.birds(i).exptnum);
    %     birdname = TrialStruct.birds(i).birdname;
    
    for ii=1:Numexpt
        %         exptname = TrialStruct.birds(i).exptnum(ii).exptname;
        
        % ---------- SKIP IF NO DATA
        if isempty(TrialStruct.birds(i).exptnum(ii).sylnum)
            disp(['[no DATA] skipping ' TrialStruct.birds(i).exptnum(ii).exptname]);
            continue
        end
        
        % --------- ignore if lMAN?
        if ignoreLMANexpt==1
            if isfield(TrialStruct.birds(i).exptnum(ii), 'LMANinactivated')
                isLMAN = TrialStruct.birds(i).exptnum(ii).LMANinactivated;
                if isLMAN==1
                    disp(['[is LMAN] skipping ' TrialStruct.birds(i).exptnum(ii).exptname]);
                    continue
                end
            end
        end
        
        % ===================== INFOR ABOUT TARGET
        %         indtarg = find([TrialStruct.birds(i).exptnum(ii).sylnum.INFO_istarget]==1);
        %         tvals_targ = TrialStruct.birds(i).exptnum(ii).sylnum(indtarg).Tvals;
        %         ffvals_targ = TrialStruct.birds(i).exptnum(ii).sylnum(indtarg).FFvals;
        
        % ############################################ SAME'
        if includeTarg==0
            sylinds_this = find([TrialStruct.birds(i).exptnum(ii).sylnum.INFO_istarget]==0);
        elseif includeTarg==1
            sylinds_this = 1:length(TrialStruct.birds(i).exptnum(ii).sylnum);
        end
        % ===================== GO THRU ALL NONTARG SYLS. FOR EACH REND COLLECT
        % DATA
        for indthis = sylinds_this
            %
            %             indthis = sylinds_this(j);
            %             sylthis = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).syl;
            samethis = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).INFO_similar;
            istargthis = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).INFO_istarget;
            
            
            
            % ################################# 1) COLLECT SINGING DENSITY
            if strcmp(densitymethod, 'refperiod')
                % ---------------------- 1) NTARGS VS. NNONTARGS IN REF WINDOW
                ntarg = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_NTargRendsInRef;
                nnontarg = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_NNonTargRendsInRef;
                
            elseif strcmp(densitymethod, 'entirebin')
                maxtime_day = maxtime_fromref/(60*24);
                mintime_day = -0.25/(60*24);
                %                 maxtime_day = (30-5)/(60*24);
                %                 mintime_day = -0.25/(60*24);
                %                 maxtime_day = 0.25/(60*24);
                
                % ---------------- COUNT NUMBER OF targ and nontarg between
                % ref syl and end of bin [NOTE: end is potentialy not exact
                % end of last bin...]
                timedevs_nontarg = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_TimeDev;
                timedevs_targ = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_TimeDevTargs;
                
                nrends = length(timedevs_nontarg);
                
                ntarg = nan(nrends,1);
                nnontarg = nan(nrends,1);
                for rr =1:nrends
                    if isempty(timedevs_nontarg{rr})
                        continue
                    end
                    %                     nnontarg(rr) = sum(timedevs_nontarg{rr}>=0 & timedevs_nontarg{rr}<maxtime_day);
                    %                     ntarg(rr) = sum(timedevs_targ{rr}>=0 & timedevs_targ{rr}<maxtime_day);
                    nnontarg(rr) = sum(timedevs_nontarg{rr}>=mintime_day & ...
                        timedevs_nontarg{rr}<maxtime_day);
                    ntarg(rr) = sum(timedevs_targ{rr}>=mintime_day & ...
                        timedevs_targ{rr}<maxtime_day);
                end
                
                % ===================== COMPARE IT TO REFPERIOD
                % ---------------------- 1) NTARGS VS. NNONTARGS IN REF WINDOW
                ntarg_REF = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_NTargRendsInRef;
                nnontarg_REF = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_NNonTargRendsInRef;
                
            elseif strcmp(densitymethod, 'beforefirstnontarg')
                
                % ------------------ COUNT RENDITIONS FROM REF (INCLUSIVE)
                % UP TO FIRST NONTARG POST REF
                timedevs_nontarg = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_TimeDev;
                timedevs_targ = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_TimeDevTargs;
                
                nrends = length(timedevs_nontarg);
                ntarg = nan(nrends,1);
                nnontarg = nan(nrends,1);
                for rr =1:nrends
                    if isempty(timedevs_nontarg{rr})
                        continue
                    end
                    
                    % ---- what is time of first nontarg post ref>?
                    timeoffirstNONTARG = timedevs_nontarg{rr}(find(timedevs_nontarg{rr}>0, 1, 'first'));
                    if isempty(timeoffirstNONTARG)
                        timeoffirstNONTARG=0;
                    end
                    
                    nnontarg(rr) = sum(timedevs_nontarg{rr}>=0 & timedevs_nontarg{rr}<timeoffirstNONTARG);
                    ntarg(rr) = sum(timedevs_targ{rr}>=0 & timedevs_targ{rr}<timeoffirstNONTARG);
                end
            end
            assert(all(isnan(ntarg)==isnan(nnontarg)));
            
            % =================== OUTPUT
            DATBYREND.Density_targ = [DATBYREND.Density_targ; ntarg];
            DATBYREND.Density_nontarg = [DATBYREND.Density_nontarg; nnontarg];
            
            
            
            
            
            % ########################### EACH REND HIGH OR LOW DENSITY?
            if strcmp(cutoffmethod, 'medianslice')
                %  for each value of nontarg, get median value of targ, and use
                %  that as threshold
                cutoffs_y_targ = grpstats(ntarg, nnontarg, {'median'});
                cutoffs_x_nontarg = unique(nnontarg);
                cutoffs_x_nontarg = unique(nnontarg(~isnan(nnontarg)));
                
                [~, indstmp] = ismember(nnontarg, cutoffs_x_nontarg);
                indstmp(indstmp==0) = 1; % temporarily make 1, since will convert to nan in a bit
                % ---- for each rend, determine whether it was high or low
                % density (in terms of targ syls)
                ishighdensity = ntarg>cutoffs_y_targ(indstmp);
                ishighdensity = single(ishighdensity);
                ishighdensity(isnan(ntarg)) = nan;
            elseif strcmp(cutoffmethod, 'medianslice_rand')
                ntarg_rand = ntarg - 0.2 + 0.4*rand(size(ntarg));
                cutoffs_y_targ = grpstats(ntarg_rand, nnontarg, {'median'});
                cutoffs_x_nontarg = unique(nnontarg);
                cutoffs_x_nontarg = unique(nnontarg(~isnan(nnontarg)));
                
                % ==== if any bin only has one y bin, then autamitcally
                % assign it to low density
                [ymin, ymax] = grpstats(ntarg, nnontarg, {'min', 'max'}); assert(length(ymin) == length(cutoffs_y_targ));
                cutoffs_y_targ((ymax-ymin)==0) = 1000;
                
                [~, indstmp] = ismember(nnontarg, cutoffs_x_nontarg);
                indstmp(indstmp==0) = 1; % temporarily make 1, since will convert to nan in a bit
                % ---- for each rend, determine whether it was high or low
                % density (in terms of targ syls)
                ishighdensity = ntarg>cutoffs_y_targ(indstmp);
                ishighdensity = single(ishighdensity);
                ishighdensity(isnan(ntarg)) = nan;
                
                % ---------------- if any x bin only has one y value, then
                % ignore if (since cannot split that bin)
                if (0)
                    [ymin, ymax] = grpstats(ntarg, nnontarg, {'min', 'max'});
                    xbintoremove = cutoffs_x_nontarg((ymax-ymin)==0);
                    indstoremove = ismember(nnontarg, xbintoremove);
                    if (0) % plots separately those with only one, or multiple y bins.
                        lt_figure; hold on;
                        plot(nnontarg(indstoremove), ntarg(indstoremove), 'ok');
                        plot(nnontarg(~indstoremove), ntarg(~indstoremove), 'ob');
                    end
                    % --- remove them
                    ishighdensity(indstoremove) = nan;
                end
                
            elseif strcmp(cutoffmethod, 'medianoverall')
                
                cutoffthis = nanmedian(ntarg);
                
                ishighdensity = ntarg>cutoffthis;
                ishighdensity = single(ishighdensity);
                ishighdensity(isnan(ntarg)) = nan;
                
            elseif strcmp(cutoffmethod, 'targffdev')
                % --- get global distribution of deviations for target
                TrialStruct.birds(i).exptnum(ii).sylnum(indthis);
                % CONTINUE HERE...
                
                % --- get median of that distribution
                
                
                
            end
            
            DATBYREND.Density_isHigh = [DATBYREND.Density_isHigh; ishighdensity];
            TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_IsHighDensTarg = ishighdensity;
            
            
            
            
            
            
            % ############################## NOTE DOWN AMOUNT OF LEARNING FOR THIS SYL
            %             learndir = TrialStruct.birds(i).exptnum(ii).targlearndir; % direction of learning at targ
            % ===== 2 ways to be called significnat learnign: i) slope, ii)
            % endpoint diff from baseline
            
            
            % ---------- 1) First measure, slope
            WNon = TrialStruct.birds(i).exptnum(ii).WNontime;
            ttmp = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).Tvals;
            indstrain = ttmp>WNon;
            
            if (1) % this version take last 1 renditions of baseline as well (appends to onset)
                ind1 = find(indstrain, 1, 'first');
                ind2 = find(indstrain, 1, 'last');
                if ind1>9
                    ind1 = ind1-10;
                else
                    ind1 = 1;
                end
                ttmp = ttmp(ind1:ind2);
                fftmp = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).FFvals(ind1:ind2);
            else % this version takes only trials during training
                ttmp = ttmp(indstrain);
                fftmp = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).FFvals(indstrain);
            end
            [b,bint] = lt_regress(fftmp, ttmp, 0); % slope of learning at nontarg
            
            
            % ------- 2nd measure, last N rends
            assert(all(diff(ttmp)>=0), 'expect to be sorted...');
            if (1) % version 1, does one test for last N renditions.
                NrendsOrig = 25;
                Nrends = NrendsOrig;
                %                 Nrends = max([NrendsOrig floor(length(ttmp)/2)]);
                if length(fftmp)<Nrends
                    p_end = 1;
                    ffendval = mean(fftmp);
                else
                    p_end = signrank(fftmp(end-Nrends+1:end));
                    ffendval = mean(fftmp(end-Nrends+1:end));
                end
            elseif (0) % version 1 - does 3 tests on 3 last 10-rend bins. this is
                % probably better as less prone to false positives.
                if length(fftmp)<3*10
                    p_end = 1;
                    ffendval = mean(fftmp);
                else
                    [~, p1] = ttest(fftmp(end-9:end));
                    [~, p2] = ttest(fftmp(end-19:end-10));
                    [~, p3] = ttest(fftmp(end-29:end-20));
                    p_end = max([p1 p2 p3]);
                    ffendval = mean(fftmp(end-19:end));
                end
            elseif (0) % take first and last half
                ntmp = round(length(fftmp)/2);
                %                 p_end = ranksum(fftmp(1:25), fftmp(end-24:end));
                p_end = ranksum(fftmp(1:ntmp), fftmp(ntmp+1:end));
                ffendval = mean(fftmp(ntmp+1:end));
            end
            
            isSigLearn = (sign(bint(2,1))==sign(bint(2,2)) | p_end<0.01) & ... % if CI of slope does not overlap 0;
                ffendval>0; % if slope is in direction of learning ...
            %             isSigLearn = (sign(bint(2,1))==sign(bint(2,2))); % if slope is in direction of learning ...
            
            DATBYREND.SigLearn  = [DATBYREND.SigLearn; isSigLearn*ones(size(ishighdensity))];
            
            % --- also save learning magnitude
            learnmag = (ttmp(end) - ttmp(1))*b(2);
            %             learnmag = fftmp(end) - fftmp(1);
            DATBYREND.LearnMag_regr  = [DATBYREND.LearnMag_regr; learnmag*ones(size(ishighdensity))];
            
            
            
            % ############################## COLLECT FF AND TIME DEVIATIONS
            timedev = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_TimeDev;
            ffdev = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_FFDev;
            timedev_TARG = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_TimeDevTargs;
            
            % ---------- ONLY KEEP FIRST RENDITION?
            % i.e. only keep 1 redition, the first one that has positive
            % timedev
            if singleRendOnly==1
                ntmp = length(timedev);
                
                % -------------- NONTARG
                for nn = 1:ntmp
                    indtmp = find(timedev{nn}>0, 1, 'first');
                    
                    if isempty(indtmp)
                        timedev{nn} = [];
                        ffdev{nn} = [];
                    else
                        timedev{nn} = timedev{nn}(indtmp);
                        ffdev{nn} = ffdev{nn}(indtmp);
                    end
                end
                
                %                 % ------------ TARGET
                %                 for nn = 1:ntmp
                %                     indtmp = find(timedev_TARG{nn}>0, 1, 'first');
                %
                %                     if isempty(indtmp)
                %                       timedev_TARG{nn} = [];
                %                     else
                %                       timedev_TARG{nn} = timedev_TARG{nn}(indtmp);
                %                     end
                %                 end
            elseif singleRendOnly==0
                % then only keep positive time deviations
                ntmp = length(timedev);
                
                % -------------- NONTARG
                for nn = 1:ntmp
                    
                    if isempty(timedev{nn})
                        continue
                    end
                    
                    indstmp = timedev{nn}>0;
                    
                    timedev{nn} = timedev{nn}(indstmp);
                    ffdev{nn} = ffdev{nn}(indstmp);
                    
                end
                
            end
            
%             if std(cell2mat(ffdev))>500
%                 disp('STOPPED');
%                 keyboard
%             end
            DATBYREND.FF_dev = [DATBYREND.FF_dev; ffdev];
            DATBYREND.Time_dev = [DATBYREND.Time_dev; timedev];
            DATBYREND.Time_dev_targ = [DATBYREND.Time_dev_targ; timedev_TARG];
            
            
            % ===================== THINGS ABOUT TRAINING FEEDBACK
            % --- catch
            if isfield(TrialStruct.birds(i).exptnum(ii).sylnum(indthis), 'isCatch')
                iscatch = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).isCatch;
            else
                iscatch = nan(size(ffdev));
            end
            
            % --- WN hit
            if isfield(TrialStruct.birds(i).exptnum(ii).sylnum(indthis), 'isWNhit')
                isWN = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).isWNhit;
            else
                isWN = nan(size(ffdev));
            end
            
            % ----- is this using song or syl as rend?
            if length(iscatch) ~= length(ffdev)
                % --- then is song as rend...
                keyboard
                
            else
                % --- is syl by rend ...
                DATBYREND.IsWN = [DATBYREND.IsWN; isWN];
                DATBYREND.IsCatch = [DATBYREND.IsCatch; iscatch];
            end
            
            % ===================== IS THIS DUERING TRAINING OR BASELINE?
            istrain = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).IsDurTrain;
            DATBYREND.IsDurTrain = [DATBYREND.IsDurTrain; istrain];
            
            
            %             if i==14 & ii==1 & indthis==7
            %                 keyboard
            %             end
            
            % ======================
%                 if i==1 & ii==4
%                     keyboard
%                 end
try
            nhits = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).WN_nhits;
            nmiss = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).WN_nmiss;
catch err
    nhits = nan;
    nmiss = nan;
end
          DATBYREND.WN_hits = [DATBYREND.WN_hits; nhits];
          DATBYREND.WN_miss = [DATBYREND.WN_miss; nmiss];

            
            % ============================ LOCAL LAERNING
            learnlocal = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).LearnLocal;
            DATBYREND.LearnLocal = [DATBYREND.LearnLocal; learnlocal];
            
            learnlocal_targ = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).LearnLocal_Targ;
            DATBYREND.LearnLocal_targ = [DATBYREND.LearnLocal_targ; learnlocal_targ];
            % --- sanity check
            if (0)
                lt_figure; hold on;
                ftmp = [TrialStruct.birds(i).exptnum(ii).sylnum(indthis).FFvals];
                [nan diff(ftmp)']
                ftmp2 = [TrialStruct.birds(i).exptnum(ii).sylnum(indthis).LearnLocal]
                figure; plot([nan diff(ftmp)'], ftmp2, 'ok')
                xlabel('local learnnig');
                ylabel('ff diffs');
            end
            
            
            % ================== NOTE DOWN WHETHER THIS EXPERIMENT IS
            % EXTRACTED FROM GENERLAIZAITON STRUCT OR SEQDEPPITCJH
            isfrom_SDP = isfield(TrialStruct.birds(i).exptnum(ii).sylnum(1), ...
                'INFO_SylDimensions'); % only SDP experiments have this. all SDP expts have this.
            
            
%             if sylcount==30
%                 keyboard
%             end
%             
            % ###################### KEEP TRACK OF EXPT/SYL
            DATBYREND.Birdnum= [DATBYREND.Birdnum; i*ones(size(istrain))];
            DATBYREND.Exptnum= [DATBYREND.Exptnum; ii*ones(size(istrain))];
            DATBYREND.Sylnum= [DATBYREND.Sylnum; indthis*ones(size(istrain))];
            DATBYREND.Sylcounter = [DATBYREND.Sylcounter; sylcount*ones(size(istrain))];
            DATBYREND.IsSame = [DATBYREND.IsSame; samethis*ones(size(istrain))];
            DATBYREND.IsTarg = [DATBYREND.IsTarg; istargthis*ones(size(istrain))];
            
            DATBYREND.Isfrom_SDP = [DATBYREND.Isfrom_SDP; isfrom_SDP*ones(size(istrain))];
            
            sylcount = sylcount+1;
            
            
        end
    end
end
