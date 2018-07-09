function DATBYREND = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S2(TrialStruct, ...
    singleRendOnly, densitymethod, cutoffmethod, mintime_fromref, ...
    maxtime_fromref, ignoreLMANexpt)
%% lt 6/13/18 - pulls out all data by rend (for all nontarget syls)

if singleRendOnly==1
    % needs to be this to make sense ...
    densitymethod = 'beforefirstnontarg';
end


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
        
        % ############################################ SAME
        sylinds_this = find([TrialStruct.birds(i).exptnum(ii).sylnum.INFO_istarget]==0);
        
        % ===================== GO THRU ALL NONTARG SYLS. FOR EACH REND COLLECT
        % DATA
        for j=1:length(sylinds_this)
            
            indthis = sylinds_this(j);
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
            WNon = TrialStruct.birds(i).exptnum(ii).WNontime;
            ttmp = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).Tvals;
            fftmp = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).FFvals;
            indstrain = ttmp>WNon;
            [b,bint] = lt_regress(fftmp(indstrain), ttmp(indstrain), 0);
            isSigLearn = sign(bint(2,1))==sign(bint(2,2)); % if CI of slope does not overlap 0;
            
            DATBYREND.SigLearn  = [DATBYREND.SigLearn; isSigLearn*ones(size(ishighdensity))];
            
            
            
            
            
            
            
            % ############################## COLLECT FF AND TIME DEVIATIONS
            timedev = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_TimeDev;
            ffdev = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_FFDev;
            timedev_TARG = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_TimeDevTargs;
            istrain = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).IsDurTrain;
            
            DATBYREND.FF_dev = [DATBYREND.FF_dev; ffdev];
            DATBYREND.Time_dev = [DATBYREND.Time_dev; timedev];
            DATBYREND.Time_dev_targ = [DATBYREND.Time_dev_targ; timedev_TARG];
            DATBYREND.IsDurTrain = [DATBYREND.IsDurTrain; istrain];

            
            
            
            
            
            % ###################### KEEP TRACK OF EXPT/SYL
            
            DATBYREND.Birdnum= [DATBYREND.Birdnum; i*ones(size(istrain))];
            DATBYREND.Exptnum= [DATBYREND.Exptnum; ii*ones(size(istrain))];
            DATBYREND.Sylnum= [DATBYREND.Sylnum; indthis*ones(size(istrain))];
            DATBYREND.Sylcounter = [DATBYREND.Sylcounter; sylcount*ones(size(istrain))];
            DATBYREND.IsSame = [DATBYREND.IsSame; samethis*ones(size(istrain))];
            DATBYREND.IsTarg = [DATBYREND.IsTarg; istargthis*ones(size(istrain))];
                        
            sylcount = sylcount+1;
            
            
        end
    end
end
