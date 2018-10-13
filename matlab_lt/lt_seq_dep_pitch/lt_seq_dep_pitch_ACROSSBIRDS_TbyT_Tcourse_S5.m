function TrialStruct = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S5(TrialStruct, ...
    ignoreLMANexpt, tflankplot, twind_plot, flanktime_targ, flipTargNontarg, ...
    useMeanFFRef, collectTarg, songasrend)
%%
% twind_plot = [-tflankplot tflankplot]; % for rend-locking, how much flanking time to plot (hrs)
if ~exist('songasrend', 'var')
    collect_targ_learn=0;
elseif songasrend==1
    % i.e. can only collect trial locked stuff at target if looking at
    % songs.
    collect_targ_learn=1;
else
    collect_targ_learn = 0;
end

%% by default gets just same and diff.

if ~exist('collectTarg', 'var')
    collectTarg = [];
end
if isempty(collectTarg)
    collectTarg = 0;
end

%% note: checked closely on 7/31, everything is correct.

Numbirds = length(TrialStruct.birds);

%% for sanity checks later
ntargs_inref_ALL = [];
n_nontargs_inref_ALL =[];


%%
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
                isLMAN = TrialStruct.birds(i).exptnum(ii).LMANinactivated;
                if isLMAN==1
                    disp(['[is LMAN] skipping ' TrialStruct.birds(i).exptnum(ii).exptname]);
                    continue
                end
            end
        end
        
        
        % =============== for a single target, go thru all the sametype syls
        Sylind_targ = find([TrialStruct.birds(i).exptnum(ii).sylnum.INFO_istarget]); assert(length(Sylind_targ)==1, 'only allowed one targ...');
        
        if collectTarg == 1
        Sylinds_toget = 1:length(TrialStruct.birds(i).exptnum(ii).sylnum); 
        elseif collectTarg == 0
        Sylinds_toget = find([TrialStruct.birds(i).exptnum(ii).sylnum.INFO_istarget]==0);
        end
        
        % ================= COLLECT TARG DATA
        %             ffdev_targ= TrialStruct.birds(i).exptnum(ii).sylnum(indtarg).ff_deviations;
        
        WNon =TrialStruct.birds(i).exptnum(ii).WNontime;
        
        
        % ===================== GO THRU ALL NONTARG SYLS. FOR EACH REND COLLECT
        % DATA
        for j=1:length(Sylinds_toget)
            indthis = Sylinds_toget(j);
            
            if flipTargNontarg==0
                % default
                t_nontarg = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).Tvals;
                ffvals_nontarg = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).FFvals;
                
                t_targ = TrialStruct.birds(i).exptnum(ii).sylnum(Sylind_targ).Tvals;
                ffvals_targ = TrialStruct.birds(i).exptnum(ii).sylnum(Sylind_targ).FFvals;
            elseif flipTargNontarg==1
                t_nontarg = TrialStruct.birds(i).exptnum(ii).sylnum(Sylind_targ).Tvals;
                ffvals_nontarg = TrialStruct.birds(i).exptnum(ii).sylnum(Sylind_targ).FFvals;
                
                t_targ = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).Tvals;
                
            end
            assert(all(diff(t_targ)>=0), 'code below assumes sorted');
            assert(all(diff(t_nontarg)>=0), 'code below assumes sorted');
           
            % -------------- for each NONTARG rendition, find flanking nontarg
            % renditions
            nrends = length(ffvals_nontarg);
            
            Tvals = cell(nrends,1);
            FFvals = cell(nrends,1);
            Ntargs = nan(nrends,1);
            Nnontargs = nan(nrends,1);
            Targtimes = cell(nrends,1);
            IsTrain = nan(nrends,1);
            LearnLocal = nan(nrends, 1);
            LearnLocal_Targ = nan(nrends, 1);
            
            % ======================= GO THRU ALL NTARG RENDS
            for r=1:nrends
                
                %                 twind = twind_plot./24 + t_targ(r);
                %                 ffdev_targ_this = ffdev_targ(r);
                %                 tthis_targ = t_targ(r);
                
                
                
                % ====================== stats for nontarg rend
                t_this = t_nontarg(r);
                ff_this = ffvals_nontarg(r);
                
                % -------- skip iftoo close to edge
                daythis = floor(t_this);
                tthisday = t_nontarg(floor(t_nontarg)==daythis);
                t_nontarg_earliest = min(tthisday);
                t_nontarg_latest = max(tthisday);
                
                tmp1 = t_this - t_nontarg_earliest; % flank time, preceding
                tmp2 = t_nontarg_latest - t_this; % ", following
                mintime = tflankplot/24;
                if tmp1<mintime | tmp2<mintime
                    continue
                end
                
                
                % ==================== COLLECT ALL NONTARG RENDS THAT FLANK
                % ------------ 1) WHAT IS FLANKING WINDOW?
                twind = twind_plot./24 + t_this;
                
                % ------------ 2) COLLECT ALL RENDS IN FLANKING IWNODW
                inds_flank= t_nontarg>twind(1) & t_nontarg<twind(2);
                
                                
                % ------------ 1) WHAT IS REF WINDOW?
                flanktime_day = flanktime_targ/(60*24); % convert from min to day
                twindthis = [t_this-flanktime_day t_this+flanktime_day];
                
                
                % -------------- 2) FIND ALL REF NONTARGS
                ind_refNtarg = t_nontarg>=twindthis(1) & t_nontarg<=twindthis(2);
                assert(any(ind_refNtarg)==1, 'must be, by def');
                n_nontargs_inref_ALL = [n_nontargs_inref_ALL; sum(ind_refNtarg)];
                
                % ---------- 3) WHAT IS FF OF REFERENCE POINT?
                if useMeanFFRef==1
                    ff_ref = nanmean(ffvals_nontarg(ind_refNtarg));
                elseif useMeanFFRef==0
                    ff_ref = ff_this;
                end
                
                
                % ======================= FINALLY: SUBTRACT REF FROM FLANK
                t_flank = t_nontarg(inds_flank) - t_this;
                ff_flank = ffvals_nontarg(inds_flank) ...
                    - ff_ref; % all ff vals, subtract reference point.
                
                
                % ====================== HOW MANY RENDITIONS OF TARGET
                % WITHIN THIS REFERENCE WINDOW?
                indstmp_targ = t_targ>=twindthis(1) & t_targ<=twindthis(2);
                ntargs_inref = sum(indstmp_targ);
                ntargs_inref_ALL = [ntargs_inref_ALL; ntargs_inref];
                
                
                % ===================== SAVE TIMING FOR ALL TARG RENDS
                % WITHIN LARGER FLANK.
                indtmp_targ_large = t_targ>twind(1) & t_targ<twind(2);
                targtimes_this = t_targ(indtmp_targ_large) - t_this;
                
                
                % =====================
                % is this during training or baseline?
                istrain = t_this>WNon;
                
                
                % ==================== STATS FOR LOCAL LEARNING (E.G.
                % RECENT DEVIATION)
                if r>1
                    locallearn = ffvals_nontarg(r) - ffvals_nontarg(r-1);
                else
                    locallearn = nan;
                end
                
                % =================== LOCAL LEARNING FOR TARGET RENDITIONS
                % AT SAME TIME AS THIS RENDITION.
                % -- can only do this if using song as rend
                if collect_targ_learn ==1
                    indtargtmp = find(t_targ == t_this);
                    assert(length(indtargtmp)<2, 'why multiple target rends for this nontarg?');
                    
                    if isempty(indtargtmp)
                        % -- ignore
                        locallearn_targ = nan;
                    elseif indtargtmp<2
                        % ignore
                        locallearn_targ = nan;
                    else
                        % --- collect
                        locallearn_targ = ffvals_targ(indtargtmp) - ffvals_targ(indtargtmp-1);
                    end
                    
                else
                    locallearn_targ = nan;
                end
                
                
                % ======================= SAVE FOR THIS REND
                Tvals{r} = t_flank;
                FFvals{r} = ff_flank;
                Ntargs(r) = ntargs_inref;
                Nnontargs(r) = sum(ind_refNtarg);
                Targtimes{r} = targtimes_this;
                IsTrain(r) = istrain;
                LearnLocal(r) = locallearn;
                LearnLocal_Targ(r) = locallearn_targ;
            end
            
            % =============== SAVE FOR THIS SYL
            if std(cell2mat(FFvals))>500
               disp('asfdas');
                keyboard
            end
            TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_TimeDev = Tvals;
            TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_FFDev = FFvals;
            TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_NTargRendsInRef = Ntargs;
            TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_NNonTargRendsInRef = Nnontargs;
            TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_TimeDevTargs = Targtimes;
            TrialStruct.birds(i).exptnum(ii).sylnum(indthis).IsDurTrain = IsTrain;
            
            TrialStruct.birds(i).exptnum(ii).sylnum(indthis).LearnLocal = LearnLocal;
            TrialStruct.birds(i).exptnum(ii).sylnum(indthis).LearnLocal_Targ = LearnLocal_Targ;
            
        end
        
        
    end
end


lt_figure; hold on;

% ==============
lt_subplot(2,2,1); hold on;
lt_plot_histogram(ntargs_inref_ALL);
title('number of targets in ref point');

% ==============
lt_subplot(2,2,2); hold on;
lt_plot_histogram(n_nontargs_inref_ALL);
title('number of NONtargets in ref point');

% ==============
lt_subplot(2,2,3); hold on;
xlabel('num nontarg in ref point');
ylabel('num TARG in ref point');
plot(n_nontargs_inref_ALL, ntargs_inref_ALL, 'x');

title('number of NONtargets in ref point');