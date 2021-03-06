function lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse(TrialStruct, ParamsTrial, ...
    ignoreLMANexpt, plotraw)
%% lt 5/28/18 - asking how learning in target context relates to nontarg?

windtouseforNontargDev = 'preceding';
% 'preceding': then "timewindow" duration preceeding rend
% 'entire': then +/- timewindow, relative to rend

% ---- Params related to when want to use "reference point" to calcualte
% deviation in pitch

useFFdevWithinSylTiming = 0;
% -1: calculates deviation from time window containint targ syl (but slides
% one flanktime over to left, so that is not in window concurrent with
% target syls
% 0: calcualtes deviation from time window containing target syl
% [windtouseforNontargDev becomes obsolete]
% 1: uses each syl renditions own deviation from its own history

% flanktime_targ = 2.5; % minutes (will do +/- this number). This only used if
% useFFdevWithinSylTiming = 0;
% flanktime_targ = 0.25; % minutes (will do +/- this number). This only used if

nonsingonly = 0; % if 1, then makes sure the targ syl is not sung again between
% rendition of target and second nontarg. NOTE; works for both if useFFdevWithinSylTiming=0 or 1
% will include a nontarg if it is at same time as next targ (i.e.
% inclusive)

% =========== FILTERING NONTARGS BY HOW MUCH LEARNING OVERALL THEY SHOWED
zscorethresh = 0; % end of learning (last half day)

% === defaults:
% 1) use window preceding syl to determing deviation
% windtouseforNontargDev = 'preceding';
% useFFdevWithinSylTiming = 1;
% flanktime_targ = [];

% 1) use deviation from nontarg at time of target
% windtouseforNontargDev = 'preceding';
% useFFdevWithinSylTiming = 0;
% flanktime_targ = 2; % 2 minutes?

% control analyses
windtouseforNontargDev = 'preceding';
useFFdevWithinSylTiming = 0;
flanktime_targ = 0.25; % minutes (will do +/- this number). This only used if
nonsingonly = 0; % if 1, then makes sure the targ syl is not sung again between

% =========== FILTERING NONTARGS BY HOW MUCH LEARNING OVERALL THEY SHOWED
zscorethresh = 0; % end of learning (last half day)


tflankplot = 0.75; % for rend-locking, how much flanking time to plot (hrs) (i.e. for each
% targ rend, collect his much nontarg data. also, ignore any targ rends
% that are this close to start or end of day.
timewindow = 0.75; % in hours, for calculating pitch deviations.
binsize = 20; % renditions for running mean

twind_plot = [-tflankplot tflankplot]; % for rend-locking, how much flanking time to plot (hrs)

% ================= CONTROL ANALYSES (2ND HALF OF CODE)

%%

Numbirds = length(TrialStruct.birds);


%% CALCULATE FF DEVIATIONS

disp('NOTE TO SELF: ffvals will be automatically flipped to be dir of learning (so no need to flip again later)');
pause;

% DATSTRUCT.AllBirdnum =[];
% DATSTRUCT.AllExptnum =[];
% DATSTRUCT.AllSylname ={};
%
% DATSTRUCT.AllIsTarg = [];
% DATSTRUCT.AllIsSame = [];
%
% DATSTRUCT.AllDayList = [];
% DATSTRUCT.AllTedges = [];
% DATSTRUCT.AllFFedges = [];
%
% DATSTRUCT.AllBasedays = [];

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

            
            % --------------- subtract mean and flip if negative larning
            ff = (ff-ffmean_base).*targlearndir;
            TrialStruct.birds(i).exptnum(ii).sylnum(ss).FFvals = ff;
            
            % ------------- for all renditions calculate deviation from
            % recent trials
            % === method1 - fit regression line to one hour of data
            % (directly preceding this rendition...) record deviation from
            % that hour's prediction
            
            % ==== method 2 - deviation from mean of last hour(or whatever).
            nrends = length(ff);
            ff_deviations = nan(size(ff));
            for j=1:nrends
                
                tthis = t(j);
                ffthis = ff(j);
                
                % ----------------- 1) don't record if this rendition is too early in the
                % day
                daythis = floor(tthis); % today
                t_earliest = min(t(floor(t) == daythis));
                
                timeFromDayOn = tthis - t_earliest; % time from first rend today
                
                if timeFromDayOn<timewindow/24
                    continue
                end
                
                
                
                % ----------------- 2) collect rends in last hour
                if strcmp(windtouseforNontargDev, 'preceding')
                    tmin = tthis - timewindow/24;
                    tmax = tthis;
                elseif strcmp(windtouseforNontargDev, 'entire')
                    tmin = tthis - timewindow/24;
                    tmax = tthis + timewindow/24;
                end
                
                indstmp = t>=tmin & t<tmax;
                ttmp = t(indstmp);
                fftmp = ff(indstmp);
                
                % ----------------- 3) calculate deviation for this trial
                ffdev = ffthis - mean(fftmp);
                
                ff_deviations(j) = ffdev;
                
            end
            
            % ============================ OUTPUT
            TrialStruct.birds(i).exptnum(ii).sylnum(ss).ff_deviations = ff_deviations;
            
        end
    end
end



%% PLOT EXAMPLE TO SEE WHAT TRANSFOMRAITON FROM FF TO FF DEVIATION LOOKS LIKE
if plotraw==1
    
    i=13;
    birdname = TrialStruct.birds(i).birdname;
    
    % Numexpt = length(TrialStruct.birds(i).exptnum);
    ii=1;
    exptname = TrialStruct.birds(i).exptnum(ii).exptname;
    
    % ---------- SKIP IF NO DATA
    % if isempty(TrialStruct.birds(i).exptnum(ii).sylnum)
    %     disp(['[no DATA] skipping ' TrialStruct.birds(i).exptnum(ii).exptname]);
    %     continue
    % end
    %
    % % --------- ignore if lMAN?
    % if ignoreLMANexpt==1
    %     isLMAN = TrialStruct.birds(i).exptnum(ii).LMANinactivated;
    %     if isLMAN==1
    %         disp(['[is LMAN] skipping ' TrialStruct.birds(i).exptnum(ii).exptname]);
    %         continue
    %     end
    % end
    
    figcount=1;
    subplotrows=4;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    Numsyls = length(TrialStruct.birds(i).exptnum(ii).sylnum);
    plotcols = lt_make_plot_colors(Numsyls, 0,0);
    hsplots = [];
    
    
    for ss=1:Numsyls
        
        t = TrialStruct.birds(i).exptnum(ii).sylnum(ss).Tvals;
        ff = TrialStruct.birds(i).exptnum(ii).sylnum(ss).FFvals;
        ffdev = TrialStruct.birds(i).exptnum(ii).sylnum(ss).ff_deviations;
        istarg = TrialStruct.birds(i).exptnum(ii).sylnum(ss).INFO_istarget;
        issame = TrialStruct.birds(i).exptnum(ii).sylnum(ss).INFO_similar;
        sylname = TrialStruct.birds(i).exptnum(ii).sylnum(ss).syl;
        
        if istarg==1
            pcol = 'k';
            titstr = [sylname '[' birdname '-' exptname ']'];
        elseif istarg==0 & issame==1
            pcol ='b';
            titstr = [sylname];
        elseif istarg==0 & issame==0
            pcol ='r';
            titstr = [sylname];
        end
        
        
        % =========== plot raw trial by trial data
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(titstr);
        plot(t, ff, 'o', 'Color', pcol);
        
        % =============== wn on time
        wnon = TrialStruct.birds(i).exptnum(ii).WNontime;
        line([wnon wnon], ylim, 'Color', 'm');
        
        % ============ plot deviations
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(titstr);
        plot(t, ffdev, 's', 'Color', pcol);
        lt_plot_zeroline;
        
        % =============== wn on time
        line([wnon wnon], ylim, 'Color', 'm');
        
        % -----------------
        linkaxes(hsplots, 'x');
        
        
    end
end


%% ============== [COLLECT] "targ locked deviation of nontarg"
% ---- for each rendition of targ, plot change in nontarg as function of
% teporal lag from the target.
% ---- make separate plots for targ ff that is posistive (i.e. laerning)
% vs. negative.


AllBirdnum = [];
AllExptnum =[];
AllSylNum = [];
AllIssame = [];
AllTvals_nontarg = {};
AllFFvals_nontarg = {};
AllFFvals_targ = {};
AllIndsBase_targ = {};

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
        
        % ############################################ SAME
        issame=1;
        
        % =============== for a single target, go thru all the sametype syls
        sylind_targ = find([TrialStruct.birds(i).exptnum(ii).sylnum.INFO_istarget]);
        sylinds_same = find([TrialStruct.birds(i).exptnum(ii).sylnum.INFO_similar]==issame ...
            & [TrialStruct.birds(i).exptnum(ii).sylnum.INFO_istarget]==0);
        
        % -------------- figures
        figcount=1;
        subplotrows=4;
        subplotcols=2;
        fignums_alreadyused=[];
        hfigs=[];
        hsplots = [];
        
        
        for j=1:length(sylinds_same)
            
            indthis = sylinds_same(j);
            indtarg = sylind_targ;
            
            t_nontarg = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).Tvals;
            t_targ = TrialStruct.birds(i).exptnum(ii).sylnum(indtarg).Tvals;
            ffdev_nontarg = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).ff_deviations;
            ffdev_targ= TrialStruct.birds(i).exptnum(ii).sylnum(indtarg).ff_deviations;
            sylname = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).syl;
            
            ffvals_nontarg = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).FFvals;
            
            assert(all(diff(t_targ)>=0), 'code below assumes sorted');
            % -------------- for each targ rendition, find flanking nontarg
            % renditions
            nrends = length(ffdev_targ);
            Tvals = cell(nrends,1);
            FFvals = cell(nrends,1);
            FFdev_targ = nan(nrends,1);
            for r =1:nrends
                
                twind = twind_plot./24 + t_targ(r);
                ffdev_targ_this = ffdev_targ(r);
                tthis_targ = t_targ(r);
                
                % -------- skip if i) targ dev is nan or ii) too close to
                % edge
                if isnan(ffdev_targ_this)
                    continue
                end
                daythis = floor(t_targ(r));
                t_nontarg_earliest = min(t_nontarg(floor(t_nontarg)==daythis));
                t_nontarg_latest = max(t_nontarg(floor(t_nontarg)==daythis));
                
                tmp1 = t_targ(r) - t_nontarg_earliest; % flank time, preceding
                tmp2 = t_nontarg_latest - t_targ(r); % ", following
                mintime = tflankplot/24;
                if tmp1<mintime | tmp2<mintime
                    continue
                end
                
                
                % -------- collect nontarg trials that are in flanking
                % winodw
                indsthis_nontarg = t_nontarg>twind(1) & t_nontarg<twind(2);
                
                tthis_nontarg = t_nontarg(indsthis_nontarg) - t_targ(r); % time relative to time of target
                
                if useFFdevWithinSylTiming==1
                    % -------- then use each syls deviation from it's own
                    % recent histroy
                    ffdevthis_nontarg = ffdev_nontarg(indsthis_nontarg);
                elseif useFFdevWithinSylTiming==0 | useFFdevWithinSylTiming==-1
                    
                    % -------- then get FF of each syl relative to the
                    % timing of the target syl
                    flanktime_day = flanktime_targ/(60*24); % convert from min to day
                    
                    
                    % ------------- WHAT IS REFERENCE WINDOW?
                    if useFFdevWithinSylTiming==0
                        % WINDOW = window concurrent with targ syl
                        twindthis = [tthis_targ-flanktime_day tthis_targ+flanktime_day];
                    elseif useFFdevWithinSylTiming==-1
                        % WINDOW = window one timebin shifted left of targ
                        % syl (not concurrent)
                        twindthis = [tthis_targ-3*flanktime_day tthis_targ-flanktime_day];
                    end
                    
                    
                    % ---- find nontarg renditions that are in this window
                    indtmpthis = t_nontarg>=twindthis(1) & t_nontarg<=twindthis(2);
                    if ~any(indtmpthis)
                        % then there is no reference point for nontarg.
                        % IGNORE THIS RENDITION
                        disp('skipping rend... lacking nontarg in reference point');
                        ffdevthis_nontarg = [];
                        tthis_nontarg = [];
                        ffdev_targ_this = nan;
                    else
                        % ----- get all ff values for nontarget as
                        % difference from reference point
                        ff_ref = mean(ffvals_nontarg(indtmpthis));
                        
                        ffdevthis_nontarg = ffvals_nontarg(indsthis_nontarg) ...
                            - ff_ref; % all ff vals, subtract reference point.
                    end
                end
                
                % ====================== WANT TO ONLY KEEP NONTARGS THAT
                % OCCUR WITHOUT EXPRESSION OF TARG?
                if nonsingonly==1
                    % then find rends that occur after next nontarg, remove
                    % them
                    % ---- what is time of next targ rend?
                    tthis_nextnontarg = t_targ(r+1);
                    
                    indstoremove = tthis_nontarg+t_targ(r)>tthis_nextnontarg;
                    
                    % -------- remove those rends too late
                    tthis_nontarg(indstoremove) = [];
                    ffdevthis_nontarg(indstoremove) = [];
                end
                
                Tvals{r} = tthis_nontarg;
                FFvals{r} = ffdevthis_nontarg;
                FFdev_targ(r) = ffdev_targ_this;
            end
            
            % ===================== WHICH IS BASELINE/TRAINING RENDITIONS
            WNontime = TrialStruct.birds(i).exptnum(ii).WNontime;
            inds_base = t_targ<WNontime;
            
            %             basedaylast = TrialStruct.birds(i).exptnum(ii).BaseDays(end);
            %             inds_base = t_targ<basedaylast+1;
            
            
            % ===============  PLOT
            if rand<0
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
                title(sylname);
                for r=1:nrends
                    t = Tvals{r};
                    ff = FFvals{r};
                    ff_targ = FFdev_targ(r);
                    
                    if ff_targ>0
                        pcol = 'k';
                    elseif ff_targ<=0
                        pcol = 'r';
                    end
                    
                    plot(t, ff, 'x', 'Color', pcol);
                    
                end
            end
            
            % ================== OUTPUT
            AllBirdnum = [AllBirdnum; i];
            AllExptnum =[AllExptnum; ii];
            AllIssame = [AllIssame; issame];
            AllSylNum = [AllSylNum; indthis];
            AllTvals_nontarg = [AllTvals_nontarg; {Tvals}];
            AllFFvals_nontarg = [AllFFvals_nontarg; {FFvals}];
            AllFFvals_targ = [AllFFvals_targ; {FFdev_targ}];
            AllIndsBase_targ = [AllIndsBase_targ; {inds_base}];
            
        end
        
        
        % ============ combined plot
        
    end
end





%% ============== [PLOT RAW] for each experiment plot some examples at level of rendition

Maxbirds = max(AllBirdnum);
Maxexpts = max(AllExptnum);



% -------------- figures
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% ------------------
nrendstoplot = 20;
if nonsingonly==1
    nrendstoplot = 10*nrendstoplot; % will be less dense, accomodate
end
for i=1:Maxbirds
    
    for ii=1:Maxexpts
        
        %  ############################ SAME TYPE
        inds = find(AllBirdnum==i & AllExptnum==ii & AllIssame==1);
        if isempty(inds)
            continue
        end
        
        tvals_all = [];
        ffvals_all = [];
        fftarg_all = [];
        for j=1:length(inds)
            indthis = inds(j);
            
            % -------------- DATA FOR THIS SYL
            tvals_ntarg = AllTvals_nontarg{indthis};
            ffvals_ntarg = AllFFvals_nontarg{indthis};
            ffvals_targ = AllFFvals_targ{indthis};
            
            % ------------- FIGURE
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            bname = TrialStruct.birds(i).birdname;
            ename = TrialStruct.birds(i).exptnum(ii).exptname;
            sylname = TrialStruct.birds(i).exptnum(ii).sylnum(AllSylNum(indthis)).syl;
            title([bname '-' ename '-' sylname]);
            ylabel('dev (b=targup; r=targdn)');
            
            % ================ GO THRU ALL RENDS (OF TARG) AND PLOT
            nrends = length(tvals_ntarg);
            
            % ==================== PLOT SUBSET OF RENDS
            if nrends>nrendstoplot
                rendslist = randperm(nrends, nrendstoplot);
            else
                rendslist = 1:nrends;
                
            end
            for r=rendslist
                t = tvals_ntarg{r}*(24*60); % convert to minutes
                ff = ffvals_ntarg{r};
                ff_targ = ffvals_targ(r);
                
                if isempty(ff)
                    continue
                end
                
                if ff_targ>0
                    pcol = 'b';
                elseif ff_targ<=0
                    pcol = 'r';
                end
                
                plot(t, ff, '-x', 'Color', pcol);
                
                
                % ======= collect
                tvals_all = [tvals_all; t];
                ffvals_all = [ffvals_all; ff];
                fftarg_all = [fftarg_all; ff_targ*ones(size(ff))];
            end
            
            
            % -----------------------
            lt_plot_zeroline;
            lt_plot_zeroline_vert;
            xlim([-60 60]);
        end
    end
end



%% ============ [PLOT] - raw and individual syls
if plotraw==1
    disp('NOTE: might be mixing baseline and WN on periods ... CHECK!!');
    pause;
    Maxbirds = max(AllBirdnum);
    Maxexpts = max(AllExptnum);
    
    
    hsplots = [];
    for i=1:Maxbirds
        
        for ii=1:Maxexpts
            
            %  ############################ SAME TYPE
            inds = find(AllBirdnum==i & AllExptnum==ii & AllIssame==1);
            if isempty(inds)
                continue
            end
            
            % -------------- figures
            figcount=1;
            subplotrows=4;
            subplotcols=2;
            fignums_alreadyused=[];
            hfigs=[];
            hsplots = [];
            
            tvals_all = [];
            ffvals_all = [];
            fftarg_all = [];
            for j=1:length(inds)
                indthis = inds(j);
                
                % -------------- DATA FOR THIS SYL
                tvals_ntarg = AllTvals_nontarg{indthis};
                ffvals_ntarg = AllFFvals_nontarg{indthis};
                ffvals_targ = AllFFvals_targ{indthis};
                
                % ------------- FIGURE
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
                %             title(sylname);
                
                % ================ GO THRU ALL RENDS (OF TARG) AND PLOT
                nrends = length(tvals_ntarg);
                for r=1:nrends
                    t = tvals_ntarg{r}*(24*60); % convert to minutes
                    ff = ffvals_ntarg{r};
                    ff_targ = ffvals_targ(r);
                    
                    if isempty(ff)
                        continue
                    end
                    
                    if ff_targ>0
                        pcol = 'k';
                    elseif ff_targ<=0
                        pcol = 'r';
                    end
                    
                    plot(t, ff, 'x', 'Color', pcol);
                    
                    % ======= collect
                    tvals_all = [tvals_all; t];
                    ffvals_all = [ffvals_all; ff];
                    fftarg_all = [fftarg_all; ff_targ*ones(size(ff))];
                end
                
            end
            
            
            % ============ combined plot
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title('all nontargs combined');
            
            % --- targ positive learn
            indtmp = fftarg_all>0;
            pcol = 'k';
            
            x = tvals_all(indtmp);
            y = ffvals_all(indtmp);
            plot(x, y, 'x', 'Color', pcol);
            
            % ---- plot running mean
            [~, indsort] = sort(x); % sort
            x = x(indsort);
            y = y(indsort);
            xrun = lt_running_stats(x, binsize);
            yrun = lt_running_stats(y, binsize);
            
            shadedErrorBar(xrun.Median, yrun.Mean, yrun.SEM, {'Color', pcol}, 1);
            
            % --- targ negative learn
            indtmp = fftarg_all<=0;
            pcol = 'r';
            
            x = tvals_all(indtmp);
            y = ffvals_all(indtmp);
            plot(x, y, 'x', 'Color', pcol);
            
            % ---- plot running mean
            [~, indsort] = sort(x); % sort
            x = x(indsort);
            y = y(indsort);
            xrun = lt_running_stats(x, binsize);
            yrun = lt_running_stats(y, binsize);
            
            shadedErrorBar(xrun.Median, yrun.Mean, yrun.SEM, {'Color', pcol}, 1);
        end
    end
    
end
%% ################# [GRAND MEAN] - across all renditions, all expts/syls


Maxbirds = max(AllBirdnum);
Maxexpts = max(AllExptnum);

hsplots = [];
tvals_all = [];
ffvals_all = [];
fftarg_all = [];
isbase_all = [];
birdnum_all =[];
exptnum_all = [];
sylnum_all = [];

for i=1:Maxbirds
    
    for ii=1:Maxexpts
        disp([num2str(i) '-' num2str(ii)]);
        
        %  ############################ SAME TYPE
        inds = find(AllBirdnum==i & AllExptnum==ii & AllIssame==1);
        if isempty(inds)
            continue
        end
        
        for j=1:length(inds)
            indthis = inds(j);
            
            % -------------- DATA FOR THIS SYL
            tvals_ntarg = AllTvals_nontarg{indthis};
            ffvals_ntarg = AllFFvals_nontarg{indthis};
            ffvals_targ = AllFFvals_targ{indthis};
            indsbase = AllIndsBase_targ{indthis};
            syl_ntarg = AllSylNum(indthis);
            
            % ================ GO THRU ALL RENDS (OF TARG) AND PLOT
            nrends = length(tvals_ntarg);
            for r=1:nrends
                t = tvals_ntarg{r}*(24*60); % convert to minutes
                ff = ffvals_ntarg{r};
                ff_targ = ffvals_targ(r);
                indbase_this = indsbase(r);
                
                % ======= collect
                tvals_all = [tvals_all; t];
                ffvals_all = [ffvals_all; ff];
                fftarg_all = [fftarg_all; ff_targ*ones(size(ff))];
                isbase_all = [isbase_all; indbase_this*ones(size(ff))];
                
                birdnum_all =[birdnum_all; i*ones(size(ff))];
                exptnum_all = [exptnum_all; ii*ones(size(ff))];
                sylnum_all = [sylnum_all; syl_ntarg*ones(size(ff))];
            end
            
        end
        
    end
end

%% [PLOT] GRAND MEAN, ALL RENDS ################################
usemedian = 0;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lt_figure; hold on;
xedges = -65:10:65; % minutes
% xedges = -62.5:5:62.5; % minutes
% xedges = -61.25:2.5:61.25; % minutes
binsize = xedges(2)-xedges(1);
xcenters = xedges(1:end-1)+binsize/2;

Ytrain = nan(1, length(xcenters)); % for combined data
Ybase = nan(1, length(xcenters));

% ################################################ BASELINE
lt_subplot(3,1,1); hold on;
plotbase =1;
title('[BASE] all renditions [GRAND MEAN]');
% binsize = round(length(fftarg_all)/100);
ylabel('ff dev (b=targpos; k=targneg)');


% =============================== 1) TARGET POSITIVE LEARINING
indtmp = fftarg_all>0 & isbase_all==plotbase;
pcol = 'b';

x = tvals_all(indtmp);
y = ffvals_all(indtmp);

% ---------------- slide into bins
xbins = discretize(x, xedges);
assert(~any(isnan(xbins)), 'edges not wide enough ...');

[ymean, ysem] = grpstats(y, xbins, {'mean', 'sem'});
X = xcenters(unique(xbins));
lt_plot(X, ymean, {'Errors', ysem, 'Color', pcol});
plot(X, ymean, '-', 'Color', pcol);


% =============================== 1) TARGET NEGATIVE LEARINING
indtmp = fftarg_all<=0 & isbase_all==plotbase;
pcol = 'r';

x = tvals_all(indtmp);
y = ffvals_all(indtmp);

% ---------------- slide into bins
xbins = discretize(x, xedges);
assert(~any(isnan(xbins)), 'edges not wide enough ...');

[ymean, ysem] = grpstats(y, xbins, {'mean', 'sem'});
X = xcenters(unique(xbins));
lt_plot(X, ymean, {'Errors', ysem, 'Color', pcol});
plot(X, ymean, '-', 'Color', pcol);


% =========================== 3) ALL COMBINED (IGNORE DIR OF TARG DEV)
indtmp = isbase_all==plotbase;
pcol = 'K';

x = tvals_all(indtmp);
y = ffvals_all(indtmp);

% ---------------- slide into bins
xbins = discretize(x, xedges);
assert(~any(isnan(xbins)), 'edges not wide enough ...');

[ymean, ysem] = grpstats(y, xbins, {'mean', 'sem'});
X = xcenters(unique(xbins));
lt_plot(X, ymean, {'Errors', ysem, 'Color', pcol});
plot(X, ymean, '-', 'Color', pcol);

Ybase(unique(xbins)) = ymean;

% ----------------
lt_plot_zeroline;
ylim([-20 20]);
lt_plot_zeroline_vert;

% ################################################# DURING TRAINING
lt_subplot(3,1,2); hold on;
plotbase =0;
title('[TRAINING] all renditions [GRAND MEAN]');
% binsize = round(length(fftarg_all)/100);
ylabel('ff dev (b=targpos; k=targneg)');


% =============================== 1) TARGET POSITIVE LEARINING
indtmp = fftarg_all>0 & isbase_all==plotbase;
pcol = 'b';

x = tvals_all(indtmp);
y = ffvals_all(indtmp);

% ---------------- slide into bins
xbins = discretize(x, xedges);
assert(~any(isnan(xbins)), 'edges not wide enough ...');

[ymean, ysem] = grpstats(y, xbins, {'mean', 'sem'});
X = xcenters(unique(xbins));
lt_plot(X, ymean, {'Errors', ysem, 'Color', pcol});
plot(X, ymean, '-', 'Color', pcol);


% =============================== 1) TARGET NEGATIVE LEARINING
indtmp = fftarg_all<=0 & isbase_all==plotbase;
pcol = 'r';

x = tvals_all(indtmp);
y = ffvals_all(indtmp);

% ---------------- slide into bins
xbins = discretize(x, xedges);
assert(~any(isnan(xbins)), 'edges not wide enough ...');

[ymean, ysem] = grpstats(y, xbins, {'mean', 'sem'});
X = xcenters(unique(xbins));
lt_plot(X, ymean, {'Errors', ysem, 'Color', pcol});
plot(X, ymean, '-', 'Color', pcol);


% =========================== 3) ALL COMBINED (IGNORE DIR OF TARG DEV)
indtmp = isbase_all==plotbase;
pcol = 'K';

x = tvals_all(indtmp);
y = ffvals_all(indtmp);

% ---------------- slide into bins
xbins = discretize(x, xedges);
assert(~any(isnan(xbins)), 'edges not wide enough ...');

[ymean, ysem] = grpstats(y, xbins, {'mean', 'sem'});
X = xcenters(unique(xbins));
lt_plot(X, ymean, {'Errors', ysem, 'Color', pcol});
plot(X, ymean, '-', 'Color', pcol);

Ytrain(unique(xbins)) = ymean;


% --------------
lt_plot_zeroline;
lt_plot_zeroline_vert;
ylim([-20 20]);

% ################################################ TRAINING MINUS BASELINE
Y = Ytrain - Ybase;
lt_subplot(3,1,3); hold on;
title('train minus base (bk lines)');
plot(xcenters, Y, '-ok');



% ----------------
lt_plot_zeroline;
ylim([-20 20]);


%% ################# [EXPTS] - across all renditions, all expts/syls
%         xedges = -115:10:115; % minutes
%         binsize = xedges(2)-xedges(1);
%         xcenters = xedges(1:end-1)+binsize/2;

for i=1:Maxbirds
    for ii=1:Maxexpts
        
        inds = birdnum_all==i & exptnum_all==ii;
        if ~any(inds)
            continue
        end
        
        bname = TrialStruct.birds(i).birdname;
        ename = TrialStruct.birds(i).exptnum(ii).exptname;
        
        hsplots = [];
        lt_figure; hold on;
        
        Ytrain = nan(1, length(xcenters)); % for combined data
        Ybase = nan(1, length(xcenters));
        
        % ################################################ BASELINE
        hsplot = lt_subplot(3,1,1); hold on;
        hsplots = [hsplots hsplot];
        plotbase =1;
        title(['[BASE], grandmean, ' bname '-' ename]);
        % binsize = round(length(fftarg_all)/100);
        ylabel('ff dev (b=targpos; k=targneg)');
        
        
        % =============================== 1) TARGET POSITIVE LEARINING
        indtmp = fftarg_all>0 & isbase_all==plotbase & birdnum_all==i & exptnum_all==ii;
        pcol = 'b';
        
        x = tvals_all(indtmp);
        y = ffvals_all(indtmp);
        
        % ---------------- slide into bins
        xbins = discretize(x, xedges);
        assert(~any(isnan(xbins)), 'edges not wide enough ...');
        
        if usemedian==1
            [ymean, ysem] = grpstats(y, xbins, {'median', 'sem'});
        else
            [ymean, ysem] = grpstats(y, xbins, {'mean', 'sem'});
        end
        X = xcenters(unique(xbins));
        lt_plot(X, ymean, {'Errors', ysem, 'Color', pcol});
        plot(X, ymean, '-', 'Color', pcol);
        
        
        % =============================== 1) TARGET NEGATIVE LEARINING
        indtmp = fftarg_all<=0 & isbase_all==plotbase & birdnum_all==i & exptnum_all==ii;
        pcol = 'r';
        
        x = tvals_all(indtmp);
        y = ffvals_all(indtmp);
        
        % ---------------- slide into bins
        xbins = discretize(x, xedges);
        assert(~any(isnan(xbins)), 'edges not wide enough ...');
        
        if usemedian==1
            [ymean, ysem] = grpstats(y, xbins, {'median', 'sem'});
        else
            [ymean, ysem] = grpstats(y, xbins, {'mean', 'sem'});
        end
        X = xcenters(unique(xbins));
        lt_plot(X, ymean, {'Errors', ysem, 'Color', pcol});
        plot(X, ymean, '-', 'Color', pcol);
        
        
        % =========================== 3) ALL COMBINED (IGNORE DIR OF TARG DEV)
        indtmp = isbase_all==plotbase & birdnum_all==i & exptnum_all==ii;
        pcol = 'K';
        
        x = tvals_all(indtmp);
        y = ffvals_all(indtmp);
        
        % ---------------- slide into bins
        xbins = discretize(x, xedges);
        assert(~any(isnan(xbins)), 'edges not wide enough ...');
        
        if usemedian==1
            [ymean, ysem] = grpstats(y, xbins, {'median', 'sem'});
        else
            [ymean, ysem] = grpstats(y, xbins, {'mean', 'sem'});
        end
        X = xcenters(unique(xbins));
        lt_plot(X, ymean, {'Errors', ysem, 'Color', pcol});
        plot(X, ymean, '-', 'Color', pcol);
        
        Ybase(unique(xbins)) = ymean;
        
        % ----------------
        lt_plot_zeroline;
        ylim([-20 20]);
        
        
        % ################################################# DURING TRAINING
        hsplot = lt_subplot(3,1,2); hold on;
        hsplots = [hsplots hsplot];
        plotbase =0;
        title('[TRAINING] all renditions [GRAND MEAN]');
        % binsize = round(length(fftarg_all)/100);
        ylabel('ff dev (b=targpos; k=targneg)');
        
        
        % =============================== 1) TARGET POSITIVE LEARINING
        indtmp = fftarg_all>0 & isbase_all==plotbase & birdnum_all==i & exptnum_all==ii;
        pcol = 'b';
        
        x = tvals_all(indtmp);
        y = ffvals_all(indtmp);
        
        % ---------------- slide into bins
        xbins = discretize(x, xedges);
        assert(~any(isnan(xbins)), 'edges not wide enough ...');
        
        if usemedian==1
            [ymean, ysem] = grpstats(y, xbins, {'median', 'sem'});
        else
            [ymean, ysem] = grpstats(y, xbins, {'mean', 'sem'});
        end
        X = xcenters(unique(xbins));
        lt_plot(X, ymean, {'Errors', ysem, 'Color', pcol});
        plot(X, ymean, '-', 'Color', pcol);
        
        
        % =============================== 1) TARGET NEGATIVE LEARINING
        indtmp = fftarg_all<=0 & isbase_all==plotbase & birdnum_all==i & exptnum_all==ii;
        pcol = 'r';
        
        x = tvals_all(indtmp);
        y = ffvals_all(indtmp);
        
        % ---------------- slide into bins
        xbins = discretize(x, xedges);
        assert(~any(isnan(xbins)), 'edges not wide enough ...');
        
        if usemedian==1
            [ymean, ysem] = grpstats(y, xbins, {'median', 'sem'});
        else
            [ymean, ysem] = grpstats(y, xbins, {'mean', 'sem'});
        end
        X = xcenters(unique(xbins));
        lt_plot(X, ymean, {'Errors', ysem, 'Color', pcol});
        plot(X, ymean, '-', 'Color', pcol);
        
        
        % =========================== 3) ALL COMBINED (IGNORE DIR OF TARG DEV)
        indtmp = isbase_all==plotbase & birdnum_all==i & exptnum_all==ii;
        pcol = 'K';
        
        x = tvals_all(indtmp);
        y = ffvals_all(indtmp);
        
        % ---------------- slide into bins
        xbins = discretize(x, xedges);
        assert(~any(isnan(xbins)), 'edges not wide enough ...');
        
        if usemedian==1
            [ymean, ysem] = grpstats(y, xbins, {'median', 'sem'});
        else
            [ymean, ysem] = grpstats(y, xbins, {'mean', 'sem'});
        end
        X = xcenters(unique(xbins));
        lt_plot(X, ymean, {'Errors', ysem, 'Color', pcol});
        plot(X, ymean, '-', 'Color', pcol);
        
        Ytrain(unique(xbins)) = ymean;
        
        % --------------
        lt_plot_zeroline;
        
        % ################################################ TRAINING MINUS BASELINE
        Y = Ytrain - Ybase;
        hsplot = lt_subplot(3,1,3); hold on;
        hsplots = [hsplots hsplot];
        title('train minus base (bk lines)');
        plot(xcenters, Y, '-ok');
        
        
        
        % ----------------
        lt_plot_zeroline;
        ylim([-20 20]);
        
        linkaxes(hsplots, 'xy');
        
    end
end

%% =================== [SUMMARY PLOT] - EACH EXPT ONE DATAPOINT

TMPSTRUCT.train.x = {};
TMPSTRUCT.train.x_inds = {};
TMPSTRUCT.train.y = {};
TMPSTRUCT.base.x = {};
TMPSTRUCT.base.x_inds = {};
TMPSTRUCT.base.y = {};


for i=1:Maxbirds
    for ii=1:Maxexpts
        
        inds = birdnum_all==i & exptnum_all==ii;
        if ~any(inds)
            continue
        end
        
        %         bname = TrialStruct.birds(i).birdname;
        %         ename = TrialStruct.birds(i).exptnum(ii).exptname;
        
        % ################################################ TRAINING
        plotbase =0;
        
        % =========================== 3) ALL COMBINED (IGNORE DIR OF TARG DEV)
        indtmp = isbase_all==plotbase & inds;
        
        x = tvals_all(indtmp);
        y = ffvals_all(indtmp);
        
        % ---------------- slide into bins
        xbins = discretize(x, xedges);
        assert(~any(isnan(xbins)), 'edges not wide enough ...');
        if usemedian==1
            [ymean, ysem] = grpstats(y, xbins, {'median', 'sem'});
        else
            [ymean, ysem] = grpstats(y, xbins, {'mean', 'sem'});
        end
        X = xcenters(unique(xbins));
        
        TMPSTRUCT.train.x = [TMPSTRUCT.train.x; X];
        TMPSTRUCT.train.x_inds = [TMPSTRUCT.train.x_inds; unique(xbins)];
        TMPSTRUCT.train.y = [TMPSTRUCT.train.y; ymean];
        
        
        % ################################################# BASELINE
        plotbase =1;
        
        % =========================== 3) ALL COMBINED (IGNORE DIR OF TARG DEV)
        indtmp = isbase_all==plotbase & birdnum_all==i & exptnum_all==ii;
        
        x = tvals_all(indtmp);
        y = ffvals_all(indtmp);
        
        % ---------------- slide into bins
        xbins = discretize(x, xedges);
        assert(~any(isnan(xbins)), 'edges not wide enough ...');
        
        if usemedian==1
            [ymean, ysem] = grpstats(y, xbins, {'median', 'sem'});
        else
            [ymean, ysem] = grpstats(y, xbins, {'mean', 'sem'});
        end
        X = xcenters(unique(xbins));
        
        TMPSTRUCT.base.x = [TMPSTRUCT.base.x; X];
        TMPSTRUCT.base.y = [TMPSTRUCT.base.y; ymean];
        TMPSTRUCT.base.x_inds = [TMPSTRUCT.base.x_inds; unique(xbins)];
    end
end


% ============= [PLOT] ----- overlay all experiments
lt_figure; hold on;

% =================== baseline
fname = 'base';
lt_subplot(3,1,1); hold on;
title(fname);
xlabel('time from targ syl rend (min)');
ylabel('ff deviation');

numexpt = length(TMPSTRUCT.base.x);
Yall = nan(numexpt, length(xcenters)); % expt x bins
for i=1:numexpt
    x = TMPSTRUCT.(fname).x{i};
    y = TMPSTRUCT.(fname).y{i};
    
    plot(x,y, '-o', 'Color', [0.7 0.7 0.7]);
    
    % -------- slide into overall mat
    xinds = TMPSTRUCT.(fname).x_inds{i};
    Yall(i, xinds) = y;
end

% ---- plot mean
ymean = nanmean(Yall,1);
plot(xcenters, ymean, '-ok', 'LineWidth', 3);
lt_plot_zeroline;




% =================== TRAINING
fname = 'train';
lt_subplot(3,1,2); hold on;
title(fname);
xlabel('time from targ syl rend (min)');
ylabel('ff deviation');

numexpt = length(TMPSTRUCT.train.x);
Yall = nan(numexpt, length(xcenters)); % expt x bins
for i=1:numexpt
    x = TMPSTRUCT.(fname).x{i};
    y = TMPSTRUCT.(fname).y{i};
    
    plot(x,y, '-o', 'Color', [0.7 0.7 0.7]);
    
    % -------- slide into overall mat
    xinds = TMPSTRUCT.(fname).x_inds{i};
    Yall(i, xinds) = y;
end

% ---- plot mean
ymean = nanmean(Yall,1);
plot(xcenters, ymean, '-ok', 'LineWidth', 3);
lt_plot_zeroline;



%% =================== [SUMMARY PLOT] - EACH SYL ONE DATAPOINT

Maxsyls = max(sylnum_all);

TMPSTRUCT.train.x = {};
TMPSTRUCT.train.x_inds = {};
TMPSTRUCT.train.y = {};
TMPSTRUCT.base.x = {};
TMPSTRUCT.base.x_inds = {};
TMPSTRUCT.base.y = {};


for i=1:Maxbirds
    for ii=1:Maxexpts
        
        for ss=1:Maxsyls
            
            inds = birdnum_all==i & exptnum_all==ii & sylnum_all==ss;
            if ~any(inds)
                continue
            end
            
            %         bname = TrialStruct.birds(i).birdname;
            %         ename = TrialStruct.birds(i).exptnum(ii).exptname;
            
            % ################################################ TRAINING
            plotbase =0;
            
            % =========================== 3) ALL COMBINED (IGNORE DIR OF TARG DEV)
            indtmp = isbase_all==plotbase & inds;
            
            x = tvals_all(indtmp);
            y = ffvals_all(indtmp);
            
            % ---------------- slide into bins
            xbins = discretize(x, xedges);
            assert(~any(isnan(xbins)), 'edges not wide enough ...');
            if usemedian==1
                [ymean, ysem] = grpstats(y, xbins, {'median', 'sem'});
            else
                [ymean, ysem] = grpstats(y, xbins, {'mean', 'sem'});
            end
            X = xcenters(unique(xbins));
            
            TMPSTRUCT.train.x = [TMPSTRUCT.train.x; X];
            TMPSTRUCT.train.x_inds = [TMPSTRUCT.train.x_inds; unique(xbins)];
            TMPSTRUCT.train.y = [TMPSTRUCT.train.y; ymean];
            
            
            % ################################################# BASELINE
            plotbase =1;
            
            % =========================== 3) ALL COMBINED (IGNORE DIR OF TARG DEV)
            indtmp = isbase_all==plotbase & inds;
            
            x = tvals_all(indtmp);
            y = ffvals_all(indtmp);
            
            % ---------------- slide into bins
            xbins = discretize(x, xedges);
            assert(~any(isnan(xbins)), 'edges not wide enough ...');
            
            if usemedian==1
                [ymean, ysem] = grpstats(y, xbins, {'median', 'sem'});
            else
                [ymean, ysem] = grpstats(y, xbins, {'mean', 'sem'});
            end
            X = xcenters(unique(xbins));
            
            TMPSTRUCT.base.x = [TMPSTRUCT.base.x; X];
            TMPSTRUCT.base.y = [TMPSTRUCT.base.y; ymean];
            TMPSTRUCT.base.x_inds = [TMPSTRUCT.base.x_inds; unique(xbins)];
        end
    end
end


% ============= [PLOT] ----- overlay all experiments
lt_figure; hold on;

% =================== baseline
fname = 'base';
lt_subplot(3,1,1); hold on;
title(fname);
xlabel('time from targ syl rend (min)');
ylabel('ff deviation');

numexpt = length(TMPSTRUCT.base.x);
Yall = nan(numexpt, length(xcenters)); % expt x bins
for i=1:numexpt
    x = TMPSTRUCT.(fname).x{i};
    y = TMPSTRUCT.(fname).y{i};
    
    plot(x,y, '-o', 'Color', [0.7 0.7 0.7]);
    
    % -------- slide into overall mat
    xinds = TMPSTRUCT.(fname).x_inds{i};
    Yall(i, xinds) = y;
end

% ---- plot mean
ymean = nanmean(Yall,1);
plot(xcenters, ymean, '-ok', 'LineWidth', 3);
lt_plot_zeroline;




% =================== TRAINING
fname = 'train';
lt_subplot(3,1,2); hold on;
title(fname);
xlabel('time from targ syl rend (min)');
ylabel('ff deviation');

numexpt = length(TMPSTRUCT.train.x);
Yall = nan(numexpt, length(xcenters)); % expt x bins
for i=1:numexpt
    x = TMPSTRUCT.(fname).x{i};
    y = TMPSTRUCT.(fname).y{i};
    
    plot(x,y, '-o', 'Color', [0.7 0.7 0.7]);
    
    % -------- slide into overall mat
    xinds = TMPSTRUCT.(fname).x_inds{i};
    Yall(i, xinds) = y;
end

% ---- plot mean
ymean = nanmean(Yall,1);
plot(xcenters, ymean, '-ok', 'LineWidth', 3);
lt_plot_zeroline;



%% ###############################################################
%% ############################### [CONTROL]
% FOR EACH RENDITION OF NONTARG, NOTE DOWN HOW MUCH TARG IS SUNG IN THAT
% WINDOW, AND GET OTHER NONTARGS AS DEVIATION.

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

TrialStruct = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S5(TrialStruct, ...
    ignoreLMANexpt, tflankplot, twind_plot, flanktime_targ, flipTargNontarg, ...
    useMeanFFRef);

%
% if flipTargNontarg==0
%     TrialStruct = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S5(TrialStruct, ...
%         ignoreLMANexpt, tflankplot, twind_plot, flanktime_targ, flipTargNontarg);
% elseif flipTargNontarg==1
%     TrialStruct = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S6(TrialStruct, ...
%         ignoreLMANexpt, tflankplot, twind_plot, flanktime_targ, flipTargNontarg);
% end

%% ====== TO COLLECT DATA ND PLOT - OLD VERSION, TOO COMPLICATED  ...
% NOTE: this only does for same type. is a mess to deal with. replaced
% entirely with code below and confirmed that output is identical.

plotEachSylRaw = 0;
% 0: none
% 1: simple stats for each syl
% 2: raw data, overlaied with extracted stats.

lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S1;






%% ====== TO COLLECT DATA TO PLOT (NEW VERSION, COLLECTS FIRST, THEN PLOTS)
% Instead of collecitng and plotting at same time

% ======== what data to take?
% dotrain = 1; % during training?
singleRendOnly=1; % only take data up to first nontarg after the referen
% if 0, then takes all rends


% ====== what data to use to compute density?
densitymethod = 'refperiod';
% refperiod: will use reference period
% entirebin: will use time from rendtion to maxtime_fromref (see below)
% beforefirstnontarg = will get number that is >= ref time and < time of first target post-reference.
% [if singleRendOnly is 1, then will by default use beforefirstnontarg]


% ============ method for decideing hi and lo density trials
cutoffmethod = 'medianoverall';
% medianslice: for each nontarg bin, finds median for targ
% medianoverall: overall median of targ
% medianslice_rand: gets median by first addaing random jitter. useful for
% i) when even num bins. if only one y bin, then forces that to be in low
% density bin. if odd num bins, then will split favoring putting
% more data into low density category. 
% targffdev: then splits based on mean FF dev of target syllables 


% ======== summary plot: what time period to plot (locked to ref period)
mintime_fromref = 5;
maxtime_fromref = 30; % minutes


% =============== DEFAULTS
% singleRendOnly=0; % only take data up to first nontarg after the referen
% densitymethod = 'refperiod';
% cutoffmethod = 'medianslice';
% mintime_fromref = 5;
% maxtime_fromref = 30; % minutes
% minrends_inbin = 20;

% ================ RUN
DATBYREND = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S2(TrialStruct, ...
    singleRendOnly, densitymethod, cutoffmethod, mintime_fromref, ...
    maxtime_fromref, ignoreLMANexpt);


assert(all(isnan(DATBYREND.IsDurTrain) == isnan(DATBYREND.Density_isHigh)), 'then nan likely matches up');
assert(all(isnan(DATBYREND.Density_targ) == isnan(DATBYREND.Density_nontarg)), 'then nan likely matches up');



%% =========== [PLOT DIAGNOSTIC] FOR INDIVIDUAL EXPERIMENTS

if (1)
    birdtoplot = 10;
    expttotplot =2;
    
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


%% ========= [COLLECT FOR PLOTTING] SEPARATE BY HIGH AND LO DENS

% ============================================ TRAIN
dotrain = 1; % 1 = train, 0 = base
dosametype = 1; % 1=same; 0: diff type;
onlyifSigLearn = 0; % 1=yes, 0=don't care

% ----------------- % method for equalizing density
densitymethod = 'default'; % then whatever was used to classify has high and low density
% FFwindow; % then from 0(inclusive) to whatever window used to get FFdev

% ----------- min num trials in this bin...
minrends_inbin = 20;
minsongs_inbin = 1;

% ============================ TRAINING
[NtargCell, NnontargCell, FFbinnedCell, TbinnedCell, FFsinglebinMat, ...
    Nratio_hilo_targ, NumDatPerRend] = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S3(DATBYREND, ...
    dotrain, dosametype, onlyifSigLearn, densitymethod, xedges, mintime_fromref, ...
    maxtime_fromref, minrends_inbin, minsongs_inbin);


% ============================================== BASE
dotrain = 0; % 1 = train, 0 = base
[NtargCell_base, NnontargCell_base, FFbinnedCell_base, TbinnedCell_base, ...
    FFsinglebinMat_base, Nratio_hilo_targ_base, NumDatPerRend_base] = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S3(DATBYREND, ...
    dotrain, dosametype, onlyifSigLearn, densitymethod, xedges, mintime_fromref, ...
    maxtime_fromref, minrends_inbin, minsongs_inbin);



%% ==================== distribution on datapoints per rendition
lt_figure; hold on;

% ======= scatter, high vs. low dens
lt_subplot(2,2,1); hold on;
title('num dat (in dat bin) per rend');
xlabel('lo dens trials');
ylabel('hi dens trials');
plot(NumDatPerRend(:,1), NumDatPerRend(:,2), 'ok');
lt_plot_makesquare_plot45line(gca, 'b', -1);

lt_subplot(2,2,2); hold on;
title('log2');
xlabel('lo dens trials');
ylabel('hi dens trials');
plot(log2(NumDatPerRend(:,1)), log2(NumDatPerRend(:,2)), 'ok');
lt_plot_makesquare_plot45line(gca, 'b');

lt_subplot(2,2,3); hold on;
y = mean(NumDatPerRend,2);
lt_plot_histogram(y(~isnan(y)));
title('mean of hi and lo dens');
xlabel('num dat per rend');

lt_subplot(2,2,4); hold on;
y = mean(NumDatPerRend,2);
lt_plot_histogram(log2(y(~isnan(y))));
title('mean of hi and lo dens');
xlabel('log2(num dat per rend)');



%% ===================== [PLOT] WHICH INDS TO PLOT?
onlyPlotIfDensDiff = 1; % makes sure that targ syl has higher N for hi density relative to lo density.
% only check targ, makes sure ratio>1;
onlyPlotIfDenseLabeling = 0; % decides if dense by looking at number of datapoints
% per rendition.
labelthresh = 1; % in log2 units, only matters if onlyPlotIfDenseLabeling=1

% ================== 1) ONLY PLOT IF HAVE DATA FOR BOTH LO AND HIGH DENS
tmp = [cellfun('isempty', FFbinnedCell(:,1)) cellfun('isempty', FFbinnedCell(:,2))];

if onlyPlotIfDensDiff==1
    indstoplot = ~any(tmp') & Nratio_hilo_targ'>1;
else
    indstoplot = ~any(tmp');
end
if onlyPlotIfDenseLabeling==1
     indstoplot = indstoplot & log2(mean(NumDatPerRend,2))'>labelthresh;
end
indstoplot = find(indstoplot);


% ====================== BASE
tmp = [cellfun('isempty', FFbinnedCell_base(:,1)) cellfun('isempty', FFbinnedCell_base(:,2))];

if onlyPlotIfDensDiff==1
    indstoplot_base = ~any(tmp') & Nratio_hilo_targ_base'>1;
else
    indstoplot_base = ~any(tmp');
end
if onlyPlotIfDenseLabeling==1
     indstoplot_base = indstoplot_base & log2(mean(NumDatPerRend_base,2))'>labelthresh;
end
indstoplot_base = find(indstoplot_base);

% 
% 
% tmp = [cellfun('isempty', FFbinnedCell_base(:,1)) cellfun('isempty', FFbinnedCell_base(:,2))];
% 
% if onlyPlotIfDensDiff==1
%     indstoplot_base = find(~any(tmp') & Nratio_hilo_targ_base'>1);
% else
%     indstoplot_base = find(~any(tmp'));
% end



%% ==================== [PLOT]
plotTrainAndBase =0;


if plotTrainAndBase == 0
    % ---------- Then only training
    lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S4(NtargCell, NnontargCell, ...
        FFbinnedCell, TbinnedCell, FFsinglebinMat, indstoplot)
elseif plotTrainAndBase==1
    % ---------- Then base and training, matched syllables
    
    % ----- 1) collect matched syllables
    indstoplot = intersect(indstoplot, indstoplot_base);
    
    
    % -------------------- 2) PLOT TRAIN
    lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S4(NtargCell, NnontargCell, ...
        FFbinnedCell, TbinnedCell, FFsinglebinMat, indstoplot)
    
    % ---------------------- 3) PLOT BASELINE
    lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S4(NtargCell_base, NnontargCell_base, ...
        FFbinnedCell_base, TbinnedCell_base, FFsinglebinMat_base, indstoplot)
    
end





%% ################################################################
%% ############################ NOT SEPARATING BY TARGET DENSITY



%% ====================== extract binned data
% ============================================ TRAIN
dosametype = 1; % 1=same; 0: diff type;
onlyifSigLearn = 0; % 1=yes, 0=don't care

% ----------------- % method for equalizing density
densitymethod = 'default'; % then whatever was used to classify has high and low density
% FFwindow; % then from 0(inclusive) to whatever window used to get FFdev

% ----------- min num trials in this bin...
minrends_inbin = 20;
minsongs_inbin = 2;

% ============================= TRAINING
dotrain = 1; % 1 = train, 0 = base
[~, ~, FFbinnedCell, TbinnedCell, ~, ...
    ~, NumDatPerRend, BirdNumMat] = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_noDens(DATBYREND, ...
    dotrain, dosametype, onlyifSigLearn, densitymethod, xedges, mintime_fromref, ...
    maxtime_fromref, minrends_inbin, minsongs_inbin);


% ============================= BASELINE
dotrain = 0;
[~, ~, FFbinnedCell_BASE, TbinnedCell_BASE, ~, ...
    ~, NumDatPerRend_BASE, BirdNumMat_BASE] = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_noDens(DATBYREND, ...
    dotrain, dosametype, onlyifSigLearn, densitymethod, xedges, mintime_fromref, ...
    maxtime_fromref, minrends_inbin, minsongs_inbin);

%% =================== distrubution of labeling densityl;
lt_figure; hold on;

% === distribution
lt_subplot(2,2,1); hold on;
indstmp = ~isnan(NumDatPerRend);
lt_plot_histogram(NumDatPerRend(indstmp));


lt_subplot(2,2,2); hold on;
xlabel('birdnum');
ylabel('renditions per locked syl');

indstmp = ~isnan(NumDatPerRend);
x = BirdNumMat(indstmp);
y = NumDatPerRend(indstmp);
plot(x,y, 'ok');


%% ===================== [PLOT] WHICH INDS TO PLOT?
onlyPlotIfDenseLabeling = 0; % decides if dense by looking at number of datapoints
% per rendition.
labelthresh = 0; % in log2 units, only matters if onlyPlotIfDenseLabeling=1
birdstokeep = [13:17]; % leave empty to take all birds

% ================== 1) ONLY PLOT IF HAVE DATA FOR BOTH LO AND HIGH DENS
% -------- 1) only plot if have data
indstoplot = ~cellfun('isempty', FFbinnedCell);

% -------- 2) only keep if dense labeling?
if onlyPlotIfDenseLabeling==1
     indstoplot = indstoplot & log2(NumDatPerRend)>labelthresh;
end

% --------- 3) birds to keep
if ~isempty(birdstokeep)
    indstoplot = indstoplot & ismember(BirdNumMat, birdstokeep);
end


% --------- FINAL, convert to numerical indeices
indstoplot = find(indstoplot)';



%% =================== [PLOT] SHOW SAMPLE SIZE

lt_figure; hold on;

% ############################ BINNED FF DEV 
lt_subplot(3,2,1); hold on;

for i=indstoplot

    x = TbinnedCell{i};
    y = FFbinnedCell{i};
    plot(x,y, '-x', 'Color', 'b');
end
lt_plot_zeroline;
lt_plot_zeroline_vert;
xlim([-30 30]);

% ############################ BINNED FF DEV (BASELINE)
lt_subplot(3,2,2); hold on;

for i=indstoplot

    x = TbinnedCell_BASE{i};
    y = FFbinnedCell_BASE{i};
    plot(x,y, '-x', 'Color', 'b');
end
lt_plot_zeroline;
lt_plot_zeroline_vert;
xlim([-30 30]);


% ############################## IN EACH BIN, COMPARE DAT VS. BASE
lt_subplot(3,2,3); hold on;

for i=indstoplot

    x = TbinnedCell{i};
    y = FFbinnedCell{i};
    
    xbase = TbinnedCell_BASE{i};
    ybase = FFbinnedCell_BASE{i};
    
%     X = [xbase'-2 x'+2];
%     Y = [ybase y];

    for cc = xcenters
       
        X = [cc-2 cc+2];
        Y = [ybase(xbase==cc) y(x==cc)];
        
        if length(Y)==2
            plot(X, Y, '-');
        end
        
    end
%     for ii=1:size(X,1)
%         plot(X(ii,:), Y(ii,:), '-');
%     end
end
lt_plot_zeroline;
lt_plot_zeroline_vert;
xlim([-30 30]);



% ############################ BINNED FF DEV [mean across syls] 
lt_subplot(3,2,4); hold on;
title('mean');
X = [TbinnedCell{indstoplot}];
Y = cell2mat(FFbinnedCell(indstoplot))';

[ymean, ysem] = grpstats(Y, X, {'mean', 'sem'});
xmean = unique(X);
lt_plot(xmean, ymean, {'Errors', ysem, 'Color', 'k'});


% % ###################################### TAKE ALL RENDS IN ONE BIN
% lt_subplot(3,2,6); hold on;
% 
% Y = FFsinglebinMat(indstoplot,:);
% X = [1 2];
% plot(X, Y', '-', 'Color', [0.7 0.7 0.7]);
% % -- means
% lt_plot(X+0.1, mean(Y), {'Errors', lt_sem(Y)});
% 
% xlim([0 3]);
% lt_plot_zeroline;
% % --- signifnicace
% [~, p] = ttest(Y(:,1), Y(:,2));
% lt_plot_pvalue(p, 'ttest', 1);
% 
% 




