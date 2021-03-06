function lt_seq_dep_pitch_ACROSSBIRDS_TbyT_TcourseBACK(TrialStruct, ParamsTrial, ...
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
flanktime_targ = 2.5; % minutes (will do +/- this number). This only used if
% useFFdevWithinSylTiming = 0; 
nonsingonly = 1; % if 1, then makes sure the targ syl is not sung again between
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


timewindow = 1; % in hours
binsize = 20; % renditions for running mean
%%

Numbirds = length(TrialStruct.birds);

%% CALCULATE FF DEVIATIONS

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
            isLMAN = TrialStruct.birds(i).exptnum(ii).LMANinactivated;
            if isLMAN==1
                disp(['[is LMAN] skipping ' TrialStruct.birds(i).exptnum(ii).exptname]);
                continue
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
            istarg = TrialStruct.birds(i).exptnum(ii).sylnum(ss).INFO_istarget;
            issame = TrialStruct.birds(i).exptnum(ii).sylnum(ss).INFO_similar;
            sylname = TrialStruct.birds(i).exptnum(ii).sylnum(ss).syl;
            
            % --------------- lines for base mean and 1std
            basedays = TrialStruct.birds(i).exptnum(ii).BaseDays;
            wndays = TrialStruct.birds(i).exptnum(ii).WNDays;
            indsbase = t<basedays(end)+1;
            ffmean_base = mean(ff(indsbase));
            ffstd_base = std(ff(indsbase));
            
            
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

i=3;
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
    
    % ============ plot deviations
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(titstr);
    plot(t, ffdev, 's', 'Color', pcol);
    lt_plot_zeroline;
    
    % -----------------
    linkaxes(hsplots, 'x');
    
    
end
end


%% ============== [COLLECT] "targ locked deviation of nontarg"
% ---- for each rendition of targ, plot change in nontarg as function of
% teporal lag from the target.
% ---- make separate plots for targ ff that is posistive (i.e. laerning)
% vs. negative.

tflankplot = 1; % for rend-locking, how much flanking time to plot (hrs)
twind_plot = [-tflankplot tflankplot]; % for rend-locking, how much flanking time to plot (hrs)

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
            isLMAN = TrialStruct.birds(i).exptnum(ii).LMANinactivated;
            if isLMAN==1
                disp(['[is LMAN] skipping ' TrialStruct.birds(i).exptnum(ii).exptname]);
                continue
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
            basedaylast = TrialStruct.birds(i).exptnum(ii).BaseDays(end);
            inds_base = t_targ<basedaylast+1;
            
            
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
        
%         
%         % ============ combined plot
%         [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%         hsplots = [hsplots hsplot];
%         title('all nontargs combined');
%         
%         % --- targ positive learn
%         indtmp = fftarg_all>0;
%         pcol = 'k';
%         
%         x = tvals_all(indtmp);
%         y = ffvals_all(indtmp);
%         plot(x, y, 'x', 'Color', pcol);
%         
%         % ---- plot running mean
%         [~, indsort] = sort(x); % sort
%         x = x(indsort);
%         y = y(indsort);
%         xrun = lt_running_stats(x, binsize);
%         yrun = lt_running_stats(y, binsize);
%         
%         shadedErrorBar(xrun.Median, yrun.Mean, yrun.SEM, {'Color', pcol}, 1);
%         
%         % --- targ negative learn
%         indtmp = fftarg_all<=0;
%         pcol = 'r';
%         
%         x = tvals_all(indtmp);
%         y = ffvals_all(indtmp);
%         plot(x, y, 'x', 'Color', pcol);
%         
%         % ---- plot running mean
%         [~, indsort] = sort(x); % sort
%         x = x(indsort);
%         y = y(indsort);
%         xrun = lt_running_stats(x, binsize);
%         yrun = lt_running_stats(y, binsize);
%         
%         shadedErrorBar(xrun.Median, yrun.Mean, yrun.SEM, {'Color', pcol}, 1);
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

for i=1:Maxbirds
    
    for ii=1:Maxexpts
        
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

            end
            
        end
        
    end
end

%% [PLOT] GRAND MEAN, ALL RENDS ################################

lt_figure; hold on;
% xedges = -65:10:65; % minutes
xedges = -62.5:5:62.5; % minutes
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
        
        [ymean, ysem] = grpstats(y, xbins, {'mean', 'sem'});
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
        
        [ymean, ysem] = grpstats(y, xbins, {'mean', 'sem'});
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
        
        [ymean, ysem] = grpstats(y, xbins, {'mean', 'sem'});
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
        
        [ymean, ysem] = grpstats(y, xbins, {'mean', 'sem'});
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
        
        [ymean, ysem] = grpstats(y, xbins, {'mean', 'sem'});
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
        
        [ymean, ysem] = grpstats(y, xbins, {'mean', 'sem'});
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
        indtmp = isbase_all==plotbase & birdnum_all==i & exptnum_all==ii;
                
        x = tvals_all(indtmp);
        y = ffvals_all(indtmp);
        
        % ---------------- slide into bins
        xbins = discretize(x, xedges);
        assert(~any(isnan(xbins)), 'edges not wide enough ...');
        
        [ymean, ysem] = grpstats(y, xbins, {'mean', 'sem'});
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
        
        [ymean, ysem] = grpstats(y, xbins, {'mean', 'sem'});
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





