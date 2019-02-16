function [BirdExptIncluded] = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Slopes(TrialStruct, ParamsTrial, ...
    ignoreLMANexpt, threshOnSametype, scalemethod, combineSylsInExpt, ...
    onlyifhaveAllSylTypes, throwoutnan, plotEachExptRaw)
%% 10/6/18 -

% BirdExptIncluded = []; 'i-ii' where i and ii are numbers and concat into
% string.

%% lt 5/25/18 -
% align all experiments by day and plot fitted slopes
normmethod = 'base_overall';
% base_edges: separately for onset and offset (accounts for circadian, using baseline regression);
% base_overall: one value, mean baseline across rends.

% scalemethod = '';
% % lastdaymean: mean of last day becomes 1

% combineSylsInExpt=1; % if 0, datapt is syls, if 1, datapt is experiments
onlyKeepIfHaveNontarg=1; % throw out experiment if it doesn't have a sametype.

useregression=0; % then gets day on and off by fitting each day separately
% otherwise get first and last N renditions
Nrend = 10;

% ======== throw out same-type syls that don't show any learning?
% threshOnSametype = 1; % threshold is 0 for mean of mrinig and night of last day.

%%


%%

Numbirds = length(TrialStruct.birds);

%% COLLECT ALL DATA

DATSTRUCT.AllBirdnum =[];
DATSTRUCT.AllExptnum =[];
DATSTRUCT.AllSylname ={};

DATSTRUCT.AllIsTarg = [];
DATSTRUCT.AllIsSame = [];

DATSTRUCT.AllDayList = [];
DATSTRUCT.AllTedges = [];
DATSTRUCT.AllFFedges = [];

DATSTRUCT.AllBasedays = [];

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
        birdname =  TrialStruct.birds(i).birdname;
        exptname = TrialStruct.birds(i).exptnum(ii).exptname;
        targlearndir = TrialStruct.birds(i).exptnum(ii).targlearndir;
        
        % =========== collect syls for this experiment
        ffedges_allsyls =[];
        tedges_allsyls = [];
        istarg_allsyls = [];
        issame_allsyls =[];
        sylnames_allsyls = {};
        for ss =1:Numsyls
            
            % ============== subplot for this syl
            t = TrialStruct.birds(i).exptnum(ii).sylnum(ss).Tvals;
            ff = TrialStruct.birds(i).exptnum(ii).sylnum(ss).FFvals;
            istarg = TrialStruct.birds(i).exptnum(ii).sylnum(ss).INFO_istarget;
            issame = TrialStruct.birds(i).exptnum(ii).sylnum(ss).INFO_similar;
            sylname = TrialStruct.birds(i).exptnum(ii).sylnum(ss).syl;
            
            if any(isnan(ff))
                keyboard
            end
            
            % --------------- lines for base mean and 1std
            basedays = TrialStruct.birds(i).exptnum(ii).BaseDays;
            wndays = TrialStruct.birds(i).exptnum(ii).WNDays;
            indsbase = t<basedays(end)+1;
            ffmean_base = mean(ff(indsbase));
            ffstd_base = std(ff(indsbase));
            
            
            % ======================== FIT SIMPLE REGRESSION LINES FOR EACH
            % DAY
            daylist = [basedays wndays];
            Tedges_byday = nan(1, 2*length(daylist)); % COLUMNS, day, split to ons and off (e.g. [day1on, day1off, day2on, day2off...];
            FFedges_byday = nan(1, 2*length(daylist));
            ff_fit = [];
            for dd = 1:length(daylist)
                daythis = daylist(dd);
                indthis = floor(t)==daythis;
                
                % ------------- skip if not eniogh rends
                if sum(indthis)<1.5*Nrend
                    disp('NOT ENOGU REND - SKIPPING');
                    continue
                end
                
                if ~any(indthis)
                    % -- then no data for this day...
                    tthis = nan;
                    ff_fit = nan;
                else
                    tthis = t(indthis);
                    ffthis = ff(indthis);
                    
                    % -- sort
                    [~, indtmp] = sort(tthis);
                    tthis = tthis(indtmp);
                    ffthis = ffthis(indtmp);
                    
                    if useregression ==1
                        % -- fit regression
                        [b] =lt_regress(ffthis, tthis, 0);
                        
                        % -- collect FF at beginning and end of each day (from
                        % regression fit)
                        ff_fit = b(1) + b(2)*(tthis);
                    else
                        ff_fit(1) = mean(ffthis(1:Nrend));
                        ff_fit(2) = mean(ffthis(end-Nrend+1:end));
                    end
                end
                
                xtmp = [dd*2-1 dd*2];
                Tedges_byday(xtmp) = [tthis(1) tthis(end)];
                FFedges_byday(xtmp) = [ff_fit(1) ff_fit(end)];
            end
            
            
            % ================= NORMALIZE
            if strcmp(normmethod, 'base_overall')
                % -------- 1) subtract one overall baseline mean
                FFedges_byday = FFedges_byday - ffmean_base;
                FFedges_byday = FFedges_byday*targlearndir;
                ff_minusbase = (ff - ffmean_base)*targlearndir;
            elseif strcmp(normmethod, 'base_edges')
                % ---------- 2) subtract baseline mean separately for start
                % and end of day (account for circadian rhythm)
                numbasedays = length(basedays);
                % --- get mean baselines for early and late day
                baseON = mean(FFedges_byday(1:2:numbasedays*2));
                baseOFF = mean(FFedges_byday(2:2:numbasedays*2));
                % --- for each day subtract baseline for early and late
                FFedges_byday(1:2:end) = FFedges_byday(1:2:end) - baseON;
                FFedges_byday(2:2:end) = FFedges_byday(2:2:end) - baseOFF;
                % --- make up direction of learning
                FFedges_byday = FFedges_byday*targlearndir;
                ff_minusbase = [];
            end
            
            
            % =========== save edge values for this expt
            TrialStruct.birds(i).exptnum(ii).sylnum(ss).FFedges_byday = FFedges_byday;
            TrialStruct.birds(i).exptnum(ii).sylnum(ss).ff_minusbase = ff_minusbase;
            
            %             % =================== SCALE
            %             if strcmp(scalemethod, 'lastdaymean')
            %                 assert(1==2, 'PRIOBLME: need to have nontarg scaled to target...');
            %                 assert(strcmp(normmethod, 'base_overall') | strcmp(normmethod, 'base_edges'), ...
            %                     'to scale must norm relative to baseline');
            %                 % --------- scale so that end of last day is 1 and baseline
            %                 % is 0
            %                 maxff = mean(FFedges_byday(end-1:end));
            %                 FFedges_byday = FFedges_byday./abs(maxff);
            %                 % --- take abs of max to make sure the sign of the last day still holds (i.e.
            %                 % did not take negative laerning and then flip it)
            %
            %             end
            %
            
            % ===================== CLEANING UP EXPT
            
            % =================== THROW OUT EXPT IF NO LEARNING?
            if threshOnSametype==1 & istarg==0 & issame==1
                if mean(FFedges_byday(end-1:end))<0
                    continue
                end
            end
            if throwoutnan==1
                if any(isnan(FFedges_byday(:)))
                    continue
                end
            end
            
            % ================= OUTPUT
            ffedges_allsyls =[ffedges_allsyls; FFedges_byday];
            tedges_allsyls = [tedges_allsyls; Tedges_byday];
            istarg_allsyls = [istarg_allsyls; istarg];
            issame_allsyls =[issame_allsyls; issame];
            sylnames_allsyls = [sylnames_allsyls; sylname];
            
        end
        
        % ==================== OUTPUTS
        indstoremove = []; % the ones that I am combining
        % ---------- combine within syl types?
        if combineSylsInExpt==1
            
            % -- same
            indstmp = find(istarg_allsyls==0 & issame_allsyls==1);
            if ~isempty(indstmp)
                ffedges_allsyls = [ffedges_allsyls; nanmean(ffedges_allsyls(indstmp,:),1)];
                tedges_allsyls = [tedges_allsyls; nanmean(tedges_allsyls(indstmp,:),1)];
                istarg_allsyls = [istarg_allsyls; 0];
                issame_allsyls = [issame_allsyls; 1];
                
                indstoremove = [indstoremove; indstmp];
            end
            
            % --- DIFF
            indstmp = find(istarg_allsyls==0 & issame_allsyls==0);
            if ~isempty(indstmp)
                
                ffedges_allsyls = [ffedges_allsyls; nanmean(ffedges_allsyls(indstmp,:),1)];
                tedges_allsyls = [tedges_allsyls; nanmean(tedges_allsyls(indstmp,:),1)];
                istarg_allsyls = [istarg_allsyls; 0];
                issame_allsyls = [issame_allsyls; 0];
                
                indstoremove = [indstoremove; indstmp];
            end
            
            % ----- remove syls
            ffedges_allsyls(indstoremove,:) = [];
            tedges_allsyls(indstoremove,:) = [];
            istarg_allsyls(indstoremove,:) = [];
            issame_allsyls(indstoremove,:) = [];
            sylnames_allsyls = [];
            
            if any(isnan(ffedges_allsyls(:)))
                keyboard
            end
        end
        
        % ======== put those into output struct
        if onlyKeepIfHaveNontarg==1
            % --- check, does sametype exist?
            indtmp = istarg_allsyls==0 & issame_allsyls ==1;
            if ~any(indtmp)
                continue
            end
            
            tmp = ~isnan(ffedges_allsyls(indtmp,:));
            if ~any(tmp(:))
                continue
            end
            
        end
        
        
        % ========= need all syl types?
        if onlyifhaveAllSylTypes==1
            if ~(any(istarg_allsyls==1) & any(istarg_allsyls==0 & issame_allsyls==1) ...
                    & any(istarg_allsyls==0 & issame_allsyls==0))
                disp('SKIP, dosnt have all syltypes');
                continue
            end
        end
        
        
        DATSTRUCT.AllBirdnum = [DATSTRUCT.AllBirdnum; i*ones(size(ffedges_allsyls,1),1)];
        DATSTRUCT.AllExptnum =[DATSTRUCT.AllExptnum; ii*ones(size(ffedges_allsyls,1),1)];
        %         DATSTRUCT.AllDayList = [DATSTRUCT.AllDayList; daylist];
        %         DATSTRUCT.AllBasedays = [DATSTRUCT.AllBasedays; basedays];
        
        DATSTRUCT.AllSylname = [DATSTRUCT.AllSylname; sylnames_allsyls];
        DATSTRUCT.AllIsTarg = [DATSTRUCT.AllIsTarg; istarg_allsyls];
        DATSTRUCT.AllIsSame = [DATSTRUCT.AllIsSame; issame_allsyls];
        
        
        DATSTRUCT.AllTedges = [DATSTRUCT.AllTedges; tedges_allsyls];
        DATSTRUCT.AllFFedges = [DATSTRUCT.AllFFedges; ffedges_allsyls];
        
        
    end
end


maxbirds = max(DATSTRUCT.AllBirdnum);
maxexpts = max(DATSTRUCT.AllExptnum);


%% ===== sample size

disp(['N(birds) = ' num2str(max(DATSTRUCT.AllBirdnum))]);
disp(['N(expt) = ' num2str(max(lt_tools_grp2idx({DATSTRUCT.AllBirdnum, DATSTRUCT.AllExptnum})))]);

% --- collect unqiue expts
ind_expts = lt_tools_grp2idx({DATSTRUCT.AllBirdnum, DATSTRUCT.AllExptnum});
col1 = grpstats(DATSTRUCT.AllBirdnum, ind_expts, {'mean'});
col2 = grpstats(DATSTRUCT.AllExptnum, ind_expts, {'mean'});
BirdExptIncluded  = [col1 col2];


%% ================= PLOT EACH EXPERIEMNT, OPVERLAYING CALCULATED EDGE VALUES.
if plotEachExptRaw==1
    disp('ONLY WORKS WELL IF USE OVERALL BASE AS NORM...');
    
    pcol = {};
    pcol{2,2} = 'k'; %(targ, smae)
    pcol{1,2} = 'b';
    pcol{1,1} = 'r';
    
    for i=1:maxbirds
        for ii=1:maxexpts
            
            indsthis = DATSTRUCT.AllBirdnum==i & DATSTRUCT.AllExptnum==ii;
            if ~any(indsthis)
                continue
            end
            
            bname = TrialStruct.birds(i).birdname;
            ename = TrialStruct.birds(i).exptnum(ii).exptname;
            
            % ======= plot this experiment
            figcount=1;
            subplotrows=6;
            subplotcols=2;
            fignums_alreadyused=[];
            hfigs=[];
            hsplots = [];
            
            % ====== go thru all syls
            numsyls = length(TrialStruct.birds(i).exptnum(ii).sylnum);
            for ss=1:numsyls
                
                sylname = TrialStruct.birds(i).exptnum(ii).sylnum(ss).syl;
                t = TrialStruct.birds(i).exptnum(ii).sylnum(ss).Tvals;
                ff_minusbase = TrialStruct.birds(i).exptnum(ii).sylnum(ss).ff_minusbase;
                ffedges = TrialStruct.birds(i).exptnum(ii).sylnum(ss).FFedges_byday;
                tedges = reshape([unique(floor(t))'+0.2917; unique(floor(t))'+0.875], length(ffedges), [])';
                assert(length(tedges)==length(ffedges));
                
                % is targ?
                istarg = TrialStruct.birds(i).exptnum(ii).sylnum(ss).INFO_istarget;
                issame = TrialStruct.birds(i).exptnum(ii).sylnum(ss).INFO_similar;
                pcthis = pcol{istarg+1, issame+1};
                
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
                title(sylname);
                if ss==1
                    ylabel([bname '-' ename]);
                end
                plot(t, ff_minusbase, 'x', 'Color', pcthis);
                lt_plot(tedges, ffedges);
                axis tight;
                lt_plot_zeroline;
            end
            
            % ====== plot summary for this expt
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            % hsplots = [hsplots hsplot];
            title('summary - actual analysis dat');
            % -- targ
            indtmp = DATSTRUCT.AllBirdnum==i & DATSTRUCT.AllExptnum==ii & DATSTRUCT.AllIsTarg==1;
            ff = DATSTRUCT.AllFFedges(indtmp,:);
            for j=1:length(ff)/2
                plot(j*2-1:j*2, ff(:, j*2-1:j*2)', '-xk');
            end
            
            % -- same
            indtmp = DATSTRUCT.AllBirdnum==i & DATSTRUCT.AllExptnum==ii & DATSTRUCT.AllIsTarg==0 & ...
                DATSTRUCT.AllIsSame==1;
            ff = DATSTRUCT.AllFFedges(indtmp,:);
            for j=1:length(ff)/2
                plot(j*2-1:j*2, ff(:, j*2-1:j*2)', '-xb');
            end
            % -- diff
            indtmp = DATSTRUCT.AllBirdnum==i & DATSTRUCT.AllExptnum==ii & DATSTRUCT.AllIsTarg==0 & ...
                DATSTRUCT.AllIsSame==0;
            ff = DATSTRUCT.AllFFedges(indtmp,:);
            for j=1:length(ff)/2
                plot(j*2-1:j*2, ff(:, j*2-1:j*2)', '-xr');
            end
            lt_plot_zeroline;
            
            % ----------------
            linkaxes(hsplots, 'xy');
        end
    end
end


%% ================= scale so that last day is 1 for all expts?

if strcmp(scalemethod, 'lastdaymean')
    
    assert(strcmp(normmethod, 'base_overall') | strcmp(normmethod, 'base_edges'), ...
        'to scale must norm relative to baseline');
    
    % ============ go thru all expts and do separately for each
    for i=1:maxbirds
        for ii=1:maxexpts
            
            inds = DATSTRUCT.AllBirdnum==i & DATSTRUCT.AllExptnum==ii;
            if ~any(inds)
                continue
            end
            
            % =========== get max ff for target
            indtarg = DATSTRUCT.AllBirdnum==i & DATSTRUCT.AllExptnum==ii & ...
                DATSTRUCT.AllIsTarg==1;
            assert(sum(indtarg)==1, 'why not exactly one targ?');
            
            
            maxff = mean(DATSTRUCT.AllFFedges(indtarg,end-1:end));
            DATSTRUCT.AllFFedges(inds, :) = DATSTRUCT.AllFFedges(inds, :)./abs(maxff);
            % --- take abs of max to make sure the sign of the last day still holds (i.e.
            % did not take negative laerning and then flip it)
            
        end
    end
    
end

%% ================= [PLOT] --- means across experiments

lt_figure; hold on;

% ================= TARGET
lt_subplot(3,1,1); hold on;
title('target');
inds = DATSTRUCT.AllIsTarg==1 & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

ffedges = DATSTRUCT.AllFFedges(inds,:);
x1 = 0.3:1:(size(ffedges,2)/2);
x2 = 0.7:1:(size(ffedges,2)/2);
x = [x1; x2];
x = x(:)';

% ==== plot, onsets and offsets combined
plot(x, ffedges, '-ok');
% -- plot means
ffmeans = mean(ffedges,1);
ffsem = lt_sem(ffedges);
lt_plot(x+0.1, ffmeans, {'Errors', ffsem, 'Color', 'r'});
% ---
lt_plot_zeroline;
lt_plot_text(0, max(ffmeans), ['N=' num2str(sum(inds))], 'b')

% ================= SAME TYPE
lt_subplot(3,1,2); hold on;
title('Same (syl = datapoint)');
inds = DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsSame==1 ...
    & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

ffedges = DATSTRUCT.AllFFedges(inds,:);
x1 = 0.3:1:(size(ffedges,2)/2);
x2 = 0.7:1:(size(ffedges,2)/2);
x = [x1; x2];
x = x(:)';

% ==== plot, onsets and offsets combined
plot(x, ffedges, '-ok');
% -- plot means
ffmeans = mean(ffedges,1);
ffsem = lt_sem(ffedges);
lt_plot(x+0.1, ffmeans, {'Errors', ffsem, 'Color', 'r'});
% ---
lt_plot_zeroline;
lt_plot_text(0, max(ffmeans), ['N=' num2str(sum(inds))], 'b')


% ================= DIFF TYPE
lt_subplot(3,1,3); hold on;
title('Diff (syl = datapoint)');
inds = DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsSame==0 ...
    & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

ffedges = DATSTRUCT.AllFFedges(inds,:);
x1 = 0.3:1:(size(ffedges,2)/2);
x2 = 0.7:1:(size(ffedges,2)/2);
x = [x1; x2];
x = x(:)';

% ==== plot, onsets and offsets combined
plot(x, ffedges, '-ok');
% -- plot means
ffmeans = mean(ffedges,1);
ffsem = lt_sem(ffedges);
lt_plot(x+0.1, ffmeans, {'Errors', ffsem, 'Color', 'r'});
% ---
lt_plot_zeroline;
lt_plot_text(0, max(ffmeans), ['N=' num2str(sum(inds))], 'b')



%% ================ [PLOT] Single metric of night and morning generalization for each nontarg
% i.e. night and the folliwng mornign.

numbasedays  = -ParamsTrial.DayWindow(1);
numdays = -ParamsTrial.DayWindow(1) + ParamsTrial.DayWindow(2);
bins_night = numbasedays*2+2:2:(numdays-1)*2;
bins_morning = bins_night+1;


disp('NOTE: this is best done WITHOUT combined across syls ...');
% FFtargAll = [];
% FFsumAll = [];
Yall_same = [];
Yall_diff = [];
Yall_targ = [];

Yall_same_gen = [];
Yall_diff_gen = [];

for i=1:maxbirds
    for ii=1:maxexpts
        
        indtarg = DATSTRUCT.AllBirdnum==i & DATSTRUCT.AllExptnum==ii & DATSTRUCT.AllIsTarg==1;
        indrest = DATSTRUCT.AllBirdnum==i & DATSTRUCT.AllExptnum==ii & DATSTRUCT.AllIsTarg==0;
        indsame = DATSTRUCT.AllBirdnum==i & DATSTRUCT.AllExptnum==ii & DATSTRUCT.AllIsTarg==0 ...
            & DATSTRUCT.AllIsSame==1;
        inddiff= DATSTRUCT.AllBirdnum==i & DATSTRUCT.AllExptnum==ii & DATSTRUCT.AllIsTarg==0 ...
            & DATSTRUCT.AllIsSame==0;
        
        if ~any(indtarg)
            continue
        end
        
        fftarg = DATSTRUCT.AllFFedges(indtarg,:);
        ffrest = DATSTRUCT.AllFFedges(indrest,:);
        ffsame = DATSTRUCT.AllFFedges(indsame,:);
        ffdiff = DATSTRUCT.AllFFedges(inddiff,:);
        
        % ---- for each syl, calculate mean learning at night and morning
        
        % NIGHT BIn
        ytarg = [mean(fftarg(:, bins_night),2) mean(fftarg(:, bins_morning),2)];
        ysame = [mean(ffsame(:, bins_night),2) mean(ffsame(:, bins_morning),2)];
        ydiff = [mean(ffdiff(:, bins_night),2) mean(ffdiff(:, bins_morning),2)];
        
        % ================= OUTPUT
        Yall_targ = [Yall_targ; ytarg];
        Yall_same = [Yall_same; ysame];
        Yall_diff = [Yall_diff; ydiff];
        
        Yall_same_gen = [Yall_same_gen; ysame./ytarg];
        Yall_diff_gen = [Yall_diff_gen; ydiff./ytarg];
        
        % === calculate specificity index
        %         ffall = sum([fftarg; ffrest],1);
        
        %         FFtargAll = [FFtargAll; fftarg];
        %         FFsumAll = [FFsumAll; ffall];
    end
end

% ==================== PLOT
lt_figure;

% ------
lt_subplot(3,1,1); hold on;
% ylabel('ff, mean across days');
title('(ff) each syl one val (mean across days)');
ylabel('NIGHT(day n) - MORNING(day n+1)');
xlabel('TARG - SAME - SAME(learn<0) - SAME(learn>0) - DIFF');

% - targ
x = [1 2];
y = Yall_targ;
plot(x, y', '-', 'Color', 'r');
lt_plot_bar(x, mean(y,1), {'Errors', lt_sem(y)});
% compare nigha nd mornig
p = signrank(y(:,1), y(:,2));
if p<0.2
    lt_plot_pvalue(p, ['vs,'], 1);
end

% - same
x = [4 5];
y = Yall_same;
plot(x, y', '-', 'Color', 'r');
lt_plot_bar(x, mean(y,1), {'Errors', lt_sem(y)});
% compare nigha nd mornig
p = signrank(y(:,1), y(:,2));
if p<0.2
    lt_plot_pvalue(p, ['vs,'], 1);
end

% - same
x = [7 8];
y = Yall_same(mean(Yall_same,2)<0,:);
plot(x, y', '-', 'Color', 'r');
lt_plot_bar(x, mean(y,1), {'Errors', lt_sem(y)});
% compare nigha nd mornig
p = signrank(y(:,1), y(:,2));
if p<0.2
    lt_plot_pvalue(p, ['vs,'], 1);
end

% - same
x = [10 11];
y = Yall_same(mean(Yall_same,2)>0,:);
plot(x, y', '-', 'Color', 'r');
lt_plot_bar(x, mean(y,1), {'Errors', lt_sem(y)});
% compare nigha nd mornig
p = signrank(y(:,1), y(:,2));
if p<0.2
    lt_plot_pvalue(p, ['vs,'], 1);
end


% - diff
x = [13 14];
y = Yall_diff;
plot(x, y', '-', 'Color', 'r');
lt_plot_bar(x, mean(y,1), {'Errors', lt_sem(y)});
% compare nigha nd mornig
p = signrank(y(:,1), y(:,2));
if p<0.2
    lt_plot_pvalue(p, ['vs,'], 1);
end



% =================== GENERALIZATION
lt_subplot(3,1,2); hold on;
title('similar, but plotting generalization');
ylabel('first mean ff across days, then norm to targ');
xlabel('same(all) -- same(learn<0) -- same(learn>0) -- diff');

% - same
x = [1 2];
y = Yall_same_gen;
plot(x, y', '-', 'Color', 'r');
lt_plot_bar(x, mean(y,1), {'Errors', lt_sem(y)});
% compare nigha nd mornig
p = signrank(y(:,1), y(:,2));
if p<0.2
    lt_plot_pvalue(p, ['vs,'], 1);
end

% - same (neg learn)
x = [4 5];
y = Yall_same_gen(mean(Yall_same,2)<0, :);
plot(x, y', '-', 'Color', 'r');
lt_plot_bar(x, mean(y,1), {'Errors', lt_sem(y)});
% compare nigha nd mornig
p = signrank(y(:,1), y(:,2));
if p<0.2
    lt_plot_pvalue(p, ['vs,'], 1);
end

% - same (pos learn)
x = [7 8];
y = Yall_same_gen(mean(Yall_same,2)>0, :);
plot(x, y', '-', 'Color', 'r');
lt_plot_bar(x, mean(y,1), {'Errors', lt_sem(y)});
% compare nigha nd mornig
p = signrank(y(:,1), y(:,2));
if p<0.2
    lt_plot_pvalue(p, ['vs,'], 1);
end

% - diff
x = [10 11];
y = Yall_diff_gen;
plot(x, y', '-', 'Color', 'r');
lt_plot_bar(x, mean(y,1), {'Errors', lt_sem(y)});
% compare nigha nd mornig
p = signrank(y(:,1), y(:,2));
if p<0.2
    lt_plot_pvalue(p, ['vs,'], 1);
end


% ======================= CORRELATION BETWEEN END OF DAY VS. REVERSION?
lt_subplot(3,3,7); hold on;
xlabel('end of day(n)')
ylabel('morning(n+1)');
title('targ');
x = Yall_targ(:,1);
y = Yall_targ(:,2);
plot(x,y, 'ok');
lt_plot_makesquare_plot45line(gca, 'k');

lt_subplot(3,3,8); hold on;
xlabel('end of day(n)')
ylabel('morning(n+1)');
title('same');
x = Yall_same(:,1);
y = Yall_same(:,2);
plot(x,y, 'ok');
lt_plot_makesquare_plot45line(gca, 'k');

lt_subplot(3,3,9); hold on;
xlabel('end of day(n)')
ylabel('morning(n+1)');
title('diff');
x = Yall_diff(:,1);
y = Yall_diff(:,2);
plot(x,y, 'ok');
lt_plot_makesquare_plot45line(gca, 'k');

%% ================ [PLOT] specificity index

disp('NOTE: this is best done WITHOUT combined across syls ...');
FFtargAll = [];
FFsumAll = [];
for i=1:maxbirds
    for ii=1:maxexpts
        
        indtarg = DATSTRUCT.AllBirdnum==i & DATSTRUCT.AllExptnum==ii & DATSTRUCT.AllIsTarg==1;
        indrest = DATSTRUCT.AllBirdnum==i & DATSTRUCT.AllExptnum==ii & DATSTRUCT.AllIsTarg==0;
        
        if ~any(indtarg)
            continue
        end
        
        fftarg = DATSTRUCT.AllFFedges(indtarg,:);
        ffrest = DATSTRUCT.AllFFedges(indrest,:);
        
        % === calculate specificity index
        ffall = sum([fftarg; ffrest],1);
        
        FFtargAll = [FFtargAll; fftarg];
        FFsumAll = [FFsumAll; ffall];
    end
end

% == one slope for each edge
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for j=1:size(FFtargAll,2)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots; hsplot];
    title(['edgenum: ' num2str(j)]);
    xlabel('sum');
    ylabel('targ');
    
    x = FFsumAll(:,j);
    y = FFtargAll(:,j);
    
    lt_regress(y, x, 1, 0, 1, 1, 'k');
    
end
linkaxes(hsplots, 'xy');
%% ================= [PLOT] --- MEAN CHANGE IN FF..

lt_figure; hold on;


% =============================== TARGET
lt_subplot(3,1,1); hold on;
title('target');
inds = DATSTRUCT.AllIsTarg==1 & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,:),1,2);
xdiff = 1:size(ffdiff,2);

% --- line connecting all
plot(xdiff, ffdiff', '-', 'Color', [0.7 0.7 0.7]);

% ---- DAYTIME
pcol ='k';
ffdifftmp = ffdiff(:,1:2:end);
xdifftmp = xdiff(1:2:end);
plot(xdifftmp, ffdifftmp', 'x', 'Color', pcol);
% -- plot bar
% lt_plot_bar(xdifftmp, mean(ffdifftmp,1), {'Errors', lt_sem(ffdifftmp), ...
%     'Color', 'none', 'BarWidth', 0.2});
lt_plot_bar(xdifftmp+0.2, mean(ffdifftmp,1), {'Errors', lt_sem(ffdifftmp), ...
    'Color', 'none', 'BarWidth', 0.2});

% ---- NIGHTTIME
pcol ='r';
ffdifftmp = ffdiff(:,2:2:end);
xdifftmp = xdiff(2:2:end);
plot(xdifftmp, ffdifftmp', 'x', 'Color', pcol);
% -- plot bar
lt_plot_bar(xdifftmp+0.2, mean(ffdifftmp,1), {'Errors', lt_sem(ffdifftmp), ...
    'Color', 'none', 'BarWidth', 0.2});

% ---- OTHER STUFF
lt_plot_zeroline;



% =============================== SAME TYPE
lt_subplot(3,1,2); hold on;
title('same type');
inds = DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsSame==1 ...
    & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,:),1,2);
xdiff = 1:size(ffdiff,2);

% --- line connecting all
plot(xdiff, ffdiff', '-', 'Color', [0.7 0.7 0.7]);

% ---- DAYTIME
pcol ='k';
ffdifftmp = ffdiff(:,1:2:end);
xdifftmp = xdiff(1:2:end);
plot(xdifftmp, ffdifftmp', 'x', 'Color', pcol);
% -- plot bar
% lt_plot_bar(xdifftmp, mean(ffdifftmp,1), {'Errors', lt_sem(ffdifftmp), ...
%     'Color', 'none', 'BarWidth', 0.2});
lt_plot_bar(xdifftmp+0.2, mean(ffdifftmp,1), {'Errors', lt_sem(ffdifftmp), ...
    'Color', 'none', 'BarWidth', 0.2});

% ---- NIGHTTIME
pcol ='r';
ffdifftmp = ffdiff(:,2:2:end);
xdifftmp = xdiff(2:2:end);
plot(xdifftmp, ffdifftmp', 'x', 'Color', pcol);
% -- plot bar
lt_plot_bar(xdifftmp+0.2, mean(ffdifftmp,1), {'Errors', lt_sem(ffdifftmp), ...
    'Color', 'none', 'BarWidth', 0.2});

% ---- OTHER STUFF
lt_plot_zeroline;




% =============================== DIFF TYPE
lt_subplot(3,1,3); hold on;
title('diff type');
inds = DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsSame==0 ...
    & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,:),1,2);
xdiff = 1:size(ffdiff,2);

% --- line connecting all
plot(xdiff, ffdiff', '-', 'Color', [0.7 0.7 0.7]);

% ---- DAYTIME
pcol ='k';
ffdifftmp = ffdiff(:,1:2:end);
xdifftmp = xdiff(1:2:end);
plot(xdifftmp, ffdifftmp', 'x', 'Color', pcol);
% -- plot bar
% lt_plot_bar(xdifftmp, mean(ffdifftmp,1), {'Errors', lt_sem(ffdifftmp), ...
%     'Color', 'none', 'BarWidth', 0.2});
lt_plot_bar(xdifftmp+0.2, mean(ffdifftmp,1), {'Errors', lt_sem(ffdifftmp), ...
    'Color', 'none', 'BarWidth', 0.2});

% ---- NIGHTTIME
pcol ='r';
ffdifftmp = ffdiff(:,2:2:end);
xdifftmp = xdiff(2:2:end);
plot(xdifftmp, ffdifftmp', 'x', 'Color', pcol);
% -- plot bar
lt_plot_bar(xdifftmp+0.2, mean(ffdifftmp,1), {'Errors', lt_sem(ffdifftmp), ...
    'Color', 'none', 'BarWidth', 0.2});

% ---- OTHER STUFF
lt_plot_zeroline;



%% ############## [SCATTER PLOT] TARG/NONTARG, OVER EDGES (MONRING, NIGHT)
if combineSylsInExpt==1 & onlyifhaveAllSylTypes==1
    numbasedays  = -ParamsTrial.DayWindow(1);
    numdays = -ParamsTrial.DayWindow(1) + ParamsTrial.DayWindow(2);
    timeinds = numbasedays*2:numdays*2; % from ngiht of last base day...
    
    % ============ PLOT EACH DAY-NIGHT PAIR INDIVIDUALLY
    % ========= COLLECT (this ensures that same and targ are matched data)
    Yall_targ = [];
    Yall_same = [];
    for i=1:maxbirds
        for ii=1:maxexpts
            
            indtarg = DATSTRUCT.AllBirdnum==i & DATSTRUCT.AllExptnum==ii & DATSTRUCT.AllIsTarg==1;
            indsame = DATSTRUCT.AllBirdnum==i & DATSTRUCT.AllExptnum==ii & DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsSame==1;
            
            if ~any(indtarg)
                continue
            end
            
            if sum(indsame)~=1
                disp('need 1:1 to do this analysis...');
                continue
            end
            
            Yall_targ = [Yall_targ; DATSTRUCT.AllFFedges(indtarg,:)];
            Yall_same = [Yall_same; DATSTRUCT.AllFFedges(indsame,:)];
        end
    end
    
    
    % ======================== PLOTS
    figcount=1;
    subplotrows=3;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    daystoget = [1:numdays-1]'; % will get the day plus following night.
    timeinds_edges = [daystoget*2-1 daystoget*2+1];
    
    % ========= COLLECT THINGS
    Learn_day_targ = []; % targid x daynum
    Learn_day_nontarg = [];
    Learn_night_targ = [];
    Learn_night_nontarg = [];
    
    for j=1:size(timeinds_edges,1)
        
        timebinsthis = timeinds_edges(j,1):timeinds_edges(j,2);
        
        % -- targs
        ytarg = Yall_targ(:, timebinsthis);
        ysame = Yall_same(:, timebinsthis);
        
        % =========== plot
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        xlabel('targ');
        ylabel('same');
        title(['timebins: ' num2str(timebinsthis)]);
        
        % ----- INDIVIDUAL LINES (entire monring-night-morning lines)
        for jj=1:size(ytarg,1)
            plot(ytarg(jj,:), ysame(jj,:), '-', 'Color', [0.7 0.7 0.7]);
        end
        
        % --- COLOR DOTS FOR MORNING AND NIGHT
        % - morning 1
        indtmp = 1;
        pcol = 'b';
        plot(ytarg(:,indtmp), ysame(:,indtmp), 'x', 'Color', pcol);
        lt_plot(mean(ytarg(:,indtmp)), mean(ysame(:, indtmp)), {'Color', pcol});
        % - night 1
        indtmp = 2;
        pcol = 'k';
        plot(ytarg(:,indtmp), ysame(:,indtmp), 'x', 'Color', pcol);
        lt_plot(mean(ytarg(:,indtmp)), mean(ysame(:, indtmp)), {'Color', pcol});
        % - morning 2
        indtmp = 3;
        pcol = 'r';
        plot(ytarg(:,indtmp), ysame(:,indtmp), 'x', 'Color', pcol);
        lt_plot(mean(ytarg(:,indtmp)), mean(ysame(:, indtmp)), {'Color', pcol});
        lt_plot_makesquare_plot45line(gca, 'k');
        
        
        % ================ COLLECT FOR SUMMARY PLOT
        Learn_day_targ = [Learn_day_targ ytarg(:,2)-ytarg(:,1)]; % targid x daynum
        Learn_day_nontarg = [Learn_day_nontarg ysame(:,2)-ysame(:,1)];
        Learn_night_targ = [Learn_night_targ ytarg(:,3)-ytarg(:,2)];
        Learn_night_nontarg = [Learn_night_nontarg ysame(:,3)-ysame(:,2)];
        
    end
    linkaxes(hsplots, 'xy');
    
    % ================== summary figures
    lt_figure; hold on;
    lt_subplot(2,2,1); hold on;
    xlabel('targ')
    ylabel('nontarg');
    title('day(b); night(r)');
    
    % -- day
    pcol = 'b';
    x = Learn_day_targ(:, numbasedays+1:end);
    y = Learn_day_nontarg(:, numbasedays+1:end);
    x = x(:); y = y(:);
    plot(x,y, 'x', 'Color', pcol);
    
    % -- night
    pcol = 'r';
    x = Learn_night_targ(:, numbasedays+1:end);
    y = Learn_night_nontarg(:, numbasedays+1:end);
    x = x(:); y = y(:);
    plot(x,y, 'x', 'Color', pcol);
    
    % --
    lt_plot_makesquare_plot45line(gca, 'k')
    
    
    % ===============================
    lt_subplot(2,2,2); hold on;
    xlabel('targ')
    ylabel('nontarg');
    title('day(b); night(r) [each expt one point');
    
    % -- day
    pcol = 'b';
    x = Learn_day_targ(:, numbasedays+1:end);
    y = Learn_day_nontarg(:, numbasedays+1:end);
    x = mean(x,2); y = mean(y,2);
    plot(x,y, 'x', 'Color', pcol);
    
    % -- night
    pcol = 'r';
    x = Learn_night_targ(:, numbasedays+1:end);
    y = Learn_night_nontarg(:, numbasedays+1:end);
    x = mean(x,2); y = mean(y,2);
    plot(x,y, 'x', 'Color', pcol);
    
    % --
    lt_plot_makesquare_plot45line(gca, 'k')
    
    
    % =========
    lt_subplot(2,2,3); hold on;
    xlabel('day');
    ylabel('night');
    title('targ(k), nontarg(m)');
    
    % -- targ
    pcol = 'k';
    x = Learn_day_targ(:, numbasedays+1:end);
    y = Learn_night_targ(:, numbasedays+1:end);
    x = x(:); y = y(:);
    plot(x,y, 'x', 'Color', pcol);
    
    % -- nontarg
    pcol = 'm';
    x = Learn_day_nontarg(:, numbasedays+1:end);
    y = Learn_night_nontarg(:, numbasedays+1:end);
    x = x(:); y = y(:);
    plot(x,y, 'x', 'Color', pcol);
    
    % --
    lt_plot_makesquare_plot45line(gca, 'k')
    
    
    % =========
    lt_subplot(2,2,4); hold on;
    xlabel('day');
    ylabel('night');
    title('targ(k), nontarg(m) [N=expts]');
    
    % -- targ
    pcol = 'k';
    x = Learn_day_targ(:, numbasedays+1:end);
    y = Learn_night_targ(:, numbasedays+1:end);
    x = mean(x,2);
    y = mean(y,2);
    plot(x,y, 'x', 'Color', pcol);
    
    % -- nontarg
    pcol = 'm';
    x = Learn_day_nontarg(:, numbasedays+1:end);
    y = Learn_night_nontarg(:, numbasedays+1:end);
    x = mean(x,2);
    y = mean(y,2);
    plot(x,y, 'x', 'Color', pcol);
    
    % --
    lt_plot_makesquare_plot45line(gca, 'k')
end
%% =============== [PLOT] OVERDAY AND OVERNIGHT, ONE VALUE FOR EACH SYL
% ==== FIGURE OUT INDS FOR DURING WN ONLY
numbasedays = -ParamsTrial.DayWindow(1);
ind_wnon = (2*numbasedays)+1; % because each day has 2 values ...
plottext=1;

% ######################################## PLOTS
lt_figure; hold on;


% =============================== TARGET
lt_subplot(2,2,1); hold on;
title('target');
inds = DATSTRUCT.AllIsTarg==1 & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end),1,2);

Y = [];
% ---- DAYTIME
ffdifftmp = ffdiff(:,1:2:end);
Y(:,1) = mean(ffdifftmp,2);
% ---- OVERNIGHT
ffdifftmp = ffdiff(:,2:2:end);
Y(:,2) = mean(ffdifftmp,2);

% ======== plot
plot([1 2], Y', '-', 'Color', [0.7 0.7 0.7]);
lt_plot([1 2]+0.1, mean(Y), {'Errors', lt_sem(Y)});
lt_plot_zeroline;
% --- stats
for j=1:2
    y = Y(:,j);
    p = signrank(y);
    %     [~, p] = ttest(y);
    if p<0.1
        lt_plot_text(j, 1.1*max(y), ['p=' num2str(p)], 'r');
    end
end
p = signrank(Y(:,1), Y(:,2));
% [~, p] = ttest(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'srank vs', 1);

% ------------ PLOT TEXT
if plottext==1
    indnums = find(inds)';
    for j=1:length(indnums)
        
        bnum = DATSTRUCT.AllBirdnum(indnums(j));
        bname = TrialStruct.birds(bnum).birdname;
        enum = DATSTRUCT.AllExptnum(indnums(j));
        ename = TrialStruct.birds(bnum).exptnum(enum).exptname;
        
        if combineSylsInExpt==0
            syl = DATSTRUCT.AllSylname{indnums(j)};
            textstring = [bname '-' ename '-' syl];
        else
            textstring = [bname '-' ename];
        end
        lt_plot_text(2.1, Y(j,2), textstring , 'b', 6);
    end
end



% =============================== SAME
lt_subplot(2,2,2); hold on;
title('SAME');
inds = DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsSame==1 ...
    & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end),1,2);

Y = [];
% ---- DAYTIME
ffdifftmp = ffdiff(:,1:2:end);
Y(:,1) = mean(ffdifftmp,2);
% ---- OVERNIGHT
ffdifftmp = ffdiff(:,2:2:end);
Y(:,2) = mean(ffdifftmp,2);

% ======== plot
plot([1 2], Y', '-', 'Color', [0.7 0.7 0.7]);
lt_plot([1 2]+0.1, mean(Y), {'Errors', lt_sem(Y)});
lt_plot_zeroline;
% --- stats
for j=1:2
    y = Y(:,j);
    p = signrank(y);
    %     [~, p] = ttest(y);
    if p<0.1
        lt_plot_text(j, 1.1*max(y), ['p=' num2str(p)], 'r');
    end
end
p = signrank(Y(:,1), Y(:,2));
% [~, p] = ttest(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'srank vs', 1);

% ------------ PLOT TEXT
if plottext==1
    indnums = find(inds)';
    for j=1:length(indnums)
        
        bnum = DATSTRUCT.AllBirdnum(indnums(j));
        bname = TrialStruct.birds(bnum).birdname;
        enum = DATSTRUCT.AllExptnum(indnums(j));
        ename = TrialStruct.birds(bnum).exptnum(enum).exptname;
        
        if combineSylsInExpt==0
            syl = DATSTRUCT.AllSylname{indnums(j)};
            textstring = [bname '-' ename '-' syl];
        else
            textstring = [bname '-' ename];
        end
        lt_plot_text(2.1, Y(j,2), textstring , 'b', 6);
    end
end



% =============================== DIFF
lt_subplot(2,2,3); hold on;
title('DIFF');
inds = DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsSame==0 ...
    & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end),1,2);

Y = [];
% ---- DAYTIME
ffdifftmp = ffdiff(:,1:2:end);
Y(:,1) = mean(ffdifftmp,2);
% ---- OVERNIGHT
ffdifftmp = ffdiff(:,2:2:end);
Y(:,2) = mean(ffdifftmp,2);

% ======== plot
plot([1 2], Y', '-', 'Color', [0.7 0.7 0.7]);
lt_plot([1 2]+0.1, mean(Y), {'Errors', lt_sem(Y)});
lt_plot_zeroline;
% --- stats
for j=1:2
    y = Y(:,j);
    p = signrank(y);
    %     [~, p] = ttest(y);
    if p<0.1
        lt_plot_text(j, 1.1*max(y), ['p=' num2str(p)], 'r');
    end
end
p = signrank(Y(:,1), Y(:,2));
% [~, p] = ttest(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'srank vs', 1);
% ------------ PLOT TEXT
if plottext==1
    indnums = find(inds)';
    for j=1:length(indnums)
        
        bnum = DATSTRUCT.AllBirdnum(indnums(j));
        bname = TrialStruct.birds(bnum).birdname;
        enum = DATSTRUCT.AllExptnum(indnums(j));
        ename = TrialStruct.birds(bnum).exptnum(enum).exptname;
        
        if combineSylsInExpt==0
            syl = DATSTRUCT.AllSylname{indnums(j)};
            textstring = [bname '-' ename '-' syl];
        else
            textstring = [bname '-' ename];
        end
        lt_plot_text(2.1, Y(j,2), textstring , 'b', 6);
    end
    
end

%% =============== [PLOT, POOLING OVER DAYS] OVERDAY AND OVERNIGHT, ONE VALUE FOR EACH SYL
% ==== FIGURE OUT INDS FOR DURING WN ONLY
numbasedays = -ParamsTrial.DayWindow(1);
ind_wnon = (2*numbasedays)+1; % because each day has 2 values ...

% ######################################## PLOTS
lt_figure; hold on;


% =============================== TARGET
lt_subplot(2,2,1); hold on;
title('target');
inds = DATSTRUCT.AllIsTarg==1 & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end-1),1,2);

Y = {};
% ---- DAYTIME
ffdifftmp = ffdiff(:,1:2:end);
Y{1} = ffdifftmp(:);
% ---- OVERNIGHT
ffdifftmp = ffdiff(:,2:2:end);
Y{2} = ffdifftmp(:);

% ======== plot
lt_plot_MultDist(Y, [1 2], 0);
lt_plot_zeroline;
% --- stats
for j=1:2
    y = Y{j};
    p = signrank(y);
    %     [~, p] = ttest(y);
    if p<0.1
        lt_plot_text(j, 1.1*max(y), ['p=' num2str(p)], 'r');
    end
end
p = ranksum(Y{1}, Y{2});
% [~, p] = ttest(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'ranksum vs', 1);


% =============================== TARGET
lt_subplot(2,2,1); hold on;
title('target');
inds = DATSTRUCT.AllIsTarg==1 & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end-1),1,2);

Y = {};
% ---- DAYTIME
ffdifftmp = ffdiff(:,1:2:end);
Y{1} = ffdifftmp(:);
% ---- OVERNIGHT
ffdifftmp = ffdiff(:,2:2:end);
Y{2} = ffdifftmp(:);

% ======== plot
lt_plot_MultDist(Y, [1 2], 0);
lt_plot_zeroline;
% --- stats
for j=1:2
    y = Y{j};
    p = signrank(y);
    %     [~, p] = ttest(y);
    if p<0.1
        lt_plot_text(j, 1.1*max(y), ['p=' num2str(p)], 'r');
    end
end
p = ranksum(Y{1}, Y{2});
% [~, p] = ttest(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'ranksum vs', 1);


% =============================== SAME TYPE
lt_subplot(2,2,2); hold on;
title('SAME');
inds = DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsSame==1 ...
    & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end-1),1,2);

Y = {};
% ---- DAYTIME
ffdifftmp = ffdiff(:,1:2:end);
Y{1} = ffdifftmp(:);
% ---- OVERNIGHT
ffdifftmp = ffdiff(:,2:2:end);
Y{2} = ffdifftmp(:);

% ======== plot
lt_plot_MultDist(Y, [1 2], 0);
lt_plot_zeroline;
% --- stats
for j=1:2
    y = Y{j};
    p = signrank(y);
    %     [~, p] = ttest(y);
    if p<0.1
        lt_plot_text(j, 1.1*max(y), ['p=' num2str(p)], 'r');
    end
end
p = ranksum(Y{1}, Y{2});
% [~, p] = ttest(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'ranksum vs', 1);



% =============================== SAME TYPE
lt_subplot(2,2,3); hold on;
title('DIFF');
inds = DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsSame==0 ...
    & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end-1),1,2);

Y = {};
% ---- DAYTIME
ffdifftmp = ffdiff(:,1:2:end);
Y{1} = ffdifftmp(:);
% ---- OVERNIGHT
ffdifftmp = ffdiff(:,2:2:end);
Y{2} = ffdifftmp(:);

% ======== plot
lt_plot_MultDist(Y, [1 2], 0);
lt_plot_zeroline;
% --- stats
for j=1:2
    y = Y{j};
    p = signrank(y);
    %     [~, p] = ttest(y);
    if p<0.1
        lt_plot_text(j, 1.1*max(y), ['p=' num2str(p)], 'r');
    end
end
p = ranksum(Y{1}, Y{2});
% [~, p] = ttest(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'ranksum vs', 1);



%% =============== [PLOT, MIXED EFFECTS] OVERDAY AND OVERNIGHT, ONE VALUE FOR EACH SYL
% ==== FIGURE OUT INDS FOR DURING WN ONLY
numbasedays = -ParamsTrial.DayWindow(1);
ind_wnon = (2*numbasedays)+1; % because each day has 2 values ...

% formula = 'Y ~ IsDay +(IsDay-1|ExptID) + (1|ExptID)';
formula = 'Y ~ IsDay +(IsDay-1|ExptID) + (1|ExptID)';

% ######################################################  TARGET
inds = DATSTRUCT.AllIsTarg==1 & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end-1),1,2);
bnums = DATSTRUCT.AllBirdnum(inds);
enums= DATSTRUCT.AllExptnum(inds);
exptIDs = grp2idx(num2str([bnums enums]));
exptIDs = repmat(exptIDs, 1, size(ffdiff,2));

% ==================== COLLECT INTO TABLE
Y = [];
IsDay = [];
ExptID = [];

% ---- DAYTIME
ffdifftmp = ffdiff(:,1:2:end);
exptidthis = exptIDs(:,1:2:end);

Y = [Y; ffdifftmp(:)];
IsDay = [IsDay; 1*ones(size(ffdifftmp(:),1),1)];
ExptID = [ExptID; exptidthis(:)];

% ------- OVERNIGHT
ffdifftmp = ffdiff(:,2:2:end);
exptidthis = exptIDs(:,2:2:end);

Y = [Y; ffdifftmp(:)];
IsDay = [IsDay; 0*ones(size(ffdifftmp(:),1),1)];
ExptID = [ExptID; exptidthis(:)];

% ======================= DO LME
tbl = table(Y, IsDay, ExptID);
lme = fitlme(tbl, formula);

% figure; hold on;
% lt_plot_MultDist({Y(IsDay==1), Y(IsDay==0)}, [1 2]);


% ######################################################  SAME
inds = DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsSame==1 ...
    & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end-1),1,2);
bnums = DATSTRUCT.AllBirdnum(inds);
enums= DATSTRUCT.AllExptnum(inds);
exptIDs = grp2idx(num2str([bnums enums]));
exptIDs = repmat(exptIDs, 1, size(ffdiff,2));

% ==================== COLLECT INTO TABLE
Y = [];
IsDay = [];
ExptID = [];

% ---- DAYTIME
ffdifftmp = ffdiff(:,1:2:end);
exptidthis = exptIDs(:,1:2:end);

Y = [Y; ffdifftmp(:)];
IsDay = [IsDay; 1*ones(size(ffdifftmp(:),1),1)];
ExptID = [ExptID; exptidthis(:)];

% ------- OVERNIGHT
ffdifftmp = ffdiff(:,2:2:end);
exptidthis = exptIDs(:,2:2:end);

Y = [Y; ffdifftmp(:)];
IsDay = [IsDay; 0*ones(size(ffdifftmp(:),1),1)];
ExptID = [ExptID; exptidthis(:)];

% ======================= DO LME
tbl = table(Y, IsDay, ExptID);
lme = fitlme(tbl, formula);

figure; hold on;
lt_plot_MultDist({Y(IsDay==1), Y(IsDay==0)}, [1 2]);

%% =============== [SCATTER PLOT] OVERDAY AND OVERNIGHT, ONE VALUE FOR EACH SYL
% ==== FIGURE OUT INDS FOR DURING WN ONLY
numbasedays = -ParamsTrial.DayWindow(1);
ind_wnon = (2*numbasedays)+1; % because each day has 2 values ...


% ######################################## PLOTS
lt_figure; hold on;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  TARGET
lt_subplot(3,2,1); hold on;
title('target (syl=datpoint)');
xlabel('over day');
ylabel('over night');
inds = DATSTRUCT.AllIsTarg==1 & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end-1),1,2);

Y = [];
% ---- DAYTIME
ffdifftmp = ffdiff(:,1:2:end);
Y(:,1) = mean(ffdifftmp,2);
% ---- OVERNIGHT
ffdifftmp = ffdiff(:,2:2:end);
Y(:,2) = mean(ffdifftmp,2);

% ======== plot
plot(Y(:,1), Y(:,2), 'ok');
lt_plot_makesquare_plot45line(gca, 'b');

% ##########################
lt_subplot(3,2,2); hold on;
title('target (day/night=datpoint)');
xlabel('over day');
ylabel('over night');
inds = DATSTRUCT.AllIsTarg==1 & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end-1),1,2);

Y = [];
% ---- DAYTIME
ffdifftmp = ffdiff(:,1:2:end);
Y(:,1) = ffdifftmp(:);
% ---- OVERNIGHT
ffdifftmp = ffdiff(:,2:2:end);
Y(:,2) = ffdifftmp(:);

% ======== plot
plot(Y(:,1), Y(:,2), 'ok');
lt_plot_makesquare_plot45line(gca, 'b');


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SAME
lt_subplot(3,2,3); hold on;
title('SAME');
inds = DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsSame==1 ...
    & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end-1),1,2);

Y = [];
% ---- DAYTIME
ffdifftmp = ffdiff(:,1:2:end);
Y(:,1) = mean(ffdifftmp,2);
% ---- OVERNIGHT
ffdifftmp = ffdiff(:,2:2:end);
Y(:,2) = mean(ffdifftmp,2);

% ======== plot

plot(Y(:,1), Y(:,2), 'ok');
lt_plot_makesquare_plot45line(gca, 'b');

% ###############
lt_subplot(3,2,4); hold on;
title('SAME');
inds = DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsSame==1 ...
    & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end-1),1,2);

Y = [];
% ---- DAYTIME
ffdifftmp = ffdiff(:,1:2:end);
Y(:,1) = ffdifftmp(:);
% ---- OVERNIGHT
ffdifftmp = ffdiff(:,2:2:end);
Y(:,2) = ffdifftmp(:);

% ======== plot

plot(Y(:,1), Y(:,2), 'ok');
lt_plot_makesquare_plot45line(gca, 'b');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DIFF
lt_subplot(3,2,5); hold on;
title('DIFF');
inds = DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsSame==0 ...
    & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end-1),1,2);

Y = [];
% ---- DAYTIME
ffdifftmp = ffdiff(:,1:2:end);
Y(:,1) = mean(ffdifftmp,2);
% ---- OVERNIGHT
ffdifftmp = ffdiff(:,2:2:end);
Y(:,2) = mean(ffdifftmp,2);

% ======== plot

plot(Y(:,1), Y(:,2), 'ok');
lt_plot_makesquare_plot45line(gca, 'b');

% ####################
lt_subplot(3,2,6); hold on;
title('DIFF');
inds = DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsSame==0 ...
    & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end-1),1,2);

Y = [];
% ---- DAYTIME
ffdifftmp = ffdiff(:,1:2:end);
Y(:,1) = ffdifftmp(:);
% ---- OVERNIGHT
ffdifftmp = ffdiff(:,2:2:end);
Y(:,2) = ffdifftmp(:);

% ======== plot

plot(Y(:,1), Y(:,2), 'ok');
lt_plot_makesquare_plot45line(gca, 'b');


%% =============== [SCATTER PLOT] OVERDAY AND OVERNIGHT, ONE VALUE FOR EACH SYL
% OVERNIGHT PREDICT NEXT DAY?
% ==== FIGURE OUT INDS FOR DURING WN ONLY
numbasedays = -ParamsTrial.DayWindow(1);
ind_wnon = (2*numbasedays)+1; % because each day has 2 values ...


% ######################################## PLOTS
lt_figure; hold on;


% ##########################
lt_subplot(3,2,1); hold on;
title('target (day/night=datpoint)');
xlabel('over night');
ylabel('over day(next)');
inds = DATSTRUCT.AllIsTarg==1 & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon+1:end),1,2);

Y = [];
% ---- NGIHT
ffdifftmp = ffdiff(:,1:2:end);
Y(:,1) = ffdifftmp(:);
% ---- NEXT DAY
ffdifftmp = ffdiff(:,2:2:end);
Y(:,2) = ffdifftmp(:);

% ======== plot
plot(Y(:,1), Y(:,2), 'ok');
lt_plot_makesquare_plot45line(gca, 'b');


% ########################## SAME
lt_subplot(3,2,2); hold on;
title('same (day/night=datpoint)');
xlabel('over night');
ylabel('over day(next)');
inds = DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsSame==1 & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon+1:end),1,2);

Y = [];
% ---- NGIHT
ffdifftmp = ffdiff(:,1:2:end);
Y(:,1) = ffdifftmp(:);
% ---- NEXT DAY
ffdifftmp = ffdiff(:,2:2:end);
Y(:,2) = ffdifftmp(:);

% ======== plot
plot(Y(:,1), Y(:,2), 'ok');
lt_plot_makesquare_plot45line(gca, 'b');


% ########################## DIFF
lt_subplot(3,2,3); hold on;
title('diff (day/night=datpoint)');
xlabel('over night');
ylabel('over day(next)');
inds = DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsSame==0 & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon+1:end),1,2);

Y = [];
% ---- NGIHT
ffdifftmp = ffdiff(:,1:2:end);
Y(:,1) = ffdifftmp(:);
% ---- NEXT DAY
ffdifftmp = ffdiff(:,2:2:end);
Y(:,2) = ffdifftmp(:);

% ======== plot
plot(Y(:,1), Y(:,2), 'ok');
lt_plot_makesquare_plot45line(gca, 'b');


% ======== plot

plot(Y(:,1), Y(:,2), 'ok');
lt_plot_makesquare_plot45line(gca, 'b');

%% IS GENERALIZATION GREATER IN MORNING OR EVENING?
% =========== FOR EVERY SYL, ASK WHETHER OVER DAY OR OVER NGIHT CONTRIBUTES
% MORE TO EXPRESSION OF LEARNING (totla learning is last morning minus
% first morning. over day and over night are sum over all days)
lt_figure; hold on;
binsize = 6; % for sommothign.
% ==== FIGURE OUT INDS FOR DURING WN ONLY
numbasedays = -ParamsTrial.DayWindow(1);
ind_wnon = (2*numbasedays)+1; % because each day has 2 values ...


% ######################## TARGET
lt_subplot(2,2,1); hold on;
title('target');
xlabel('morning-morning');
ylabel('learn(day=r; night=b)');
inds = DATSTRUCT.AllIsTarg==1 & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- TOTAL LEARNING (all morning-morning pairs)
mornings = DATSTRUCT.AllFFedges(inds,ind_wnon:2:end);
diff_mornings = diff(mornings,1,2);
diff_mornings = diff_mornings(:);

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end),1,2);

Y = [];
% ---- DAYTIME
ffdifftmp = ffdiff(:,1:2:end-1);
Y = [Y ffdifftmp(:)];

% ---- NIGHTTIME
ffdifftmp = ffdiff(:,2:2:end-1);
Y = [Y ffdifftmp(:)];

% -- sort by diffmorings
[~, indsort] = sort(diff_mornings);
diff_mornings = diff_mornings(indsort);
Y = Y(indsort,:);
% =================== PLOT
xtmp = lt_running_stats(diff_mornings, binsize);
% ---- overday
plot(diff_mornings, Y(:,1), 'or');
ytmp = lt_running_stats(Y(:,1), binsize);
shadedErrorBar(xtmp.Mean, ytmp.Mean, ytmp.SEM, {'Color', 'r'}, 1);
% ---- overnight
plot(diff_mornings, Y(:,2), 'ob');
ytmp = lt_running_stats(Y(:,2), binsize);
shadedErrorBar(xtmp.Mean, ytmp.Mean, ytmp.SEM, {'Color', 'b'}, 1);
% --
lt_plot_makesquare_plot45line(gca, 'k');

% -------------------- fit linear model
LearnTot = [diff_mornings; diff_mornings];
LearnSub = [Y(:,1); Y(:,2)];
IsDay = [ones(size(Y,1),1); 0*ones(size(Y,1),1)];

tbl = table(LearnTot, LearnSub, IsDay);
formula = 'LearnSub ~ LearnTot + LearnTot:IsDay';
lme = fitlme(tbl, formula)




% ######################## SAME
lt_subplot(2,2,2); hold on;
title('SAME');
xlabel('total learn');
ylabel('learn(day=r; night=b)');
inds = DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsSame==1 ...
    & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- TOTAL LEARNING (all morning-morning pairs)
mornings = DATSTRUCT.AllFFedges(inds,ind_wnon:2:end);
diff_mornings = diff(mornings,1,2);
diff_mornings = diff_mornings(:);

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end),1,2);

Y = [];
% ---- DAYTIME
ffdifftmp = ffdiff(:,1:2:end-1);
Y = [Y ffdifftmp(:)];

% ---- NIGHTTIME
ffdifftmp = ffdiff(:,2:2:end-1);
Y = [Y ffdifftmp(:)];


% -- sort by diffmorings
[~, indsort] = sort(diff_mornings);
diff_mornings = diff_mornings(indsort);
Y = Y(indsort,:);

% =================== PLOT
xtmp = lt_running_stats(diff_mornings, binsize);
% ---- overday
plot(diff_mornings, Y(:,1), 'or');
ytmp = lt_running_stats(Y(:,1), binsize);
shadedErrorBar(xtmp.Mean, ytmp.Mean, ytmp.SEM, {'Color', 'r'}, 1);
% ---- overnight
plot(diff_mornings, Y(:,2), 'ob');
ytmp = lt_running_stats(Y(:,2), binsize);
shadedErrorBar(xtmp.Mean, ytmp.Mean, ytmp.SEM, {'Color', 'b'}, 1);

% --
lt_plot_makesquare_plot45line(gca, 'k');

% -------------------- fit linear model
LearnTot = [diff_mornings; diff_mornings];
LearnSub = [Y(:,1); Y(:,2)];
IsDay = [ones(size(Y,1),1); 0*ones(size(Y,1),1)];

tbl = table(LearnTot, LearnSub, IsDay);
formula = 'LearnSub ~ LearnTot + LearnTot:IsDay';
lme = fitlme(tbl, formula)



% ######################## DIFF
lt_subplot(2,2,3); hold on;
title('DIFF');
xlabel('total learn');
ylabel('learn(day=r; night=b)');
inds = DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsSame==0 ...
    & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- TOTAL LEARNING (all morning-morning pairs)
mornings = DATSTRUCT.AllFFedges(inds,ind_wnon:2:end);
diff_mornings = diff(mornings,1,2);
diff_mornings = diff_mornings(:);

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end),1,2);

Y = [];
% ---- DAYTIME
ffdifftmp = ffdiff(:,1:2:end-1);
Y = [Y ffdifftmp(:)];

% ---- NIGHTTIME
ffdifftmp = ffdiff(:,2:2:end-1);
Y = [Y ffdifftmp(:)];


% -- sort by diffmorings
[~, indsort] = sort(diff_mornings);
diff_mornings = diff_mornings(indsort);
Y = Y(indsort,:);

% =================== PLOT
xtmp = lt_running_stats(diff_mornings, binsize);
% ---- overday
plot(diff_mornings, Y(:,1), 'or');
ytmp = lt_running_stats(Y(:,1), binsize);
shadedErrorBar(xtmp.Mean, ytmp.Mean, ytmp.SEM, {'Color', 'r'}, 1);
% ---- overnight
plot(diff_mornings, Y(:,2), 'ob');
ytmp = lt_running_stats(Y(:,2), binsize);
shadedErrorBar(xtmp.Mean, ytmp.Mean, ytmp.SEM, {'Color', 'b'}, 1);
% --
lt_plot_makesquare_plot45line(gca, 'k');

% -------------------- fit linear model
LearnTot = [diff_mornings; diff_mornings];
LearnSub = [Y(:,1); Y(:,2)];
IsDay = [ones(size(Y,1),1); 0*ones(size(Y,1),1)];

tbl = table(LearnTot, LearnSub, IsDay);
formula = 'LearnSub ~ LearnTot + LearnTot:IsDay';
lme = fitlme(tbl, formula)



%% IS GENERALIZATION GREATER IN MORNING OR EVENING?
% [ SIMIALR TO ABOVE, but each day one value (ie.. each morning-->mroning
% is one value)
lt_figure; hold on;

% ==== FIGURE OUT INDS FOR DURING WN ONLY
numbasedays = -ParamsTrial.DayWindow(1);
ind_wnon = (2*numbasedays)+1; % because each day has 2 values ...


% ######################## TARGET
lt_subplot(2,2,1); hold on;
title('target');
xlabel('total learn');
ylabel('learn(day=r; night=b)');
inds = DATSTRUCT.AllIsTarg==1 & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- TOTAL LEARNING (first morning to last morning)
tmp = DATSTRUCT.AllFFedges(inds,ind_wnon:end-1);
learnTot = tmp(:,end) - tmp(:,1);

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end),1,2);

Y = [];
% ---- DAYTIME
ffdifftmp = sum(ffdiff(:,1:2:end),2);
Y = [Y ffdifftmp];

% ---- OVERNIGHT
ffdifftmp = sum(ffdiff(:,2:2:end),2);
Y = [Y ffdifftmp];

% =================== PLOT
% ---- overday
plot(learnTot, Y(:,1), 'or');
% ---- overnight
plot(learnTot, Y(:,2), 'ob');
% --
lt_plot_makesquare_plot45line(gca, 'k');

% -------------------- fit linear model
LearnTot = [learnTot; learnTot];
LearnSub = [Y(:,1); Y(:,2)];
IsDay = [ones(size(Y,1),1); 0*ones(size(Y,1),1)];

tbl = table(LearnTot, LearnSub, IsDay);
formula = 'LearnSub ~ LearnTot + LearnTot:IsDay';
lme = fitlme(tbl, formula)





% ######################## SAME
lt_subplot(2,2,2); hold on;
title('SAME');
xlabel('total learn');
ylabel('learn(day=r; night=b)');
inds = DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsSame==1 ...
    & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- TOTAL LEARNING (first morning to last morning)
tmp = DATSTRUCT.AllFFedges(inds,ind_wnon:end-1);
learnTot = tmp(:,end) - tmp(:,1);

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end),1,2);

Y = [];
% ---- DAYTIME
ffdifftmp = sum(ffdiff(:,1:2:end),2);
Y = [Y ffdifftmp];

% ---- OVERNIGHT
ffdifftmp = sum(ffdiff(:,2:2:end),2);
Y = [Y ffdifftmp];

% =================== PLOT
% ---- overday
plot(learnTot, Y(:,1), 'or');
% ---- overnight
plot(learnTot, Y(:,2), 'ob');

% --
lt_plot_makesquare_plot45line(gca, 'k');

% -------------------- fit linear model
LearnTot = [learnTot; learnTot];
LearnSub = [Y(:,1); Y(:,2)];
IsDay = [ones(size(Y,1),1); 0*ones(size(Y,1),1)];

tbl = table(LearnTot, LearnSub, IsDay);
formula = 'LearnSub ~ LearnTot + LearnTot:IsDay';
lme = fitlme(tbl, formula)


% ######################## DIFF
lt_subplot(2,2,3); hold on;
title('DIFF');
xlabel('total learn');
ylabel('learn(day=r; night=b)');
inds = DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsSame==0 ...
    & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- TOTAL LEARNING (first morning to last morning)
tmp = DATSTRUCT.AllFFedges(inds,ind_wnon:end-1);
learnTot = tmp(:,end) - tmp(:,1);

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end),1,2);

Y = [];
% ---- DAYTIME
ffdifftmp = sum(ffdiff(:,1:2:end),2);
Y = [Y ffdifftmp];

% ---- OVERNIGHT
ffdifftmp = sum(ffdiff(:,2:2:end),2);
Y = [Y ffdifftmp];

% =================== PLOT
% ---- overday
plot(learnTot, Y(:,1), 'or');
% ---- overnight
plot(learnTot, Y(:,2), 'ob');
% --
lt_plot_makesquare_plot45line(gca, 'k');

% -------------------- fit linear model
LearnTot = [learnTot; learnTot];
LearnSub = [Y(:,1); Y(:,2)];
IsDay = [ones(size(Y,1),1); 0*ones(size(Y,1),1)];

tbl = table(LearnTot, LearnSub, IsDay);
formula = 'LearnSub ~ LearnTot + LearnTot:IsDay';
lme = fitlme(tbl, formula)


%% =================== NONTARG LEARNING (DAY/NIGHT SEPARATE) VS. TARG LEARNING

% ==== FIGURE OUT INDS FOR DURING WN ONLY
numbasedays = -ParamsTrial.DayWindow(1);
ind_wnon = (2*numbasedays)+1; % because each day has 2 values ...

TMPSTRUCT = struct; % to collect target, nontargs

% ######################## TARGET
fieldthis = 'target';
inds = DATSTRUCT.AllIsTarg==1 & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- TOTAL LEARNING (first morning to last morning)
tmp = DATSTRUCT.AllFFedges(inds,ind_wnon:end-1);
learnTot = tmp(:,end) - tmp(:,1);

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end-1),1,2);

%  DAYTIME
learnDay_allday = ffdiff(:,1:2:end);
learnDay = sum(ffdiff(:,1:2:end),2);

%  OVERNIGHT
learnNight_allday = ffdiff(:,2:2:end);
learnNight = sum(ffdiff(:,2:2:end),2);

% ---- COLLECT BIRD AND EXPT NUMBERS, TO MATCH TO NONTARG
bnums = DATSTRUCT.AllBirdnum(inds);
enums = DATSTRUCT.AllExptnum(inds);

% ============================= OUTPUT
TMPSTRUCT.(fieldthis).learnTot = learnTot;
TMPSTRUCT.(fieldthis).learnDay_allday = learnDay_allday;
TMPSTRUCT.(fieldthis).learnNight_allday = learnNight_allday;
TMPSTRUCT.(fieldthis).learnDay = learnDay;
TMPSTRUCT.(fieldthis).learnNight = learnNight;
TMPSTRUCT.(fieldthis).bnums = bnums;
TMPSTRUCT.(fieldthis).enums = enums;


% ######################## SAME
fieldthis = 'same';
inds = DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsSame==1 ...
    & ~any(isnan(DATSTRUCT.AllFFedges),2); % makes sure each datapoint has all days

% ----------- TOTAL LEARNING (first morning to last morning)
tmp = DATSTRUCT.AllFFedges(inds,ind_wnon:end-1);
learnTot = tmp(:,end) - tmp(:,1);

% ----------- all differences
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end-1),1,2);

%  DAYTIME
learnDay_allday = ffdiff(:,1:2:end);
learnDay = sum(ffdiff(:,1:2:end),2);

%  OVERNIGHT
learnNight_allday = ffdiff(:,2:2:end);
learnNight = sum(ffdiff(:,2:2:end),2);

% ---- COLLECT BIRD AND EXPT NUMBERS, TO MATCH TO NONTARG
bnums = DATSTRUCT.AllBirdnum(inds);
enums = DATSTRUCT.AllExptnum(inds);

% ============================= OUTPUT
TMPSTRUCT.(fieldthis).learnTot = learnTot;
TMPSTRUCT.(fieldthis).learnDay_allday = learnDay_allday;
TMPSTRUCT.(fieldthis).learnDay = learnDay;
TMPSTRUCT.(fieldthis).learnNight_allday = learnNight_allday;
TMPSTRUCT.(fieldthis).learnNight = learnNight;
TMPSTRUCT.(fieldthis).bnums = bnums;
TMPSTRUCT.(fieldthis).enums = enums;



% ######################################################## PLOT
lt_figure; hold on;

% ========================= SAME
fieldthis = 'same';
lt_subplot(2,2,1); hold on;

% ----------
xlabel('target (total learn, morning-miornig');
ylabel('nontarg (day=r; night = b)');
title(fieldthis);

nsyls = length(TMPSTRUCT.(fieldthis).learnDay);
x_targ = [];
y_day = [];
y_night = [];
for j=1:nsyls
    
    % ---- skip if not positive learning at nontarg (i.e generalization)
    %     if TMPSTRUCT.(fieldthis).learnTot(j)<=0
    %         disp('/adadasdd')
    %         continue
    %     end
    
    % ---------------------- find corresponding target data
    bnumthis = TMPSTRUCT.(fieldthis).bnums(j);
    enumthis = TMPSTRUCT.(fieldthis).enums(j);
    
    indtarg = TMPSTRUCT.target.bnums==bnumthis & TMPSTRUCT.target.enums == enumthis;
    assert(sum(indtarg)==1, 'asdfasd');
    
    % --- extract total learning for corresponding targ
    x_targ = [x_targ; TMPSTRUCT.target.learnTot(indtarg)];
    
    % --- extract day and night learning for nontarg
    y_day = [y_day; TMPSTRUCT.(fieldthis).learnDay(j)];
    y_night = [y_night; TMPSTRUCT.(fieldthis).learnNight(j)];
end

plot(x_targ, y_day, 'or');
lt_regress(y_day, x_targ, 1, 0, 1, 1, 'r');
plot(x_targ, y_night, 'ob');
lt_regress(y_night, x_targ, 1, 0, 1, 1, 'b');

lt_plot_makesquare_plot45line(gca, 'k')

% ---------
lt_subplot(2,2,2); hold on;
xlabel('day-night(bars), LOW-HIGH TARG LEARN[SETS OF BARS]');
ylabel('hz(total over days)');
title('effect of total learn at target on day/night for nontarg (each syl one dat');
% --- low targ learning
indsthis = x_targ<median(x_targ);
x = [1 2];
y = [y_day(indsthis) y_night(indsthis)];
plot(x,y, '-k');
lt_plot_bar(x, mean(y), {'Errors', lt_sem(y)});
% -- pvals
p = signrank(y(:,1));
lt_plot_text(x(1), max(y(:,1)), ['p=' num2str(p)], 'r');
p = signrank(y(:,2));
lt_plot_text(x(2), max(y(:,2)), ['p=' num2str(p)], 'r');
p = signrank(y(:,2), y(:,1));
lt_plot_text(x(1)+0.5, max(y(:,1))+10, ['(vs)p=' num2str(p)], 'r');

% --- high targ learning
indsthis = x_targ>=median(x_targ);
x = [4 5];
y = [y_day(indsthis) y_night(indsthis)];
plot(x,y, '-k');
lt_plot_bar(x, mean(y), {'Errors', lt_sem(y)})
% -- pvals
p = signrank(y(:,1));
lt_plot_text(x(1), max(y(:,1)), ['p=' num2str(p)], 'r');
p = signrank(y(:,2));
lt_plot_text(x(2), max(y(:,2)), ['p=' num2str(p)], 'r');
p = signrank(y(:,2), y(:,1));
lt_plot_text(x(1)+0.5, max(y(:,1))+10, ['(vs)p=' num2str(p)], 'r');


% ------- stats, compare low and high learn expts
p = ranksum(y_day(x_targ<median(x_targ)), y_day(x_targ>=median(x_targ)));
lt_plot_pvalue(p, 'day(low vs hi)',1);
p = ranksum(y_night(x_targ<median(x_targ)), y_night(x_targ>=median(x_targ)));
lt_plot_pvalue(p, 'night(low vs hi)',1);


% ========================= SAME [using days as dat]
fieldthis = 'same';
lt_subplot(2,2,3); hold on;

% ----------
xlabel('target (total learn, morning-miornig');
ylabel('nontarg (day=r; night = b)');
title('each day = dat');
title(fieldthis);

nsyls = length(TMPSTRUCT.(fieldthis).learnDay);
x_targ = [];
y_day = [];
y_night = [];
ndays = size(TMPSTRUCT.(fieldthis).learnDay_allday,2); % i.e. num indiv days.
for j=1:nsyls
    
    % ---- skip if not positive learning at nontarg (i.e generalization)
    %     if TMPSTRUCT.(fieldthis).learnTot(j)<=0
    %         disp('/adadasdd')
    %         continue
    %     end
    
    % ---------------------- find corresponding target data
    bnumthis = TMPSTRUCT.(fieldthis).bnums(j);
    enumthis = TMPSTRUCT.(fieldthis).enums(j);
    
    indtarg = TMPSTRUCT.target.bnums==bnumthis & TMPSTRUCT.target.enums == enumthis;
    assert(sum(indtarg)==1, 'asdfasd');
    
    % --- extract total learning for corresponding targ
    x_targ = [x_targ; ones(ndays,1)*TMPSTRUCT.target.learnTot(indtarg)];
    
    % --- extract day and night learning for nontarg
    y_day = [y_day; TMPSTRUCT.(fieldthis).learnDay_allday(j,:)'];
    y_night = [y_night; TMPSTRUCT.(fieldthis).learnNight_allday(j,:)'];
end

plot(x_targ, y_day, 'or');
lt_regress(y_day, x_targ, 1, 0, 1, 1, 'r');
plot(x_targ, y_night, 'ob');
lt_regress(y_night, x_targ, 1, 0, 1, 1, 'b');

lt_plot_makesquare_plot45line(gca, 'k')

% ---------
lt_subplot(2,2,4); hold on;
xlabel('day-night(bars), LOWTARG LEARN[SETS OF BARS]');
ylabel('hz');
title('effect of total learn at target on day/night for nontarg (each syl one dat');
% --- low targ learning
indsthis = x_targ<median(x_targ);
x = [1 2];
y = [y_day(indsthis) y_night(indsthis)];
plot(x,y, '-k');
lt_plot_bar(x, mean(y), {'Errors', lt_sem(y)})
% -- pvals
p = signrank(y(:,1));
lt_plot_text(x(1), max(y(:,1)), ['p=' num2str(p)], 'r');
p = signrank(y(:,2));
lt_plot_text(x(2), max(y(:,2)), ['p=' num2str(p)], 'r');
p = signrank(y(:,2), y(:,1));
lt_plot_text(x(1)+0.5, max(y(:,1))+10, ['(vs)p=' num2str(p)], 'r');

% --- high targ learning
indsthis = x_targ>=median(x_targ);
x = [4 5];
y = [y_day(indsthis) y_night(indsthis)];
plot(x,y, '-k');
lt_plot_bar(x, mean(y), {'Errors', lt_sem(y)})
% -- pvals
p = signrank(y(:,1));
lt_plot_text(x(1), max(y(:,1)), ['p=' num2str(p)], 'r');
p = signrank(y(:,2));
lt_plot_text(x(2), max(y(:,2)), ['p=' num2str(p)], 'r');
p = signrank(y(:,2), y(:,1));
lt_plot_text(x(1)+0.5, max(y(:,1))+10, ['(vs)p=' num2str(p)], 'r');


% ------- stats, compare low and high learn expts
p = ranksum(y_day(x_targ<median(x_targ)), y_day(x_targ>=median(x_targ)));
lt_plot_pvalue(p, 'day(low vs hi)',1);
p = ranksum(y_night(x_targ<median(x_targ)), y_night(x_targ>=median(x_targ)));
lt_plot_pvalue(p, 'night(low vs hi)',1);




