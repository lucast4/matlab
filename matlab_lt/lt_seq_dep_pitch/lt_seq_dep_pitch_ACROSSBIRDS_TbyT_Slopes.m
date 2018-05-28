function lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Slopes(TrialStruct, ParamsTrial, ...
    ignoreLMANexpt)
%% lt 5/25/18 -
% align all experiments by day and plot fitted slopes
normmethod = 'base_edges';
% base_edges: separately for onset and offset (accounts for circadian, using baseline regression);
% base_overall: one value, mean baseline across rends.

scalemethod = '';
% lastdaymean: mean of last day becomes 1

combineSylsInExpt=0; % if 0, datapt is syls, if 1, datapt is experiments
onlyKeepIfHaveNontarg=1; % throw out experiment if it doesn't have a sametype.

useregression=1; % then gets day on and off by fitting each day separately
% otherwise get first and last N renditions
Nrend = 8;

%%


%%

Numbirds = length(TrialStruct.birds);

%%

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
                
            end
            
            
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
            
            ffedges_allsyls = [ffedges_allsyls; mean(ffedges_allsyls(indstmp,:),1)];
            tedges_allsyls = [tedges_allsyls; mean(tedges_allsyls(indstmp,:),1)];
            istarg_allsyls = [istarg_allsyls; 0];
            issame_allsyls = [issame_allsyls; 1];
            
            indstoremove = [indstoremove; indstmp];
            
            % --- DIFF
            indstmp = find(istarg_allsyls==0 & issame_allsyls==0);
            
            ffedges_allsyls = [ffedges_allsyls; mean(ffedges_allsyls(indstmp,:),1)];
            tedges_allsyls = [tedges_allsyls; mean(tedges_allsyls(indstmp,:),1)];
            istarg_allsyls = [istarg_allsyls; 0];
            issame_allsyls = [issame_allsyls; 0];
            
            indstoremove = [indstoremove; indstmp];
            
            
            % ----- remove syls
            ffedges_allsyls(indstoremove,:) = [];
            tedges_allsyls(indstoremove,:) = [];
            istarg_allsyls(indstoremove,:) = [];
            issame_allsyls(indstoremove,:) = [];
            sylnames_allsyls = [];
        else
            
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



%%
%
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
%
% for i=1:Numbirds
%     Numexpt = length(TrialStruct.birds(i).exptnum);
%
%     for ii=1:Numexpt
%
%         % --------- ignore if lMAN?
%         if ignoreLMANexpt==1
%             isLMAN = TrialStruct.birds(i).exptnum(ii).LMANinactivated;
%             if isLMAN==1
%                 disp(['[is LMAN] skipping ' TrialStruct.birds(i).exptnum(ii).exptname]);
%                 continue
%             end
%         end
%
%         Numsyls = length(TrialStruct.birds(i).exptnum(ii).sylnum);
%         birdname =  TrialStruct.birds(i).birdname;
%         exptname = TrialStruct.birds(i).exptnum(ii).exptname;
%         targlearndir = TrialStruct.birds(i).exptnum(ii).targlearndir;
%
%
%         for ss =1:Numsyls
%
%             % ============== subplot for this syl
%             t = TrialStruct.birds(i).exptnum(ii).sylnum(ss).Tvals;
%             ff = TrialStruct.birds(i).exptnum(ii).sylnum(ss).FFvals;
%             istarg = TrialStruct.birds(i).exptnum(ii).sylnum(ss).INFO_istarget;
%             issame = TrialStruct.birds(i).exptnum(ii).sylnum(ss).INFO_similar;
%             sylname = TrialStruct.birds(i).exptnum(ii).sylnum(ss).syl;
%
%             % --------------- lines for base mean and 1std
%             basedays = TrialStruct.birds(i).exptnum(ii).BaseDays;
%             wndays = TrialStruct.birds(i).exptnum(ii).WNDays;
%             indsbase = t<basedays(end)+1;
%             ffmean_base = mean(ff(indsbase));
%             ffstd_base = std(ff(indsbase));
%
%
%             % ======================== FIT SIMPLE REGRESSION LINES FOR EACH
%             % DAY
%             daylist = [basedays wndays];
%             Tedges_byday = nan(1, 2*length(daylist)); % COLUMNS, day, split to ons and off (e.g. [day1on, day1off, day2on, day2off...];
%             FFedges_byday = nan(1, 2*length(daylist));
%             for dd = 1:length(daylist)
%                 daythis = daylist(dd);
%                 indthis = floor(t)==daythis;
%                 if ~any(indthis)
%                     % -- then no data for this day...
%                     tthis = nan;
%                     ff_fit = nan;
%                 else
%                     tthis = t(indthis);
%                     ffthis = ff(indthis);
%                     % -- sort
%                     [~, indtmp] = sort(tthis);
%                     tthis = tthis(indtmp);
%                     ffthis = ffthis(indtmp);
%
%                     % -- fit regression
%                     [b] =lt_regress(ffthis, tthis, 0);
%
%                     % -- collect FF at beginning and end of each day (from
%                     % regression fit)
%                     ff_fit = b(1) + b(2)*(tthis);
%                 end
%
%                 xtmp = [dd*2-1 dd*2];
%                 Tedges_byday(xtmp) = [tthis(1) tthis(end)];
%                 FFedges_byday(xtmp) = [ff_fit(1) ff_fit(end)];
%             end
%
%
%             % ================= NORMALIZE
%             if strcmp(normmethod, 'base_overall')
%                 % -------- 1) subtract one overall baseline mean
%                 FFedges_byday = FFedges_byday - ffmean_base;
%                 FFedges_byday = FFedges_byday*targlearndir;
%             elseif strcmp(normmethod, 'base_edges')
%                 % ---------- 2) subtract baseline mean separately for start
%                 % and end of day (account for circadian rhythm)
%                 numbasedays = length(basedays);
%                 % --- get mean baselines for early and late day
%                 baseON = mean(FFedges_byday(1:2:numbasedays*2));
%                 baseOFF = mean(FFedges_byday(2:2:numbasedays*2));
%                 % --- for each day subtract baseline for early and late
%                 FFedges_byday(1:2:end) = FFedges_byday(1:2:end) - baseON;
%                 FFedges_byday(2:2:end) = FFedges_byday(2:2:end) - baseOFF;
%                 % --- make up direction of learning
%                 FFedges_byday = FFedges_byday*targlearndir;
%
%             end
%
%
%             %             % =================== SCALE
%             %             if strcmp(scalemethod, 'lastdaymean')
%             %                 assert(1==2, 'PRIOBLME: need to have nontarg scaled to target...');
%             %                 assert(strcmp(normmethod, 'base_overall') | strcmp(normmethod, 'base_edges'), ...
%             %                     'to scale must norm relative to baseline');
%             %                 % --------- scale so that end of last day is 1 and baseline
%             %                 % is 0
%             %                 maxff = mean(FFedges_byday(end-1:end));
%             %                 FFedges_byday = FFedges_byday./abs(maxff);
%             %                 % --- take abs of max to make sure the sign of the last day still holds (i.e.
%             %                 % did not take negative laerning and then flip it)
%             %
%             %             end
%             %
%
%             % ===================== CLEANING UP EXPT
%
%
%          % ==================== OUTPUTS
%         DATSTRUCT.AllBirdnum = [DATSTRUCT.AllBirdnum; i];
%         DATSTRUCT.AllExptnum =[DATSTRUCT.AllExptnum; ii];
%         DATSTRUCT.AllDayList = [DATSTRUCT.AllDayList; daylist];
%         DATSTRUCT.AllBasedays = [DATSTRUCT.AllBasedays; basedays];
%
%         DATSTRUCT.AllSylname = [DATSTRUCT.AllSylname; sylname];
%         DATSTRUCT.AllIsTarg = [DATSTRUCT.AllIsTarg; istarg];
%         DATSTRUCT.AllIsSame = [DATSTRUCT.AllIsSame; issame];
%
%
%         DATSTRUCT.AllTedges = [DATSTRUCT.AllTedges; Tedges_byday];
%         DATSTRUCT.AllFFedges = [DATSTRUCT.AllFFedges; FFedges_byday];
%
%         end
%
%
%
%     end
% end
%
%
% maxbirds = max(DATSTRUCT.AllBirdnum);
% maxexpts = max(DATSTRUCT.AllExptnum);

%% ==================== COMBINE ALL SYLS IN AN EXPERIMENT
% if combineSylsInExpt==1
%
%     DATSTRUCT2.AllBirdnum = [DATSTRUCT.AllBirdnum; i];
%     DATSTRUCT2.AllExptnum =[DATSTRUCT.AllExptnum; ii];
%     DATSTRUCT2.AllSylname = [DATSTRUCT.AllSylname; sylname];
%
%     DATSTRUCT.AllIsTarg = [DATSTRUCT.AllIsTarg; istarg];
%     DATSTRUCT.AllIsSame = [DATSTRUCT.AllIsSame; issame];
%
%     DATSTRUCT.AllDayList = [DATSTRUCT.AllDayList; daylist];
%     DATSTRUCT.AllTedges = [DATSTRUCT.AllTedges; Tedges_byday];
%     DATSTRUCT.AllFFedges = [DATSTRUCT.AllFFedges; FFedges_byday];
%
%     DATSTRUCT.AllBasedays = [DATSTRUCT.AllBasedays; basedays];
%
%     % ========== go thru all expts
%     for i=1:maxbirds
%         for ii=1:maxexpts
%             inds = DATSTRUCT.AllBirdnum==i & DATSTRUCT.AllExptnum==ii;
%             if ~any(inds)
%                 continue
%             end
%
%             % =========== TARG
%
%
%             % =========== same
%             inds = DATSTRUCT.AllBirdnum==i & DATSTRUCT.AllExptnum==ii & ...
%                 DATSTRUCT.AllIsTarg==0 & DATSTRUCT.AllIsSame==1;
%
%             ffedges = mean(DATSTRUCT.AllFFedges(inds,:),1);
%             tedges = mean(DATSTRUCT.AllTedges(inds,:),1);
%
%             DATSTRUCT.AllBirdnum = [DATSTRUCT.AllBirdnum; i];
%             DATSTRUCT.AllExptnum =[DATSTRUCT.AllExptnum; ii];
%             DATSTRUCT.AllSylname = [DATSTRUCT.AllSylname; sylname];
%
%             DATSTRUCT.AllIsTarg = [DATSTRUCT.AllIsTarg; istarg];
%             DATSTRUCT.AllIsSame = [DATSTRUCT.AllIsSame; issame];
%
%             DATSTRUCT.AllDayList = [DATSTRUCT.AllDayList; daylist];
%             DATSTRUCT.AllTedges = [DATSTRUCT.AllTedges; Tedges_byday];
%             DATSTRUCT.AllFFedges = [DATSTRUCT.AllFFedges; FFedges_byday];
%
%             DATSTRUCT.AllBasedays = [DATSTRUCT.AllBasedays; basedays];
%
%             % ========== DIFF
%
%         end
%     end
%
%
%     % ==================== OUTPUTS
%     %             DATSTRUCT.AllBirdnum = [DATSTRUCT.AllBirdnum; i];
%     %             DATSTRUCT.AllExptnum =[DATSTRUCT.AllExptnum; ii];
%     %             DATSTRUCT.AllSylname = [DATSTRUCT.AllSylname; sylname];
%     %
%     %             DATSTRUCT.AllIsTarg = [DATSTRUCT.AllIsTarg; istarg];
%     %             DATSTRUCT.AllIsSame = [DATSTRUCT.AllIsSame; issame];
%     %
%     %             DATSTRUCT.AllDayList = [DATSTRUCT.AllDayList; daylist];
%     %             DATSTRUCT.AllTedges = [DATSTRUCT.AllTedges; Tedges_byday];
%     %             DATSTRUCT.AllFFedges = [DATSTRUCT.AllFFedges; FFedges_byday];
%     %
%     %             DATSTRUCT.AllBasedays = [DATSTRUCT.AllBasedays; basedays];
% end

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
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end),1,2);

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
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end),1,2);

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
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end),1,2);

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
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end),1,2);
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
ffdiff = diff(DATSTRUCT.AllFFedges(inds,ind_wnon:end),1,2);
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
lt_subplot(2,2,1); hold on;
title('target');
xlabel('over day');
ylabel('over night');
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
plot(Y(:,1), Y(:,2), 'ok');
lt_plot_makesquare_plot45line(gca, 'b');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SAME
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

plot(Y(:,1), Y(:,2), 'ok');
lt_plot_makesquare_plot45line(gca, 'b');


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DIFF
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

plot(Y(:,1), Y(:,2), 'ok');
lt_plot_makesquare_plot45line(gca, 'b');

%% IS GENERALIZATION GREATER IN MORNING OR EVENING?
% =========== FOR EVERY SYL, ASK WHETHER OVER DAY OR OVER NGIHT CONTRIBUTES
% MORE TO EXPRESSION OF LEARNING (totla learning is last morning minus
% first morning. over day and over night are sum over all days)
lt_figure; hold on;

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


% =================== PLOT
% ---- overday
plot(diff_mornings, Y(:,1), 'or');
% ---- overnight
plot(diff_mornings, Y(:,2), 'ob');
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


% =================== PLOT
% ---- overday
plot(diff_mornings, Y(:,1), 'or');
% ---- overnight
plot(diff_mornings, Y(:,2), 'ob');

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


% =================== PLOT
% ---- overday
plot(diff_mornings, Y(:,1), 'or');
% ---- overnight
plot(diff_mornings, Y(:,2), 'ob');
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


