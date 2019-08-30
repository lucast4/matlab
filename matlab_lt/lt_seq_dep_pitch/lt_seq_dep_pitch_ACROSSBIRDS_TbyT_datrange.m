function lt_seq_dep_pitch_ACROSSBIRDS_TbyT_datrange(TrialStruct, ParamsTrial, ...
    BirdExptIncluded)
% ======= 
if ~exist('BirdExptIncluded', 'var') % this is list of birds to plot (2 cols,[bird expt])(get from later function)
    BirdExptIncluded = [];
end

%% plots time range for data across all syls, days, and expts


Numbirds = length(TrialStruct.birds);

%%
count = 1;
subplotrows=5;
subplotcols=1;
figcount=1;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

All_bnum = [];
All_enum = [];
All_sylnum =[];
All_day = [];
All_tminmax = [];


for i=1:Numbirds
    Numexpt = length(TrialStruct.birds(i).exptnum);
    birdname = TrialStruct.birds(i).birdname;
    for ii=1:Numexpt
        exptname = TrialStruct.birds(i).exptnum(ii).exptname;
        
        if ~isempty(BirdExptIncluded)
            
            if ~any(BirdExptIncluded(:,1)==i & BirdExptIncluded(:,2)==ii)
            disp('skipping bird, not in desired list BirdExptIncluded..');
                % thne skip this bird
                continue
            end
                
        end
        % ========================figure
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title([birdname '-' exptname]);
        xlabel('syl');
        ylabel('day (up later');
        
        Numsyls = length(TrialStruct.birds(i).exptnum(ii).sylnum);
        DayList = unique(floor(TrialStruct.birds(i).exptnum(ii).sylnum(1).Tvals));
        Numdays = length(DayList);
        pcolsyls = lt_make_plot_colors(Numsyls, 0,0);
        for ss =1:Numsyls
            
            tvalminmax = nan(Numdays, 2);
            
            for dd=1:Numdays
                daythis = DayList(dd);
                indtmp = floor(TrialStruct.birds(i).exptnum(ii).sylnum(ss).Tvals)==daythis;
                
                if ~any(indtmp)
                    continue
                end
                tvals = TrialStruct.birds(i).exptnum(ii).sylnum(ss).Tvals(indtmp);
                
                tvalminmax(dd,:) = [min(tvals) max(tvals)];
                
                % ========== plot
                
                trange = [min(tvals) max(tvals)]-floor(tvals);
                x = [ss + trange]; % sylklable location
                y = count + (1/(Numdays+1))*dd; % expt counter + daynum
                line(x, [y y], 'LineWidth', 3, 'Color', pcolsyls{ss});
                
                % ======================= SAVE FOR OTUPUT
                All_bnum = [All_bnum; i];
                All_enum = [All_enum; ii];
                All_sylnum =[All_sylnum; ss];
                All_day = [All_day; dd];
                All_tminmax = [All_tminmax; [min(tvals) max(tvals)]];

            end
            
            % ---------- plot lines for 7am and 9pm
            x = [ss+7/24 ss+21/24];
            line([x(1) x(1)], ylim, 'Color', 'k', 'LineStyle', '--');
            line([x(2) x(2)], ylim, 'Color', 'k', 'LineStyle', '--');
            
            
            % ======= plot
            %         x = ss; % syllable
            %         incr = []; % experiment counter + day
            %         y = count+incr;
            
            % =============== COLLECT ACROSS ALL DATA
            
        end
        set(gca, 'XTickLabel', {TrialStruct.birds(i).exptnum(ii).sylnum.syl});
        count = count+1;
    end
end

All_tminmax_orig = All_tminmax;

%% ====== PLOT DISTRUBTIONS OF ON/OFF OVER ALL DAYS/SYLS/EXPTS

lt_figure;

assert(all(floor(All_tminmax(:,1))==floor(All_tminmax(:,2))));

% ==== 1) reformat time
% get time from start of day
All_tminmax = All_tminmax - floor(All_tminmax(:,1));
% convert to hours
All_tminmax = All_tminmax.*24;

% === plot
lt_subplot(2,3,1); hold on;
title('all data');
[~, x] = lt_plot_histogram(All_tminmax(:,1), 6:0.5:23, 1, 1, 0.5, 1, 'b');
lt_plot_histogram(All_tminmax(:,2), x, 1, 1, 0.5, 1, 'r');

% === put lines for 7am and 9pm
line([7 7], ylim, 'Color', 'k');
line([22 22], ylim, 'Color', 'k');


% #################### ONE PT FOR EACH EXPT
[indsgrp, indsgrpU] = lt_tools_grp2idx({All_bnum, All_enum});

Y = [];
for i=1:length(indsgrpU)
   indsthis = indsgrp == indsgrpU(i);
   
   Y = [Y; mean(All_tminmax(indsthis,:),1)];
end

lt_subplot(2,3,2); hold on;
title('dat = expt');
[~, x] = lt_plot_histogram(Y(:,1), 6:0.5:23, 1, 1, 0.5, 1, 'b');
lt_plot_histogram(Y(:,2), x, 1, 1, 0.5, 1, 'r');

% === put lines for 7am and 9pm
line([7 7], ylim, 'Color', 'k');
line([22 22], ylim, 'Color', 'k');







