function lt_seq_dep_pitch_ACROSSBIRDS_TbyT_datrange(TrialStruct, ParamsTrial)
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


for i=1:Numbirds
    Numexpt = length(TrialStruct.birds(i).exptnum);
    birdname = TrialStruct.birds(i).birdname;
    for ii=1:Numexpt
        exptname = TrialStruct.birds(i).exptnum(ii).exptname;
        
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
            end
            
            % ---------- plot lines for 7am and 9pm
            x = [ss+7/24 ss+21/24];
            line([x(1) x(1)], ylim, 'Color', 'k', 'LineStyle', '--');
            line([x(2) x(2)], ylim, 'Color', 'k', 'LineStyle', '--');
            
            
            % ======= plot
            %         x = ss; % syllable
            %         incr = []; % experiment counter + day
            %         y = count+incr;
        end
        set(gca, 'XTickLabel', {TrialStruct.birds(i).exptnum(ii).sylnum.syl});
        count = count+1;
    end
end