function lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Raw(TrialStruct, ParamsTrial, ...
    ignoreDiffType)
%% plots all experiments, raw FF

%%

Numbirds = length(TrialStruct.birds);

%%
for i=1:Numbirds
   Numexpt = length(TrialStruct.birds(i).exptnum);
   
   for ii=1:Numexpt
      
       % ========================figure
       if ignoreDiffType==1
        subplotrows=4;
        subplotcols=1;
       else
       subplotrows=5;
        subplotcols=2;
       end
              figcount=1;
        fignums_alreadyused=[];
        hfigs=[];
        hsplots = [];

        Numsyls = length(TrialStruct.birds(i).exptnum(ii).sylnum);
        
        for ss =1:Numsyls
           
            % ============== subplot for this syl
            t = TrialStruct.birds(i).exptnum(ii).sylnum(ss).Tvals;
            ff = TrialStruct.birds(i).exptnum(ii).sylnum(ss).FFvals;
            istarg = TrialStruct.birds(i).exptnum(ii).sylnum(ss).INFO_istarget;
            issame = TrialStruct.birds(i).exptnum(ii).sylnum(ss).INFO_similar;
            sylname = TrialStruct.birds(i).exptnum(ii).sylnum(ss).syl;
            birdname =  TrialStruct.birds(i).birdname;
            exptname = TrialStruct.birds(i).exptnum(ii).exptname;
            
            % ------------ ignore diff type?
            if ignoreDiffType==1
                if issame==0
                    continue
                end
            end                    
                
            
            if istarg==1 
                pcol ='k';
            elseif istarg==0 & issame==1
                pcol = 'b';
            elseif istarg==0 & issame==0
                pcol = 'r';
            end
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            if istarg==1
            title([sylname ',' birdname '-' exptname]);    
            else
            title([sylname]);
            end
            ylabel('ff');
            
            
            % --------------- lines for base mean and 1std
            basedays = TrialStruct.birds(i).exptnum(ii).BaseDays;
            indsbase = t<basedays(end)+1;
            ffmean_base = mean(ff(indsbase));
            ffstd_base = std(ff(indsbase));
            shadedErrorBar([t(1) t(end)], [ffmean_base ffmean_base], ...
                [ffstd_base ffstd_base], {'Color', pcol},1);
            
            % ----------------- PLOT DATAPOINTS
            plot(t, ff, 'x', 'Color', pcol);
            
            
            % ======================== FIT SIMPLE REGRESSION LINES FOR EACH
            % DAY
            daylist = unique(floor(t));
            for daythis = daylist'
                indthis = floor(t)==daythis;
                tthis = t(indthis);
                ffthis = ff(indthis);
                % -- sort
                [~, indtmp] = sort(tthis);
               tthis = tthis(indtmp);
               ffthis = ffthis(indtmp);
                % -- fit regression
            [b] =lt_regress(ffthis, tthis, 0);
            ff_fit = b(1) + b(2)*(tthis);
            line([tthis(1) tthis(end)], [ff_fit(1) ff_fit(end)], 'Color', pcol, ...
                'LineWidth', 2);
            end
            
            % ------------------- put lines for expt onset
            baseend = TrialStruct.birds(i).exptnum(ii).BaseDays(end);
            line([baseend+1 baseend+1], ylim, 'Color', 'k','LineStyle', '--');            
                
            
        end
        linkaxes(hsplots, 'x');
   end
end
